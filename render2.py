import networkx as nx
import numpy as np
from collections import defaultdict
from itertools import combinations, pairwise, product
import math
from enum import IntEnum
import drawsvg as svg

from util.collections import merge_alternating
from util.geometry import (
    get_side,
    dist_euclidean,
    get_closest_point,
    centroid,
    logical_coords_to_physical,
    are_faces_adjacent,
    interpolate_biarcs,
    get_angle,
    offset_edge,
    is_point_inside_circle,
    get_segment_circle_intersection,
)
from util.graph import (
    incident_edges,
    get_longest_simple_paths,
    get_closest_pair,
    get_shortest_path_between,
)


class BendPenalty(IntEnum):
    ONE_EIGHTY = 1
    ONE_THIRTY_FIVE = 2
    NINETY = 4
    FORTY_FIVE = 8


class EdgeType(IntEnum):
    PHYSICAL = -1
    CENTER = 0  # center to center edge
    DRAW = 1  # physical edge that is actually drawn


class NodeType(IntEnum):
    CENTER = 0
    PORT = 1
    CORNER = 2


DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 100,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
}


def add_ports_to_hex_node(G, node, data, side_length=0.25):
    hex_corners = [
        (0, side_length),  # N
        (side_length * math.sqrt(3) / 2, side_length / 2),  # NE
        (side_length * math.sqrt(3) / 2, -side_length / 2),  # SE
        (0, -side_length),  # S
        (-side_length * math.sqrt(3) / 2, -side_length / 2),  # SW
        (-side_length * math.sqrt(3) / 2, side_length / 2),  # NW
    ]
    pos = data["pos"]
    hex_corners = [(x + pos[0], y + pos[1]) for x, y in hex_corners]
    dirs = ["n", "ne", "se", "s", "sw", "nw"]
    hex_corners = list(zip(dirs, hex_corners))
    for dir, corner in hex_corners:
        G.add_node(corner, node=NodeType.CORNER, belongs_to=node, pos=corner, port=dir)

    sides = pairwise(hex_corners + [hex_corners[0]])
    hex_sides = []
    for c1, c2 in sides:
        d1, p1 = c1
        d2, p2 = c2
        n = centroid([p1, p2])
        hex_sides.append(n)
        G.add_node(n, node=NodeType.PORT, belongs_to=node, pos=n, port=(d1, d2))

    penalties = [
        BendPenalty.NINETY,
        BendPenalty.ONE_THIRTY_FIVE,
        BendPenalty.ONE_EIGHTY,
        BendPenalty.ONE_THIRTY_FIVE,
        BendPenalty.NINETY,
    ]
    for i in range(len(hex_sides)):
        p = 0
        n1 = hex_sides[i]
        for j in range(i + 1, len(hex_sides)):
            n2 = hex_sides[j]
            G.add_edge(
                n1, n2, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL, weight=penalties[p]
            )
            p += 1

    # connect center to sides
    for n in hex_sides:
        G.add_edge(
            node,
            n,
            EdgeType.PHYSICAL,
            edge=EdgeType.PHYSICAL,
            weight=0,
        )

    # find neighbors and replace center edges with physical edges to ports
    neighbors = list(nx.neighbors(G, node))
    for neighbor in neighbors:
        G.remove_edge(node, neighbor)
        # use physically closest port for any neighbor
        neighbor_pos = G.nodes[neighbor]["pos"]
        port = get_closest_point(neighbor_pos, hex_sides)

        weight = dist_euclidean(port, neighbor_pos)
        G.add_edge(
            neighbor, port, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL, weight=weight
        )
    return G


def make_hex_graph(m, n):
    # each hex is 17 nodes and...  136 edges :zipface:
    # = 1 center, 6 corners and 6 sides.
    # fully connect corners, sides and center
    # neighbor edges run between adjacent sides of hexes

    # start with square graph
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.CENTER
    for logpos, n in G.nodes(data=True):
        x, y = logpos
        n["node"] = NodeType.CENTER
        n["belongs_to"] = logpos
        n["logpos"] = logpos
        n["pos"] = logical_coords_to_physical(x, y, "hex")
        # TODO support glyph spacing in x and y separately

    G_ = nx.MultiGraph(incoming_graph_data=G)

    for logpos, hex_center in G.nodes(data=True):
        x, y = logpos

        # add diagonal edges
        if y % 2 != 0 and y > 0 and x < m - 1:
            # tilting to right in odd rows
            G_.add_edge((x, y), (x + 1, y - 1), EdgeType.CENTER, edge=EdgeType.CENTER)
        if y % 2 == 0 and y > 0 and x > 0:
            # tilting to left in even rows
            G_.add_edge((x, y), (x - 1, y - 1), EdgeType.CENTER, edge=EdgeType.CENTER)

    for logpos, hex_center in G.nodes(data=True):
        G_ = add_ports_to_hex_node(G_, logpos, hex_center, side_length=0.25)

    for logpos, n in G.nodes(data=True):
        x, y = logpos
        if "pos" not in n:
            n["pos"] = logpos

    return G_


def get_routing_graph(lattice_type, lattice_size):
    """Returns a networkx graph. It contains two types of nodes and two types of edges.
    Nodes can be 'glyph' nodes, which correspond to spots where glyphs will be placed.
    Nodes can also be 'anchor' nodes, which corresponds to anchor points in the margins of the layout along which we trace lines and anchor polygons. We place them in the center of 'glyph' node faces.
    Edges can be 'neighbor' edges, which corresponds to direct neighbor relations of glyphs, e.g., in a sqr grid most nodes have 4 neighbors.
    Edges can also be 'anchor' edges. They connect both 'glyph' and 'anchor' nodes. In a hex lattice, faces (polygons between 'glyph' nodes) are triangles. So each 'anchor' node has 6 'anchor' edges incident: 3 for neighboring anchors and 3 for glyphs on the boundary of the face.
    """
    m, n = lattice_size
    G = None
    match lattice_type:
        case "hex":
            G = make_hex_graph(m, n)
        case "sqr":
            # G = make_sqr_graph(m, n)
            pass
        case _:
            raise Exception(f"unknown lattice type {lattice_type}")

    return G


def geometrize(instance, M):
    one_unit_px = DEFAULT_PARAMS["unit_size_in_px"]
    margins = np.array((0.5, 0.5)) * one_unit_px
    factor = one_unit_px
    geometries = []
    mx, my = margins

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        M.nodes[i]["pos"] = (x * factor + mx, -y * factor + my)

    # glyph nodes
    for i, n in M.nodes(data=True):
        if n["node"] == NodeType.CENTER:
            corners = [
                m
                for m in nx.subgraph_view(
                    M,
                    filter_node=lambda p: G.nodes[p]["node"] == NodeType.CORNER
                    and G.nodes[p]["belongs_to"] == i,
                )
            ]
            for corner in corners:
                cx, cy = M.nodes[corner]["pos"]
                """
                color = {
                    "n": "red",
                    "ne": "orange",
                    "se": "green",
                    "s": "blue",
                    "sw": "magenta",
                    "nw": "pink",
                }"""
                geometries.append(svg.Circle(cx=cx, cy=cy, r=2, fill="black", z=1000))
            xs, ys = zip(*map(lambda m: M.nodes[m]["pos"], corners))
            # so if we wanted to render to something else than svg, we could
            # replace the drawsvg elements with some dicts and make drawsvg
            # elements in draw_svg function. convert to something else in
            # draw_smth_else.
            # for now it saves time to not do that
            glyph = svg.Lines(
                *merge_alternating(xs, ys),
                close=True,
                fill="transparent",
                stroke="black",
            )
            # geometries.append(glyph)

    for u, v, k in M.edges(keys=True):
        x1, y1 = M.nodes[u]["pos"]
        x2, y2 = M.nodes[v]["pos"]
        if M.edges[(u, v, k)]["edge"] == EdgeType.PHYSICAL:
            color = {
                BendPenalty.NINETY: "red",
                BendPenalty.ONE_THIRTY_FIVE: "orange",
                BendPenalty.ONE_EIGHTY: "green",
            }
            w = M.edges[(u, v, k)].get("weight", 0)
            geometries.append(
                svg.Line(
                    x1,
                    y1,
                    x2,
                    y2,
                    z=w,
                    stroke=color.get(w, "gray"),
                    stroke_width=1,
                )
            )

    # uncomment to draw edge-bundled lines
    """
    set_colors = ["red", "blue", "orange", "green", "magenta"]
    uniq_edges = list(set(M.edges()))
    for u, v in uniq_edges:
        ks = [k for w, x, k in M.edges(keys=True) if (w, x) == (u, v)]
        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        # print normal of edge vector just for debugging
        # cx1, cy1 = centroid([src, tgt])
        # cx2, cy2 = offset_point((cx1, cy1), edge_angle + math.pi / 2, 10)
        # geometries.append(
        #    svg.Line(cx1, cy1, cx2, cy2, stroke="lightgray", stroke_width=1)
        # )

        paths_at_edge = M.edges[(u, v, ks[0])]["oeb_order"][(u, v)]
        offset = factor / 30  # TODO idk some pixels
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -offset
        for i, path_id in enumerate(paths_at_edge):
            set_id, _ = path_id.split("-")  # TODO meh
            set_idx = instance["set_ftb_order"].index(set_id)
            offset_dir = math.pi / 2
            offset_length = centering_offset + i * offset
            o_u, o_v = offset_edge((src, tgt), edge_angle - offset_dir, offset_length)

            M.edges[(u, v, set_idx + 1)]["edge_pos"] = {
                (u, v): (o_u, o_v),
                (v, u): (o_v, o_u),
            }

    paths = list(set([d["path_id"] for u, v, k, d in M.edges(keys=True, data=True)]))
    for path_id in paths:
        set_id, _ = path_id.split("-")  # TODO meh
        set_idx = instance["set_ftb_order"].index(set_id)
        M_ = nx.subgraph_view(
            M, filter_edge=lambda u, v, k: M.edges[(u, v, k)]["path_id"] == path_id
        )
        endpoints = [n for n in M_.nodes() if M_.degree(n) == 1]
        path = nx.shortest_path(M_, endpoints[0], endpoints[1])
        # idk, like keyframe... find better word
        keypoints = [
            # (type of connection to following point, keypoint coordinates)
        ]

        for i, pair in enumerate(pairwise(path)):
            u, v = pair
            ux, uy = M_.nodes[u]["pos"]
            vx, vy = M_.nodes[v]["pos"]

            # at each node, we find segments from the last key point to the next node's hub circle
            hub_radius = 1 / 16 * factor
            uhub_center = (ux, uy)
            uhub_circle = (uhub_center, hub_radius)

            vhub_center = (vx, vy)
            vhub_circle = (vhub_center, hub_radius)

            edge = M_.edges[(u, v, set_idx + 1)]["edge_pos"][(u, v)]
            a, b = edge

            if is_point_inside_circle(a, uhub_circle):
                # if start point of edge is inside hub circle,
                # we want a straight line from the hub circle circumference to the following point
                keypoints.append(
                    ("L", get_segment_circle_intersection((a, b), uhub_circle))
                )
            else:
                # start point of edge is not inside the hub circle
                # so we need an arc from the hub circle circumference to the edge start
                keypoints.append(
                    (
                        "B",
                        get_segment_circle_intersection((uhub_center, a), uhub_circle),
                    ),
                )
                # and then a straight line
                keypoints.append(("L", a))

            if is_point_inside_circle(b, vhub_circle):
                keypoints.append(
                    ("B", get_segment_circle_intersection((a, b), vhub_circle))
                )
            else:
                # edge to hub
                keypoints.append(("B", b))
                keypoints.append(
                    (
                        "B" if i + 1 < len(path) - 1 else "L",
                        get_segment_circle_intersection((vhub_center, b), vhub_circle),
                    ),
                )

        kwargs = {
            "close": False,
            "stroke_width": 2,
            "fill": "none",
            "stroke": set_colors[set_idx],
            "z": set_idx,
        }
        line = interpolate_biarcs(keypoints, **kwargs)
        geometries.append(line)
    """

    if False:
        hubs = [n for n in M.nodes() if M.degree[n] > 0]
        for hub in hubs:
            cx, cy = M.nodes[hub]["pos"]
            r = 1 / 16 * factor
            geometries.append(
                svg.Circle(cx, cy, r, fill="none", stroke="gray", stroke_width=1)
            )
    return geometries


def draw_svg(geometries):
    d = svg.Drawing(2000, 2000, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


INSTANCE = {
    "lattice_type": "hex",
    "glyph_ids": [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
    ],  # list of strings (could be file urls)
    # A B C
    #  D E F
    # G H I
    "glyph_positions": [
        (0, 0),
        (1, 0),
        (2, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (0, 2),
        (1, 2),
        (2, 2),
    ],  # a panda df with two columns (x and y) corresponding to logical grid position
    "set_system": {
        "set0": ["A", "B", "H", "D", "E"],
        "set1": ["A", "I"],
        "set2": ["G", "C", "D", "E"],
        "set3": ["A", "F", "D", "E"],
        "set4": ["A", "I", "D", "E"],
    },
    "set_ftb_order": [
        "set4",
        "set3",
        "set2",
        "set1",
        "set0",
    ],  # a list of set ids that defines front to back ordering (index = 0 is most front)
}


if __name__ == "__main__":
    m = 5
    n = 5
    lattice_type = INSTANCE["lattice_type"]
    G = get_routing_graph(lattice_type, (m, n))
    # G = embed_to_routing_graph(INSTANCE, G)
    # G = render_line(INSTANCE, G, DEFAULT_PARAMS)
    geometries = geometrize(INSTANCE, G)
    img = draw_svg(geometries)

    with open("drawing.svg", "w") as f:
        f.write(img)
        f.flush()
