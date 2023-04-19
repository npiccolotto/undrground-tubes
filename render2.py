import math
import time
from collections import defaultdict
from enum import IntEnum
from itertools import combinations, pairwise, product

import drawsvg as svg
import networkx as nx
import numpy as np

from util.perf import timing
from util.collections import (
    get_elements_in_same_lists,
    merge_alternating,
    invert_dict_of_lists,
)
from util.geometry import (
    are_faces_adjacent,
    centroid,
    dist_euclidean,
    get_angle,
    get_closest_point,
    get_segment_circle_intersection,
    get_side,
    interpolate_biarcs,
    is_point_inside_circle,
    logical_coords_to_physical,
    offset_edge,
)
from util.graph import (
    approximate_steiner_tree,
    are_node_sets_connected,
    get_closest_pair,
    get_longest_simple_paths,
    get_shortest_path_between,
    incident_edges,
)


class BendPenalty(IntEnum):
    ONE_EIGHTY = 0.5
    ONE_THIRTY_FIVE = 1
    NINETY = 1.5
    FORTY_FIVE = 2


class EdgeType(IntEnum):
    SUPPORT = -2  # edge in hypergraph support
    PHYSICAL = -1  # physical edge that is actually drawn
    CENTER = 0  # center to center edge


class NodeType(IntEnum):
    CENTER = 0
    PORT = 1


DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 100,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
}


def add_ports_to_hex_node(G, node, data, side_length=0.25):
    hex_ports = [
        (0, side_length),  # N
        (side_length * math.sqrt(3) / 2, side_length / 2),  # NE
        (side_length * math.sqrt(3) / 2, -side_length / 2),  # SE
        (0, -side_length),  # S
        (-side_length * math.sqrt(3) / 2, -side_length / 2),  # SW
        (-side_length * math.sqrt(3) / 2, side_length / 2),  # NW
    ]
    pos = data["pos"]
    hex_ports = [(x + pos[0], y + pos[1]) for x, y in hex_ports]
    dirs = ["n", "ne", "se", "s", "sw", "nw"]
    hex_ports = list(zip(dirs, hex_ports))
    for dir, corner in hex_ports:
        G.add_node(corner, node=NodeType.PORT, belongs_to=node, pos=corner, port=dir)

    sides = pairwise(hex_ports + [hex_ports[0]])
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

    for logpos, _ in G.nodes(data=True):
        x, y = logpos

        # add diagonal edges
        if y % 2 != 0 and y > 0 and x < m - 1:
            # tilting to right in odd rows
            right_diagonal = (x + 1, y - 1)
            weight = dist_euclidean(
                G_.nodes[logpos]["pos"], G.nodes[right_diagonal]["pos"]
            )
            G_.add_edge(
                logpos,
                right_diagonal,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
                weight=weight,
            )
        if y % 2 == 0 and y > 0 and x > 0:
            # tilting to left in even rows
            left_diagonal = (x - 1, y - 1)
            weight = dist_euclidean(
                G_.nodes[logpos]["pos"], G.nodes[left_diagonal]["pos"]
            )
            G_.add_edge(
                logpos,
                left_diagonal,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
                weight=weight,
            )

    for logpos, hex_center in G.nodes(data=True):
        G_ = add_ports_to_hex_node(G_, logpos, hex_center, side_length=0.25)

    return G_


def add_ports_to_sqr_node(G, node, data, side_length=0.25):
    sqr_corners = [
        (0, side_length / 2),  # N
        (side_length / 2, side_length / 2),  # NE
        (side_length / 2, 0),  # E
        (side_length / 2, -side_length / 2),  # SE
        (0, -side_length / 2),  # S
        (-side_length / 2, -side_length / 2),  # SW
        (-side_length / 2, 0),  # W
        (-side_length / 2, side_length / 2),  # NW
    ]
    pos = data["pos"]
    sqr_corners = [(x + pos[0], y + pos[1]) for x, y in sqr_corners]
    dirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"]
    sqr_corners_dirs = list(zip(dirs, sqr_corners))
    for dir, corner in sqr_corners_dirs:
        G.add_node(corner, node=NodeType.PORT, belongs_to=node, pos=corner, port=dir)

    penalties = [
        BendPenalty.FORTY_FIVE,
        BendPenalty.NINETY,
        BendPenalty.ONE_THIRTY_FIVE,
        BendPenalty.ONE_EIGHTY,
        BendPenalty.ONE_THIRTY_FIVE,
        BendPenalty.NINETY,
        BendPenalty.FORTY_FIVE,
    ]
    for i in range(len(sqr_corners)):
        p = 0
        n1 = sqr_corners[i]
        for j in range(i + 1, len(sqr_corners)):
            n2 = sqr_corners[j]

            G.add_edge(
                n1, n2, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL, weight=penalties[p]
            )
            p += 1

    # connect center to ports
    for n in sqr_corners:
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
        # use physically closest port for any neighbor
        neighbor_pos = G.nodes[neighbor]["pos"]
        port = get_closest_point(neighbor_pos, sqr_corners)

        weight = dist_euclidean(port, neighbor_pos)
        G.add_edge(neighbor, port, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL, weight=1)
    return G


def make_sqr_graph(m, n):
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.CENTER
    for logpos, n in G.nodes(data=True):
        x, y = logpos
        n["node"] = NodeType.CENTER
        n["belongs_to"] = logpos
        n["logpos"] = logpos
        n["pos"] = logical_coords_to_physical(x, y, "sqr")

    G_ = nx.MultiGraph(incoming_graph_data=G)

    for logpos, _ in G.nodes(data=True):
        x, y = logpos
        if y > 0 and x < m - 1:
            # tilt to right
            right_diagonal = (x + 1, y - 1)
            weight = dist_euclidean(
                G_.nodes[logpos]["pos"], G.nodes[right_diagonal]["pos"]
            )
            G_.add_edge(
                logpos,
                right_diagonal,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
                weight=1,
            )
        if y > 0 and x > 0:
            # tilt to left
            left_diagonal = (x - 1, y - 1)
            weight = dist_euclidean(
                G_.nodes[logpos]["pos"], G.nodes[left_diagonal]["pos"]
            )
            G_.add_edge(
                logpos,
                left_diagonal,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
                weight=1,
            )

    for logpos, center in G.nodes(data=True):
        G_ = add_ports_to_sqr_node(G_, logpos, center, side_length=0.25)

    return G_


def get_routing_graph(lattice_type, lattice_size):
    m, n = lattice_size
    G = None
    match lattice_type:
        case "hex":
            G = make_hex_graph(m, n)
        case "sqr":
            G = make_sqr_graph(m, n)
        case _:
            raise Exception(f"unknown lattice type {lattice_type}")
    return G


def add_glyphs_to_nodes(instance, G):
    for n in G.nodes():
        G.nodes[n]["occupied"] = False
    for i, glyph in enumerate(instance["glyph_ids"]):
        logpos = instance["glyph_positions"][i]
        if G.nodes[logpos]["node"] != NodeType.CENTER:
            raise Exception("node to position glyph on is somehow not a glyph center")
        G.nodes[logpos]["occupied"] = True
        G.nodes[logpos]["glyph"] = glyph
    return G


def block_center_edges_except_to(G, open_nodes):
    print("open nodes", open_nodes)
    for u, v in G.edges():
        un = G.nodes[u]["node"]
        vn = G.nodes[v]["node"]
        if un == NodeType.CENTER or vn == NodeType.CENTER:
            center = u if un == NodeType.CENTER else v
            if G.nodes[center]["occupied"] and center not in open_nodes:
                print(f"blocking {u}-{v} ({un},{vn})")
                G.edges[(u, v)]["weight"] = math.inf
    return G


def unblock_center_edges(G):
    for u, v in G.edges():
        un = G.nodes[u]["node"]
        vn = G.nodes[v]["node"]
        if un == NodeType.CENTER or vn == NodeType.CENTER:
            G.edges[(u, v)]["weight"] = 0
    return G


def route_set_lines(instance, G):
    # steps:
    # 1. identify intersection groups: elements that belong to the same sets
    intersection_groups = get_elements_in_same_lists(instance["set_system"])
    sets_per_element = invert_dict_of_lists(instance["set_system"])

    sorted_intersection_groups = sorted(
        zip(intersection_groups.keys(), intersection_groups.values()),
        key=lambda x: len(x[1]),
        reverse=True,
    )
    processed_elements = set()

    G_ = nx.Graph()
    G_.add_nodes_from([(n, d) for n, d in G.nodes(data=True)])
    G_.add_edges_from(
        [
            (u, v, d)
            for u, v, k, d in G.edges(keys=True, data=True)
            if d["edge"] == EdgeType.PHYSICAL
        ]
    )

    # 2. then process from biggest to smallest group. for each:
    for elements, sets in sorted_intersection_groups:
        print(sets, elements)
        S = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER
            and d["occupied"] == True
            and d["glyph"] in elements
        ]
        # determine union of all processed elements sharing any set with the current set of sets o_O
        shared_processed_elements = [
            el
            for el in processed_elements
            if len(set(sets_per_element[el]).intersection(sets)) > 0
        ]

        # close edges to all center nodes not belonging to elements
        G_ = block_center_edges_except_to(G_, S)
        # 3. determine approximated steiner tree acc. to kou1981 / wu1986.
        with timing("steiner tree"):
            steiner = approximate_steiner_tree(G_, S)

        # add those edges to support
        for u, v in steiner.edges():
            if not G.has_edge(u, v, EdgeType.SUPPORT):
                G.add_edge(
                    u, v, EdgeType.SUPPORT, edge=EdgeType.SUPPORT, sets=set(sets)
                )
            else:
                G.edges[(u, v, EdgeType.SUPPORT)]["sets"] = set(sets).union(
                    G.edges[(u, v, EdgeType.SUPPORT)]["sets"]
                )

        support = nx.subgraph_view(G, filter_edge=lambda u, v, k: k == EdgeType.SUPPORT)

        # if support is not connected but should be, find cut and fix it by connecting via shortest path
        if len(shared_processed_elements) > 0:
            nodes_of_shared_processed_elements = [
                n
                for n, d in G_.nodes(data=True)
                if d["node"] == NodeType.CENTER
                and d["occupied"]
                and d["glyph"] in shared_processed_elements
            ]
            if not are_node_sets_connected(
                support, S, nodes_of_shared_processed_elements
            ):
                _, shortest_path_edgelist = get_shortest_path_between(
                    G, S, nodes_of_shared_processed_elements
                )
                for u, v in shortest_path_edgelist:
                    G.add_edge(
                        u, v, EdgeType.SUPPORT, edge=EdgeType.SUPPORT, sets=set(sets)
                    )

        processed_elements = processed_elements.union(set(S))

        G_ = unblock_center_edges(G_)

    # 4. update host graph: lighten weight on edges in F as suggested in castermans2019, penalize crossing dual edges as suggested by bast2020
    # 5. stop when all elements have been processed
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
        nx, ny = (x * factor + mx, -y * factor + my)
        M.nodes[i]["pos"] = (nx, ny)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            geometries.append(svg.Circle(cx=nx, cy=ny, r=one_unit_px / 4))

    for u, v, k in M.edges(keys=True):
        x1, y1 = M.nodes[u]["pos"]
        x2, y2 = M.nodes[v]["pos"]
        if M.edges[(u, v, k)]["edge"] == EdgeType.PHYSICAL:
            w = M.edges[(u, v, k)].get("weight", 1)
            geometries.append(
                svg.Line(
                    x1,
                    y1,
                    x2,
                    y2,
                    stroke="gray",
                    stroke_width=1,
                )
            )
        if M.edges[(u, v, k)]["edge"] == EdgeType.SUPPORT:
            geometries.append(
                svg.Line(
                    x1,
                    y1,
                    x2,
                    y2,
                    z=1000,
                    stroke="black",
                    stroke_width=4,
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
    d = svg.Drawing(3000, 3000, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


INSTANCE = {
    "lattice_type": "sqr",
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
        (10, 0),
        (20, 0),
        (0, 10),
        (10, 10),
        (20, 10),
        (0, 20),
        (10, 20),
        (20, 20),
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
    m = 30
    n = 30
    lattice_type = INSTANCE["lattice_type"]
    G = get_routing_graph(lattice_type, (m, n))

    G = add_glyphs_to_nodes(INSTANCE, G)

    G = route_set_lines(INSTANCE, G)
    # G = embed_to_routing_graph(INSTANCE, G)
    # G = render_line(INSTANCE, G, DEFAULT_PARAMS)
    geometries = geometrize(INSTANCE, G)
    img = draw_svg(geometries)

    with open("drawing.svg", "w") as f:
        f.write(img)
        f.flush()
