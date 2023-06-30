import json
import math
import time
from util.enums import NodeType, EdgeType, EdgePenalty
from itertools import combinations, pairwise, product, chain

import drawsvg as svg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from util.collections import (
    get_elements_in_same_lists,
    group_by_intersection_group,
    group_by_set,
    invert_dict_of_lists,
    list_of_lists_to_set_system_dict,
    merge_alternating,
)
from util.geometry import (
    are_faces_adjacent,
    centroid,
    dist_euclidean,
    do_lines_intersect,
    get_angle,
    get_closest_point,
    get_segment_circle_intersection,
    get_side,
    interpolate_biarcs,
    is_point_inside_circle,
    logical_coords_to_physical,
    offset_edge,
    get_linear_order,
)
from util.graph import (
    approximate_tsp_tour,
    approximate_steiner_tree,
    are_node_sets_connected,
    get_closest_pair,
    get_longest_simple_paths,
    get_shortest_path_between_sets,
    incident_edges,
    path_to_edges,
)
from util.layout import layout_qsap, layout_dr
from util.bundle import order_bundles, convert_to_bundleable
from util.perf import timing


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
        EdgePenalty.NINETY,
        EdgePenalty.ONE_THIRTY_FIVE,
        EdgePenalty.ONE_EIGHTY,
        EdgePenalty.ONE_THIRTY_FIVE,
        EdgePenalty.NINETY,
    ]
    for i in range(len(hex_sides)):
        p = 0
        n1 = hex_sides[i]
        for j in range(i + 1, len(hex_sides)):
            n2 = hex_sides[j]
            G.add_edge(
                n1,
                n2,
                EdgeType.PHYSICAL,
                edge=EdgeType.PHYSICAL,
                weight=EdgePenalty.HOP + penalties[p],
            )
            p += 1

    # connect center to sides
    for n in hex_sides:
        G.add_edge(
            node,
            n,
            EdgeType.PHYSICAL,
            edge=EdgeType.PHYSICAL,
            weight=EdgePenalty.HOP + EdgePenalty.TO_CENTER,
        )

    # find neighbors and replace center edges with physical edges to ports
    neighbors = list(nx.neighbors(G, node))
    for neighbor in neighbors:
        # use physically closest port for any neighbor
        neighbor_pos = G.nodes[neighbor]["pos"]
        port = get_closest_point(neighbor_pos, hex_sides)

        weight = dist_euclidean(port, neighbor_pos)
        G.add_edge(
            neighbor,
            port,
            EdgeType.PHYSICAL,
            edge=EdgeType.PHYSICAL,
            weight=EdgePenalty.HOP,
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
                weight=EdgePenalty.HOP,
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
                weight=EdgePenalty.HOP,
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

    penalties_cw = [
        EdgePenalty.FORTY_FIVE,
        EdgePenalty.NINETY,
        EdgePenalty.ONE_THIRTY_FIVE,
        EdgePenalty.ONE_EIGHTY,
        EdgePenalty.ONE_THIRTY_FIVE,
        EdgePenalty.NINETY,
        EdgePenalty.FORTY_FIVE,
    ]
    for i in range(0, len(sqr_corners) - 1):
        p = 0
        for j in range(i + 1, len(sqr_corners)):
            G.add_edge(
                sqr_corners[i],
                sqr_corners[j],
                EdgeType.PHYSICAL,
                edge=EdgeType.PHYSICAL,
                weight=EdgePenalty.HOP + penalties_cw[p],
            )
            p += 1

    # connect center to ports
    for n in sqr_corners:
        G.add_edge(
            node,
            n,
            EdgeType.PHYSICAL,
            edge=EdgeType.PHYSICAL,
            weight=EdgePenalty.HOP + EdgePenalty.TO_CENTER,
        )

    return G


def make_sqr_graph(m, n):
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.CENTER
    for node, n in G.nodes(data=True):
        x, y = node
        n["node"] = NodeType.CENTER
        n["logpos"] = node
        n["pos"] = logical_coords_to_physical(x, y, "sqr")

    G_ = nx.MultiGraph(incoming_graph_data=G)

    for node, _ in G.nodes(data=True):
        x, y = node
        can_tilt_right = y > 0 and x < m - 1
        if can_tilt_right:
            neighbor_nw = (x + 1, y - 1)
            G_.add_edge(
                node,
                neighbor_nw,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
            )
        can_tilt_left = y > 0 and x > 0
        if can_tilt_left:
            neighbor_ne = (x - 1, y - 1)
            G_.add_edge(
                node,
                neighbor_ne,
                EdgeType.CENTER,
                edge=EdgeType.CENTER,
            )

    for node, center in G.nodes(data=True):
        G_ = add_ports_to_sqr_node(G_, node, center, side_length=0.25)

    G_logical = nx.subgraph_view(G_, filter_edge=lambda u, v, k: k == EdgeType.CENTER)
    G_physical = nx.subgraph_view(
        G_, filter_edge=lambda u, v, k: k == EdgeType.PHYSICAL
    )

    for node in list(G_.nodes()):
        if G_.nodes[node]["node"] != NodeType.CENTER:
            continue

        ports = list(nx.neighbors(G_physical, node))

        # find neighbors and add physical edges to ports
        neighbors = list(
            nx.neighbors(
                G_logical,
                node,
            )
        )
        for neighbor in neighbors:
            ports_nb = list(nx.neighbors(G_physical, neighbor))
            # use physically closest port for any neighbor
            port_nb = get_closest_point(G_.nodes[node]["pos"], ports_nb)
            port_self = get_closest_point(G_.nodes[neighbor]["pos"], ports)

            # weight = dist_euclidean(port_nb, port_self)
            G_.add_edge(
                port_self,
                port_nb,
                EdgeType.PHYSICAL,
                edge=EdgeType.PHYSICAL,
                weight=EdgePenalty.HOP,
            )

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


# TODO could be a context manager
def block_edges_using(G, closed_nodes):
    for node in closed_nodes:
        ports = G.neighbors(node)
        for p in ports:
            for e in G.edges(p):
                u, v = e
                if (
                    G.nodes[u]["node"] == NodeType.PORT
                    and G.nodes[v]["node"] == NodeType.PORT
                ):
                    G.edges[e]["base_weight"] = G.edges[e]["weight"]
                    G.edges[e]["weight"] = float(math.inf)

    # for u, v in G.edges():
    #    if G.nodes[u]["node"] == NodeType.PORT and G.nodes[v]["node"] == NodeType.PORT:
    #        parents = set([G.nodes[u]["belongs_to"], G.nodes[v]["belongs_to"]])
    #        if len(parents.intersection(set(closed_nodes))) > 0:
    #            G.edges[(u, v)]["base_weight"] = G.edges[(u, v)]["weight"]
    #            G.edges[(u, v)]["weight"] = float(math.inf)
    return G


def unblock_edges(G, closed_nodes):
    for node in closed_nodes:
        ports = G.neighbors(node)
        for p in ports:
            for e in G.edges(p):
                if "base_weight" in G.edges[e]:
                    G.edges[e]["weight"] = G.edges[e]["base_weight"]

    # for u, v in G.edges():
    #    if G.nodes[u]["node"] == NodeType.PORT and G.nodes[v]["node"] == NodeType.PORT:
    #        parents = set([G.nodes[u]["belongs_to"], G.nodes[v]["belongs_to"]])
    #        if len(parents.intersection(set(closed_nodes))) > 0:
    #            G.edges[(u, v)]["weight"] = G.edges[(u, v)]["base_weight"]
    return G


def get_port_edge_between_centers(G, u, v):
    uports = G.neighbors(u)

    for p in uports:
        for w, x in G.edges(p):
            if G.nodes[x]["node"] == NodeType.PORT and G.nodes[x]["belongs_to"] == v:
                return (w, x)
    raise Exception("no port edge between centers")


def update_weights_for_support_edge(G, edge):
    u, v = edge
    un = G.nodes[u]["node"]
    vn = G.nodes[v]["node"]

    # only care about ports
    if un != NodeType.PORT or vn != NodeType.PORT:
        return G

    # make this edge cheaper so future routes will use it more
    G.edges[edge]["weight"] += EdgePenalty.IN_SUPPORT

    # penalize other edges crossing this one
    # is it an edge between nodes?
    uparent = G.nodes[u]["belongs_to"]
    vparent = G.nodes[v]["belongs_to"]

    if uparent == vparent and not G.nodes[uparent]["occupied"]:
        # nope, within a node
        # get all port edges and find those crossing this one
        ports = G.neighbors(uparent)
        port_edges = [G.edges(p) for p in ports]
        port_edges = list(set([item for sub_list in port_edges for item in sub_list]))
        port_edges = [
            (u, v)
            for u, v in port_edges
            if G.nodes[u]["node"] == NodeType.PORT
            and G.nodes[v]["node"] == NodeType.PORT
        ]
        crossing_edges = [
            (w, x) for w, x in port_edges if do_lines_intersect(w, x, u, v)
        ]
        for w, x in crossing_edges:
            G.edges[(w, x)]["weight"] = G.edges[(w, x)]["weight"] + EdgePenalty.CROSSING
    elif uparent != vparent:
        # yep, between nodes
        uparent, vparent = list(sorted([uparent, vparent], key=lambda x: x[1]))
        # now u is up
        ux, uy = uparent
        vx, vy = vparent
        is_left_tilt = ux < vx and uy < vy
        is_right_tilt = ux > vx and uy < vy

        if is_left_tilt:
            e = get_port_edge_between_centers(G, (ux, uy + 1), (ux + 1, uy))
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING
        if is_right_tilt:
            e = get_port_edge_between_centers(G, (ux - 1, uy), (ux, uy + 1))
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING

    return G


def route_set_lines(instance, G, element_set_partition, support_type="steiner-tree"):
    # steps:
    processed_elements = set()

    # TODO 100ms spent here
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
    for elements, sets in element_set_partition:
        print(elements, sets)
        new_edges = []
        S = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER and d["occupied"] and d["glyph"] in elements
        ]

        # close edges to all other occupied nodes so as to not route "over" them
        S_minus = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER
            and d["occupied"]
            and d["glyph"] not in elements
        ]
        G_ = block_edges_using(G_, S_minus)

        # 3. determine approximated steiner tree acc. to kou1981 / wu1986.
        # or TSP path
        if len(S) > 1:
            set_support = (
                approximate_steiner_tree(G_, S)
                if support_type == "steiner-tree"
                else approximate_tsp_tour(G_, S)
            )
            set_support_edges = (
                set_support.edges()
                if support_type == "steiner-tree"
                else path_to_edges(set_support)
            )
        else:
            set_support_edges = []
        # print([(G.nodes[u], G.nodes[v]) for u,v in set_support_edges])

        # add those edges to support
        for u, v in set_support_edges:
            if not G.has_edge(u, v, EdgeType.SUPPORT):
                new_edges.append((u, v))
                G.add_edge(
                    u, v, EdgeType.SUPPORT, edge=EdgeType.SUPPORT, sets=set(sets)
                )
            else:
                G.edges[(u, v, EdgeType.SUPPORT)]["sets"] = set(sets).union(
                    G.edges[(u, v, EdgeType.SUPPORT)]["sets"]
                )

        G_ = unblock_edges(G_, S_minus)

        for s in sets:
            processed_elements_for_s = set(instance["set_system"][s]).intersection(
                processed_elements
            )
            if len(processed_elements_for_s) > 0:
                nodes_of_processed_elements_for_s = [
                    n
                    for n, d in G_.nodes(data=True)
                    if d["node"] == NodeType.CENTER
                    and d["occupied"]
                    and d["glyph"] in processed_elements_for_s
                ]

                # if support is not connected but should be, fix it by connecting via shortest path
                support = nx.subgraph_view(
                    G, filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
                )
                is_connected = are_node_sets_connected(
                    support, S, nodes_of_processed_elements_for_s
                )

                if not is_connected:
                    S_minus = list(
                        set(instance["glyph_ids"]).difference(
                            set(elements).union(set(processed_elements_for_s))
                        )
                    )
                    S_minus = [
                        n
                        for n, d in G_.nodes(data=True)
                        if d["node"] == NodeType.CENTER
                        and d["occupied"]
                        and d["glyph"] in S_minus
                    ]
                    G_ = block_edges_using(
                        G_,
                        S_minus,
                    )
                    shortest_path_edgelist = get_shortest_path_between_sets(
                        G_, S, nodes_of_processed_elements_for_s
                    )
                    for u, v in shortest_path_edgelist:
                        if not G.has_edge(u, v, EdgeType.SUPPORT):
                            new_edges.append((u, v))
                            G.add_edge(
                                u,
                                v,
                                EdgeType.SUPPORT,
                                edge=EdgeType.SUPPORT,
                                sets=set([s]),
                            )
                        else:
                            G.edges[(u, v, EdgeType.SUPPORT)]["sets"] = set([s]).union(
                                G.edges[(u, v, EdgeType.SUPPORT)]["sets"]
                            )

                    G_ = unblock_edges(G_, S_minus)
        processed_elements = processed_elements.union(set(elements))

        # 4. update host graph: lighten weight on edges in F as suggested in castermans2019,
        # TODO only applicable if EdgePenalty.HOP is greater than zero
        # penalize crossing dual edges as suggested by bast2020
        for e in new_edges:
            u, v = e
            update_weights_for_support_edge(G_, e)
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
        px, py = (x * factor + mx, -y * factor + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            c = svg.Circle(cx=px, cy=py, r=one_unit_px / 4)
            c.append_title(M.nodes[i]["glyph"])
            geometries.append(c)

    for u, v, k in M.edges(keys=True):
        x1, y1 = M.nodes[u]["pos"]
        x2, y2 = M.nodes[v]["pos"]

        if M.edges[(u, v, k)]["edge"] == EdgeType.PHYSICAL:
            if (
                M.nodes[u]["node"] != NodeType.PORT
                or M.nodes[v]["node"] != NodeType.PORT
            ):
                continue
            w = M.edges[(u, v, k)].get("weight", -1)
            geometries.append(
                svg.Line(
                    x1,
                    y1,
                    x2,
                    y2,
                    data_weight=w,
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
                    # data_sets=M.edges[(u, v, k)]["sets"],
                    stroke="black",
                    stroke_width=4,
                )
            )

    # uncomment to draw edge-bundled lines

    set_colors = [
        "#1f78b4",
        "#33a02c",
        "#e31a1c",
        "#ff7f00",
        "#6a3d9a",
        "#ffff99",
        "#b15928",
        "#a6cee3",
        "#b2df8a",
        "#fb9a99",
    ]
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
        for i, k in enumerate(paths_at_edge):
            offset_dir = math.pi / 2
            offset_length = centering_offset + i * offset
            o_u, o_v = offset_edge((src, tgt), edge_angle - offset_dir, offset_length)

            M.edges[(u, v, k)]["edge_pos"] = {
                (u, v): (o_u, o_v),
                (v, u): (o_v, o_u),
            }

    sets = list(map(lambda x: instance["set_ftb_order"].index(x) + 1, instance["sets"]))

    for k in sets:
        M_ = nx.subgraph_view(M, filter_edge=lambda u, v, k_: k == k_)
        endpoints = [n for n in M_.nodes() if M_.degree(n) == 1]
        if len(endpoints) < 1:
            print("wtf", instance["sets"][k - 1])
            continue
        start = endpoints[0]
        start_nb = list(nx.neighbors(M_, start))
        if len(start_nb) < 1:
            continue
        paths = get_longest_simple_paths(M_, endpoints[0], start_nb[0])
        for path in paths:
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

                edge = M_.edges[(u, v, k)]["edge_pos"][(u, v)]
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
                            get_segment_circle_intersection(
                                (uhub_center, a), uhub_circle
                            ),
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
                            get_segment_circle_intersection(
                                (vhub_center, b), vhub_circle
                            ),
                        ),
                    )

            kwargs = {
                "close": False,
                "stroke_width": 2,
                "fill": "none",
                "stroke": set_colors[k - 1],
                "z": k,
            }
            line = interpolate_biarcs(keypoints, **kwargs)
            geometries.append(line)

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


def read_instance(name):
    with open(f"data/{name}.json") as f:
        data = json.load(f)
    elements = data["E"]
    sets = data["S"]
    return {
        "glyph_ids": elements,
        "sets": sets,
        "set_system": list_of_lists_to_set_system_dict(elements, data["SR"]),
        "D_EA": data["EA"],
        "D_SR": data["SA"],
        "set_ftb_order": list(sorted(sets)),
        # pipeline config
        "dr_method": "mds",
        "dr_gridification": "hagrid",  #  'hagrid' or 'dgrid'
        "support_type": "path",  #  'path' or 'steiner-tree'
        "support_partition_by": "set",  #  'set' or 'intersection-group'
    }


if __name__ == "__main__":
    m = 10
    n = 10
    instance = read_instance("wienerlinien/wienerlinien_sm")
    lattice_type = "sqr"

    with timing("layout"):
        instance["glyph_positions"] = layout_dr(
            instance["glyph_ids"],
            instance["D_EA"],
            instance["D_SR"],
            m=m,
            n=n,
            weight=0,
        )

    with timing("routing"):
        G = get_routing_graph(lattice_type, (m, n))
        G = add_glyphs_to_nodes(instance, G)
        element_set_partition = (
            group_by_intersection_group(instance["set_system"])
            if instance["support_partition_by"] == "intersection-group"
            else group_by_set(instance["set_system"])
        )
        element_set_partition = sorted(
            element_set_partition, key=lambda x: len(x[1]), reverse=True
        )
        M = route_set_lines(
            instance, G, element_set_partition, support_type=instance["support_type"]
        )

    G = convert_to_bundleable(instance, G)
    # after this we have a multigraph, with `set_id` property
    G = order_bundles(instance, G)
    # M.add_edges_from(list(G.edges(keys=True,data=True)))
    # print(list([(u,v,k,d) for u,v,k,d in M.edges(keys=True,data=True) if 'edge' not in d]))

    geometries = geometrize(instance, G)

    with timing("draw svg"):
        img = draw_svg(geometries)
    with timing("write svg"):
        with open("drawing.svg", "w") as f:
            f.write(img)
            f.flush()

    if False:
        G = nx.subgraph_view(
            G,
            filter_edge=lambda u, v, k: k == EdgeType.PHYSICAL
            # and G.nodes[u]["node"] == NodeType.PORT
            # and G.nodes[v]["node"] == NodeType.PORT
            # and G.nodes[u]["belongs_to"] != G.nodes[v]["belongs_to"],
        )
        nx.draw(G, pos=nx.get_node_attributes(G, "pos"), node_size=50)
        plt.show()
