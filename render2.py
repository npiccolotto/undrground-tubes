from collections import defaultdict
import json
import math
import subprocess
import sys
from itertools import chain, combinations, pairwise, product

import drawsvg as svg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from util.bundle import bundle_lines
from util.collections import (
    get_elements_in_same_lists,
    group_by_intersection_group,
    group_by_set,
    invert_dict_of_lists,
    list_of_lists_to_set_system_dict,
    merge_alternating,
)
from util.enums import EdgePenalty, EdgeType, NodeType
from util.geometry import (
    are_faces_adjacent,
    biarc,
    centroid,
    dist_euclidean,
    do_lines_intersect,
    draw_biarc,
    get_angle,
    get_closest_point,
    get_linear_order,
    get_segment_circle_intersection,
    get_side,
    interpolate_biarcs,
    is_point_inside_circle,
    logical_coords_to_physical,
    offset_edge,
    offset_point,
)
from util.graph import (
    approximate_steiner_tree,
    approximate_steiner_tree_nx,
    approximate_tsp_tour,
    are_node_sets_connected,
    calculate_path_length,
    get_closest_pair,
    get_longest_simple_paths,
    get_node_with_degree,
    get_shortest_path_between_sets,
    incident_edges,
    path_to_edges,
    visit_edge_pairs_starting_at_node,
)
from util.layout import layout_dr, layout_qsap
from util.perf import timing

DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 100,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
}


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
                efrom=dirs[i],
                eto=dirs[j],
                epenalty=penalties_cw[p],
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
            efrom=n,
            eto="center",
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
            length_penalty = dist_euclidean(port_nb, port_self)

            G_.add_edge(
                port_self,
                port_nb,
                EdgeType.PHYSICAL,
                edge=EdgeType.PHYSICAL,
                weight=EdgePenalty.HOP + length_penalty,
            )

    return G_


def get_routing_graph(lattice_type, lattice_size):
    m, n = lattice_size
    G = None
    match lattice_type:
        case "sqr":
            G = make_sqr_graph(m, n)
        case _:
            raise Exception(f"unknown lattice type {lattice_type}")
    return G


def add_glyphs_to_nodes(instance, G):
    for n in G.nodes():
        G.nodes[n]["occupied"] = False
    for i, element in enumerate(instance["elements"]):
        logpos = instance["glyph_positions"][i]
        if G.nodes[logpos]["node"] != NodeType.CENTER:
            raise Exception("node to position glyph on is somehow not a glyph center")
        G.nodes[logpos]["occupied"] = True
        G.nodes[logpos]["label"] = element
        if "glyph_ids" in instance:
            G.nodes[logpos]["glyph"] = instance["glyph_ids"][i]
    return G


# TODO could be a context manager
def block_edges_using(G, closed_nodes):
    for node in closed_nodes:
        if G.nodes[node]["node"] != NodeType.CENTER:
            raise BaseException("trying to block non-center node!")
        ports = G.neighbors(node)
        for p in ports:
            for e in G.edges(p):
                u, v = e
                if (
                    G.nodes[u]["node"] == NodeType.PORT
                    and G.nodes[v]["node"] == NodeType.PORT
                    and ("blocked" not in G.edges[e] or not G.edges[e]["blocked"])
                ):
                    G.edges[e]["blocked"] = True
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
                if "blocked" in G.edges[e] and G.edges[e]["blocked"]:
                    G.edges[e]["blocked"] = False
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
    """
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
    """
    if uparent != vparent:
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
            if k == EdgeType.PHYSICAL
        ]
    )

    # 2. then process from biggest to smallest group. for each:
    for elements, sets in element_set_partition:
        new_edges = []
        S = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER and d["occupied"] and d["label"] in elements
        ]

        # close edges to all other occupied nodes so as to not route "over" them
        S_minus = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER
            and d["occupied"]
            and d["label"] not in elements
        ]
        G_ = block_edges_using(G_, S_minus)

        # 3. determine approximated steiner tree acc. to kou1981 / wu1986.
        # or TSP path
        if len(S) > 1:
            with timing("sub support"):
                set_support = (
                    approximate_steiner_tree_nx(G_, S)
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
                    and d["label"] in processed_elements_for_s
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
                        set(instance["elements"]).difference(
                            set(elements).union(set(processed_elements_for_s))
                        )
                    )
                    S_minus = [
                        n
                        for n, d in G_.nodes(data=True)
                        if d["node"] == NodeType.CENTER
                        and d["occupied"]
                        and d["label"] in S_minus
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


def edge_filter_ports(G, u, v, same_centers=False):
    uparent = (
        None if G.nodes[u]["node"] == NodeType.CENTER else G.nodes[u]["belongs_to"]
    )
    vparent = (
        None if G.nodes[v]["node"] == NodeType.CENTER else G.nodes[v]["belongs_to"]
    )

    match (uparent, vparent):
        case (None, None):
            return False
        case (None, _):
            return False
        case (_, None):
            return False
        case (_, _):
            match same_centers:
                case True:
                    return uparent == vparent
                case False:
                    return uparent != vparent
                case None:
                    return True


def get_port_edges(M, node):
    all_ports = [
        p for p in nx.neighbors(M, node) if M.nodes[p]["node"] == NodeType.PORT
    ]
    all_edges_at_ports = set()
    port_dirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"]
    for port in all_ports:
        edges_at_port = [
            (a, b)
            if port_dirs.index(M.nodes[a]["port"]) < port_dirs.index(M.nodes[b]["port"])
            else (b, a)
            for a, b, k in M.edges(nbunch=port, keys=True)
            if k == EdgeType.PHYSICAL
            and M.nodes[a]["node"] == NodeType.PORT
            and M.nodes[b]["node"] == NodeType.PORT
            and M.nodes[a]["belongs_to"] == node
            and M.nodes[b]["belongs_to"] == node
        ]

        all_edges_at_ports = all_edges_at_ports.union(set(edges_at_port))
    return list(all_edges_at_ports)


def are_port_edges_crossing(us, ut, vs, vt):
    cw_dirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"] + [
        "n",
        "ne",
        "e",
        "se",
        "s",
        "sw",
        "w",
        "nw",
    ]
    ccw_dirs = ["n", "nw", "w", "sw", "s", "se", "e", "ne"] + [
        "n",
        "nw",
        "w",
        "sw",
        "s",
        "se",
        "e",
        "ne",
    ]

    # start at us and go clockwise to ut, collect all ports in between
    cw_set = []
    i = cw_dirs.index(us["port"]) + 1
    j = cw_dirs.index(ut["port"])
    if j > i:
        cw_set = cw_dirs[i:j]

    # then start at us and go ccw to ut, collect all ports in between
    ccw_set = []
    i = ccw_dirs.index(ut["port"]) + 1
    j = ccw_dirs.index(us["port"])
    if j > i:
        ccw_set = ccw_dirs[i:j]

    # edges cross iff vs in former and vt in latter set (or vice versa)
    port1 = vs["port"]
    port2 = vt["port"]

    are_crossing = True
    if port1 in cw_set and port2 in cw_set and port1 in ccw_set and port2 in ccw_set:
        are_crossing = False
    if port1 in [ut["port"], us["port"]] or port2 in [ut["port"], us["port"]]:
        are_crossing = False
    return are_crossing


def get_crossing_port_edges(G):
    # assuming G consinsts only of port edges of one node
    port_dirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"]

    edges = list(
        [
            (u, v)
            if port_dirs.index(M.nodes[u]["port"]) < port_dirs.index(M.nodes[v]["port"])
            else (v, u)
            for (u, v) in G.edges()
        ]
    )
    edges = list(sorted(edges, key=lambda e: port_dirs.index(G.nodes[e[0]]["port"])))
    crossings = []

    for i, e1 in enumerate(edges):
        u, v = e1
        for j, e2 in enumerate(edges):
            if j <= i:
                continue
            w, x = e2
            if are_port_edges_crossing(G.nodes[u], G.nodes[v], G.nodes[w], G.nodes[x]):
                crossings.append((e1, e2))
    return crossings


def geometrize(instance, M):
    geometries = []
    one_unit_px = DEFAULT_PARAMS["unit_size_in_px"]
    margins = np.array((0.5, 0.5)) * one_unit_px
    factor = one_unit_px
    mx, my = margins

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (x * factor + mx, -y * factor + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            if "glyph" not in M.nodes[i]:
                c = svg.Circle(cx=px, cy=py, r=one_unit_px / 4)
                c.append_title(M.nodes[i]["label"])
                geometries.append(c)
            else:
                w = 60
                h = 90
                img = svg.Image(
                    px - w / 2, py - h / 2, width=w, height=h, path=M.nodes[i]["glyph"]
                )
                img.append_title(M.nodes[i]["label"])
                geometries.append(img)

    set_colors = [
        "#1f78b4",
        "#33a02c",
        "#e31a1c",
        "#ff7f00",
        "#6a3d9a",
        "#ffff33",
        "#b15928",
        "#a6cee3",
        "#b2df8a",
        "#fb9a99",
        "#666666",
    ]

    for u, v, k in M.edges(keys=True):
        if k != EdgeType.SUPPORT or not edge_filter_ports(M, u, v, same_centers=False):
            continue

        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        # print normal of edge vector just for debugging
        # cx1, cy1 = centroid([src, tgt])
        # cx2, cy2 = offset_point((cx1, cy1), edge_angle + math.pi / 2, 10)
        # geometries.append(
        #    svg.Line(cx1, cy1, cx2, cy2, stroke="lightgray", stroke_width=1)
        # )

        paths_at_edge = M.edges[(u, v, EdgeType.SUPPORT)]["oeb_order"][(u, v)]
        offset = factor / 30  # TODO idk some pixels
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -offset
        M.edges[(u, v, EdgeType.SUPPORT)]["edge_pos"] = {}
        for i, set_id in enumerate(paths_at_edge):
            offset_dir = 3 * math.pi / 2
            offset_length = centering_offset + i * offset
            o_u, o_v = offset_edge((src, tgt), edge_angle - offset_dir, offset_length)

            M.edges[(u, v, EdgeType.SUPPORT)]["edge_pos"][set_id] = {
                (u, v): (o_u, o_v),
                (v, u): (o_v, o_u),
            }

    print(instance["set_ftb_order"])
    for set_id in instance["sets"]:
        print(set_id)
        G_ = nx.subgraph_view(
            M,
            filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
            and edge_filter_ports(M, u, v, same_centers=False)
            and set_id in M.edges[(u, v, EdgeType.SUPPORT)]["sets"],
        )
        # this draws straight lines between nodes
        for u, v in G_.edges():
            circle_r = one_unit_px * 0.25
            upos, vpos = G_.edges[(u, v, EdgeType.SUPPORT)]["edge_pos"][set_id][(u, v)]

            u_intersect = get_segment_circle_intersection(
                (upos, vpos), (M.nodes[M.nodes[u]["belongs_to"]]["pos"], circle_r)
            )
            v_intersect = get_segment_circle_intersection(
                (upos, vpos), (M.nodes[M.nodes[v]["belongs_to"]]["pos"], circle_r)
            )

            line = svg.Path(
                **{
                    "close": False,
                    "stroke_width": 2,
                    "fill": "none",
                    "stroke": set_colors[instance["set_ftb_order"].index(set_id)],
                }
            )
            line.M(*u_intersect)
            line.L(*v_intersect)
            geometries.append(line)

        # this draws all the connections at non-occupied nodes
        for u, v, k, d in M.edges(data=True, keys=True):
            if k != EdgeType.SUPPORT or not edge_filter_ports(
                G_, u, v, same_centers=True
            ):
                continue
            uparent = M.nodes[u]["belongs_to"]

            u_adjacent = list(
                [
                    x if w == u else w
                    for w, x, kk in M.edges(nbunch=[u], keys=True)
                    if kk == EdgeType.SUPPORT
                ]
            )
            uu = None
            for x in u_adjacent:
                xparent = M.nodes[x].get("belongs_to")
                if xparent is not None and xparent != uparent:
                    uu = x
                    break

            v_adjacent = list(
                [
                    x if w == v else w
                    for w, x, kk in M.edges(nbunch=[v], keys=True)
                    if kk == EdgeType.SUPPORT
                ]
            )
            vv = None
            for x in v_adjacent:
                xparent = M.nodes[x].get("belongs_to")
                if xparent is not None and xparent != uparent:
                    vv = x
                    break

            if uu is None or vv is None:
                print("howw")
                continue

            if (
                set_id not in M.edges[(uu, u, EdgeType.SUPPORT)]["edge_pos"]
                or set_id not in M.edges[(v, vv, EdgeType.SUPPORT)]["edge_pos"]
            ):
                continue

            uupos, upos = M.edges[(uu, u, EdgeType.SUPPORT)]["edge_pos"][set_id][
                (uu, u)
            ]
            vpos, vvpos = M.edges[(v, vv, EdgeType.SUPPORT)]["edge_pos"][set_id][
                (v, vv)
            ]

            u_intersect = get_segment_circle_intersection(
                (uupos, upos), (M.nodes[uparent]["pos"], circle_r)
            )
            v_intersect = get_segment_circle_intersection(
                (vpos, vvpos), (M.nodes[uparent]["pos"], circle_r)
            )

            uu_u_center = centroid([uupos, upos])
            vv_v_center = centroid([vvpos, vpos])

            line = svg.Path(
                **{
                    "close": False,
                    "stroke_width": 2,
                    "fill": "none",
                    "stroke": set_colors[instance["set_ftb_order"].index(set_id)],
                }
            )
            barc = biarc(uu_u_center, u_intersect, v_intersect, vv_v_center)
            draw_biarc(line, barc)
            # line.M(*uu_u_center)
            # line.L(*vv_v_center)
            geometries.append(line)

        # so this then draws connections within occupied nodes
        for node in [
            n
            for n, d in M.nodes(data=True)
            if d["node"] == NodeType.CENTER and d["occupied"]
        ]:
            # strategy
            # find out if this node is used for this set and by which ports
            # this we do by checking outgoing edges from each port and looking at `sets` property
            # when we identified the used ports by the current set
            # we take the subset of all port-port edges that are between used ports
            # we statically identify and penalize (possibly via block function?) crossing edges
            # and on this graph we do a minimum spanning tree, which should result in the optimal line-routing=

            # all_used_ports = set(
            #    [
            #        b
            #        for a, b, k in M.edges(nbunch=node, keys=True)
            #        if k == EdgeType.SUPPORT
            #    ]
            # )
            all_used_ports = []
            for p in [
                p
                for p in nx.neighbors(M, node)
                if M.nodes[p]["node"] == NodeType.PORT
                and M.nodes[p]["belongs_to"] == node
            ]:
                edges = [
                    (a, b)
                    for a, b, k, d in M.edges(nbunch=p, keys=True, data=True)
                    if k == EdgeType.SUPPORT and "sets" in d
                ]
                if len(edges) > 0:
                    all_used_ports.append(p)

            used_ports = set()
            outward_edge_at_port = dict()
            for port in all_used_ports:
                p_adjacent = list(
                    [
                        x
                        for w, x, k in M.edges(nbunch=[port], keys=True)
                        if k == EdgeType.SUPPORT and M.nodes[x]["node"] == NodeType.PORT
                    ]
                )
                for x in p_adjacent:
                    edge = M.edges[(port, x, EdgeType.SUPPORT)]
                    if "sets" in edge and set_id in edge["sets"]:
                        used_ports = used_ports.union(set([port]))
                        outward_edge_at_port[port] = (port, x)

            if len(used_ports) < 1:
                # cannot happen actually
                continue
            if len(used_ports) == 1:
                # this is a deg 1 node for this set
                # idk, could connect to center
                # or maybe draw a small mark?
                a = used_ports.pop()
                _, b = outward_edge_at_port[a]
                apos, bpos = M.edges[(a, b, EdgeType.SUPPORT)]["edge_pos"][set_id][
                    (a, b)
                ]
                cx, cy = get_segment_circle_intersection(
                    (apos, bpos), (M.nodes[node]["pos"], circle_r)
                )
                circle = svg.Circle(
                    cx=cx,
                    cy=cy,
                    r=circle_r / 8,
                    fill=set_colors[instance["set_ftb_order"].index(set_id)],
                )
                circle.append_title(set_id)
                geometries.append(circle)
                continue

            all_edges_at_ports = get_port_edges(M, node)
            port_port_edges = [
                (a, b, M.edges[(a, b, EdgeType.PHYSICAL)]["weight"])
                for a, b in all_edges_at_ports
                if a in used_ports and b in used_ports
            ]

            G_node = nx.Graph()
            G_node.add_weighted_edges_from(port_port_edges)
            for a, b, w in port_port_edges:
                G_node.add_node(a, **M.nodes[a])
                G_node.add_node(b, **M.nodes[b])

            # identify crossing edges and penalize the one with more weight
            for e1, e2 in get_crossing_port_edges(G_node):
                u, v = e1
                w, x = e2
                weight_uv = G_node.edges[u, v]["weight"]
                weight_wx = G_node.edges[w, x]["weight"]
                edge_to_penalize = e1
                if weight_wx > weight_uv:
                    edge_to_penalize = e2
                a, b = edge_to_penalize
                G_node.edges[a, b]["weight"] = float("inf")

            within_node_connections = nx.minimum_spanning_tree(G_node, weight="weight")
            for a, b in within_node_connections.edges():
                line = svg.Path(
                    **{
                        "close": False,
                        "stroke_width": 2,
                        "fill": "none",
                        "data_weight": G_node.edges[a, b]["weight"],
                        "stroke": set_colors[instance["set_ftb_order"].index(set_id)],
                    }
                )
                _, v = outward_edge_at_port[a]
                _, x = outward_edge_at_port[b]

                apos, vpos = M.edges[(a, v, EdgeType.SUPPORT)]["edge_pos"][set_id][
                    (a, v)
                ]
                bpos, xpos = M.edges[(b, x, EdgeType.SUPPORT)]["edge_pos"][set_id][
                    (b, x)
                ]

                a_intersect = get_segment_circle_intersection(
                    (vpos, apos), (M.nodes[node]["pos"], circle_r)
                )
                b_intersect = get_segment_circle_intersection(
                    (bpos, xpos), (M.nodes[node]["pos"], circle_r)
                )

                av_center = centroid([apos, vpos])
                bx_center = centroid([bpos, xpos])

                barc = biarc(av_center, a_intersect, b_intersect, bx_center)
                draw_biarc(line, barc)

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


def draw_svg(geometries, width, height):
    d = svg.Drawing(width, height, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


def read_instance(name):
    with open(f"data/{name}.json") as f:
        data = json.load(f)
    elements = data["E"]
    sets = data["S"]
    inst =  {
        "elements": elements,
        "sets": sets,
        "set_system": list_of_lists_to_set_system_dict(elements, data["SR"]),
        "D_EA": data["EA"],
        "D_SR": data["SA"],
        "set_ftb_order": list(sorted(sets)),
        # pipeline config
        "dr_method": "mds",
        "dr_gridification": "hagrid",  #  'hagrid' or 'dgrid'
        "support_type": "steiner-tree",  #  'path' or 'steiner-tree'
        "support_partition_by": "set",  #  'set' or 'intersection-group'
    }
    if "glyph_ids" in data:
        inst["glyph_ids"]= data["glyph_ids"]
    return inst


if __name__ == "__main__":
    m = 10
    n = 10
    instance = read_instance("imdb/imdb_10")
    lattice_type = "sqr"

    with timing("layout"):
        instance["glyph_positions"] = layout_dr(
            instance["elements"],
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

    with timing("bundle lines"):
        M = bundle_lines(instance, M)

    with timing("draw svg"):
        geometries = geometrize(instance, M)
        img = draw_svg(geometries, m * 100, n * 100)
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
