import json
import math
import subprocess
import sys
from collections import defaultdict
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
    invert_list,
    list_of_lists_to_set_system_dict,
    merge_alternating,
)
from util.draw import (
    CELL_SIZE_PX,
    DRAW_DEG1_MARK_SIZE_PX,
    DRAW_DEG1_MARKS,
    DRAW_GLYPHS,
    DRAW_GLYPHS_OVER_LINES,
    DRAW_HUBS,
    GLYPH_SIZE_PX,
    LINE_GAP,
    MARGINS,
    NODE_CIRCLE_RADIUS,
    SET_COLORS,
    STROKE_WIDTH,
    draw_support,
    draw_svg,
    edge_filter_ports,
)
from util.enums import EdgePenalty, EdgeType, NodeType, PortDirs
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
from util.graph import path_to_edges, get_port_edges, extract_support_layer, get_ports
from util.layout import layout_dr, layout_dr_multiple, layout_qsap
from util.perf import timing
from util.route import route_multilayer_ilp, route_multilayer_heuristic


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
    sqr_corners_dirs = list(zip(PortDirs, sqr_corners))
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
                efrom=PortDirs[i],
                eto=PortDirs[j],
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


def make_sqr_graph(m, n, with_ports=True):
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.CENTER
    for node, d in G.nodes(data=True):
        x, y = node
        d["node"] = NodeType.CENTER
        d["logpos"] = node
        d["pos"] = logical_coords_to_physical(x, y, "sqr")

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

    for u, v in list(G_.edges()):
        ux, uy = u
        vx, vy = v
        dx = vx - ux
        dy = vy - uy

        if (
            dx != 0
            and dy != 0
            and ux + dx >= 0
            and ux + dx < m
            and uy + dy >= 0
            and uy + dy < n
        ):
            t = ((ux + dx, uy), (ux, uy + dy))
            G_.edges[(u, v, EdgeType.CENTER)]["crossing"] = t

    if with_ports:
        for node, center in G.nodes(data=True):
            G_ = add_ports_to_sqr_node(G_, node, center, side_length=0.25)

        G_logical = nx.subgraph_view(
            G_, filter_edge=lambda u, v, k: k == EdgeType.CENTER
        )
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


def get_routing_graph(lattice_type, lattice_size, with_ports=True):
    m, n = lattice_size
    G = None
    match lattice_type:
        case "sqr":
            G = make_sqr_graph(m, n, with_ports=with_ports)
        case _:
            raise Exception(f"unknown lattice type {lattice_type}")
    return G


def add_glyphs_to_nodes(instance, G):
    num_layers = instance["num_layers"]

    for n in [n for n in G.nodes() if G.nodes[n]["node"] == NodeType.CENTER]:
        # G.nodes[n]["occupied"] = False
        layers = []

        for k in range(num_layers):
            layer_info = {"occupied": False}
            for i, element in enumerate(instance["elements"]):
                logpos = instance["glyph_positions"][k][i]
                if logpos == n:
                    layer_info = {
                        "occupied": True,
                        "label": element,
                    }
                    if "glyph_ids" in instance:
                        layer_info["glyph"] = instance["glyph_ids"][i]
            layers.append(layer_info)
        G.nodes[n]["layers"] = layers
    return G

    # for i, element in enumerate(instance["elements"]):
    #    logpos = instance["glyph_positions"][i]
    #    if G.nodes[logpos]["node"] != NodeType.CENTER:
    #        raise Exception("node to position glyph on is somehow not a glyph center")
    #    G.nodes[logpos]["occupied"] = True
    #    G.nodes[logpos]["label"] = element
    #    if "glyph_ids" in instance:
    #        G.nodes[logpos]["glyph"] = instance["glyph_ids"][i]
    # return G

    # for u, v in G.edges():
    #    if G.nodes[u]["node"] == NodeType.PORT and G.nodes[v]["node"] == NodeType.PORT:
    #        parents = set([G.nodes[u]["belongs_to"], G.nodes[v]["belongs_to"]])
    #        if len(parents.intersection(set(closed_nodes))) > 0:
    #            G.edges[(u, v)]["base_weight"] = G.edges[(u, v)]["weight"]
    #            G.edges[(u, v)]["weight"] = float(math.inf)
    return G


def geometrize(instance, M):
    geometries = []
    mx, my = MARGINS

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (x * CELL_SIZE_PX + mx, -y * CELL_SIZE_PX + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            if DRAW_GLYPHS and "glyph" in M.nodes[i]:
                w, h = (GLYPH_SIZE_PX, GLYPH_SIZE_PX)
                img = svg.Image(
                    px - w / 2, py - h / 2, width=w, height=h, path=M.nodes[i]["glyph"]
                )
                img.append_title(M.nodes[i]["label"])
                geometries.append(img)
            else:
                c = svg.Circle(cx=px, cy=py, r=NODE_CIRCLE_RADIUS)
                c.append_title(M.nodes[i]["label"])
                geometries.append(c)

    for u, v in M.edges():
        if not edge_filter_ports(M, u, v, same_centers=False):
            continue

        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        paths_at_edge = M.edges[(u, v)]["oeb_order"][(u, v)]
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -LINE_GAP
        M.edges[(u, v)]["edge_pos"] = {}
        for i, set_id in enumerate(paths_at_edge):
            offset_dir = 3 * math.pi / 2
            offset_length = centering_offset + i * LINE_GAP
            o_u, o_v = offset_edge((src, tgt), edge_angle - offset_dir, offset_length)

            M.edges[(u, v)]["edge_pos"][set_id] = {
                (u, v): (o_u, o_v),
                (v, u): (o_v, o_u),
            }

    for set_id in instance["sets"]:
        G_ = nx.subgraph_view(
            M,
            filter_edge=lambda u, v: edge_filter_ports(M, u, v, same_centers=False)
            and set_id in M.edges[(u, v)]["sets"],
        )
        # this draws straight lines between nodes
        for u, v in G_.edges():
            upos, vpos = G_.edges[(u, v)]["edge_pos"][set_id][(u, v)]

            u_intersect = get_segment_circle_intersection(
                (upos, vpos),
                (M.nodes[M.nodes[u]["belongs_to"]]["pos"], NODE_CIRCLE_RADIUS),
            )
            v_intersect = get_segment_circle_intersection(
                (upos, vpos),
                (M.nodes[M.nodes[v]["belongs_to"]]["pos"], NODE_CIRCLE_RADIUS),
            )

            line = svg.Path(
                **{
                    "close": False,
                    "stroke_width": STROKE_WIDTH,
                    "fill": "none",
                    "stroke": SET_COLORS[instance["set_ftb_order"].index(set_id)],
                }
            )
            line.M(*u_intersect)
            line.L(*v_intersect)
            geometries.append(line)

        # this draws all the connections at non-occupied nodes
        if False:
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
                    (uupos, upos), (M.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS)
                )
                v_intersect = get_segment_circle_intersection(
                    (vpos, vvpos), (M.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS)
                )

                uu_u_center = centroid([uupos, upos])
                vv_v_center = centroid([vvpos, vpos])

                line = svg.Path(
                    **{
                        "close": False,
                        "stroke_width": STROKE_WIDTH,
                        "fill": "none",
                        "stroke": SET_COLORS[instance["set_ftb_order"].index(set_id)],
                    }
                )
                barc = biarc(uu_u_center, u_intersect, v_intersect, vv_v_center)
                draw_biarc(line, barc)
                geometries.append(line)

        if False:
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
                            if k == EdgeType.SUPPORT
                            and M.nodes[x]["node"] == NodeType.PORT
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
                if len(used_ports) == 1 and DRAW_DEG1_MARKS:
                    # this is a deg 1 node for this set
                    # idk, could connect to center
                    # or maybe draw a small mark?
                    a = used_ports.pop()
                    _, b = outward_edge_at_port[a]
                    apos, bpos = M.edges[(a, b, EdgeType.SUPPORT)]["edge_pos"][set_id][
                        (a, b)
                    ]
                    cx, cy = get_segment_circle_intersection(
                        (apos, bpos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS)
                    )
                    circle = svg.Circle(
                        cx=cx,
                        cy=cy,
                        r=DRAW_DEG1_MARK_SIZE_PX,
                        fill=SET_COLORS[instance["set_ftb_order"].index(set_id)],
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

                within_node_connections = nx.minimum_spanning_tree(
                    G_node, weight="weight"
                )
                for a, b in within_node_connections.edges():
                    line = svg.Path(
                        **{
                            "close": False,
                            "stroke_width": STROKE_WIDTH,
                            "fill": "none",
                            "data_weight": G_node.edges[a, b]["weight"],
                            "stroke": SET_COLORS[
                                instance["set_ftb_order"].index(set_id)
                            ],
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
                        (vpos, apos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS)
                    )
                    b_intersect = get_segment_circle_intersection(
                        (bpos, xpos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS)
                    )

                    av_center = centroid([apos, vpos])
                    bx_center = centroid([bpos, xpos])

                    barc = biarc(av_center, a_intersect, b_intersect, bx_center)
                    draw_biarc(line, barc)

                    geometries.append(line)

    if DRAW_HUBS:
        hubs = [n for n in M.nodes() if M.degree[n] > 0]
        for hub in hubs:
            cx, cy = M.nodes[hub]["pos"]
            r = 1 / 16 * CELL_SIZE_PX
            geometries.append(
                svg.Circle(cx, cy, r, fill="none", stroke="gray", stroke_width=1)
            )

    if DRAW_GLYPHS_OVER_LINES:
        geometries = list(reversed(geometries))

    return geometries

def read_instance(name):
    with open(f"data/{name}.json") as f:
        data = json.load(f)
    elements = data["E"]
    sets = data["S"]
    inst = {
        "grid_x": 10,
        "grid_y": 10,
        "elements": elements,
        "elements_inv": invert_list(elements),
        "sets": sets,
        "set_system": list_of_lists_to_set_system_dict(elements, data["SR"]),
        "D_EA": data["EA"],
        "D_SR": data["SA"],
        "set_ftb_order": list(sorted(sets)),
        # pipeline config
        "strategy": "heuristic",  # 'opt' or 'heuristic'
        "dr_method": "mds",
        "dr_gridification": "hagrid",  #  'hagrid' or 'dgrid'
        "support_type": "steiner-tree",  #  'path' or 'steiner-tree'
        "support_partition_by": "set",  #  'set' or 'intersection-group'
    }
    if "glyph_ids" in data:
        inst["glyph_ids"] = data["glyph_ids"]
    return inst


if __name__ == "__main__":
    instance = read_instance("imdb/imdb_10")
    lattice_type = "sqr"
    m = instance["grid_x"]
    n = instance["grid_y"]
    num_layers = 2
    instance["num_layers"] = num_layers

    with timing("layout"):
        instance["glyph_positions"] = layout_dr_multiple(
            instance["D_EA"],
            instance["D_SR"],
            m=m,
            n=n,
            num_samples=num_layers,
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

        if instance["strategy"] == "opt":
            L = route_multilayer_ilp(
                instance,
                nx.subgraph_view(
                    G,
                    filter_edge=lambda u, v, k: k == EdgeType.CENTER,
                    filter_node=lambda n: G.nodes[n]["node"] == NodeType.CENTER,
                ),
                element_set_partition,
                support_type=instance["support_type"],
            )
        else:
            L = route_multilayer_heuristic(
                instance,
                G,
                element_set_partition,
                support_type=instance["support_type"],
            )

    with timing("bundle lines"):
        L = bundle_lines(instance, L)


    if instance["strategy"] == "opt":
        # merge grid graph data (ports at nodes) with line graph data (routing and bundling info)
        # we have
        # G = multigraph with center and physical nodes/edges
        # L = a multigraph with (layer, support) edges and center nodes
        for layer in range(num_layers):
            for i, esp in enumerate(element_set_partition):
                elements, sets = esp
                root_pos = instance["glyph_positions"][layer][
                    instance["elements_inv"][elements[0]]
                ]
                for j, el in enumerate(elements):
                    if j == 0:
                        continue
                    j_pos = instance["glyph_positions"][layer][
                        instance["elements_inv"][el]
                    ]
                    P = nx.subgraph_view(
                        L,
                        filter_edge=lambda u, v, k: k == (layer, EdgeType.SUPPORT)
                        and i in L.edges[u, v, k]["partitions"],
                    )
                    path_j = path_to_edges(
                        nx.shortest_path(P, root_pos, j_pos)
                    )  # P should actually just be a path already but we do this to order edges

                    edge_pairs_path_j = list(pairwise(path_j))
                    for l, edge_pair in enumerate(edge_pairs_path_j):
                        e1, e2 = edge_pair
                        u, v = e1
                        v, x = e2

                        port_u, port_vu = get_ports(G, u, v)
                        port_vx, port_x = get_ports(G, v, x)

                        if l == 0:
                            # add edge from center to first port
                            G.add_edge(
                                u,
                                port_u,
                                (layer, EdgeType.SUPPORT),
                                edge=EdgeType.SUPPORT,
                                sets=set(
                                    L.edges[u, v, (layer, EdgeType.SUPPORT)]["sets"]
                                ),
                            )
                        if l == len(edge_pairs_path_j) - 1:
                            # add edge from last port to center
                            G.add_edge(
                                port_x,
                                x,
                                (layer, EdgeType.SUPPORT),
                                edge=EdgeType.SUPPORT,
                                sets=set(
                                    L.edges[v, x, (layer, EdgeType.SUPPORT)]["sets"]
                                ),
                            )

                        # TODO remove MST heurisitic in drawing code

                        oeb_order_u_v = {
                            (port_u, port_vu): L.edges[u, v, (layer, EdgeType.SUPPORT)][
                                "oeb_order"
                            ][(u, v)],
                            (port_vu, port_u): L.edges[u, v, (layer, EdgeType.SUPPORT)][
                                "oeb_order"
                            ][(v, u)],
                        }
                        G.add_edge(
                            port_u,
                            port_vu,
                            (layer, EdgeType.SUPPORT),
                            **{
                                **L.edges[u, v, (layer, EdgeType.SUPPORT)],
                                "oeb_order": oeb_order_u_v,
                            },
                        )

                        oeb_order_v_x = {
                            (port_vx, port_x): L.edges[v, x, (layer, EdgeType.SUPPORT)][
                                "oeb_order"
                            ][(v, x)],
                            (port_x, port_vx): L.edges[v, x, (layer, EdgeType.SUPPORT)][
                                "oeb_order"
                            ][(x, v)],
                        }
                        G.add_edge(
                            port_vx,
                            port_x,
                            (layer, EdgeType.SUPPORT),
                            **{
                                **L.edges[v, x, (layer, EdgeType.SUPPORT)],
                                "oeb_order": oeb_order_v_x,
                            },
                        )

                        G.add_edge(
                            port_vu,
                            port_vx,
                            (layer, EdgeType.SUPPORT),
                            edge=EdgeType.SUPPORT,
                            sets=set(
                                L.edges[v, x, (layer, EdgeType.SUPPORT)]["sets"]
                            ).intersection(
                                set(L.edges[u, v, (layer, EdgeType.SUPPORT)]["sets"])
                            ),
                        )
        L = G

    # draw_support(instance, M.copy())

    with timing("draw+write svg"):
        for layer in range(num_layers):
            M_ = extract_support_layer(L, layer)
            print('layer', layer, M_)
            geometries = geometrize(instance, M_)
            img = draw_svg(geometries, m * CELL_SIZE_PX, n * CELL_SIZE_PX)
            with open(f"drawing_{layer}.svg", "w") as f:
                f.write(img)
                f.flush()

    print("Done.")
