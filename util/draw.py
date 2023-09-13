import drawsvg as svg
import numpy as np
import math
import os
import networkx as nx
from util.enums import NodeType, EdgeType
from util.graph import (
    extract_support_layer,
    get_crossing_port_edges,
    get_port_edges,
    edge_filter_ports,
)
from util.geometry import (
    get_angle,
    offset_edge,
    get_segment_circle_intersection,
    centroid,
    biarc,
)
from util.config import (
    DARK_MODE,
    DRAW_GLYPHS,
    DRAW_GLYPHS_OVER_LINES,
    STROKE_WIDTH,
    LINE_GAP,
    CELL_SIZE_PX,
    GLYPH_SIZE_PX,
    MARGINS,
    NODE_CIRCLE_RADIUS,
    DRAW_DEG1_MARKS,
    DRAW_DEG1_MARK_SIZE_PX,
    DRAW_HUBS,
    DRAW_GRID_CELLS,
    SET_COLORS,
    GRID_WIDTH,
    GRID_HEIGHT,
    WRITE_DIR,
)


def draw_biarc(line, barc):
    p0 = barc["p"][1]
    line.M(p0[0], p0[1])

    if barc["type"] == "C":
        line.C(*barc["args"])
    else:
        arc1 = barc["arc1"]
        arc2 = barc["arc2"]
        arc1_x, arc1_y, arc1_r, arc1_sweep, arc1_large = arc1
        arc2_x, arc2_y, arc2_r, arc2_sweep, arc2_large = arc2

        line.A(arc1_r, arc1_r, 0, arc1_large, arc1_sweep, arc1_x, arc1_y)
        line.M(arc1_x, arc1_y)
        line.A(arc2_r, arc2_r, 0, arc2_large, arc2_sweep, arc2_x, arc2_y)
        line.M(arc2_x, arc2_y)


def draw_embedding(X, filename, **kwargs):
    r = kwargs.get("r", 0.5)
    w = kwargs.get("width", 100)
    h = kwargs.get("height", 100)
    mx = kwargs.get("margin_x", 5)
    my = kwargs.get("margin_y", 5)

    Y = (mx, my) + ((X - np.min(X)) / np.ptp(X)) * (w - mx, h - my)

    geometries = []
    for i in range(len(Y)):
        x, y = Y[i, :]
        geometries.append(
            svg.Circle(cx=x, cy=y, r=r, fill="white" if DARK_MODE.get() else "black")
        )
    with open(os.path.join(WRITE_DIR.get(), filename), "w") as f:
        f.write(draw_svg(geometries, w + 2 * mx, h + 2 * my))


def draw_svg(geometries, width, height):
    d = svg.Drawing(width, height, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


def draw_support(instance, M, layer=0):
    dir = WRITE_DIR.get()
    path = os.path.join(dir, f"support_{layer}.svg")

    geometries = []
    # project nodes
    mx, my = MARGINS.get()
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (x * CELL_SIZE_PX.get() + mx, -y * CELL_SIZE_PX.get() + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            c = svg.Circle(
                cx=px,
                cy=py,
                r=NODE_CIRCLE_RADIUS.get(),
                fill="white" if DARK_MODE.get() else "black",
            )
            c.append_title(M.nodes[i]["label"])
            geometries.append(c)

    M_ = nx.subgraph_view(
        M,
        filter_edge=lambda u, v: edge_filter_ports(M, u, v, same_centers=False),
    )
    for u, v in M_.edges():
        src = M.nodes[M.nodes[u]["belongs_to"]]["pos"]
        tgt = M.nodes[M.nodes[v]["belongs_to"]]["pos"]

        ux, uv = src
        vx, vy = tgt

        geometries.append(
            svg.Line(
                sx=ux,
                sy=uv,
                ex=vx,
                ey=vy,
                data_sets=str(M.edges[u,v]['sets']),
                stroke_width=STROKE_WIDTH.get(),
                stroke="white" if DARK_MODE.get() else "black",
            )
        )

    w = GRID_WIDTH.get() * CELL_SIZE_PX.get()
    h = GRID_HEIGHT.get() * CELL_SIZE_PX.get()
    with open(path, "w") as f:
        f.write(draw_svg(geometries, w, h))


def geometrize(instance, L, element_set_partition, layer=0):
    geometries = []
    mx, my = MARGINS.get()

    M = extract_support_layer(L, layer)

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (x * CELL_SIZE_PX.get() + mx, -y * CELL_SIZE_PX.get() + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER:
            if M.nodes[i]["occupied"]:
                if DRAW_GLYPHS.get() and "glyph" in M.nodes[i]:
                    w, h = (GLYPH_SIZE_PX.get(), GLYPH_SIZE_PX.get())
                    img = svg.Image(
                        px - w / 2,
                        py - h / 2,
                        width=w,
                        height=h,
                        path=M.nodes[i]["glyph"],
                    )
                    img.append_title(M.nodes[i]["label"])
                    geometries.append(img)
                else:
                    c = svg.Circle(cx=px, cy=py, r=NODE_CIRCLE_RADIUS.get())
                    c.append_title(M.nodes[i]["label"])
                    geometries.append(c)
            elif DRAW_GRID_CELLS.get():
                c = svg.Circle(
                    cx=px, cy=py, r=NODE_CIRCLE_RADIUS.get() / 8, fill="gray"
                )
                geometries.append(c)

    for u, v in M.edges():
        if not edge_filter_ports(M, u, v, same_centers=False):
            continue

        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        paths_at_edge = M.edges[(u, v)]["oeb_order"][(u, v)]
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -LINE_GAP.get()
        M.edges[(u, v)]["edge_pos"] = {}
        for i, set_id in enumerate(paths_at_edge):
            offset_dir = 3 * math.pi / 2
            offset_length = centering_offset + i * LINE_GAP.get()
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
                (M.nodes[M.nodes[u]["belongs_to"]]["pos"], NODE_CIRCLE_RADIUS.get()),
            )
            v_intersect = get_segment_circle_intersection(
                (upos, vpos),
                (M.nodes[M.nodes[v]["belongs_to"]]["pos"], NODE_CIRCLE_RADIUS.get()),
            )

            line = svg.Path(
                **{
                    "close": False,
                    "stroke_width": STROKE_WIDTH.get(),
                    "fill": "none",
                    "stroke": SET_COLORS.get()[instance["set_ftb_order"].index(set_id)],
                }
            )
            line.M(*u_intersect)
            line.L(*v_intersect)
            geometries.append(line)

        # this draws all the connections at non-occupied nodes
        # because at occupied nodes the lines go through the center node
        for u, v, d in M.edges(data=True):
            if (
                not edge_filter_ports(M, u, v, same_centers=True)
                or set_id not in M.edges[(u, v)]["sets"]
            ):
                continue
            uparent = M.nodes[u]["belongs_to"]

            u_adjacent = list([x if w == u else w for w, x in M.edges(nbunch=[u])])
            uu = None
            for x in u_adjacent:
                xparent = M.nodes[x].get("belongs_to")
                if xparent is not None and xparent != uparent:
                    uu = x
                    break

            v_adjacent = list([x if w == v else w for w, x in M.edges(nbunch=[v])])
            vv = None
            for x in v_adjacent:
                xparent = M.nodes[x].get("belongs_to")
                if xparent is not None and xparent != uparent:
                    vv = x
                    break

            if uu is None or vv is None:
                print("howw")
                continue

            uupos, upos = M.edges[(uu, u)]["edge_pos"][set_id][(uu, u)]
            vpos, vvpos = M.edges[(v, vv)]["edge_pos"][set_id][(v, vv)]

            u_intersect = get_segment_circle_intersection(
                (uupos, upos), (M.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS.get())
            )
            v_intersect = get_segment_circle_intersection(
                (vpos, vvpos), (M.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS.get())
            )

            uu_u_center = centroid([uupos, upos])
            vv_v_center = centroid([vvpos, vpos])

            line = svg.Path(
                **{
                    "close": False,
                    "stroke_width": STROKE_WIDTH.get(),
                    "fill": "none",
                    "stroke": SET_COLORS.get()[instance["set_ftb_order"].index(set_id)],
                }
            )
            barc = biarc(uu_u_center, u_intersect, v_intersect, vv_v_center)
            draw_biarc(line, barc)
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

            all_used_ports = []
            for p in [
                p
                for p in nx.neighbors(M, node)
                if M.nodes[p]["node"] == NodeType.PORT
                and M.nodes[p]["belongs_to"] == node
            ]:
                edges = [
                    (a, b) for a, b, d in M.edges(nbunch=p, data=True) if "sets" in d
                ]
                if len(edges) > 0:
                    all_used_ports.append(p)

            used_ports = set()
            outward_edge_at_port = dict()
            for port in all_used_ports:
                p_adjacent = list(
                    [
                        x
                        for w, x in M.edges(nbunch=[port])
                        if M.nodes[x]["node"] == NodeType.PORT
                    ]
                )
                for x in p_adjacent:
                    edge = M.edges[(port, x)]
                    if "sets" in edge and set_id in edge["sets"]:
                        used_ports = used_ports.union(set([port]))
                        outward_edge_at_port[port] = (port, x)

            if len(used_ports) < 1:
                # cannot happen actually
                continue
            if len(used_ports) == 1 and DRAW_DEG1_MARKS.get():
                # this is a deg 1 node for this set
                # idk, could connect to center
                # or maybe draw a small mark?
                a = used_ports.pop()
                _, b = outward_edge_at_port[a]
                apos, bpos = M.edges[(a, b)]["edge_pos"][set_id][(a, b)]
                cx, cy = get_segment_circle_intersection(
                    (apos, bpos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS.get())
                )
                circle = svg.Circle(
                    cx=cx,
                    cy=cy,
                    r=DRAW_DEG1_MARK_SIZE_PX.get(),
                    fill=SET_COLORS.get()[instance["set_ftb_order"].index(set_id)],
                )
                circle.append_title(set_id)
                geometries.append(circle)
                continue

            all_edges_at_ports = get_port_edges(L, node)
            port_port_edges = [
                (a, b, L.edges[(a, b, EdgeType.PHYSICAL)]["weight"])
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
                        "stroke_width": STROKE_WIDTH.get(),
                        "fill": "none",
                        "data_weight": G_node.edges[a, b]["weight"],
                        "stroke": SET_COLORS.get()[
                            instance["set_ftb_order"].index(set_id)
                        ],
                    }
                )
                _, v = outward_edge_at_port[a]
                _, x = outward_edge_at_port[b]

                apos, vpos = M.edges[(a, v)]["edge_pos"][set_id][(a, v)]
                bpos, xpos = M.edges[(b, x)]["edge_pos"][set_id][(b, x)]

                a_intersect = get_segment_circle_intersection(
                    (vpos, apos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS.get())
                )
                b_intersect = get_segment_circle_intersection(
                    (bpos, xpos), (M.nodes[node]["pos"], NODE_CIRCLE_RADIUS.get())
                )

                av_center = centroid([apos, vpos])
                bx_center = centroid([bpos, xpos])

                barc = biarc(av_center, a_intersect, b_intersect, bx_center)
                draw_biarc(line, barc)

                geometries.append(line)

    if False:
        # TODO new approach
        # find the root (first element in partition)
        # iterate over other elements in partition
        # shortest path
        # at each occupied node draw the necessary segment
        for elements, sets in element_set_partition:
            for set_id in sets:
                G_ = nx.subgraph_view(
                    M,
                    filter_edge=lambda u, v: set_id in M.edges[(u, v)]["sets"],
                    filter_node=lambda n: True,
                )

                root = elements[0]
                rootpos = instance["glyph_positions"][layer][
                    instance["elements_inv"][root]
                ]

                for element in elements[1:]:
                    elpos = instance["glyph_positions"][layer][
                        instance["elements_inv"][element]
                    ]

                    shortest_path = path_to_edges(nx.shortest_path(G_, rootpos, elpos))
                    edge_pairs = [
                        (e1, e2)
                        for e1, e2 in pairwise(shortest_path)
                        if G_.nodes[e1[1]]["node"] == NodeType.CENTER
                    ]

                    for e1, e2 in edge_pairs:
                        u, v = e1
                        v, w = e2

                        uparent = v
                        u_adjacent = list(
                            [x if w == u else w for w, x in M.edges(nbunch=[u])]
                        )
                        uu = None
                        for x in u_adjacent:
                            xparent = G_.nodes[x].get("belongs_to")
                            if xparent is not None and xparent != uparent:
                                uu = x
                                break

                        w_adjacent = list(
                            [x if y == w else y for y, x in G_.edges(nbunch=[w])]
                        )
                        ww = None
                        for x in w_adjacent:
                            xparent = G_.nodes[x].get("belongs_to")
                            if xparent is not None and xparent != uparent:
                                ww = x
                                break

                        if uu is None or ww is None:
                            print("howw")
                            continue

                        if (
                            set_id not in G_.edges[(uu, u)]["edge_pos"]
                            or set_id not in G_.edges[(w, ww)]["edge_pos"]
                        ):
                            continue

                        uupos, upos = G_.edges[(uu, u)]["edge_pos"][set_id][(uu, u)]
                        wpos, wwpos = G_.edges[(w, ww)]["edge_pos"][set_id][(w, ww)]

                        u_intersect = get_segment_circle_intersection(
                            (uupos, upos),
                            (G_.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS.get()),
                        )
                        v_intersect = get_segment_circle_intersection(
                            (wpos, wwpos),
                            (G_.nodes[uparent]["pos"], NODE_CIRCLE_RADIUS.get()),
                        )

                        uu_u_center = centroid([uupos, upos])
                        vv_v_center = centroid([wwpos, wpos])

                        line = svg.Path(
                            **{
                                "close": False,
                                "stroke_width": STROKE_WIDTH.get(),
                                "fill": "none",
                                "stroke": SET_COLORS.get()[
                                    instance["set_ftb_order"].index(set_id)
                                ],
                            }
                        )
                        barc = biarc(uu_u_center, u_intersect, v_intersect, vv_v_center)
                        draw_biarc(line, barc)
                        geometries.append(line)

    if DRAW_HUBS.get():
        hubs = [n for n in M.nodes() if M.degree[n] > 0]
        for hub in hubs:
            cx, cy = M.nodes[hub]["pos"]
            r = 1 / 16 * CELL_SIZE_PX.get()
            geometries.append(
                svg.Circle(cx, cy, r, fill="none", stroke="gray", stroke_width=1)
            )

    if DRAW_GLYPHS_OVER_LINES.get():
        geometries = list(reversed(geometries))

    return geometries
