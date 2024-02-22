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
    orient_edge_node_inside,
    get_outgoing_edge_to_other_center_from_port,
)
from util.geometry import (
    get_angle,
    offset_edge,
    get_segment_circle_intersection,
    centroid,
    biarc,
)
from util.route import determine_router
from util.config import config_vars, get_grid


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
            svg.Circle(
                cx=x,
                cy=y,
                r=r,
                fill="white" if config_vars["draw.darkmode"].get() else "black",
            )
        )
    with open(os.path.join(config_vars["general.writedir"].get(), filename), "w") as f:
        f.write(draw_svg(geometries, w + 2 * mx, h + 2 * my))


def draw_svg(geometries, width, height):
    d = svg.Drawing(width, height, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


def draw_support(instance, M, layer=0):
    dir = config_vars["general.writedir"].get()
    path = os.path.join(dir, f"support_{layer}.svg")

    geometries = []
    # project nodes
    mx, my = (
        config_vars["draw.marginhorizontal"].get(),
        config_vars["draw.marginvertical"].get(),
    )

    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (
            x * config_vars["draw.cellsizepx"].get() + mx,
            -y * config_vars["draw.cellsizepx"].get() + my,
        )
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            c = svg.Circle(
                cx=px,
                cy=py,
                r=config_vars["draw.nodecircleradius"].get(),
                fill="white" if config_vars["draw.darkmode"].get() else "black",
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
                data_sets=str(M.edges[u, v]["sets"]),
                stroke_width=config_vars["draw.strokewidth"].get(),
                stroke="white" if config_vars["draw.darkmode"].get() else "black",
            )
        )

    w, h = get_grid()
    w = w * config_vars["draw.cellsizepx"].get()
    h = h * config_vars["draw.cellsizepx"].get()
    with open(path, "w") as f:
        f.write(draw_svg(geometries, w, h))


def geometrize(instance, L, element_set_partition, layer=0):
    geometries = []
    mx, my = (
        config_vars["draw.marginhorizontal"].get(),
        config_vars["draw.marginvertical"].get(),
    )
    set_colors = instance["set_colors"]
    M = extract_support_layer(L, layer)

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (
            x * config_vars["draw.cellsizepx"].get() + mx,
            -y * config_vars["draw.cellsizepx"].get() + my,
        )
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER:
            if M.nodes[i]["occupied"]:
                if config_vars["draw.drawglyphs"].get() and "glyph" in M.nodes[i]:
                    w, h = (
                        config_vars["draw.glyphsizepx"].get(),
                        config_vars["draw.glyphsizepx"].get(),
                    )
                    img = svg.Image(
                        px - w / 2,
                        py - h / 2,
                        width=w,
                        height=h,
                        data_node=True,
                        data_id=M.nodes[i]['label'],
                        path=M.nodes[i]["glyph"],
                    )
                    if config_vars["draw.glyphtitle"].get():
                        img.append_title(M.nodes[i]["label"])
                    geometries.append(img)
                else:
                    c = svg.Circle(
                        cx=px, cy=py, r=config_vars["draw.nodecircleradius"].get()
                    )
                    if config_vars["draw.glyphtitle"].get():
                        c.append_title(M.nodes[i]["label"])
                    geometries.append(c)
            elif config_vars["draw.drawgridcells"].get():
                c = svg.Circle(
                    cx=px,
                    cy=py,
                    r=config_vars["draw.nodecircleradius"].get() / 8,
                    fill="gray",
                )
                geometries.append(c)

    for u, v in M.edges():
        if not edge_filter_ports(M, u, v, same_centers=False):
            continue

        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        paths_at_edge = M.edges[(u, v)]["oeb_order"][(u, v)]
        linegap = config_vars["draw.linegap"].get()
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -linegap
        M.edges[(u, v)]["edge_pos"] = {}
        for i, set_id in enumerate(paths_at_edge):
            offset_dir = 3 * math.pi / 2
            offset_length = centering_offset + i * linegap
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
                (
                    M.nodes[M.nodes[u]["belongs_to"]]["pos"],
                    config_vars["draw.nodecircleradius"].get(),
                ),
            )
            v_intersect = get_segment_circle_intersection(
                (upos, vpos),
                (
                    M.nodes[M.nodes[v]["belongs_to"]]["pos"],
                    config_vars["draw.nodecircleradius"].get(),
                ),
            )

            line = svg.Path(
                **{
                    "close": False,
                    "data-weight": L.edges[u, v, EdgeType.PHYSICAL]["weight"],
                    "stroke_width": config_vars["draw.strokewidth"].get(),
                    "fill": "none",
                    "stroke": set_colors[set_id],
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

            try:
                _,uu = get_outgoing_edge_to_other_center_from_port(M, u)
                _,vv = get_outgoing_edge_to_other_center_from_port(M, v)
            except:
                # ?? sometimes there are too many edges idk
                continue

            if (
                set_id not in M.edges[(uu, u)]["sets"]
                or set_id not in M.edges[(v, vv)]["sets"]
            ):
                print('whyy')
                continue
            uupos, upos = M.edges[(uu, u)]["edge_pos"][set_id][(uu, u)]
            vpos, vvpos = M.edges[(v, vv)]["edge_pos"][set_id][(v, vv)]

            u_intersect = get_segment_circle_intersection(
                (uupos, upos),
                (M.nodes[uparent]["pos"], config_vars["draw.nodecircleradius"].get()),
            )
            v_intersect = get_segment_circle_intersection(
                (vpos, vvpos),
                (M.nodes[uparent]["pos"], config_vars["draw.nodecircleradius"].get()),
            )

            uu_u_center = centroid([uupos, upos])
            vv_v_center = centroid([vvpos, vpos])

            line = svg.Path(
                **{
                    "close": False,
                    "data-weight": L.edges[u, v, EdgeType.PHYSICAL]["weight"],
                    "stroke_width": config_vars["draw.strokewidth"].get(),
                    "fill": "none",
                    "stroke": set_colors[set_id],
                }
            )
            barc = biarc(uu_u_center, u_intersect, v_intersect, vv_v_center)
            draw_biarc(line, barc)
            geometries.append(line)

        # so this then draws connections within occupied nodes
        if False:
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
                        for a, b, d in M.edges(nbunch=p, data=True)
                        if "sets" in d
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
                            and M.nodes[x]["belongs_to"] != node
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
                if len(used_ports) == 1 and config_vars["draw.drawdeg1marks"].get():
                    # this is a deg 1 node for this set
                    # idk, could connect to center
                    # or maybe draw a small mark?
                    a = used_ports.pop()
                    _, b = orient_edge_node_inside(outward_edge_at_port[a], a)
                    apos, bpos = M.edges[(a, b)]["edge_pos"][set_id][(a, b)]
                    cx, cy = get_segment_circle_intersection(
                        (apos, bpos),
                        (
                            M.nodes[node]["pos"],
                            config_vars["draw.nodecircleradius"].get(),
                        ),
                    )
                    circle = svg.Circle(
                        cx=cx,
                        cy=cy,
                        r=config_vars["draw.deg1marksize"].get(),
                        fill=set_colors[set_id],
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

                within_node_connections = nx.minimum_spanning_tree(
                    G_node, weight="weight"
                )
                for a, b in within_node_connections.edges():
                    line = svg.Path(
                        **{
                            "close": False,
                            "data-weight": L.edges[a, b, EdgeType.PHYSICAL]["weight"],
                            "stroke_width": config_vars["draw.strokewidth"].get(),
                            "fill": "none",
                            "stroke": set_colors[set_id],
                        }
                    )
                    _, v = orient_edge_node_inside(outward_edge_at_port[a], a)
                    _, x = orient_edge_node_inside(outward_edge_at_port[b], b)

                    apos, vpos = M.edges[(a, v)]["edge_pos"][set_id][(a, v)]
                    bpos, xpos = M.edges[(b, x)]["edge_pos"][set_id][(b, x)]

                    a_intersect = get_segment_circle_intersection(
                        (vpos, apos),
                        (
                            M.nodes[node]["pos"],
                            config_vars["draw.nodecircleradius"].get(),
                        ),
                    )
                    b_intersect = get_segment_circle_intersection(
                        (bpos, xpos),
                        (
                            M.nodes[node]["pos"],
                            config_vars["draw.nodecircleradius"].get(),
                        ),
                    )

                    av_center = centroid([apos, vpos])
                    bx_center = centroid([bpos, xpos])

                    barc = biarc(av_center, a_intersect, b_intersect, bx_center)
                    draw_biarc(line, barc)

                    geometries.append(line)

    if config_vars["draw.drawhubs"].get():
        hubs = [n for n in M.nodes() if M.degree[n] > 0]
        for hub in hubs:
            cx, cy = M.nodes[hub]["pos"]
            r = 1 / 16 * config_vars["draw.cellsizepx"].get()
            geometries.append(
                svg.Circle(cx, cy, r, fill="none", stroke="gray", stroke_width=1)
            )

    if config_vars["draw.drawglyphsoverlines"].get():
        geometries = list(reversed(geometries))

    return geometries


def cosmetic_post_processing(instance, G):
    return G
