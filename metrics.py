import click
import json
import re
import os
import networkx as nx
import math
from itertools import combinations
from util.geometry import get_angle
from util.enums import EdgeType, NodeType, PortDirs
from util.graph import edge_filter_ports, extract_support_layer, are_port_edges_crossing

DOUBLE_TUPLE_REGEX = re.compile(
    "\\((?P<u>\\(-?\\d+(?:\\.\\d+)?, -?\\d+(?:\\.\\d+)?\\)), (?P<v>\\(-?\\d+(?:\\.\\d+)?, -?\\d+(?:\\.\\d+)?\\))\\)"
)


def figure_out_num_layers(G):
    max_layer = 0
    for u, v, k in G.edges(keys=True):
        layer, _ = k
        if layer > max_layer:
            max_layer = layer
    return max_layer


def figure_out_size(G):
    max_x = 0
    max_y = 0
    for n in G.nodes():
        x, y = n
        if x > max_x:
            max_x = x
        if y > max_y:
            max_y = y

    return (max_x, max_y)


def compute_crossings_outside(G, size=(10, 10), what="edges"):
    # between diagional nodes
    result = 0
    m, n = size

    if what not in ["edges", "lines"]:
        raise BaseException(f"unsure what to count: {what}")

    for u, v, d in G.edges(data=True):
        if not edge_filter_ports(
            G, u, v, same_centers=False, possibly_with_center=False
        ):
            continue

        uparent = G.nodes[u]["belongs_to"]
        vparent = G.nodes[v]["belongs_to"]

        ux, uy = uparent
        vx, vy = vparent
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
            crossing_edge = ((ux + dx, uy), (ux, uy + dy))

            # found crossing edge with centers as nodes
            # that won't be in the support graph so we have to look for
            # the same edge using the centers' repsective ports

            # so let the new edge be (w,x)
            w, x = crossing_edge
            # then get the angle from w to x, select appropriate port
            angle = get_angle(w, x) * 180 / math.pi

            wport = "ne"

            if int(angle) == 45:
                # w is southwest of x, use ne port and sw port
                pass
            elif int(angle) == 315:
                # w is northwest of x, use se port and nw port
                wport = "se"
            else:
                raise BaseException(f"did not expect angle {angle}")

            xport = "nw" if wport == "se" else "sw"

            xp_node = list(
                [
                    p
                    for p in G.nodes()
                    if G.nodes[p]["node"] == NodeType.PORT
                    and G.nodes[p]["belongs_to"] == x
                    and G.nodes[p]["port"] == xport
                ]
            )[0]
            wp_node = list(
                [
                    p
                    for p in G.nodes()
                    if G.nodes[p]["node"] == NodeType.PORT
                    and G.nodes[p]["belongs_to"] == w
                    and G.nodes[p]["port"] == wport
                ]
            )[0]

            if (wp_node, xp_node) in G.edges:
                result += (
                    1
                    if what == "edges"
                    else len(G.edges[uparent, vparent]["sets"])
                    * len(G.edges[wp_node, xp_node]["sets"])
                )

    return result


def compute_crossings_inside(G, what="edges"):
    # in unoccupied nodes
    result = 0

    if what not in ["edges", "lines"]:
        raise BaseException(f"unsure what to count: {what}")

    for n, d in G.nodes(data=True):
        if d["node"] != NodeType.CENTER or d["occupied"]:
            continue

        all_ports = set(
            [
                p
                for p in G.nodes()
                if G.nodes[p]["node"] == NodeType.PORT and G.nodes[p]["belongs_to"] == n
            ]
        )

        if len(all_ports) < 2:
            continue

        all_edges_at_ports = set()
        for p1, p2 in combinations(all_ports, 2):
            if (p1, p2) in G.edges:
                all_edges_at_ports.add((p1, p2))

        for uv, wx in combinations(all_edges_at_ports, 2):
            u, v = uv
            w, x = wx
            if are_port_edges_crossing(
                G.nodes[u],
                G.nodes[v],
                G.nodes[w],
                G.nodes[x],
                cross_when_node_shared=False,
            ):
                result += (
                    1
                    if what == "edges"
                    else len(G.edges[u, v]["sets"]) * len(G.edges[w, x]["sets"])
                )
    return result


def compute_total_length(G, what="edges"):
    result = 0

    if what not in ["edges", "lines"]:
        raise BaseException(f"unsure what to count: {what}")

    for u, v, d in G.edges(data=True):
        if edge_filter_ports(G, u, v, possibly_with_center=False, same_centers=False):
            result += 1 if what == "edges" else len(d["sets"])

    return result


def compute_metrics(G):
    result = []
    layers = figure_out_num_layers(G)
    for layer in range(layers + 1):
        G_ = extract_support_layer(G, layer)
        # print("layer", layer)
        result.append(
            {
                "total_lines": compute_total_length(G_, what="lines"),
                "total_edges": compute_total_length(G_, what="edges"),
                "total_line_crossings": compute_crossings_inside(G_, what="lines")
                + compute_crossings_outside(G_, what="lines", size=figure_out_size(G)),
                "total_edge_crossings": compute_crossings_inside(G_, what="edges")
                + compute_crossings_outside(G_, what="edges", size=figure_out_size(G)),
            }
        )
    return result


@click.command()
@click.option(
    "--graph-path", default="serialized.json", help="path to serialized json graph"
)
@click.option("--write-dir", default="./", help="where to write the results to")
def main(graph_path, write_dir):
    with open(graph_path, "r") as f:
        jayson = json.load(f)

        # convert stuff back that python doesn't like
        # mostly lists as dict keys, so restore all the tuples
        for node in jayson["nodes"]:
            node["id"] = tuple(node["id"])
            node["pos"] = tuple(node["pos"])
            if node["node"] == 0:
                node["logpos"] = tuple(node["logpos"])
            if node["node"] == 1:
                node["belongs_to"] = tuple(node["belongs_to"])

        for edge in jayson["links"]:
            edge["source"] = tuple(edge["source"])
            edge["target"] = tuple(edge["target"])
            edge["key"] = tuple(edge["key"])

        G = nx.node_link_graph(jayson)

        # now that everything is python again, restore the datatypes that
        # json didn't like
        # sets, enums, and more tuples!
        # note: enum in key not handled, i assume since it's an intenum comparisons should be fine everywhere anyhow?
        for u, v, k, d in G.edges(keys=True, data=True):
            d["edge"] = EdgeType.SUPPORT
            if "sets" in d:
                d["sets"] = set(d["sets"])
            if "oeb_order" in d:
                keys = list(d["oeb_order"])
                for key in keys:
                    match = DOUBLE_TUPLE_REGEX.match(key)
                    ux = float(match.group("u")[1:-1].split(",")[0])
                    uy = float(match.group("u")[1:-1].split(",")[1])
                    vx = float(match.group("v")[1:-1].split(",")[0])
                    vy = float(match.group("v")[1:-1].split(",")[1])
                    d["oeb_order"][((ux, uy), (vx, vy))] = d["oeb_order"][key]
                    del d["oeb_order"][key]

        for n, d in G.nodes(data=True):
            d["node"] = NodeType.CENTER if d["node"] == 0 else NodeType.PORT

        metrics = compute_metrics(G)
        print(metrics)
        with open(os.path.join(write_dir, "metrics.json"), "w") as f_out:
            json.dump(metrics, f_out)


if __name__ == "__main__":
    main()
