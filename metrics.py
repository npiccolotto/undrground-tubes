import click
import json
import re
import networkx as nx
from itertools import combinations
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


def compute_crossings_outside(G):
    # between diagional nodes
    pass


def compute_crossings_inside(G):
    # in unoccupied nodes
    result = 0
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
                result += 1
    return result


def compute_total_edges_used(G):
    result = 0

    for u, v, d in G.edges(data=True):
        if edge_filter_ports(G, u, v, possibly_with_center=False, same_centers=False):
            result += 1

    return result


def compute_total_line_length(G):
    result = 0

    for u, v, d in G.edges(data=True):
        num_sets_at_edge = len(d["sets"])
        if edge_filter_ports(G, u, v, possibly_with_center=False, same_centers=False):
            result += num_sets_at_edge

    return result


def compute_metrics(G):
    result = []
    layers = figure_out_num_layers(G)
    for layer in range(layers + 1):
        G_ = extract_support_layer(G, layer)
        #print("layer", layer)
        result.append(
            {
                "total_line_length": compute_total_line_length(G_),
                "total_edges": compute_total_edges_used(G_),
                "total_crossing_inside": compute_crossings_inside(G_),
            }
        )
    return result


@click.command()
@click.option(
    "--graph-path", default="serialized.json", help="path to serialized json graph"
)
def main(graph_path):
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

        print(compute_metrics(G))


if __name__ == "__main__":
    main()
