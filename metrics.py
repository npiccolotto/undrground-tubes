import click
import json
import re
import networkx as nx
from util.enums import EdgeType, NodeType
from util.graph import edge_filter_ports

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


def compute_total_line_length(G):
    result = {}

    for u, v, k, d in G.edges(keys=True, data=True):
        layer, _ = k
        if layer not in result:
            result[layer] = {"within-node": 0, "between-node": 0}
        num_sets_at_edge = len(d["sets"])
        if edge_filter_ports(G, u, v, possibly_with_center=False, same_centers=True):
            result[layer]["within-node"] += num_sets_at_edge
        elif edge_filter_ports(G, u, v, possibly_with_center=True, same_centers=True):
            pass
        else:
            result[layer]["between-node"] += num_sets_at_edge

    return result


def compute_metrics(G):
    total_line_length = compute_total_line_length(G)

    return {"total_line_length": total_line_length}


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
