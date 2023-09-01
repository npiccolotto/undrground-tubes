import json
import math
import subprocess
import sys
import copy
import click
from collections import defaultdict
from itertools import chain, combinations, pairwise, product
import contextvars
from functools import partial

import drawsvg as svg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from util.config import (
    DRAW_GLYPHS,
    CELL_SIZE_PX,
    SUB_SUPPORT_GROUPING,
    SUB_SUPPORT_TYPE,
    STRATEGY
)
from util.bundle import bundle_lines
from util.collections import (
    group_by_intersection_group,
    group_by_set,
    invert_list,
    list_of_lists_to_set_system_dict,
)
from util.draw import (
    draw_svg,
    geometrize,
    draw_support,
    draw_embedding,
)
from util.enums import EdgePenalty, EdgeType, NodeType, PortDirs
from util.geometry import (
    dist_euclidean,
    get_closest_point,
    logical_coords_to_physical,
)
from util.graph import (
    path_to_edges,
    extract_support_layer,
    get_ports,
)
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
                    weight=EdgePenalty.HOP - 1 + length_penalty,
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


def read_instance(directory, name):
    with open(f"{directory}/{name}.json") as f:
        data = json.load(f)
    elements = data["E"]
    sets = data["S"]
    inst = {
        "elements": elements,
        "elements_inv": invert_list(elements),
        "sets": sets,
        "set_system": list_of_lists_to_set_system_dict(elements, data["SR"]),
        "D_EA": data["EA"],
        "D_SR": data["SA"],
        "set_ftb_order": list(sorted(sets)),
        # pipeline config
        "dr_method": "mds",
        "dr_gridification": "hagrid",  #  'hagrid' or 'dgrid'
    }
    if "glyph_ids" in data:
        inst["glyph_ids"] = data["glyph_ids"]
    print(list_of_lists_to_set_system_dict(elements, data["SR"]))
    return inst


def render(
    read_dir,
    write_dir,
    dataset,
    num_weights,
    grid_width,
    grid_height,
):
    instance = read_instance(read_dir, dataset)
    lattice_type = "sqr"
    instance["grid_x"] = grid_width
    instance["grid_y"] = grid_height
    instance["num_layers"] = num_weights

    with timing("layout"):
        instance["glyph_positions"] = layout_dr_multiple(
            instance["D_EA"],
            instance["D_SR"],
            m=grid_width,
            n=grid_height,
            num_samples=num_weights,
        )

    with timing("routing"):
        G = get_routing_graph(lattice_type, (grid_width, grid_height))
        G = add_glyphs_to_nodes(instance, G)

        element_set_partition = (
            group_by_intersection_group(instance["set_system"])
            if SUB_SUPPORT_GROUPING.get() == "intersection-group"
            else group_by_set(instance["set_system"])
        )
        element_set_partition = sorted(
            element_set_partition, key=lambda x: len(x[1]), reverse=True
        )

        print("element partition", element_set_partition)

        if STRATEGY.get() == "opt":
            L = route_multilayer_ilp(
                instance,
                nx.subgraph_view(
                    G,
                    filter_edge=lambda u, v, k: k == EdgeType.CENTER,
                    filter_node=lambda n: G.nodes[n]["node"] == NodeType.CENTER,
                ),
                element_set_partition,
            )
        else:
            L = route_multilayer_heuristic(
                instance,
                G,
                element_set_partition,
            )

    for layer in range(num_weights):
        M = extract_support_layer(L, layer)
        draw_support(instance, M, write_dir, layer)

    with timing("bundle lines"):
        L = bundle_lines(instance, L)

    if STRATEGY.get() == "opt":
        # merge grid graph data (ports at nodes) with line graph data (routing and bundling info)
        # we have
        # G = multigraph with center and physical nodes/edges
        # L = a multigraph with (layer, support) edges and center nodes
        for layer in range(num_weights):
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

    with timing("serialize graph"):
        # avoid referencing issues
        L_ = copy.deepcopy(L)

        for u, v, d in L_.edges(data=True):
            # convert stuff that json doesn't like
            # sets
            if "sets" in d:
                d["sets"] = list(d["sets"])
            # tuple keys
            if "oeb_order" in d:
                keys = list(d["oeb_order"])
                for key in keys:
                    d["oeb_order"][str(key)] = [*d["oeb_order"][key]]
                    del d["oeb_order"][key]

        with open(f"{write_dir}serialized.json", "w") as f:
            json.dump(nx.node_link_data(L_), f)

    with timing("draw+write svg"):
        for layer in range(num_weights):
            # TODO dump L in some format we can read again
            L.add_edges_from(
                [
                    (u, v, k, d)
                    for u, v, k, d in G.edges(keys=True, data=True)
                    if k == EdgeType.PHYSICAL
                ]
            )
            geometries = geometrize(instance, L, element_set_partition, layer=layer)
            img = draw_svg(
                geometries,
                grid_width * CELL_SIZE_PX.get(),
                grid_height * CELL_SIZE_PX.get(),
            )
            with open(f"{write_dir}drawing_{layer}.svg", "w") as f:
                f.write(img)
                f.flush()

    print("Done.")


@click.command()
@click.option(
    "--read-dir", default="./data", help="directory to read the datasets from"
)
@click.option("--write-dir", default="./", help="directory to write the output to")
@click.option("--dataset", default="imdb/imdb_10", help="dataset to load")
@click.option("--strategy", default='heuristic', type=click.Choice(["opt", "heuristic"], case_sensitive=False))
@click.option("--grid-width", "-w", default=10, help="grid width as # cols")
@click.option("--grid-height", "-h", default=10, help="grid height as # rows")
@click.option(
    "--num-weights", default=2, help="how many samples between 0..1 to use for weights"
)
@click.option(
    "--support-type",
    type=click.Choice(["path", "steiner-tree"], case_sensitive=False),
    default="path",
    help="the support type",
)
@click.option(
    "--support-partition",
    type=click.Choice(["set", "intersection-group"], case_sensitive=False),
    default="intersection-group",
    help="the partition type",
)
@click.option(
    "--draw-glyphs/--no-draw-glyphs",
    type=bool,
    default=True,
    help="draw glyphs (true) or circles (false)",
)
def vis(
    read_dir,
    write_dir,
    dataset,
    strategy,
    num_weights,
    grid_width,
    grid_height,
    support_type,
    support_partition,
    draw_glyphs,
):

    DRAW_GLYPHS.set(draw_glyphs)
    SUB_SUPPORT_GROUPING.set(support_partition)
    SUB_SUPPORT_TYPE.set(support_type)
    STRATEGY.set(strategy)

    fun = partial(
        render,
        read_dir,
        write_dir,
        dataset,
        num_weights,
        grid_width,
        grid_height,
    )
    ctx = contextvars.copy_context()
    ctx.run(fun)


if __name__ == "__main__":
    vis()
