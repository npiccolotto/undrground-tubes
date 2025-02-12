import json
import time
import os
import math
import traceback
import copy
import click
from itertools import chain, combinations, pairwise, product
import contextvars
from functools import partial

import drawsvg as svg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from util.metrics import compute_metrics as calc_metrics
from util.config import config_vars, get_grid
from util.bundle import bundle_lines_gg, bundle_lines_lg
from util.collections import (
    group_by_intersection_group,
    group_by_set,
    invert_list,
    list_of_lists_to_set_system_dict,
    select_sets,
)
from util.draw import (
    draw_svg,
    geometrize,
    draw_support,
    draw_embedding,
    cosmetic_post_processing,
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
    convert_line_graph_to_grid_graph,
)
from util.layout import layout
from util.perf import timing
from util.route import route, determine_router


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
                weight=penalties_cw[p],
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
            weight=EdgePenalty.TO_CENTER,
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
                penalty = dist_euclidean(port_nb, port_self) - 1
                # print(length_penalty)

                G_.add_edge(
                    port_self,
                    port_nb,
                    EdgeType.PHYSICAL,
                    edge=EdgeType.PHYSICAL,
                    center_edge=(node, neighbor),
                    weight=EdgePenalty.HOP + penalty,
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
    with open(os.path.join(directory, f"{name}.json")) as f:
        data = json.load(f)
    elements = data["E"]
    sets = data["S"]
    inst = {
        "elements": elements,
        "elements_inv": invert_list(elements),
        "sets": sets,
        "set_system": list_of_lists_to_set_system_dict(elements, data["SR"]),
        "D_EA": data["EA"],
        #"D_SR": data["SA"],
        "set_colors": dict(
            zip(
                sets,
                data["SC"] if "SC" in data else config_vars["draw.setcolors"].get(),
            )
        ),
        "set_ftb_order": sets,
    }
    if "glyph_ids" in data:
        inst["glyph_ids"] = data["glyph_ids"]
    selection = config_vars["general.sets"].get()
    return select_sets(inst, inst['sets'] if selection == "all" else selection.split(","))


def render():
    dataset = config_vars["general.dataset"].get()
    instance = read_instance(config_vars["general.readdir"].get(), dataset)
    lattice_type = "sqr"
    num_weights = config_vars["general.numweights"].get()
    instance["num_layers"] = num_weights

    with timing("layout"):
        instance["glyph_positions"] = layout(
            instance["elements"],
            instance["D_EA"],
            instance["D_SR"],
            m=config_vars["general.gridwidth"].get(),
            n=config_vars["general.gridheight"].get(),
            pad=config_vars["general.gridpad"].get(),
        )

    with timing("routing"):
        G = get_routing_graph(
            lattice_type,
            get_grid(),
        )
        G = add_glyphs_to_nodes(instance, G)

        element_set_partition = group_by_set(instance["set_system"])
        element_set_partition = sorted(
            element_set_partition, key=lambda x: len(x[1]), reverse=True
        )

        L = route(instance, G, element_set_partition)

    for layer in range(num_weights):
        M = extract_support_layer(L, layer)
        draw_support(instance, M, layer)

    router = determine_router()
    with timing("bundle lines"):
        L = bundle_lines_gg(instance, L)

    # if router == "opt":
    #    # the heuristic routing uses the grid graph (ie., with port nodes), but the optimal routing does not
    #    # this is fine for the bundle step, which also uses the line graph (i.e., without port nodes)
    #    # but for drawing we need the port nodes, so a postprocessing step is required that inserts them.
    #    # we have
    #    # G = multigraph with center and physical nodes/edges
    #    # L = a multigraph with (layer, support) edges and center nodes
    #    G = convert_line_graph_to_grid_graph(instance, L, G, element_set_partition)
    #    L = G

    L = cosmetic_post_processing(instance, L)

    if config_vars['general.serializegraph'].get():
        with timing("serialize graph"):
            # avoid referencing issues
            L_ = copy.deepcopy(L)
            L_.remove_edges_from(
                [
                    (u, v, k)
                    for u, v, k in L_.edges(keys=True)
                    if not (isinstance(k, tuple) and k[1] == EdgeType.SUPPORT)
                ]
            )
            for u, v, d in L_.edges(data=True):
                # convert stuff that json doesn't like
                # sets
                if "sets" in d:
                    d["sets"] = list(d["sets"])
                if "partitions" in d:
                    d["partitions"] = list(d["partitions"])
                # tuple keys
                if "oeb_order" in d:
                    keys = list(d["oeb_order"])
                    for key in keys:
                        d["oeb_order"][str(key)] = [*d["oeb_order"][key]]
                        del d["oeb_order"][key]

            with open(
                os.path.join(config_vars["general.writedir"].get(), "serialized.json"), "w"
            ) as f:
                json.dump(nx.node_link_data(L_), f)
    R = copy.deepcopy(L)
    with timing("draw+write svg"):
        L.add_edges_from(
            [
                (u, v, k, d)
                for u, v, k, d in G.edges(keys=True, data=True)
                if k == EdgeType.PHYSICAL
            ]
        )
        for layer in range(num_weights):
            geometries = geometrize(instance, L, element_set_partition, layer=layer)
            w, h = get_grid()
            img = draw_svg(
                geometries,
                w * config_vars["draw.cellsizepx"].get(),
                h * config_vars["draw.cellsizepx"].get(),
            )
            with open(
                os.path.join(
                    config_vars["general.writedir"].get(), f"drawing_{layer}.svg"
                ),
                "w",
            ) as f:
                f.write(img)
                f.flush()
    R.remove_edges_from(
        [
            (u, v, k)
            for u, v, k in R.edges(keys=True)
            if not (isinstance(k, tuple) and k[1] == EdgeType.SUPPORT)
        ]
    )
    return R


def autogridsize(nrow):
    side = 1
    margin = 1
    base = math.ceil(math.log2(math.sqrt(nrow)))
    # if the layout is too dense then it may get difficult to connect everything properly
    # the failing instances were often around a factor elements/cells = 0.15
    # there may be better metrics to use but this is what i can do now
    while nrow / (side**2) > 0.1:
        exp = base + margin
        side = 2**exp
        margin += 1
    return (side, side)


@click.command()
@click.option("--dataset", help="dataset to load")
@click.option("--read-dir", help="directory to read the datasets from")
@click.option("--write-dir", help="directory to write the output to")
@click.option("--sets", type=str, help="selection of sets")
@click.option(
    "--strategy",
    type=click.Choice(["opt", "heuristic"], case_sensitive=False),
)
@click.option("--grid-width", "-w", type=int, help="grid width as # cols")
@click.option("--grid-height", "-h", type=int, help="grid height as # rows")
@click.option(
    "--num-weights",
    type=int,
    help="how many samples between (weight)..1 to use for weights",
)
@click.option("--weight", type=float, help="start offset for num-weights")
@click.option(
    "--support-type",
    type=click.Choice(["path", "steiner-tree"], case_sensitive=False),
    help="the support type",
)
@click.option(
    "--support-partition",
    type=click.Choice(["set", "intersection-group"], case_sensitive=False),
    help="the partition type",
)
@click.option(
    "--overlap-remover",
    type=click.Choice(["dgrid", "hagrid"], case_sensitive=False),
    help="the overlap removal algorithm",
)
@click.option(
    "--layouter",
    type=click.Choice(["auto", "mds", "qsap", 'umap'], case_sensitive=False),
    help="the layout algorithm",
)
@click.option(
    "--router",
    type=click.Choice(["auto", "opt", "heuristic"], case_sensitive=False),
    help="optimal or heuristic routing",
)
@click.option(
    "--connecter",
    type=click.Choice(["auto", "opt", "heuristic"], case_sensitive=False),
    help="optimal or heuristic connectivity",
)
@click.option(
    "--connect-objective",
    type=click.Choice(["joint", "separate"], case_sensitive=False),
    help="whether connections between nodes should be jointly optimized (minimizing total edge count) or separately.",
)
@click.option(
    "--compute-metrics", type=bool, help="whether or not to compute metrics about graph"
)
@click.option(
    "--serialize-graph", type=bool, help="whether or not to write graph to a file readable by networkx"
)
def vis(
    read_dir,
    write_dir,
    sets,
    dataset,
    strategy,
    num_weights,
    weight,
    grid_width,
    grid_height,
    support_type,
    support_partition,
    overlap_remover,
    layouter,
    router,
    connecter,
    connect_objective,
    compute_metrics,
    serialize_graph,
):
    if dataset is not None:
        config_vars["general.dataset"].set(dataset)
    if support_partition is not None:
        config_vars["general.subsupportgrouping"].set(support_partition)
    if support_type is not None:
        config_vars["general.subsupporttype"].set(support_type)
    if strategy is not None:
        config_vars["general.strategy"].set(strategy)
    if read_dir is not None:
        config_vars["general.readdir"].set(read_dir)
    if write_dir is not None:
        config_vars["general.writedir"].set(write_dir)
    if num_weights is not None:
        config_vars["general.numweights"].set(num_weights)
    if weight is not None:
        config_vars["general.weight"].set(weight)
    if overlap_remover is not None:
        config_vars["layout.overlapremover"].set(overlap_remover)
    if layouter is not None:
        config_vars["layout.layouter"].set(layouter)
    if router is not None:
        config_vars["route.router"].set(router)
    if connecter is not None:
        config_vars["connect.connecter"].set(connecter)
    if grid_width is not None:
        config_vars["general.gridwidth"].set(grid_width)
    if grid_height is not None:
        config_vars["general.gridheight"].set(grid_height)
    if connect_objective is not None:
        config_vars["connect.objective"].set(connect_objective)
    if sets is not None:
        config_vars["general.sets"].set(sets)
    if compute_metrics is not None:
        config_vars["general.computemetrics"].set(compute_metrics)
    if serialize_graph is not None:
        config_vars['general.serializegraph'].set(serialize_graph)

    grid_width, grid_height = get_grid(include_pad=False)
    print(f"Grid size is {grid_width} x {grid_height}")
    print(f"Set selection: {sets}")

    inst = read_instance(
        config_vars["general.readdir"].get(), config_vars["general.dataset"].get()
    )
    if grid_width == 0 or grid_height == 0:
        nrow = len(inst["elements"])
        grid_width, grid_height = autogridsize(nrow)
        config_vars["general.gridwidth"].set(grid_width)
        config_vars["general.gridheight"].set(grid_height)
        print(f"Automatically changed grid size to {grid_width} x {grid_height}")

    os.makedirs(config_vars["general.writedir"].get(), exist_ok=True)

    start = time.time()
    try:
        ctx = contextvars.copy_context()
        G = ctx.run(render)
        end = time.time()

        duration = end - start

        metrics = {}
        with open(
            os.path.join(config_vars["general.writedir"].get(), "call.json"), "w"
        ) as f:
            if config_vars["general.computemetrics"].get():
                metrics = calc_metrics(G)
            json.dump(
                {
                    "success": True,
                    "duration_ms": duration * 1000,
                    "metrics": metrics,
                    "ctx": {key: value.get() for key, value in config_vars.items()},
                },
                f,
            )
            print("SUCCESS")

    except BaseException as e:
        end = time.time()
        duration = end - start
        with open(
            os.path.join(config_vars["general.writedir"].get(), "call.json"), "w"
        ) as f:
            json.dump(
                {
                    "success": False,
                    "duration_ms": duration * 1000,
                    "traceback": traceback.format_exc(),
                    "ctx": {key: value.get() for key, value in config_vars.items()},
                },
                f,
            )
        print("ERROR")
        raise e


if __name__ == "__main__":
    vis()
