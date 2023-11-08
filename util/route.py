import math
import sys
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import product, combinations
from collections import defaultdict
from util.enums import EdgeType, EdgePenalty, NodeType, PortDirs
from util.geometry import get_angle
from util.perf import timing
from util.graph import (
    update_weights_for_support_edge,
    get_shortest_path_between_sets,
    are_node_sets_connected,
    block_edges_using,
    unblock_edges,
    updated_edge_weights,
    updated_port_node_edge_weights_incident_at,
    approximate_steiner_tree_nx,
    approximate_steiner_tree,
    approximate_tsp_tour,
    path_to_edges,
)
from util.collections import set_contains, flatten
from util.config import config_vars
from util.mip import write_status, write_fake_status
import scipy.spatial.distance as sp_dist
import networkx.algorithms.approximation.traveling_salesman as tsp


def get_bend(e1, e2):
    u, v = e1
    v, x = e2

    base = get_angle(u, v) * 180 / math.pi
    turn = get_angle(v, x) * 180 / math.pi

    bend = abs(180 - abs(turn - base))

    match int(bend):
        case 180:
            return 0
        case 0:
            return 0
        case 135:
            return 1
        case 90:
            return 2
        case 45:
            return 3

    raise BaseException(f"invalid bend: {bend}")


def route_single_layer_heuristic(instance, G, element_set_partition, layer=0):
    support_type = config_vars["route.subsupporttype"].get()
    # TODO 100ms spent here
    G_ = nx.Graph()
    G_.add_nodes_from(
        [
            (n, d | d["layers"][layer] if d["node"] == NodeType.CENTER else d)
            for n, d in G.nodes(data=True)
        ]
    )

    # TODO could do this on outside
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

        # find elements to connect
        S = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER and d["occupied"] and d["label"] in elements
        ]

        # lava nodes - can't touch those
        S_minus = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER
            and d["occupied"]
            and d["label"] not in elements
        ]

        current_support_for_sets = nx.subgraph_view(
            G,
            filter_node=lambda n: len(
                list(
                    [
                        (u, v, k)
                        for u, v, k in G.edges(nbunch=n, keys=True)
                        if k == EdgeType.SUPPORT
                        and set_contains(G.edges[u, v, k]["sets"], set(sets))
                    ]
                )
            )
            > 0,
            filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
            and set_contains(G.edges[u, v, k]["sets"], set(sets)),
        )

        current_support_for_sets_exists = (
            len(list(current_support_for_sets.edges())) > 0
        )
        current_support_for_sets_is_connected = (
            current_support_for_sets_exists
            and nx.is_connected(current_support_for_sets)
        )
        current_support_for_sets_contains_S = set_contains(
            set(list(current_support_for_sets.nodes())), set(S)
        )

        if (
            current_support_for_sets_exists
            and current_support_for_sets_is_connected
            and current_support_for_sets_contains_S
        ):
            # we already have a working support for this combination of sets
            print("skipping", sets, ": already supported")
            continue

        # there is a partial support, so everything gets trickier
        # because we have to be smart about incorporating the existing support
        # if a path is desired:
        #  find connected components in partial support
        #  for each component, identify deg-1 nodes
        #  make shortest path graph of all deg-1 and unconnected nodes with weight = path length (no paths within a component)
        #  add virtual start node with edge cost 0 connected to all components and unconnected nodes
        #  find a tour starting at start node

        components = [set(c) for c in nx.connected_components(current_support_for_sets)]
        nodes_in_components = components
        if support_type == "path":
            nodes_in_components = []
            for c in components:
                cs = current_support_for_sets.subgraph(c).copy()
                deg1_nodes = [n for n in c if cs.degree[n] == 1]
                if len(deg1_nodes) > 0:
                    nodes_in_components.append(set(deg1_nodes))

        unconnected_nodes = set(S).difference(
            set(list(current_support_for_sets.nodes()))
        )
        nodes_in_components.extend(list(map(lambda x: set([x]), unconnected_nodes)))

        set_support_edges = []
        if len(S) > 1:
            if support_type == "path":
                with updated_port_node_edge_weights_incident_at(G_, S_minus, math.inf):
                    with updated_edge_weights(
                        G_,
                        [(u, v) for u, v in current_support_for_sets.edges()],
                        math.inf,
                    ):
                        deg2_nodes = [
                            c
                            for i, n in enumerate(components)
                            for c in n
                            if c not in nodes_in_components[i]
                            and G_.nodes[c]["node"] == NodeType.CENTER
                        ]
                        with updated_port_node_edge_weights_incident_at(
                            G_, deg2_nodes, math.inf
                        ):
                            set_support = approximate_tsp_tour(
                                G_, nodes_in_components, current_support_for_sets
                            )
                            set_support_edges = path_to_edges(set_support)
            else:
                # steiner treee
                with updated_edge_weights(
                    G_, list(current_support_for_sets.edges()), 0
                ):
                    with updated_port_node_edge_weights_incident_at(
                        G_, S_minus, math.inf
                    ):
                        nodes_in_components = list(
                            map(
                                lambda g: [
                                    n
                                    for n in g
                                    if G_.nodes[n]["node"] == NodeType.CENTER
                                ],
                                nodes_in_components,
                            )
                        )
                        set_support = approximate_steiner_tree(
                            G_, flatten(nodes_in_components)
                        )
                        set_support_edges = set_support.edges()

        # print('in support edges?', ((0, 5) ,(2, 6)) in set_support_edges)
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

        for e in new_edges:
            u, v = e
            update_weights_for_support_edge(G_, e)
    # 5. stop when all elements have been processed

    support_graph = nx.subgraph_view(
        G, filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
    )
    SG = nx.MultiGraph()
    SG.add_nodes_from(support_graph.nodes(data=True))
    for u, v in support_graph.edges():
        SG.add_edge(
            u,
            v,
            (layer, EdgeType.SUPPORT),
            **support_graph.edges[u, v, EdgeType.SUPPORT],
        )

    return SG


def route_single_layer_ilp(instance, G, element_set_partition, layer):
    support_type = config_vars["route.subsupporttype"].get()
    num_layers = config_vars["general.numlayers"].get()
    el_idx_lookup = instance["elements_inv"]

    MAX_OUT_TERMINAL = 4
    MAX_OUT_ROOT = 4
    if support_type == "path":
        MAX_OUT_TERMINAL = 1
        MAX_OUT_ROOT = 2

    # convert to arcs
    M = nx.DiGraph(incoming_graph_data=G)

    model = gp.Model("multilayer-route")
    if config_vars["route.ilptimeoutsecs"].get() > 0:
        model.params.timeLimit = config_vars["route.ilptimeoutsecs"].get()

    model.params.MIPGap = config_vars["route.ilpmipgap"].get()

    arcs = list(M.edges())
    edges = list(G.edges())

    x = model.addVars(
        [(i, e) for i in range(len(element_set_partition)) for e in arcs],
        vtype=GRB.BINARY,
        name="x_uv_i",
    )

    model._x = x

    # bend penalties
    b = model.addVars(
        [
            (i, n, t)
            for i in range(len(element_set_partition))
            for n in G.nodes()
            for t in range(4)
        ],
        vtype=GRB.BINARY,
        name="bends",
    )

    for n in G.nodes():
        out_arcs_at_n = M.edges(nbunch=n)
        in_arcs_at_n = [(v, u) for u, v in out_arcs_at_n]
        for i in range(len(element_set_partition)):
            for e1, e2 in product(in_arcs_at_n, out_arcs_at_n):
                both_on = model.addVar(vtype=GRB.BINARY)
                model.addConstr(both_on == gp.and_(x[(i, e1)], x[(i, e2)]))
                model.addConstr(both_on <= b[(i, n, get_bend(e1, e2))])

    # connectivity
    for i, esp in enumerate(element_set_partition):
        for u, v in edges:
            # connection flows only in one direction
            model.addConstr(x[(i, (u, v))] + x[(i, (v, u))] <= 1)

        terminals, _ = esp

        # rest is steiner nodes

        for n in G.nodes():
            is_occupied = G.nodes[n]["layers"][layer]["occupied"]
            is_terminal = (
                is_occupied and G.nodes[n]["layers"][layer]["label"] in terminals
            )
            is_non_root_terminal = (
                is_terminal and G.nodes[n]["layers"][layer]["label"] != terminals[0]
            )
            is_root_terminal = (
                is_terminal and G.nodes[n]["layers"][layer]["label"] == terminals[0]
            )
            is_lava = (
                is_occupied and G.nodes[n]["layers"][layer]["label"] not in terminals
            )
            is_steiner = not is_occupied

            out_arcs_at_n = M.edges(nbunch=n)
            in_arcs_at_n = [(v, u) for u, v in out_arcs_at_n]

            sum_in = gp.quicksum([x[(i, (u, v))] for u, v in in_arcs_at_n])
            sum_out = gp.quicksum([x[(i, (u, v))] for u, v in out_arcs_at_n])

            model.addConstr(sum_in <= 1)

            if is_terminal:
                if is_non_root_terminal:
                    model.addConstr(sum_in == 1)
                    model.addConstr(sum_out <= MAX_OUT_TERMINAL)
                if is_root_terminal:
                    model.addConstr(sum_in == 0)
                    model.addConstr(sum_out >= 1)
                    model.addConstr(sum_out <= MAX_OUT_ROOT)

            if is_lava:
                model.addConstr(sum_in == 0)
                model.addConstr(sum_out == 0)

            if is_steiner:
                model.addConstr(sum_in == sum_out)

    obj = config_vars["route.bendlayerfactor"].get() * gp.quicksum(
        [
            b[(i, n, 1)] + b[(i, n, 2)] * 2 + b[(i, n, 3)] * 3
            for i in range(len(element_set_partition))
            for n in G.nodes()
        ]
    ) + config_vars["route.edgelayerlengthfactor"].get() * gp.quicksum(x)

    def addDynamicConstraints(m, xVals):
        # DCC: for each layer and partition: check if it's connected, else add constraint to connect them
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp
            root = instance["glyph_positions"][layer][el_idx_lookup[elements[0]]]

            terminals = set(
                [
                    n
                    for n in G.nodes()
                    if G.nodes[n]["layers"][layer]["occupied"]
                    and G.nodes[n]["layers"][layer]["label"] in elements
                    and n != root
                ]
            )

            # build intermediate graph
            G_ki = nx.DiGraph()
            for u, v in arcs:
                if xVals[(i, (u, v))] > 0:
                    G_ki.add_edge(u, v)
            z = set(nx.depth_first_search.dfs_preorder_nodes(G_ki, root))
            not_z = set([n for n in G.nodes()]).difference(z)
            intersect = terminals.intersection(z)
            if len(intersect) < len(terminals):
                # there are unconnected terminals
                # add a constraint that at least one arc must go from z to !z
                m.cbLazy(
                    gp.quicksum(
                        [m._x[(i, (u, v))] for u, v in arcs if u in z and v in not_z]
                    )
                    >= 1
                )

    # Callback - add lazy constraints to eliminate sub-tours for integer and fractional solutions
    def callback(model, where):
        # check integer solutions for feasibility
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            xValues = model.cbGetSolution(model._x)
            addDynamicConstraints(model, xValues)
        # check fractional solutions to find violated CECs/DCCs to strengthen the bound
        elif (
            where == GRB.Callback.MIPNODE
            and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            xValues = model.cbGetNodeRel(model._x)
            addDynamicConstraints(model, xValues)

    model.update()
    model.setObjective(obj, sense=GRB.MINIMIZE)
    model.Params.LazyConstraints = 1
    model.optimize(callback)

    write_status(f"route{layer}", model)

    MM = nx.MultiGraph(incoming_graph_data=G)

    for i, esp in enumerate(element_set_partition):
        elements, sets = esp

        for e in arcs:
            u, v = e
            if x[(i, e)].x > 0:
                if (u, v, (layer, EdgeType.SUPPORT)) not in MM.edges:
                    MM.add_edge(
                        u,
                        v,
                        (layer, EdgeType.SUPPORT),
                        layer=layer,
                        sets=sets,
                        partitions=[i],
                    )
                else:
                    merged_sets = set(
                        MM.edges[u, v, (layer, EdgeType.SUPPORT)]["sets"]
                    ).union(set(sets))
                    merged_partitions = set(
                        MM.edges[u, v, (layer, EdgeType.SUPPORT)]["partitions"]
                    ).union(set([i]))
                    MM.add_edge(
                        u,
                        v,
                        (layer, EdgeType.SUPPORT),
                        layer=layer,
                        edge=EdgeType.SUPPORT,
                        sets=merged_sets,
                        partitions=merged_partitions,
                    )

    return MM


def determine_router():
    router = config_vars["route.router"].get()
    if router == "auto":
        router = (
            "opt" if config_vars["general.strategy"].get() == "opt" else "heuristic"
        )
    return router


def get_spanning_tree_approx(D, nodeset):
    G = nx.Graph()
    for i, j in combinations(nodeset, r=2):
        G.add_edge(i, j, weight=D[i, j])
    S = nx.minimum_spanning_tree(G)
    return S.edges()


def get_tour_approx(D, nodeset):
    G = nx.Graph()
    for i, j in combinations(nodeset, r=2):
        G.add_edge(i, j, weight=D[i, j])
    S = tsp.traveling_salesman_problem(G, cycle=False)
    return path_to_edges(S)


def route_brosi_ilp(instance, C, G, esp, layer):
    G_ = nx.Graph()
    G_.add_nodes_from(
        [
            (n, d | d["layers"][layer] if d["node"] == NodeType.CENTER else d)
            for n, d in G.nodes(data=True)
        ]
    )

    # TODO could do this on outside
    G_.add_edges_from(
        [
            (u, v, d)
            for u, v, k, d in G.edges(keys=True, data=True)
            if k == EdgeType.PHYSICAL
        ]
    )

    M = nx.DiGraph(incoming_graph_data=G_)
    grid_arcs = list(M.edges())
    grid_edges = list(G_.edges())
    input_edges = list(C.edges())

    model = gp.Model("route-brosi")
    if config_vars["route.ilptimeoutsecs"].get() > 0:
        model.params.timeLimit = config_vars["route.ilptimeoutsecs"].get()

    model.params.MIPGap = config_vars["route.ilpmipgap"].get()

    # decision variable x_ew: is grid arc w used for routing connectivity graph edge e?
    x_ew = model.addVars(
        [(e, w) for e in input_edges for w in grid_arcs], vtype=GRB.BINARY, name="x_ew"
    )

    # use only one arc per grid edge
    for u, v in grid_edges:
        for e in input_edges:
            model.addConstr(x_ew[e, (u, v)] + x_ew[e, (v, u)] <= 1)

    # here the idea is to avoid an input edge carrying the same sets re-using an arc
    # while allowing it for input edges with differing sets
    unprocessed_ie = input_edges.copy()
    while len(unprocessed_ie):
        edge = unprocessed_ie.pop(0)
        edgelist = [edge]
        for i, other in list(enumerate(unprocessed_ie)):
            if C.edges[other]["sets"] == C.edges[edge]["sets"]:
                edgelist.append(other)
                unprocessed_ie.remove(other)
        print(C.edges[edge]["sets"], edgelist)
        for u, v in grid_edges:
            model.addConstr(
                gp.quicksum([x_ew[e, (u, v)] + x_ew[e, (v, u)] for e in edgelist]) <= 1
            )

    # s-t shortest path formulation
    for e in input_edges:
        s, t = e
        pos_s = instance["glyph_positions"][layer][s]
        pos_t = instance["glyph_positions"][layer][t]

        for p in M.nodes():
            out_arcs_at_p = M.edges(nbunch=p)
            in_arcs_at_p = [(v, u) for u, v in out_arcs_at_p]
            sum_out = gp.quicksum([x_ew[e, w] for w in out_arcs_at_p])
            sum_in = gp.quicksum([x_ew[e, w] for w in in_arcs_at_p])

            if M.nodes[p]["node"] == NodeType.PORT:
                parent = M.nodes[p]["belongs_to"]
                is_parent_occupied = M.nodes[parent]["occupied"]

                model.addConstr(sum_out - sum_in == 0)

                if is_parent_occupied and parent not in [pos_s, pos_t]:
                    model.addConstr(sum_in == 0)
                else:
                    model.addConstr(sum_in <= 1)

            else:
                if p == pos_s:
                    model.addConstr(sum_out == 1)
                    model.addConstr(sum_in == 0)
                elif p == pos_t:
                    model.addConstr(sum_in == 1)
                    model.addConstr(sum_out == 0)
                else:
                    model.addConstr(sum_out == 0)
                    model.addConstr(sum_in == 0)

    # TODO? used in paper as well:
    # 10) prevent grid nodes be used for more than one edge
    # 14) pick only one of the crossing twin edges
    # 15-18) ensure edge order at input/grid nodes
    # 19-21) make line bends nice at input nodes

    obj = gp.quicksum(x_ew) * 100 + gp.quicksum(
        [x_ew[(e, w)] * M.edges[w]["weight"] for e in input_edges for w in grid_arcs]
    )
    model.update()
    model.setObjective(obj, sense=GRB.MINIMIZE)
    model.optimize()

    write_status(f"route_{layer}", model)

    MM = nx.MultiGraph()
    for w in grid_arcs:
        u, v = w

        for e in input_edges:
            if x_ew[(e, w)].x > 0:
                # print('layer', layer, 'at grid edge', w, 'sets', C.edges[e]["sets"])
                existing_sets = (
                    MM.edges[(u, v)]["sets"] if (u, v) in MM.edges else set()
                )

                MM.add_edge(
                    u,
                    v,
                    (layer, EdgeType.SUPPORT),
                    edge=EdgeType.SUPPORT,
                    layer=layer,
                    sets=existing_sets.union(set(C.edges[e]["sets"])),
                )

    return MM


def route(instance, G, element_set_partition):
    router = determine_router()

    grid_graph = G
    line_graph = nx.subgraph_view(
        G,
        filter_edge=lambda u, v, k: k == EdgeType.CENTER,
        filter_node=lambda n: G.nodes[n]["node"] == NodeType.CENTER,
    )

    num_layers = config_vars["general.numlayers"].get()
    support_type = config_vars["route.subsupporttype"].get()

    layer_solutions = []
    for layer in range(num_layers):
        pos = instance["glyph_positions"][layer]
        D = sp_dist.squareform(sp_dist.pdist(pos, metric="euclidean"))

        C = nx.Graph()  # connectivity graph
        for elements, sets in element_set_partition:
            nodeset = list(map(lambda e: instance["elements_inv"][e], elements))

            # TODO exact TSP/MST
            connectivity_edges = (
                get_tour_approx(D, nodeset)
                if support_type == "path"
                else get_spanning_tree_approx(D, nodeset)
            )

            for u, v in connectivity_edges:
                existing_sets_at_edge = (
                    C.edges[(u, v)]["sets"] if (u, v) in C.edges else set()
                )
                existing_sets_at_edge = existing_sets_at_edge.union(set(sets))
                C.add_edge(u, v, sets=existing_sets_at_edge)
        print(list(C.edges(data=True)))

        L = (
            route_brosi_ilp(
                instance,
                C.copy(),
                grid_graph.copy(),
                element_set_partition,
                layer=layer,
            )
            if router == "opt"
            else route_single_layer_heuristic(
                instance, grid_graph.copy(), element_set_partition, layer=layer
            )
        )
        layer_solutions.append(L)

    MG = nx.MultiGraph()
    if router == "opt":
        MG.add_nodes_from(grid_graph.nodes(data=True))
    else:
        MG.add_nodes_from(grid_graph.nodes(data=True))

    for lsolution in layer_solutions:
        MG.add_edges_from(lsolution.edges(data=True, keys=True))

    return MG
