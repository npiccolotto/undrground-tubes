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
    approximate_steiner_tree_nx,
    approximate_tsp_tour,
    path_to_edges,
)
from util.collections import set_contains

# factor for edges on all layers
EDGE_SOFT_CONSTRAINT_WEIGHT = 1
# factor for edges of one line
EDGE_LAYER_SOFT_CONSTRAINT_WEIGHT = 1
# factor for bends
BEND_SOFT_CONSTRAINT_WEIGHT = 10


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


def route_single_layer_heuristic(
    instance, G, element_set_partition, support_type="steiner-tree", layer=0
):
    # steps:
    processed_elements = set()

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

        # check if the elements already are connected
        # could have happened e.g. from a previous, larger intersection group
        C = nx.subgraph_view(
            G,
            filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
            and set_contains(G.edges[u, v, k]["sets"], set(sets))
        )
        print(set(sets), C)
        if nx.is_connected(C):
            print(elements, "already connected")
            continue

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
                    G,
                    filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
                    and s in G.edges[u, v, k]["sets"],
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


def route_multilayer_heuristic(
    instance,
    G,
    element_set_partition,
    support_type="steiner-tree",
    multilayer_strategy=("k-of-n", 1),  # 'k-of-n' or 'prev-k'
):
    num_layers = instance.get("num_layers", 2)
    G_ = nx.MultiGraph()
    G_.add_nodes_from(list(G.nodes(data=True)))

    edge_used_in_layers = defaultdict(list)

    for layer in range(num_layers):
        L = route_single_layer_heuristic(
            instance,
            G.copy(),
            element_set_partition,
            support_type=support_type,
            layer=layer,
        )
        for u, v, k in L.edges(keys=True):
            if k != EdgeType.SUPPORT:
                continue
            edge_used_in_layers[(u, v)].append(k)
            G_.add_edge(u, v, (layer, k), **L.edges[u, v, k])

    # down-weight edges used in many layers
    # either if used in all of the previous k layers (at current layer)
    # or if used in k out of n tota layers
    # the choice should match the multilayer layout strategy
    # TODO down-weight by constant or by fraction?

    G_ = nx.MultiGraph()
    G_.add_nodes_from(list(G.nodes(data=True)))

    for layer in range(num_layers):
        L2 = G.copy()
        for edgelayers in edge_used_in_layers.items():
            edge, layers = edgelayers
            strat, k = multilayer_strategy
            should_downweight = False
            match strat:
                case "prev-k":
                    possible_layers = set(range(layer))
                    sought_layers = set(range(max(0, layer - k), layer))
                    should_downweight = (
                        len(sought_layers.intersection(possible_layers))
                        == len(sought_layers)
                        if layer > 0
                        else False
                    )
                case "k-of-n":
                    should_downweight = len(layers) >= k

            if should_downweight:
                u, v = edge
                L2.edges[(u, v, EdgeType.PHYSICAL)]["weight"] = (
                    L.edges[(u, v, EdgeType.PHYSICAL)]["weight"]
                    + EdgePenalty.COMMON_MULTILAYER
                )
        L = route_single_layer_heuristic(
            instance,
            G.copy(),
            element_set_partition,
            support_type=support_type,
            layer=layer,
        )
        for u, v, k in L.edges(keys=True):
            if k != EdgeType.SUPPORT:
                continue
            G_.add_edge(u, v, (layer, k), **L.edges[u, v, k])

    return G_


def route_multilayer_ilp(
    instance, G, element_set_partition, support_type="steiner-tree"
):
    num_layers = instance.get("num_layers", 2)
    el_idx_lookup = instance["elements_inv"]

    match support_type:
        case "steiner-tree":
            MAX_OUT_TERMINAL = 4
            MAX_OUT_ROOT = 4
        case "path":
            MAX_OUT_TERMINAL = 1
            MAX_OUT_ROOT = 2

    # convert to arcs
    M = nx.DiGraph(incoming_graph_data=G)

    model = gp.Model("multilayer-route")
    # model.params.timeLimit = 10
    model.params.MIPGap = 0

    arcs = list(M.edges())
    edges = list(G.edges())

    x = model.addVars(
        [
            (k, i, e)
            for k in range(num_layers)
            for i in range(len(element_set_partition))
            for e in arcs
        ],
        vtype=GRB.BINARY,
        name="x_k_uv_i",
    )

    model._x = x

    x_all = model.addVars([e for e in edges], vtype=GRB.BINARY, name="x_uv")
    for u, v in edges:
        for k in range(num_layers):
            for i in range(len(element_set_partition)):
                model.addConstr(x[(k, i, (u, v))] <= x_all[(u, v)])
                model.addConstr(x[(k, i, (v, u))] <= x_all[(u, v)])

    model._xa = x_all

    # bend penalties
    b = model.addVars(
        [
            (k, i, n, t)
            for k in range(num_layers)
            for i in range(len(element_set_partition))
            for n in G.nodes()
            for t in range(4)
        ],
        vtype=GRB.BINARY,
        name="bends",
    )

    for k in range(num_layers):
        for i in range(len(element_set_partition)):
            for n in G.nodes():
                model.addConstr(gp.quicksum([b[(k, i, n, t)] for t in range(4)]) <= 1)

    for n in G.nodes():
        out_arcs_at_n = M.edges(nbunch=n)
        in_arcs_at_n = [(v, u) for u, v in out_arcs_at_n]
        for k in range(num_layers):
            for i in range(len(element_set_partition)):
                for e1, e2 in product(in_arcs_at_n, out_arcs_at_n):
                    both_on = model.addVar(vtype=GRB.BINARY)
                    model.addConstr(both_on == gp.and_(x[(k, i, e1)], x[(k, i, e2)]))
                    model.addConstr(both_on <= b[(k, i, n, get_bend(e1, e2))])

    # connectivity
    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            for u, v in edges:
                # connection flows only in one direction
                model.addConstr(x[(k, i, (u, v))] + x[(k, i, (v, u))] <= 1)

            terminals, _ = esp

            # rest is steiner nodes

            for n in G.nodes():
                is_occupied = G.nodes[n]["layers"][k]["occupied"]
                is_terminal = (
                    is_occupied and G.nodes[n]["layers"][k]["label"] in terminals
                )
                is_non_root_terminal = (
                    is_terminal and G.nodes[n]["layers"][k]["label"] != terminals[0]
                )
                is_root_terminal = (
                    is_terminal and G.nodes[n]["layers"][k]["label"] == terminals[0]
                )
                is_lava = (
                    is_occupied and G.nodes[n]["layers"][k]["label"] not in terminals
                )
                is_steiner = not is_occupied

                # print(
                #    "layer",
                #    k,
                #    "partition",
                #    i,
                #    "node",
                #    n,
                #    "| = ",
                #    "[occupied]" if is_occupied else "",
                #    "[steiner]" if is_steiner else "",
                #    "[lava]" if is_lava else "",
                #    "[terminal]" if is_terminal else "",
                #    "[root]" if is_root_terminal else "",
                # )

                out_arcs_at_n = M.edges(nbunch=n)
                in_arcs_at_n = [(v, u) for u, v in out_arcs_at_n]

                sum_in = gp.quicksum([x[(k, i, (u, v))] for u, v in in_arcs_at_n])
                sum_out = gp.quicksum([x[(k, i, (u, v))] for u, v in out_arcs_at_n])

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

    obj = (
        BEND_SOFT_CONSTRAINT_WEIGHT
        * gp.quicksum(
            [
                b[(k, i, n, 1)] + b[(k, i, n, 2)] * 2 + b[(k, i, n, 3)] * 3
                for k in range(num_layers)
                for i in range(len(element_set_partition))
                for n in G.nodes()
            ]
        )
        + EDGE_SOFT_CONSTRAINT_WEIGHT * gp.quicksum(x_all)
        + EDGE_LAYER_SOFT_CONSTRAINT_WEIGHT * gp.quicksum(x)
    )

    crossings = set()
    for u, v, k in G.edges(keys=True):
        if k == EdgeType.CENTER and "crossing" in G.edges[u, v, k]:
            crossing_edge = G.edges[u, v, k]["crossing"]
            if crossing_edge is not None:
                crossing = ((u, v), crossing_edge)
                dual_crossing = (crossing_edge, (u, v))
                if not dual_crossing in crossings:
                    crossings.add(crossing)

    # alternative soft constraint
    # cross = model.addVars(crossings, vtype=GRB.BINARY)
    # model.addConstr(cross == gp.and_(x[e1],x[e2]))
    # obj += gp.quicksum(cross) * EdgePenalty.CROSS

    def addDynamicConstraints(m, xVals, xaVals):
        # CROSS: for each diagonal edge: check if both itself and its crossing twin are used, if so disallow it
        for e1, e2 in crossings:
            e1 = tuple(reversed(e1)) if e1 not in xaVals else e1
            e2 = tuple(reversed(e2)) if e2 not in xaVals else e2
            if xaVals[e1] > 0 and xaVals[e2] > 0:
                m.cbLazy(m._xa[e1] + m._xa[e2] <= 1)

        # DCC: for each layer and partition: check if it's connected, else add constraint to connect them
        for k in range(num_layers):
            for i, esp in enumerate(element_set_partition):
                elements, sets = esp
                root = instance["glyph_positions"][k][el_idx_lookup[elements[0]]]

                terminals = set(
                    [
                        n
                        for n in G.nodes()
                        if G.nodes[n]["layers"][k]["occupied"]
                        and G.nodes[n]["layers"][k]["label"] in elements
                        and n != root
                    ]
                )

                # build intermediate graph
                G_ki = nx.DiGraph()
                for u, v in arcs:
                    if xVals[(k, i, (u, v))] > 0:
                        G_ki.add_edge(u, v)
                z = set(nx.depth_first_search.dfs_preorder_nodes(G_ki, root))
                not_z = set([n for n in G.nodes()]).difference(z)
                intersect = terminals.intersection(z)
                if len(intersect) < len(terminals):
                    # there are unconnected terminals
                    # add a constraint that at least one arc must go from z to !z
                    m.cbLazy(
                        gp.quicksum(
                            [
                                m._x[(k, i, (u, v))]
                                for u, v in arcs
                                if u in z and v in not_z
                            ]
                        )
                        >= 1
                    )

    # Callback - add lazy constraints to eliminate sub-tours for integer and fractional solutions
    def callback(model, where):
        # check integer solutions for feasibility
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            xValues = model.cbGetSolution(model._x)
            xaValues = model.cbGetSolution(model._xa)
            addDynamicConstraints(model, xValues, xaValues)
        # check fractional solutions to find violated CECs/DCCs to strengthen the bound
        elif (
            where == GRB.Callback.MIPNODE
            and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            xValues = model.cbGetNodeRel(model._x)
            xaValues = model.cbGetNodeRel(model._xa)
            addDynamicConstraints(model, xValues, xaValues)

    model.update()
    model.setObjective(obj, sense=GRB.MINIMIZE)
    model.Params.LazyConstraints = 1
    model.optimize(callback)

    MM = nx.MultiGraph(incoming_graph_data=G)

    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp

            for e in arcs:
                u, v = e
                if x[(k, i, e)].x > 0:
                    if (u, v, (k, EdgeType.SUPPORT)) not in MM.edges:
                        MM.add_edge(
                            u,
                            v,
                            (k, EdgeType.SUPPORT),
                            layer=k,
                            sets=sets,
                            partitions=[i],
                        )
                    else:
                        merged_sets = set(
                            MM.edges[u, v, (k, EdgeType.SUPPORT)]["sets"]
                        ).union(set(sets))
                        merged_partitions = set(
                            MM.edges[u, v, (k, EdgeType.SUPPORT)]["partitions"]
                        ).union(set([i]))
                        MM.add_edge(
                            u,
                            v,
                            (k, EdgeType.SUPPORT),
                            layer=k,
                            edge=EdgeType.SUPPORT,
                            sets=merged_sets,
                            partitions=merged_partitions,
                        )

    return MM


def route_multilayer_ilp_gg(
    instance, G, element_set_partition, support_type="steiner-tree"
):
    """Same function as the other but on a grid graph, i.e., the bend penalties are
    not modeled as variables but as edge weights. Seems worse than doing it on the
    line graph."""
    num_layers = instance.get("num_layers", 2)
    el_idx_lookup = instance["elements_inv"]

    match support_type:
        case "steiner-tree":
            MAX_OUT_TERMINAL = 4
            MAX_OUT_ROOT = 4
        case "path":
            MAX_OUT_TERMINAL = 1
            MAX_OUT_ROOT = 2

    # convert to arcs
    M = nx.DiGraph(incoming_graph_data=G)

    model = gp.Model("multilayer-route-grid-graph")
    # model.params.timeLimit = 10
    model.params.MIPGap = 0

    arcs = list(M.edges())
    arc_weights = list([M.edges[a]["weight"] for a in arcs])
    edges = list(G.edges())
    edge_weights = list(
        [G.edges[(u, v, EdgeType.PHYSICAL)]["weight"] for u, v in edges]
    )

    x = model.addVars(
        [
            (k, i, e)
            for k in range(num_layers)
            for i in range(len(element_set_partition))
            for e in arcs
        ],
        vtype=GRB.BINARY,
        name="x_k_uv_i",
    )

    model._x = x

    x_all = model.addVars([e for e in edges], vtype=GRB.BINARY, name="x_uv")
    for u, v in edges:
        for k in range(num_layers):
            for i in range(len(element_set_partition)):
                model.addConstr(x[(k, i, (u, v))] <= x_all[(u, v)])
                model.addConstr(x[(k, i, (v, u))] <= x_all[(u, v)])

    model._xa = x_all

    # connectivity
    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            for u, v in edges:
                # connection flows only in one direction
                model.addConstr(x[(k, i, (u, v))] + x[(k, i, (v, u))] <= 1)

            terminals, _ = esp

            for n in G.nodes():
                is_center = G.nodes[n]["node"] == NodeType.CENTER
                is_port = G.nodes[n]["node"] == NodeType.PORT

                is_occupied = is_center and G.nodes[n]["layers"][k]["occupied"]
                is_terminal = (
                    is_occupied and G.nodes[n]["layers"][k]["label"] in terminals
                )
                is_non_root_terminal = (
                    is_terminal and G.nodes[n]["layers"][k]["label"] != terminals[0]
                )
                is_root_terminal = (
                    is_terminal and G.nodes[n]["layers"][k]["label"] == terminals[0]
                )
                is_lava = (
                    is_occupied and G.nodes[n]["layers"][k]["label"] not in terminals
                )
                is_steiner = not is_occupied

                # print(
                #    "layer",
                #    k,
                #    "partition",
                #    i,
                #    "node",
                #    n,
                #    "| = ",
                #    "[occupied]" if is_occupied else "",
                #    "[steiner]" if is_steiner else "",
                #    "[lava]" if is_lava else "",
                #    "[terminal]" if is_terminal else "",
                #    "[root]" if is_root_terminal else "",
                # )

                out_arcs_at_n = M.edges(nbunch=n)
                in_arcs_at_n = [(v, u) for u, v in out_arcs_at_n]

                sum_in = gp.quicksum([x[(k, i, (u, v))] for u, v in in_arcs_at_n])
                sum_out = gp.quicksum([x[(k, i, (u, v))] for u, v in out_arcs_at_n])

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

    obj = EDGE_SOFT_CONSTRAINT_WEIGHT * gp.quicksum(
        [x_all[e] * edge_weights[i] for i, e in enumerate(edges)]
    ) + EDGE_LAYER_SOFT_CONSTRAINT_WEIGHT * gp.quicksum(
        [
            x[(k, i, a)] * arc_weights[j]
            for k in range(num_layers)
            for i in range(len(element_set_partition))
            for j, a in enumerate(arcs)
        ]
    )

    # crossings = set()
    # for u, v, k in G.edges(keys=True):
    #    if k == EdgeType.CENTER and "crossing" in G.edges[u, v, k]:
    #        crossing_edge = G.edges[u, v, k]["crossing"]
    #        if crossing_edge is not None:
    #            crossing = ((u, v), crossing_edge)
    #            dual_crossing = (crossing_edge, (u, v))
    #            if not dual_crossing in crossings:
    #                crossings.add(crossing)

    # alternative soft constraint
    # cross = model.addVars(crossings, vtype=GRB.BINARY)
    # model.addConstr(cross == gp.and_(x[e1],x[e2]))
    # obj += gp.quicksum(cross) * EdgePenalty.CROSS

    def addDynamicConstraints(m, xVals, xaVals):
        # CROSS: for each diagonal edge: check if both itself and its crossing twin are used, if so disallow it
        # for e1, e2 in crossings:
        #    e1 = tuple(reversed(e1)) if e1 not in xaVals else e1
        #    e2 = tuple(reversed(e2)) if e2 not in xaVals else e2
        #    if xaVals[e1] > 0 and xaVals[e2] > 0:
        #        m.cbLazy(m._xa[e1] + m._xa[e2] <= 1)

        # DCC: for each layer and partition: check if it's connected, else add constraint to connect them
        for k in range(num_layers):
            for i, esp in enumerate(element_set_partition):
                elements, sets = esp
                root = instance["glyph_positions"][k][el_idx_lookup[elements[0]]]

                terminals = set(
                    [
                        n
                        for n in G.nodes()
                        if G.nodes[n]["node"] == NodeType.CENTER
                        and G.nodes[n]["layers"][k]["occupied"]
                        and G.nodes[n]["layers"][k]["label"] in elements
                        and n != root
                    ]
                )

                # build intermediate graph
                G_ki = nx.DiGraph()
                for u, v in arcs:
                    if xVals[(k, i, (u, v))] > 0:
                        G_ki.add_edge(u, v)
                z = set(nx.depth_first_search.dfs_preorder_nodes(G_ki, root))
                not_z = set([n for n in G.nodes()]).difference(z)
                intersect = terminals.intersection(z)
                if len(intersect) < len(terminals):
                    # there are unconnected terminals
                    # add a constraint that at least one arc must go from z to !z
                    m.cbLazy(
                        gp.quicksum(
                            [
                                m._x[(k, i, (u, v))]
                                for u, v in arcs
                                if u in z and v in not_z
                            ]
                        )
                        >= 1
                    )

    # Callback - add lazy constraints to eliminate sub-tours for integer and fractional solutions
    def callback(model, where):
        # check integer solutions for feasibility
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            xValues = model.cbGetSolution(model._x)
            xaValues = model.cbGetSolution(model._xa)
            addDynamicConstraints(model, xValues, xaValues)
        # check fractional solutions to find violated CECs/DCCs to strengthen the bound
        elif (
            where == GRB.Callback.MIPNODE
            and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            xValues = model.cbGetNodeRel(model._x)
            xaValues = model.cbGetNodeRel(model._xa)
            addDynamicConstraints(model, xValues, xaValues)

    model.update()
    model.setObjective(obj, sense=GRB.MINIMIZE)
    model.Params.LazyConstraints = 1
    model.optimize(callback)

    MM = nx.MultiGraph(incoming_graph_data=G)

    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp

            for arc in arcs:
                u, v = arc
                if x[(k, i, arc)].x > 0:
                    if (u, v, (k, EdgeType.SUPPORT)) not in MM.edges:
                        MM.add_edge(
                            u,
                            v,
                            (k, EdgeType.SUPPORT),
                            layer=k,
                            sets=sets,
                            partitions=[i],
                        )
                    else:
                        merged_sets = set(
                            MM.edges[u, v, (k, EdgeType.SUPPORT)]["sets"]
                        ).union(set(sets))
                        merged_partitions = set(
                            MM.edges[u, v, (k, EdgeType.SUPPORT)]["partitions"]
                        ).union(set([i]))
                        MM.add_edge(
                            u,
                            v,
                            (k, EdgeType.SUPPORT),
                            layer=k,
                            edge=EdgeType.SUPPORT,
                            sets=merged_sets,
                            partitions=merged_partitions,
                        )

    return MM
