import math
import sys
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import product, combinations, pairwise
from collections import defaultdict
from util.enums import EdgeType, EdgePenalty, NodeType, PortDirs
from util.geometry import get_angle, do_lines_intersect, get_segment_circle_intersection
from util.perf import timing
from util.graph import (
    update_weights_for_support_edge,
    get_shortest_path_between_sets,
    are_node_sets_connected,
    block_edges_using,
    unblock_edges,
    get_port_edge_between_centers,
    updated_edge_weights,
    updated_port_node_edge_weights_incident_at,
    approximate_steiner_tree_nx,
    approximate_steiner_tree,
    approximate_tsp_tour,
    path_to_edges,
    sort_ports,
    get_port_edges,
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


def line_deg(G, n):
    nbs = nx.neighbors(G, n)
    ldeg = 0
    for nb in nbs:
        ldeg += len(G.edges[(n, nb)]["sets"])
    return ldeg


def heuristic_edge_order(C):
    unprocessed_nodes = list(map(lambda n: (n, line_deg(C, n)), C.nodes()))
    unprocessed_nodes = sorted(unprocessed_nodes, key=lambda un: un[1], reverse=True)

    edge_order = []

    while len(unprocessed_nodes):
        dangling = [unprocessed_nodes[0]]

        while len(dangling):
            dangling = sorted(dangling, key=lambda un: un[1], reverse=True)
            n, n_deg = dangling.pop(0)

            nbs = list(map(lambda nb: (nb, line_deg(C, nb)), nx.neighbors(C, n)))
            unprocessed = list(map(lambda un: un[0], unprocessed_nodes))
            nbs = filter(lambda nb: nb[0] in unprocessed, nbs)
            nbs = sorted(nbs, key=lambda nb: nb[1], reverse=True)

            for nb, nb_deg in nbs:
                edge_order.append((n, nb))
                if (nb, nb_deg) not in dangling:
                    dangling.append((nb, nb_deg))

            unprocessed_nodes.remove((n, n_deg))

    return edge_order


def route_single_layer_heuristic2(instance, C, G, element_set_partition, layer=0):
    G_ = nx.Graph()
    G_.add_nodes_from(
        [
            (n, d | d["layers"][layer] if d["node"] == NodeType.CENTER else d)
            for n, d in G.nodes(data=True)
        ]
    )

    G_.add_edges_from(
        [
            (u, v, d)
            for u, v, k, d in G.edges(keys=True, data=True)
            if k == EdgeType.PHYSICAL
        ]
    )

    for u, v in [
        (u, v)
        for u, v in sorted(
            C.edges(), reverse=True, key=lambda uv: len(C.edges[uv]["sets"])
        )
    ]:
        sets_at_uv = C.edges[u, v]["sets"]
        elements = [instance["elements"][u], instance["elements"][v]]

        # lava nodes - can't touch those along the route
        S_minus = [
            n
            for n, d in G_.nodes(data=True)
            if d["node"] == NodeType.CENTER
            and d["occupied"]
            and d["label"] not in elements
        ]

        edges_carrying_uv_sets = [
            (w, x)
            for w, x, k in G.edges(keys=True)
            if k == EdgeType.SUPPORT
            and len(G.edges[w, x, k]["sets"].intersection(sets_at_uv)) > 0
        ]

        crossings_of_edges_carrying_uv_sets = []
        nodes_at_edges_carrying_uv_sets = []
        for w, x in edges_carrying_uv_sets:
            w_is_port = G_.nodes[w]["node"] == NodeType.PORT
            x_is_port = G_.nodes[x]["node"] == NodeType.PORT

            if w_is_port and x_is_port:
                are_from_different_centers = (
                    G_.nodes[w]["belongs_to"] != G_.nodes[x]["belongs_to"]
                )
                if are_from_different_centers:
                    # center edge belonging to physical edge
                    cw, cx = G.edges[(w, x, EdgeType.PHYSICAL)]["center_edge"]
                    if "crossing" in G.edges[(cw, cx, EdgeType.CENTER)]:
                        # crossing center edge
                        crw, crx = G.edges[(cw, cx, EdgeType.CENTER)]["crossing"]
                        # phsical edge belonging to crossing center
                        crossings_of_edges_carrying_uv_sets.append(
                            get_port_edge_between_centers(G_, crw, crx)
                        )
                else:
                    nodes_at_edges_carrying_uv_sets.append(G_.nodes[w]["belongs_to"])

        # avoid lava nodes
        with updated_port_node_edge_weights_incident_at(G_, S_minus, math.inf):
            # avoid forks/merges (=reusing edges)
            with updated_edge_weights(G_, edges_carrying_uv_sets, math.inf):
                # avoid self-crossings on diagonals
                with updated_edge_weights(
                    G_, crossings_of_edges_carrying_uv_sets, math.inf
                ):
                    # avoid self-crossings at grid cells
                    with updated_port_node_edge_weights_incident_at(
                        G_, nodes_at_edges_carrying_uv_sets, math.inf
                    ):
                        new_edges = []
                        pos_u = instance["glyph_positions"][layer][u]
                        pos_v = instance["glyph_positions"][layer][v]

                        set_support_edges = path_to_edges(
                            nx.shortest_path(G_, pos_u, pos_v, weight="weight")
                        )

                        for w, x in set_support_edges:
                            if not G.has_edge(w, x, EdgeType.SUPPORT):
                                new_edges.append((w, x))
                                G.add_edge(
                                    w,
                                    x,
                                    EdgeType.SUPPORT,
                                    edge=EdgeType.SUPPORT,
                                    sets=set(sets_at_uv),
                                )
                            else:
                                G.edges[(w, x, EdgeType.SUPPORT)]["sets"] = set(
                                    sets_at_uv
                                ).union(G.edges[(w, x, EdgeType.SUPPORT)]["sets"])

                        for e in new_edges:
                            update_weights_for_support_edge(G_, e)

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


def route_single_layer_heuristic(instance, G, element_set_partition, layer=0):
    support_type = config_vars["general.subsupporttype"].get()
    G_ = nx.Graph()
    G_.add_nodes_from(
        [
            (n, d | d["layers"][layer] if d["node"] == NodeType.CENTER else d)
            for n, d in G.nodes(data=True)
        ]
    )

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
    support_type = config_vars["general.subsupporttype"].get()
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

    write_status(f"route_{layer}", model)

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


def determine_connecter():
    connecter = config_vars["connect.connecter"].get()
    if connecter == "auto":
        connecter = (
            "opt" if config_vars["general.strategy"].get() == "opt" else "heuristic"
        )
    return connecter


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
    S = tsp.traveling_salesman_problem(G, cycle=False, method=tsp.greedy_tsp)
    return path_to_edges(S)


def get_dir(p1, p2):
    x1, y1 = p1
    x2, y2 = p2

    dx = x2 - x1
    dy = y2 - y1

    if dx > 0 and dy == 0:
        return 0
    elif dx > 0 and dy < 0:
        return 1
    elif dx == 0 and dy < 0:
        return 2
    elif dx < 0 and dy < 0:
        return 3
    elif dx < 0 and dy == 0:
        return 4
    elif dx < 0 and dy > 0:
        return 5
    elif dx == 0 and dy > 0:
        return 6
    elif dx > 0 and dy > 0:
        return 7


def dir_to_penalty(dir):
    if dir not in range(8):
        raise BaseException(f"invalid dir {dir}")
    if dir in (1, 7):
        return 3
    elif dir in (2, 6):
        return 2
    elif dir in (3, 5):
        return 1
    return 0


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
    n_nodes = len(list(C.nodes()))

    model = gp.Model("route-brosi")
    if config_vars["route.ilptimeoutsecs"].get() > 0:
        model.params.timeLimit = config_vars["route.ilptimeoutsecs"].get()

    model.params.MIPGap = config_vars["route.ilpmipgap"].get()

    # decision variable x_ew: is grid arc w used for routing input graph edge e?
    x_ew = model.addVars(
        [(e, w) for e in input_edges for w in grid_arcs], vtype=GRB.BINARY, name="x_ew"
    )
    model._x_ew = x_ew

    # convenience variables that tell whether a grid edge is used by any input edge
    x = model.addVars([(u, v) for u, v in grid_edges], vtype=GRB.BINARY)
    for u, v in grid_edges:
        for e in input_edges:
            model.addConstr(x[(u, v)] >= x_ew[(e, (u, v))])
            model.addConstr(x[(u, v)] >= x_ew[(e, (v, u))])
    model._x = x

    obj = gp.quicksum(
        [x_ew[(e, w)] * M.edges[w]["weight"] for e in input_edges for w in grid_arcs]
    ) + config_vars["route.ilpsnugfactor"].get() * gp.quicksum(x[a] for a in grid_edges)

    # use only one arc per grid edge
    for u, v in grid_edges:
        for e in input_edges:
            model.addConstr(x_ew[e, (u, v)] + x_ew[e, (v, u)] <= 1)

    # here the idea is to avoid an input edge carrying the same sets re-using an arc
    # while allowing it for input edges with differing sets
    # this is to avoid forks/merges
    sets_at_input_edges = set(map(lambda e: frozenset(C.edges[e]["sets"]), input_edges))
    while len(sets_at_input_edges):
        s = sets_at_input_edges.pop()
        edgelist = [
            e for e in input_edges if len(C.edges[e]["sets"].intersection(s)) > 0
        ]
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

    # each input edge may use only one arc at a grid node
    # that SHOULD be obvious to the solver given the edge weights but better safe than sorry
    for n in M.nodes():
        if M.nodes[n]["node"] == NodeType.CENTER:
            ports = nx.neighbors(M, n)
            port_edges = list(combinations(ports, r=2))
            port_arcs = port_edges + list(map(lambda e: tuple(reversed(e)), port_edges))
            center_edges = list(map(lambda p: (n, p), ports))
            center_arcs = center_edges + list(
                map(lambda e: tuple(reversed(e)), center_edges)
            )
            for e in input_edges:
                model.addConstr(
                    gp.quicksum([x_ew[(e, w)] for w in port_arcs + center_arcs]) <= 1
                )

    if True:
        # 15-18) ensure edge order at input/grid nodes <- this one can prevent some crossings
        delta = model.addVars(
            [(v, e) for v in range(n_nodes) for e in input_edges if v in e],
            vtype=GRB.INTEGER,
            name="delta",
            lb=0,
            ub=7,
        )
        for e in input_edges:
            s, t = e
            dir_sum_s = gp.LinExpr(0)
            dir_sum_t = gp.LinExpr(0)

            for u in M.nodes():
                if M.nodes[u]["node"] == NodeType.CENTER:
                    ports = sort_ports(M, nx.neighbors(M, u), origin="e")
                    dir_sum_s += gp.quicksum(
                        [x_ew[(e, (u, p))] * (i) for i, p in enumerate(ports)]
                    )
                    dir_sum_t += gp.quicksum(
                        [x_ew[(e, (p, u))] * (i) for i, p in enumerate(ports)]
                    )

            model.addConstr(delta[(s, e)] - dir_sum_s == 0)
            model.addConstr(delta[(t, e)] - dir_sum_t == 0)

        beta = model.addVars(
            [
                (p, u)
                for u in range(n_nodes)
                for p in range(C.degree[u])
                if C.degree[u] > 2
            ],
            vtype=GRB.BINARY,
            name="beta",
        )

        for u in range(n_nodes):
            if C.degree[u] > 2:
                ps = range(C.degree[u])
                neighbors = nx.neighbors(C, u)
                neighbors_cw = sorted(
                    neighbors,
                    reverse=False,
                    key=lambda v: get_angle(
                        instance["glyph_positions"][layer][u],
                        instance["glyph_positions"][layer][v],
                    ),
                )
                for p, p1 in pairwise(ps):
                    e1 = (
                        (u, neighbors_cw[p1])
                        if (u, neighbors_cw[p1]) in input_edges
                        else (neighbors_cw[p1], u)
                    )
                    e2 = (
                        (u, neighbors_cw[p])
                        if (u, neighbors_cw[p]) in input_edges
                        else (neighbors_cw[p], u)
                    )
                    assert e1 in input_edges
                    assert e2 in input_edges
                    model.addConstr(
                        delta[(u, e1)] - delta[(u, e2)] + 8 * beta[p, u] >= 1
                    )
                model.addConstr(gp.quicksum([beta[p, u] for p in ps]) <= 1)

    # 19-21) make line bends nice at input nodes

    # 10) prevent grid nodes be used for more than one edge
    # 14) pick only one of the crossing twin edges <- not super promising because our graph is not necessarily planar

    def addDynamicConstraints(m, x, x_ew):
        # check if any two selected edges overlap
        # if they overlap, check if they carry the same sets
        # if they do, add constrain flow of those sets to not use both input edges
        for n in M.nodes():
            if M.nodes[n]["node"] == NodeType.CENTER:
                port_edges = get_port_edges(M, n)
                port_arcs = port_edges + list(
                    map(lambda e: tuple(reversed(e)), port_edges)
                )
                for p1, p2 in combinations(port_arcs, r=2):
                    e_p1 = [e for e in input_edges if x_ew[e, p1] > 0]
                    e_p2 = [e for e in input_edges if x_ew[e, p2] > 0]
                    if len(e_p1) > 0 and len(e_p2) > 0:
                        l_p1 = set.union(*[C.edges[e]["sets"] for e in e_p1])
                        l_p2 = set.union(*[C.edges[e]["sets"] for e in e_p2])
                        common_l = l_p1.intersection(l_p2)

                        if len(common_l) > 0:
                            m.cbLazy(
                                gp.quicksum(
                                    [
                                        m._x_ew[e, p]
                                        for e in e_p1 + e_p2
                                        for p in port_arcs
                                    ]
                                )
                                <= 1
                            )

    def callback(model, where):
        # check integer solutions for feasibility
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            xValues = model.cbGetSolution(model._x)
            xwValues = model.cbGetSolution(model._x_ew)
            addDynamicConstraints(model, xValues, xwValues)
        elif (
            where == GRB.Callback.MIPNODE
            and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            xValues = model.cbGetNodeRel(model._x)
            xwValues = model.cbGetNodeRel(model._x_ew)
            addDynamicConstraints(model, xValues, xwValues)

    model.update()
    model.setObjective(obj, sense=GRB.MINIMIZE)
    if True:
        model.Params.LazyConstraints = 1
        model.optimize(callback)
    else:
        model.optimize()

    write_status(f"route_{layer}", model)

    MM = nx.MultiGraph()
    for w in grid_arcs:
        u, v = w

        for e in input_edges:
            if x_ew[(e, w)].x > 0:
                existing_sets = (
                    MM.edges[(u, v, (layer, EdgeType.SUPPORT))]["sets"]
                    if (u, v, (layer, EdgeType.SUPPORT)) in MM.edges
                    else set()
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


def get_optimal_connectivity(instance, D, element_set_partition, layer=0, tour=False):
    # another ILP but let's do a MCF formulation with a spanning tree

    model = gp.Model("connect")
    if config_vars["connect.ilptimeoutsecs"].get() > 0:
        model.params.timeLimit = config_vars["connect.ilptimeoutsecs"].get()

    model.params.MIPGap = config_vars["connect.ilpmipgap"].get()

    n_nodes = len(instance["elements"])
    pos = instance["glyph_positions"][layer]
    labels = []
    for i in range(n_nodes):
        i_labels = set()
        for elements, sets in element_set_partition:
            if instance["elements"][i] in elements:
                i_labels.add(frozenset(sets))
        labels.append(frozenset(i_labels))

    G = nx.Graph()
    edges = list(combinations(range(n_nodes), r=2))
    for i, j in edges:
        G.add_node(i, labels=labels[i])
        G.add_node(j, labels=labels[j])

        G.add_edge(i, j, weight=D[i, j])

    G_d = nx.DiGraph(incoming_graph_data=G)
    G_d.add_node("root")
    for i in range(n_nodes):
        G_d.add_edge("root", i, weight=0)

    arcs = list(G_d.edges())

    comm_nodes = []
    elements, sets = zip(*element_set_partition)
    all_labels = list(map(lambda s: frozenset(s), sets))

    for i, l in enumerate(all_labels):
        for e in elements[i]:
            comm_nodes.append((l, instance["elements_inv"][e]))

    # flow tell if an arc transports a commodity
    flow = model.addVars(
        [(a, ce) for a in arcs for ce in comm_nodes], vtype=GRB.BINARY, name="flow"
    )

    # only in one direction at an edge
    for i, j in edges:
        for ce in comm_nodes:
            model.addConstr(flow[((i, j), ce)] + flow[((j, i), ce)] <= 1)

    # x tell if an edge is used to transport any commodity
    x = model.addVars(edges, vtype=GRB.BINARY, name="x")
    for e in edges:
        for ce in comm_nodes:
            i, j = e
            model.addConstr(x[(i, j)] >= flow[((i, j), ce)])
            model.addConstr(x[(i, j)] >= flow[((j, i), ce)])

    # x_l tell if an arc is used to transport commodities of a given label
    x_l = model.addVars([(a, l) for a in arcs for l in all_labels], vtype=GRB.BINARY)

    model._x = x
    model._xl = x_l
    model._f = flow

    # outgoing arcs of root transport all commodities
    root_arcs = [("root", n) for n in range(n_nodes)]
    for ce in comm_nodes:
        model.addConstr(gp.quicksum([flow[(a, ce)] for a in root_arcs]) == 1)

    for l in all_labels:
        for e in edges:
            model.addConstr(x[e] >= x_l[e, l])

    # but all commodities of the same label must go out on the same arc
    for l in all_labels:
        model.addConstr(gp.quicksum([x_l[(a, l)] for a in root_arcs]) <= 1)

        commodities = [ce for ce in comm_nodes if ce[0] == l]
        # link x_l and flow
        for ce in commodities:
            for a in arcs:
                model.addConstr(x_l[(a, l)] >= flow[(a, ce)])

    # allow transport of commodities only between nodes that also belong to that label
    for i, j in arcs:
        if i == "root":
            continue

        for ce in comm_nodes:
            if ce[0] not in labels[i].intersection(labels[j]):
                model.addConstr(flow[((i, j), ce)] == 0)

    for i in range(n_nodes):
        in_arcs = [a for a in arcs if a[1] == i]
        out_arcs = [a for a in arcs if a[0] == i]

        # a node must give out all commodities it receives that doesn't belong to it
        not_my_stuff = [ce for ce in comm_nodes if ce[1] != i]
        for ce in not_my_stuff:
            model.addConstr(
                gp.quicksum([flow[a, ce] for a in in_arcs])
                - gp.quicksum([flow[a, ce] for a in out_arcs])
                == 0
            )

        # a node must receive and keep all commodities that belong to it
        my_stuff = [ce for ce in comm_nodes if ce[1] == i]
        for ce in my_stuff:
            model.addConstr(gp.quicksum([flow[a, ce] for a in in_arcs]) == 1)
            model.addConstr(gp.quicksum([flow[a, ce] for a in out_arcs]) == 0)

        # limit edges incident at each node
        # because we can't plot more than 8
        edges_at_i = [e for e in edges if i in e]
        model.addConstr(gp.quicksum([x[e] for e in edges_at_i]) <= 8)

        # commodities of the same label must enter at one edge
        for l in all_labels:
            model.addConstr(gp.quicksum([x_l[a, l] for a in in_arcs]) <= 1)

            # when we look for a tour for each label, they must also exit at one edge
            if tour:
                model.addConstr(gp.quicksum([x_l[a, l] for a in out_arcs]) <= 1)

    # minimize the edge length in the drawing
    obj = gp.quicksum([x[e] * G_d.edges[e]["weight"] for e in edges])
    # obj = gp.quicksum([x_l[tuple(reversed(e)),l]*x_l[e,l] * G_d.edges[e]["weight"] for e in edges for l in all_labels])

    def addDynamicConstraints(m, x, xl):
        # check if any two selected edges overlap
        # if they overlap, check if they carry the same sets
        # if they do, add constrain flow of those sets to not use both input edges
        selected_edges = set([e for e in edges if x[e] > 0])
        # with timing(f'{len(selected_edges)} selected edges'):
        for e1, e2 in combinations(selected_edges, r=2):
            u, v = e1
            w, z = e2
            common_labels = [
                l
                for l in all_labels
                if (xl[((u, v), l)] > 0 or xl[((v, u), l)] > 0)
                and (xl[((w, z), l)] > 0 or xl[((z, w), l)] > 0)
            ]

            if len(common_labels) > 0 and do_lines_intersect(
                pos[u], pos[v], pos[w], pos[z]
            ):
                for l in common_labels:
                    # of the four possible arcs, use only one
                    m.cbLazy(
                        m._xl[((u, v), l)]
                        + m._xl[((v, u), l)]
                        + m._xl[((w, z), l)]
                        + m._xl[((z, w), l)]
                        <= 1
                    )

    def callback(model, where):
        # check integer solutions for feasibility
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            xValues = model.cbGetSolution(model._x)
            xlValues = model.cbGetSolution(model._xl)
            addDynamicConstraints(model, xValues, xlValues)
        elif (
            where == GRB.Callback.MIPNODE
            and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            xValues = model.cbGetNodeRel(model._x)
            xlValues = model.cbGetNodeRel(model._xl)
            addDynamicConstraints(model, xValues, xlValues)

    model.update()
    # print(flow)
    model.setObjective(obj, sense=GRB.MINIMIZE)
    if True:
        model.Params.LazyConstraints = 1
        model.optimize(callback)
    else:
        model.optimize()

    write_status(f"connect_{layer}", model)

    MM = nx.Graph()
    for i, j in arcs:
        if i == "root":
            continue

        for ce in comm_nodes:
            if flow[((i, j), ce)].x > 0:
                sets_at_ij = MM.edges[i, j]["sets"] if (i, j) in MM.edges else set()
                l = set([s for s in ce[0]])
                MM.add_edge(i, j, sets=sets_at_ij.union(l))

        # if x[(i, j)].x > 0:
        #    labels = [l for l in all_labels if x_l[((i, j), l)].x > 0]
        #    MM.add_edge(i, j, sets=set(flatten(labels)))

    return MM


def connect(instance, G, element_set_partition, layer=0):
    support_type = config_vars["general.subsupporttype"].get()
    pos = instance["glyph_positions"][layer]
    D = sp_dist.squareform(sp_dist.pdist(pos, metric="euclidean"))
    connecter = determine_connecter()

    C = nx.Graph()  # connectivity graph
    if connecter == "opt":
        C = (
            get_optimal_connectivity(
                instance, D, element_set_partition, layer=layer, tour=True
            )
            if support_type == "path"
            else get_optimal_connectivity(
                instance, D, element_set_partition, layer=layer
            )
        )
    else:
        for elements, sets in element_set_partition:
            nodeset = list(map(lambda e: instance["elements_inv"][e], elements))
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

        write_fake_status(f"connect_{layer}")

    return C


def route(instance, grid_graph, element_set_partition):
    router = determine_router()

    num_layers = config_vars["general.numlayers"].get()

    layer_solutions = []
    for layer in range(num_layers):
        C = connect(instance, grid_graph, element_set_partition, layer=layer)

        L = (
            route_brosi_ilp(
                instance,
                C,
                grid_graph,
                element_set_partition,
                layer=layer,
            )
            if router == "opt"
            else route_single_layer_heuristic2(
                instance, C, grid_graph.copy(), element_set_partition, layer=layer
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
