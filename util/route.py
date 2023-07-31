import math
import sys
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import product, combinations
from util.enums import EdgeType, EdgePenalty, NodeType
from util.geometry import get_angle

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


def route_multilayer_ilp(
    instance, G, element_set_partition, support_type="steiner-tree"
):
    num_layers = instance.get("num_layers", 2)
    el_idx_lookup = instance["elements_inv"]

    # theoretically, if we set these to 1 then the individual supports should be paths
    # this is how we adapt to the desired support type
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
