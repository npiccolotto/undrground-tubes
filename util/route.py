import math
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import product, combinations
from util.enums import EdgeType, EdgePenalty, NodeType
from util.geometry import get_angle

USE_MAX_DEG_CONSTR = False
PREVENT_FORKS_AT_GLYPHS = False
EDGE_SOFT_CONSTRAINT_WEIGHT = 2
BEND_SOFT_CONSTRAINT_WEIGHT = 1
CROSS_SOFT_CONSTRAINT_WEIGHT = 1


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


def route_multilayer(instance, G, element_set_partition, support_type="tree"):
    num_layers = instance.get("num_layers", 2)
    el_idx_lookup = instance["elements_inv"]

    # convert to arcs
    M = nx.DiGraph(incoming_graph_data=G)

    model = gp.Model("multilayer-route")
    # model.params.timeLimit = 10
    model.params.MIPGap = 0.2

    edges = list(M.edges())

    f = model.addVars(
        [
            (k, i, j, e)
            for k in range(num_layers)
            for i in range(len(element_set_partition))
            for j in range(1, len(element_set_partition[i][0]))
            for e in edges
        ],
        vtype=GRB.BINARY,
        name="flow",
    )

    x = model.addVars(edges, vtype=GRB.BINARY, name="x_uv")

    x_k = model.addVars(
        [(k, e) for k in range(num_layers) for e in edges], vtype=GRB.BINARY
    )

    for u, v in G.edges():
        for k in range(num_layers):
            model.addConstr(x[(u, v)] + x[(v, u)] <= 1)

    obj = gp.quicksum(x) * EdgePenalty.HOP * EDGE_SOFT_CONSTRAINT_WEIGHT

    # use only one of the two arcs on each layer
    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp
            for j in range(1, len(elements)):
                for u, v in G.edges():
                    model.addConstr(f[k, i, j, (u, v)] + f[k, i, j, (v, u)] <= 1)

    # linking constraints flow to x_uv
    for e in edges:
        for k in range(num_layers):
            for i, esp in enumerate(element_set_partition):
                elements, sets = esp
                for j in range(1, len(elements)):
                    model.addConstr(f[(k, i, j, e)] <= x_k[(k, e)])
                    model.addConstr(f[(k, i, j, e)] <= x[e])

    if BEND_SOFT_CONSTRAINT_WEIGHT > 0:
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
                    model.addConstr(
                        gp.quicksum([b[(k, i, n, t)] for t in range(4)]) <= 1
                    )

        x_ki = model.addVars(
            [
                (k, i, e)
                for k in range(num_layers)
                for i in range(len(element_set_partition))
                for e in edges
            ],
            vtype=GRB.BINARY,
            name="x_ki",
        )
        # linking constraints flow to x_ki
        for e in edges:
            for k in range(num_layers):
                for i, esp in enumerate(element_set_partition):
                    elements, sets = esp
                    for j in range(1, len(elements)):
                        model.addConstr(f[(k, i, j, e)] <= x_ki[(k, i, e)])

        for n in G.nodes():
            out_edges_at_n = M.edges(nbunch=n)
            in_edges_at_n = [(v, u) for u, v in out_edges_at_n]
            for k in range(num_layers):
                for i in range(len(element_set_partition)):
                    for e1, e2 in product(in_edges_at_n, out_edges_at_n):
                        both_on = model.addVar(vtype=GRB.BINARY)
                        model.addConstr(
                            both_on == gp.and_(x_ki[(k, i, e1)], x_ki[(k, i, e2)])
                        )
                        model.addConstr(both_on <= b[(k, i, n, get_bend(e1, e2))])

        obj += gp.quicksum(
            [
                b[(k, i, n, 0)] * EdgePenalty.ONE_EIGHTY
                + b[(k, i, n, 1)] * EdgePenalty.ONE_THIRTY_FIVE
                + b[(k, i, n, 2)] * EdgePenalty.NINETY
                + b[(k, i, n, 3)] * EdgePenalty.FORTY_FIVE
                for k in range(num_layers)
                for i in range(len(element_set_partition))
                for n in G.nodes()
            ]
        )

    if CROSS_SOFT_CONSTRAINT_WEIGHT > 0:
        # collect all crossing edges between nodes
        crossings = set(
            [
                set(list((u, v), G.edges[u, v, EdgeType.CENTER]["crossing"]))
                for u, v in G.edges()
                if k == EdgeType.CENTER
                and G.nodes[u]["node"] == NodeType.CENTER
                and G.nodes[v]["node"] == NodeType.CENTER
                and G.edges[u, v]["crossing"] is not None
            ]
        )

        # alternative soft constraint
        # cross = model.addVars(crossings, vtype=GRB.BINARY)
        # model.addConstr(cross == gp.and_(x[e1],x[e2]))
        # obj += gp.quicksum(cross) * EdgePenalty.CROSS

        # avoid them
        for e1, e2 in crossings:
            for k in range(num_layers):
                model.addConstr(x_k[(k, e1)] + x_k[(k, e2)] <= 1)

        # TODO avoid also crossings within nodes

    model.setObjective(obj)

    # max-degree constraints for unoccupied nodes (prevent forks at those)
    if USE_MAX_DEG_CONSTR:
        # TODO alternatively we could compute the degree and put it into objective with penalty
        for n in G.nodes():
            out_edges_at_n = nx.edges(M, nbunch=n)
            in_edges_at_n = [(v, u) for u, v in out_edges_at_n]

            is_occupied = G.nodes[n]["layers"][k]["occupied"]

            if (PREVENT_FORKS_AT_GLYPHS and is_occupied) or (
                not PREVENT_FORKS_AT_GLYPHS and not is_occupied
            ):
                for k in range(num_layers):
                    for i in range(len(element_set_partition)):
                        num_js = len(element_set_partition[i][0])
                        for j1, j2 in combinations(range(1, num_js), 2):
                            for e1, e2 in product(in_edges_at_n, repeat=2):
                                if e1 == e2:
                                    continue
                                model.addConstr(
                                    f[(k, i, j1, e1)] + f[(k, i, j2, e2)] <= 1
                                )
                            for e1, e2 in product(out_edges_at_n, repeat=2):
                                if e1 == e2:
                                    continue
                                model.addConstr(
                                    f[(k, i, j1, e1)] + f[(k, i, j2, e2)] <= 1
                                )

    # MCF constraints

    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp
            root = instance["glyph_positions"][k][el_idx_lookup[elements[0]]]
            edges_at_root = list(M.edges(nbunch=root))

            for j, el in enumerate(elements):
                if j == 0:
                    continue

                # root sends out commodities
                model.addConstr(
                    gp.quicksum([f[(k, i, j, e)] for e in edges_at_root]) == 1
                )

                # target node must consume the commodity
                target = instance["glyph_positions"][k][el_idx_lookup[el]]
                edges_at_target = M.edges(nbunch=target)
                model.addConstr(
                    gp.quicksum([f[(k, i, j, (v, u))] for u, v in edges_at_target])
                    - gp.quicksum([f[(k, i, j, (u, v))] for u, v in edges_at_target])
                    == 1
                )

    # non-target nodes must pass commodity on
    # commodity cannot flow over occupied nodes
    for k in range(num_layers):
        for n, d in G.nodes(data=True):
            layer_info = d["layers"][k]

            out_edges_at_n = [(u, v) for u, v in M.edges(nbunch=n)]
            in_edges_at_n = [(v, u) for u, v in out_edges_at_n]

            if layer_info["occupied"]:
                element = layer_info["label"]

                node_subsets_with_el = [
                    i
                    for i, esp in enumerate(element_set_partition)
                    if element in esp[0]
                ]
                node_subsets_without_el = set(
                    range(len(element_set_partition))
                ).difference(set(node_subsets_with_el))

                # if the node is occupied, then
                # - commodities of the same node subset (except itself) on that layer must flow through: inflow = outflow
                # - commodities of other node subsets on that layer must not flow through: inflow = 0
                for i in node_subsets_with_el:
                    if element_set_partition[i][0].index(element) == 0:
                        # this is the root node for subset i
                        # do nothing
                        continue
                    for l, p in enumerate(element_set_partition[i][0]):
                        if p == element or l == 0:
                            continue
                        model.addConstr(
                            gp.quicksum([f[(k, i, l, e)] for e in in_edges_at_n])
                            - gp.quicksum([f[(k, i, l, e)] for e in out_edges_at_n])
                            == 0
                        )
                for i in node_subsets_without_el:
                    for l in range(1, len(element_set_partition[i][0])):
                        model.addConstr(
                            gp.quicksum([f[(k, i, l, e)] for e in in_edges_at_n]) == 0
                        )
                        model.addConstr(
                            gp.quicksum([f[(k, i, l, e)] for e in out_edges_at_n]) == 0
                        )
            else:
                # if the node is not occupied, then all commodities must flow through
                for i, esp in enumerate(element_set_partition):
                    elements, sets = esp
                    for j in range(1, len(elements)):
                        model.addConstr(
                            gp.quicksum([f[(k, i, j, e)] for e in in_edges_at_n])
                            - gp.quicksum([f[(k, i, j, e)] for e in out_edges_at_n])
                            == 0
                        )

    # TODO bend penalties
    # TODO crossing penalties

    model.update()
    model.optimize()

    # create support graph from model output
    # ignore root edges
    # make it a multigraph, i guess
    MM = nx.MultiGraph(incoming_graph_data=G)

    for k in range(num_layers):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp

            for j, el in enumerate(elements):
                # print((k, i, j), instance["glyph_positions"][k][el_idx_lookup[el]])
                for e in edges:
                    u, v = e
                    if j > 0 and f[(k, i, j, e)].x > 0:
                        # ((k, i, j, e))
                        if (u, v, (k, EdgeType.SUPPORT)) not in MM.edges:
                            MM.add_edge(
                                u,
                                v,
                                (k, EdgeType.SUPPORT),
                                layer=k,
                                sets=sets,
                                partitions=[(i, j)],
                            )
                        else:
                            merged_sets = set(
                                MM.edges[u, v, (k, EdgeType.SUPPORT)]["sets"]
                            ).union(set(sets))
                            merged_partitions = set(
                                MM.edges[u, v, (k, EdgeType.SUPPORT)]["partitions"]
                            ).union(set([(i, j)]))
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
