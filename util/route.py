import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from util.enums import EdgeType
from itertools import product


def route_multilayer(instance, G, element_set_partition, support_type="tree"):
    num_layers = instance.get("num_layers", 2)
    el_idx_lookup = instance["elements_inv"]

    # convert to arcs
    M = nx.DiGraph(incoming_graph_data=G)

    model = gp.Model("multilayer-route")
    # model.params.timeLimit = 10
    model.params.MIPGap = 0.05

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
    )

    x = model.addVars(edges, vtype=GRB.BINARY, name="x_uv")

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
                    model.addConstr(f[(k, i, j, e)] <= x[e])

    model.setObjective(gp.quicksum(x))

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
                print((k, i, j), instance["glyph_positions"][k][el_idx_lookup[el]])
                for e in edges:
                    u, v = e
                    if j > 0 and f[(k, i, j, e)].x > 0:
                        print((k, i, j, e))
                        if (u, v, EdgeType.SUPPORT) not in MM.edges:
                            MM.add_edge(u, v, EdgeType.SUPPORT, layer=k, sets=sets)
                        else:
                            existing_sets = set(
                                MM.edges[u, v, EdgeType.SUPPORT]["sets"]
                            )
                            merged_sets = existing_sets.union(set(sets))
                            MM.add_edge(
                                u, v, (k, EdgeType.SUPPORT), layer=k, edge=EdgeType.SUPPORT, sets=merged_sets
                            )

    return MM
