import numpy as np
import networkx as nx
import math
import random
import networkx.algorithms.approximation.traveling_salesman as tsp
import networkx.algorithms.approximation.steinertree as steinertree
from contextlib import ContextDecorator
from copy import deepcopy
from itertools import product, pairwise, combinations
from functools import partial

from util.geometry import dist_euclidean
from util.perf import timing
from util.collections import set_contains, flatten
from util.enums import NodeType, EdgePenalty, PortDirs, EdgeType
from pathos.multiprocessing import ProcessingPool as Pool


def path_to_edges(path):
    # path is a list of nodes
    # we want that as edges
    return list(pairwise(path))


def incident_edges(G, start, avoid_node=None):
    possible_edges = G.edges(nbunch=start, keys=True)
    edges = [(u, v, k) for u, v, k in possible_edges if v != avoid_node]
    return edges


def get_longest_simple_paths(G, src, tgt, visited=[], current_path=[]):
    # follow the path as long as nodes have deg = 2
    # when paths fork (deg>2), add path up until here to list
    # restart procedure along the forks
    # eventually we have visited all nodes

    next_visited = visited + [src]
    if G.degree[tgt] > 2:
        # fork, split here
        next_tgts = [i for i in G.neighbors(tgt) if i not in next_visited]
        next_paths = [
            get_longest_simple_paths(G, tgt, next_tgt, visited=[], current_path=[])
            for next_tgt in next_tgts
        ]
        current_path = current_path + [src, tgt]
        ret = [current_path]
        for p in next_paths:
            ret += p
        return ret

    if G.degree[tgt] == 1:
        # end of the simple path
        current_path = current_path + [src, tgt]
        return [current_path]

    # regular case
    next_src = tgt
    next_tgt = [i for i in G.neighbors(tgt) if i not in next_visited]
    if len(next_tgt) > 1:
        raise Exception("uhm what")
    return get_longest_simple_paths(
        G,
        next_src,
        next_tgt[0],
        visited=next_visited,
        current_path=current_path + [src],
    )


def calculate_path_length(G, path, weight=None):
    # path is a list of edges
    if len(path) == 0:
        return float("inf")

    total = 0
    for e in path:
        length = 0
        if weight:
            length = G.edges[e][weight]
        else:
            a, b = e
            length = dist_euclidean(G.nodes[a]["pos"], G.nodes[b]["pos"])

        total += length
    return total


def distance_matrix(list1, list2, dist):
    n = len(list1)
    m = len(list2)
    D = np.ndarray(shape=(n, m))
    for i, j in product(range(n), range(m)):
        D[i, j] = dist(list1[i], list2[j])
    return D


def get_closest_pair(list1, list2):
    D = distance_matrix(list1, list2, dist_euclidean)
    m = np.where(D == np.min(D))
    m = (m[0][0], m[1][0])
    return m


def get_shortest_path_between(G, list1, list2, weight=None):
    # TODO rewrite to use virtual start and end nodes
    # so that we only use one shortest path call
    shortest = float("inf"), []
    for i in list1:
        for j in list2:
            p = nx.shortest_path(G, i, j, weight=weight)
            p = path_to_edges(p)
            l = calculate_path_length(G, p, weight=weight)
            if l < shortest[0]:
                shortest = (l, p)
    return shortest


def get_shortest_path_between_sets(G, S, T):
    # idea: add temp nodes s (connected with cost 0 to all S) and t (analog)
    # then run dijkstra and remove s and t from shortest path and graph
    s = "virtual_source"
    t = "virtual_target"
    temp_edges = []

    G.add_node(s)
    for i in S:
        G.add_edge(s, i, weight=0)
        temp_edges.append((s, i))
    G.add_node(t)
    for j in T:
        G.add_edge(t, j, weight=0)
        temp_edges.append((t, j))

    def cleanup():
        for u, v in temp_edges:
            G.remove_edge(u, v)
        G.remove_node(s)
        G.remove_node(t)

    try:
        sp = nx.shortest_path(G, s, t, weight="weight")
    except nx.exception.NetworkXNoPath as e:
        cleanup()
        raise e

    cleanup()
    return [node for node in sp if node != s and node != t]


def shortest_path_graph(G, t=1):
    """Meulemans et al., 2013: Kelpfusion: A Hybrid Set Visualization Technique"""
    S = nx.Graph(incoming_graph_data=G)
    S.remove_edges_from(list(G.edges()))

    # sort by weight ascending
    edges = sorted(list(G.edges(data=True)), key=lambda u, v, d: d["weight"])

    for u, v, d in edges:
        w = d["weight"]
        try:
            p = nx.shortest_path_length(S, u, v)
        except nx.NetworkXNoPath:
            p = -1

        if p < 0 or p >= w**t:
            S.add_edge(u, v, weight=w**t)

    return S


def are_node_sets_connected(G, S, T):
    try:
        get_shortest_path_between_sets(G, S, T)
        return True
    except nx.exception.NetworkXNoPath:
        return False


def count_node_occurrence_in_path(G, S, path):
    # S are center nodes
    # count how often we use nodes belonging to these centers in path
    count = 0
    for node in path:
        if G.nodes[node]["node"] == NodeType.CENTER and node in S:
            count += 1
        if G.nodes[node]["node"] == NodeType.PORT and G.nodes[node]["belongs_to"] in S:
            count += 1
    return count


def approximate_steiner_tree_nx(G, S):
    return steinertree.steiner_tree(G, S, method="mehlhorn")


def approximate_steiner_tree(G, S, C=None):
    """Steiner tree but works with groups of nodes too"""
    # S is a group of nodes
    # C is the current support in G

    unconnected_nodes = list(filter(lambda s: len(s) == 1, S))
    groups = list(filter(lambda s: len(s) > 1, S))

    G_ = G.copy()
    C = nx.Graph() if C is None else C
    R = nx.Graph()

    if len(unconnected_nodes) > 1:
        groups.append(flatten(unconnected_nodes))
        ST = approximate_steiner_tree_nx(G_, flatten(unconnected_nodes))
        # R.add_edges_from(ST.edges())
        R.add_weighted_edges_from(
            [(u, v, G_.edges[u, v]["weight"]) for u, v in ST.edges()]
        )

    R.add_weighted_edges_from([(u, v, G_.edges[u, v]["weight"]) for u, v in C.edges()])

    if len(groups) > 1:
        while not nx.is_connected(R):
            SP = nx.Graph()
            SP.add_nodes_from(map(lambda g: frozenset(g), groups))
            missing_connections = []
            connections = []
            for g1, g2 in combinations(groups, 2):
                try:
                    sp = get_shortest_path_between_sets(R, g1, g2)
                    SP.add_edge(
                        frozenset(g1),
                        frozenset(g2),
                        path=sp,
                        weight=calculate_path_length(
                            R, path_to_edges(sp), weight="weight"
                        ),
                    )
                    connections.append((frozenset(g1), frozenset(g2)))
                except nx.exception.NetworkXNoPath:
                    missing_connections.append((frozenset(g1), frozenset(g2)))

            if not nx.is_connected(SP):
                for g1, g2 in missing_connections:
                    sp = get_shortest_path_between_sets(G_, g1, g2)
                    SP.add_edge(
                        g1,
                        g2,
                        path=sp,
                        weight=calculate_path_length(
                            G_, path_to_edges(sp), weight="weight"
                        ),
                    )
                MST = nx.minimum_spanning_tree(SP)

                # find cheapest connection
                min_segment = (math.inf, None)
                for conn in missing_connections:
                    g1, g2 = conn
                    cost, _ = min_segment
                    if MST.has_edge(g1, g2) and cost > MST.edges[g1, g2]["weight"]:
                        min_segment = (MST.edges[g1, g2]["weight"], conn)

                _, conn = min_segment
                if conn is None:
                    break

                u, v = conn
                sp = MST.edges[u, v]["path"]

                R.add_weighted_edges_from(
                    [(u, v, G_.edges[u, v]["weight"]) for u, v in path_to_edges(sp)]
                )
                update_edge_weights(G_, path_to_edges(sp), math.inf)
    assert nx.is_connected(R)
    return R


def update_edge_weights(G, edges, weight):
    for u, v in edges:
        if (u, v) in G.edges:
            G.edges[u, v]["weight"] = (
                weight[(u, v)]
                if isinstance(weight, dict)
                else (weight(G, u, v) if callable(weight) else weight)
            )


class updated_edge_weights(ContextDecorator):
    def __init__(self, G, edges, weight=0):
        self.G = G
        self.edges = edges
        self.weight = weight
        self.initial_edge_weights = {}
        for u, v in edges:
            self.initial_edge_weights[(u, v)] = G.edges[u, v]["weight"]

    def __enter__(self):
        update_edge_weights(self.G, self.edges, self.weight)

    def __exit__(self, *exc):
        update_edge_weights(self.G, self.edges, self.initial_edge_weights)


class updated_port_node_edge_weights_incident_at(ContextDecorator):
    def __init__(self, G, closed_nodes, weight=math.inf):
        self.G = G
        self.nodes = closed_nodes
        self.weight = weight
        self.edges = [
            (u, v)
            for node in closed_nodes
            for p in G.neighbors(node)
            for u, v in G.edges(p)
            if (
                G.nodes[u]["node"] == NodeType.PORT
                and G.nodes[v]["node"] == NodeType.PORT
            )
        ]
        self.initial_edge_weights = {}
        for u, v in self.edges:
            self.initial_edge_weights[(u, v)] = G.edges[u, v]["weight"]

    def __enter__(self):
        update_edge_weights(self.G, self.edges, self.weight)

    def __exit__(self, *exc):
        update_edge_weights(self.G, self.edges, self.initial_edge_weights)


# TODO could be a context manager
def block_edges_using(G, closed_nodes):
    for node in closed_nodes:
        if G.nodes[node]["node"] != NodeType.CENTER:
            raise BaseException("trying to block non-center node!")
        ports = G.neighbors(node)
        for p in ports:
            for e in G.edges(p):
                u, v = e
                if (
                    G.nodes[u]["node"] == NodeType.PORT
                    and G.nodes[v]["node"] == NodeType.PORT
                    and ("blocked" not in G.edges[e] or not G.edges[e]["blocked"])
                ):
                    G.edges[e]["blocked"] = True
                    G.edges[e]["base_weight"] = G.edges[e]["weight"]
                    G.edges[e]["weight"] = float(math.inf)
    return G


def unblock_edges(G, closed_nodes):
    for node in closed_nodes:
        ports = G.neighbors(node)
        for p in ports:
            for e in G.edges(p):
                if "blocked" in G.edges[e] and G.edges[e]["blocked"]:
                    G.edges[e]["blocked"] = False
                    G.edges[e]["weight"] = G.edges[e]["base_weight"]

    # for u, v in G.edges():
    #    if G.nodes[u]["node"] == NodeType.PORT and G.nodes[v]["node"] == NodeType.PORT:
    #        parents = set([G.nodes[u]["belongs_to"], G.nodes[v]["belongs_to"]])
    #        if len(parents.intersection(set(closed_nodes))) > 0:
    #            G.edges[(u, v)]["weight"] = G.edges[(u, v)]["base_weight"]
    return G


def calc_path_matrix(graph, S=None, heuristic=None, weight="weight", threads=1):
    """Calculates shortest paths between specified nodes in the graph."""

    nodes = S if S is not None else list(map(lambda x: set([x]), graph.nodes()))

    def _find_path():
        def path(node_pair):
            source, target = node_pair
            return (
                nx.astar_path(
                    graph,
                    source,
                    target,
                    heuristic=heuristic,
                    weight=weight,
                )
                if heuristic is not None
                else nx.shortest_path(graph, source, target, weight=weight)
            )

        return path

    find_path = _find_path()
    node_pairs = [
        (source, target)
        for s1, s2 in combinations(nodes, 2)
        for source, target in product(s1, s2)
    ]
    with Pool(nodes=threads) as pool:
        paths = pool.map(find_path, node_pairs)
    return {frozenset(node_pair): path for node_pair, path in zip(node_pairs, paths)}


def calc_distance_matrix(graph, pathmatrix, weight="weight"):
    return {
        node_pair: calculate_path_length(graph, path_to_edges(path), weight)
        for node_pair, path in pathmatrix.items()
    }


def remove_segment_in_circle(circle, s, e):
    pos1 = circle.index(s) + 1
    pos2 = circle[pos1:].index(e) + pos1
    return circle[pos2:] + circle[1:pos1]


def group_aware_greedy_tsp(G, weight="weight", groups=None, source=None):
    """Follows implementation according to [1]. Returns a tour in G.

    Properties (desired):
    - Uses an edge in G at most once -> prevents side-steps in a path that don't look like paths in the drawing.
    - Connects at least one element in each group
    - Connects each element at most once

    [1] https://networkx.org/documentation/stable/_modules/networkx/algorithms/approximation/traveling_salesman.html#greedy_tsp
    """
    if groups is None:
        groups = list(map(lambda n: set([n]), G.nodes()))

    if len(groups) == 2:
        sp = get_shortest_path_between_sets(G, groups[0], groups[1])
        return sp

    if source is None:
        # pick any but deterministically
        source = list(groups[0])[0]

    G_ = G.copy()

    nodeset = set.union(*groups)
    node_conflicts = dict()
    for n in nodeset:
        for g in groups:
            if n in g:
                gs = set(g)
                gs.remove(n)
                node_conflicts[n] = gs
                continue

    nodeset.remove(source)
    cycle = [source]
    cycle_dists = []
    next_node = source
    paths = []


    while nodeset:
        nodelist = [n for n in nodeset if n not in node_conflicts[cycle[-1]]]
        if not nodelist:
            # there are no non-conflicting nodes we could connect to
            # assume that we're done? because it means the only other nodes
            # left are from the same group, which by assumption is connected already
            break
        shortest_paths = [
            nx.shortest_path(G_, cycle[-1], n, weight=weight) for n in nodelist
        ]
        dists = [
            calculate_path_length(G_, path_to_edges(sp), weight=weight)
            for sp in shortest_paths
        ]
        argmin = np.argmin(dists)
        min_idx = argmin[0] if isinstance(argmin, list) else argmin
        next_node = nodelist[min_idx]
        cycle_dists.append((cycle[-1], next_node, dists[min_idx]))
        cycle.append(next_node)
        paths.extend(shortest_paths[min_idx][:-1])
        nodeset.remove(next_node)

        # update G
        # block used edges, must not be used again
        update_edge_weights(G_, path_to_edges(shortest_paths[min_idx]), math.inf)
        # block also the crossing twins?

    sp = nx.shortest_path(G_, cycle[-1], cycle[0], weight=weight)
    dist = calculate_path_length(G_, path_to_edges(sp), weight=weight)
    cycle_dists.append((cycle[-1], cycle[0], dist))
    cycle.append(cycle[0])
    paths.extend(sp)

    # remove biggest step in cycle
    dists = list(map(lambda x: x[2], cycle_dists))
    longest_step = np.argmax(dists)
    a, b, _ = cycle_dists[longest_step]
    paths = remove_segment_in_circle(paths, a, b)

    assert all(
        [(u, v) in G_.edges for u, v in path_to_edges(paths)]
    ), "at least one edge not in graph"
    assert all(
        [
            (G_.nodes[u]["node"], G_.nodes[v]["node"])
            != (NodeType.CENTER, NodeType.CENTER)
            for u, v in path_to_edges(paths)
        ],
    ), "at least one edge connects two centers"
    assert G_.nodes[paths[0]]["node"] == NodeType.CENTER, "first node not a center"
    assert G_.nodes[paths[-1]]["node"] == NodeType.CENTER, "last node not a center"

    return paths


def approximate_tsp_tour(G, S):
    """Returns an approximation of the shortest tour in G visiting all nodes S exactly once while using each edge at most once.
    Which is almost but not really a TSP tour (e.g., G is not completely connected), but let's continue calling it that.
    """
    tour = group_aware_greedy_tsp(G, weight="weight", groups=S)
    return tour


def get_node_with_degree(G, deg):
    nodes_and_deg = [(n, G.degree(n)) for n in G.nodes()]
    deg_nodes = [n for n, d in nodes_and_deg if d == deg]
    return deg_nodes[0] if len(deg_nodes) > 0 else None


def visit_edge_pairs_starting_at_node(G, start):
    # so here the idea is that for drawing we need edge pairs
    # we find edge pairs by starting at a node
    edge_pairs = []
    stack = []
    current_edge = None

    incident_edges = list(G.edges(start))
    if len(incident_edges) == 0:
        return []
    current_edge = incident_edges.pop()
    if len(incident_edges) > 0:
        stack = incident_edges

    while True:
        # not currently looking at an edge
        if current_edge is None:
            # is something on the stack? then continue with that
            if len(stack) > 0:
                current_edge = stack.pop()
            # nope. did we find any edge pairs already? if so we assume we're done
            elif len(edge_pairs) > 0:
                return edge_pairs

        u, v = current_edge
        next_edges = list([(w, x) for w, x in G.edges(v) if (w, x) != (v, u)])

        match len(next_edges):
            case 0:
                # there is no next edge
                # we're done for now, let's hope there's something on the stack in the next iteration
                current_edge = None
                continue
            case 1:
                # there is 1 next edge
                next_edge = next_edges[0]
                edge_pairs.append([current_edge, next_edge])
                current_edge = next_edge
                continue
            case _:
                next_edge = next_edges[0]
                pairs = product([current_edge], next_edges)
                edge_pairs = edge_pairs + list(pairs)
                stack = stack + next_edges[1:]
                current_edge = next_edge
                continue


def get_port_edge_between_centers(G, u, v):
    uports = G.neighbors(u)

    for p in uports:
        for w, x in G.edges(p):
            if G.nodes[x]["node"] == NodeType.PORT and G.nodes[x]["belongs_to"] == v:
                return (w, x)
    raise Exception("no port edge between centers")


def update_weights_for_support_edge(G, edge):
    u, v = edge
    un = G.nodes[u]["node"]
    vn = G.nodes[v]["node"]

    # only care about ports
    if un != NodeType.PORT or vn != NodeType.PORT:
        return G

    # make this edge cheaper so future routes will use it more
    G.edges[edge]["weight"] = max(0, G.edges[edge]["weight"] + EdgePenalty.IN_SUPPORT)

    # penalize other edges crossing this one
    # is it an edge between nodes?
    uparent = G.nodes[u]["belongs_to"]
    vparent = G.nodes[v]["belongs_to"]
    """
    if uparent == vparent and not G.nodes[uparent]["occupied"]:
        # nope, within a node
        # get all port edges and find those crossing this one
        ports = G.neighbors(uparent)
        port_edges = [G.edges(p) for p in ports]
        port_edges = list(set([item for sub_list in port_edges for item in sub_list]))
        port_edges = [
            (u, v)
            for u, v in port_edges
            if G.nodes[u]["node"] == NodeType.PORT
            and G.nodes[v]["node"] == NodeType.PORT
        ]
        crossing_edges = [
            (w, x) for w, x in port_edges if do_lines_intersect(w, x, u, v)
        ]
        for w, x in crossing_edges:
            G.edges[(w, x)]["weight"] = G.edges[(w, x)]["weight"] + EdgePenalty.CROSSING
    """
    if uparent != vparent:
        # yep, between nodes
        uparent, vparent = list(sorted([uparent, vparent], key=lambda x: x[1]))
        # now u is up
        ux, uy = uparent
        vx, vy = vparent
        is_left_tilt = ux < vx and uy < vy
        is_right_tilt = ux > vx and uy < vy

        if is_left_tilt:
            e = get_port_edge_between_centers(G, (ux, uy + 1), (ux + 1, uy))
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING_OUTSIDE
        if is_right_tilt:
            e = get_port_edge_between_centers(G, (ux - 1, uy), (ux, uy + 1))
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING_OUTSIDE
    else:
        # nope, within a node
        all_ports = set(
            [p for p in nx.neighbors(G, uparent) if G.nodes[p]["node"] == NodeType.PORT]
        )
        all_edges_at_ports = set()
        for port in all_ports:
            edges_at_port = [
                (a, b)
                if PortDirs.index(G.nodes[a]["port"])
                < PortDirs.index(G.nodes[b]["port"])
                else (b, a)
                for a, b in G.edges(nbunch=port)
                if G.nodes[a]["node"] == NodeType.PORT
                and G.nodes[b]["node"] == NodeType.PORT
                and G.nodes[a]["belongs_to"] == uparent
                and G.nodes[b]["belongs_to"] == uparent
            ]

            all_edges_at_ports = all_edges_at_ports.union(set(edges_at_port))

        for w, x in all_edges_at_ports:
            if are_port_edges_crossing(G.nodes[u], G.nodes[v], G.nodes[w], G.nodes[x]):
                penalty = (
                    EdgePenalty.CROSSING_INSIDE_GLYPH
                    if G.nodes[uparent]["occupied"]
                    else EdgePenalty.CROSSING_INSIDE_CELL
                )
                G.edges[(w, x)]["weight"] = G.edges[(w, x)]["weight"] + penalty

    return G


def are_port_edges_crossing(us, ut, vs, vt, cross_when_node_shared=True):
    if frozenset([ut["port"], us["port"]]) == frozenset([vs["port"], vt["port"]]):
        return False

    cw_dirs = PortDirs + PortDirs

    # order cw
    u1 = us if PortDirs.index(us["port"]) < PortDirs.index(ut["port"]) else ut
    u2 = ut if PortDirs.index(us["port"]) < PortDirs.index(ut["port"]) else us

    # start at us and go clockwise to ut, collect all ports in between
    cw_set = []
    i = cw_dirs.index(u1["port"])
    j = cw_dirs.index(u2["port"])
    cw_set = cw_dirs[(i + 1) : j]

    # then start at ut and go cw to us, collect all ports in between
    ccw_set = []
    i = cw_dirs.index(u2["port"])
    offset_j = cw_dirs[(i + 1) :].index(u1["port"])
    ccw_set = cw_dirs[(i + 1) : (i + offset_j + 1)]

    assert len(cw_set) > 0 or len(ccw_set) > 0, (
        us["port"],
        ut["port"],
        vs["port"],
        vt["port"],
    )

    # edges cross iff vs in former and vt in latter set (or vice versa)
    port1 = vs["port"]
    port2 = vt["port"]

    if (port1 in cw_set and port2 in cw_set) or (port1 in ccw_set and port2 in ccw_set):
        return False

    if not cross_when_node_shared:
        if u1["port"] in [port1, port2] or u2["port"] in [port1, port2]:
            return False

    return True


def get_crossing_port_edges(G):
    # assuming G consinsts only of port edges of one node
    edges = list(
        [
            (u, v)
            if PortDirs.index(G.nodes[u]["port"]) < PortDirs.index(G.nodes[v]["port"])
            else (v, u)
            for (u, v) in G.edges()
        ]
    )
    edges = list(sorted(edges, key=lambda e: PortDirs.index(G.nodes[e[0]]["port"])))
    crossings = []

    for i, e1 in enumerate(edges):
        u, v = e1
        for j, e2 in enumerate(edges):
            if j <= i:
                continue
            w, x = e2
            if are_port_edges_crossing(G.nodes[u], G.nodes[v], G.nodes[w], G.nodes[x]):
                crossings.append((e1, e2))
    return crossings


def get_port_edges(M, node):
    all_ports = [
        p for p in nx.neighbors(M, node) if M.nodes[p]["node"] == NodeType.PORT
    ]
    all_edges_at_ports = set()
    for port in all_ports:
        edges_at_port = [
            (a, b)
            if PortDirs.index(M.nodes[a]["port"]) < PortDirs.index(M.nodes[b]["port"])
            else (b, a)
            for a, b, k in M.edges(nbunch=port, keys=True)
            if k == EdgeType.PHYSICAL
            and M.nodes[a]["node"] == NodeType.PORT
            and M.nodes[b]["node"] == NodeType.PORT
            and M.nodes[a]["belongs_to"] == node
            and M.nodes[b]["belongs_to"] == node
        ]

        all_edges_at_ports = all_edges_at_ports.union(set(edges_at_port))
    return list(all_edges_at_ports)


def extract_support_layer(M, layer):
    G = nx.Graph()
    G.add_nodes_from(
        [
            (n, d | d["layers"][layer] if d["node"] == NodeType.CENTER else d)
            for n, d in M.nodes(data=True)
        ]
    )
    G.add_edges_from(
        [
            (u, v, d)
            for u, v, k, d in M.edges(keys=True, data=True)
            if k == (layer, EdgeType.SUPPORT)
        ]
    )
    return G


def get_port(G, parent, side):
    all_ports = nx.neighbors(G, parent)
    for p in all_ports:
        if (
            G.nodes[p]["node"] == NodeType.PORT
            and G.nodes[p]["belongs_to"] == parent
            and G.nodes[p]["port"] == side
        ):
            return p
    return None


def get_relative_ports(u, v):
    ux, uy = u
    vx, vy = v

    dx = vx - ux
    dy = vy - uy

    northsouth_u = ""
    match dy:
        case 1:
            northsouth_u = "s"
        case -1:
            northsouth_u = "n"
        case 0:
            northsouth_u = ""

    eastwest_u = ""
    match dx:
        case 1:
            eastwest_u = "e"
        case -1:
            eastwest_u = "w"
        case 0:
            eastwest_u = ""

    v_relative_to_u = f"{northsouth_u}{eastwest_u}"

    northsouth = ""
    match -dy:
        case 1:
            northsouth = "s"
        case -1:
            northsouth = "n"
        case 0:
            northsouth = ""

    eastwest = ""
    match -dx:
        case 1:
            eastwest = "e"
        case -1:
            eastwest = "w"
        case 0:
            eastwest = ""

    u_relative_to_v = f"{northsouth}{eastwest}"

    return (u_relative_to_v, v_relative_to_u)


def get_ports(G, u, v):
    """for two center nodes u and v in G, returns their respective ports that an edge would use"""
    u_relative_to_v, v_relative_to_u = get_relative_ports(u, v)

    return (get_port(G, u, v_relative_to_u), get_port(G, v, u_relative_to_v))


def orient_edge_node_inside(edge, node):
    u, v = edge
    return (u, v) if u == node else (v, u)


def edge_filter_ports(G, u, v, same_centers=False, possibly_with_center=False):
    uparent = (
        None if G.nodes[u]["node"] == NodeType.CENTER else G.nodes[u]["belongs_to"]
    )
    vparent = (
        None if G.nodes[v]["node"] == NodeType.CENTER else G.nodes[v]["belongs_to"]
    )

    match (uparent, vparent):
        case (None, None):
            return False
        case (None, _):
            return possibly_with_center
        case (_, None):
            return possibly_with_center
        case (_, _):
            match same_centers:
                case True:
                    return uparent == vparent
                case False:
                    return uparent != vparent
                case None:
                    return True


def convert_line_graph_to_grid_graph(instance, L, G, element_set_partition):
    num_weights = instance["num_layers"]
    for layer in range(num_weights):
        for i, esp in enumerate(element_set_partition):
            elements, sets = esp
            root_pos = instance["glyph_positions"][layer][
                instance["elements_inv"][elements[0]]
            ]
            for j, el in enumerate(elements):
                if j == 0:
                    continue
                j_pos = instance["glyph_positions"][layer][instance["elements_inv"][el]]
                P = nx.subgraph_view(
                    L,
                    filter_edge=lambda u, v, k: k == (layer, EdgeType.SUPPORT)
                    and i in L.edges[u, v, k]["partitions"],
                )
                path_j = path_to_edges(
                    nx.shortest_path(P, root_pos, j_pos)
                )  # P should actually just be a path already but we do this to order edges

                if len(path_j) == 1:
                    u, v = path_j[0]
                    port_u, port_v = get_ports(G, u, v)

                    # needs connection u_center -> port_u -> port_v -> v_center

                    # u_center -> port_u
                    G.add_edge(
                        u,
                        port_u,
                        (layer, EdgeType.SUPPORT),
                        edge=EdgeType.SUPPORT,
                        sets=set(L.edges[u, v, (layer, EdgeType.SUPPORT)]["sets"]),
                    )

                    # port_u -> port_v
                    oeb_order_pu_pv = {
                        (port_u, port_v): L.edges[u, v, (layer, EdgeType.SUPPORT)][
                            "oeb_order"
                        ][(u, v)],
                        (port_v, port_u): L.edges[u, v, (layer, EdgeType.SUPPORT)][
                            "oeb_order"
                        ][(v, u)],
                    }
                    G.add_edge(
                        port_u,
                        port_v,
                        (layer, EdgeType.SUPPORT),
                        **{
                            **L.edges[u, v, (layer, EdgeType.SUPPORT)],
                            "oeb_order": oeb_order_pu_pv,
                        },
                    )

                    # port_v -> v_center
                    G.add_edge(
                        port_v,
                        v,
                        (layer, EdgeType.SUPPORT),
                        edge=EdgeType.SUPPORT,
                        sets=set(L.edges[u, v, (layer, EdgeType.SUPPORT)]["sets"]),
                    )
                else:
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

    G.remove_edges_from(
        [(u, v, k) for u, v, k in G.edges(keys=True) if k == EdgeType.CENTER]
    )
    return G
