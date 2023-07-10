import numpy as np
import networkx as nx
import networkx.algorithms.approximation.traveling_salesman as tsp
import networkx.algorithms.approximation.steinertree as steinertree
from util.geometry import dist_euclidean
from util.perf import timing
from itertools import product, pairwise, combinations
from util.enums import NodeType


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
    s = "s"
    t = "t"
    temp_edges = []
    G.add_node(s)
    for i in S:
        G.add_edge(s, i, weight=0)
        temp_edges.append((s, i))
    for j in T:
        G.add_edge(t, j, weight=0)
        temp_edges.append((t, j))
    sp = nx.shortest_path(G, s, t, weight="weight")
    sp = [node for node in sp if node != s and node != t]
    # cleanup
    for u, v in temp_edges:
        G.remove_edge(u, v)
    G.remove_node(s)
    G.remove_node(t)
    return path_to_edges(sp)


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
    # take any node as source
    src = S[0]
    # keep track of nodes we must eventually find
    # assume that S is already connected
    N = set(T)
    for n in nx.dfs_preorder_nodes(G, src):
        if n in N:
            N = N.difference([n])
    return len(N) == 0


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


def approximate_steiner_tree(G, S):
    """G is a weighted graph, S is the set of terminal nodes. Returns a tree subgraph of G that is a Steiner tree.

    Wu et al., 1986: A faster approximation algorithm for the Steiner problem in graphs
    """
    # step 1: make G1, a complete graph of all the terminals where edge weight is shortest path length. pick any shortest path if it's not unique
    shortest_paths_dict = {}
    G1 = nx.Graph()
    G1.add_nodes_from(S)
    for n1, n2 in combinations(S, 2):
        sp = nx.shortest_path(G, n1, n2, weight="weight")
        sp_edgelist = path_to_edges(sp)
        spl = calculate_path_length(G, sp_edgelist, weight="weight")
        shortest_paths_dict[(n1, n2)] = (spl, sp_edgelist)
        shortest_paths_dict[(n2, n1)] = (spl, sp_edgelist)
        G1.add_edge(n1, n2, weight=spl)

    # step 2: make G2, a MST on G1
    G2 = nx.minimum_spanning_tree(G1, weight="weight", algorithm="prim")

    # step 3: make G3 by starting with G nodes and no edges. for every edge in G2 add a shortest path between the endpoints in G
    G3 = nx.Graph()
    G3.add_nodes_from(G)
    for u, v in G2.edges():
        spl, sp = shortest_paths_dict[(u, v)]
        for w, x in sp:
            G3.add_edge(w, x, weight=G.edges[(w, x)]["weight"])

    # step 4: make G4, a MST on G3
    G4 = nx.minimum_spanning_tree(G3, weight="weight")

    # step 5: make G5 by removing any non-terminal leaves from G4
    while (
        len(
            non_terminal_leaves := [
                n for n in G4.nodes() if G4.degree[n] == 1 and n not in S
            ]
        )
        > 0
    ):
        G4.remove_nodes_from(non_terminal_leaves)

    return G4


def approximate_tsp_tour(G, S):
    """Returns an approximation of the shortest tour in G visiting all nodes S exactly once"""
    path = tsp.traveling_salesman_problem(
        G, nodes=S, cycle=False, weight="weight", method=tsp.greedy_tsp
    )
    return path

def get_node_with_degree(G, deg):
    nodes_and_deg = [(n, G.degree(n)) for n in G.nodes()]
    deg_nodes = [n for n, d in nodes_and_deg if d == deg]
    return deg_nodes[0] if len(deg_nodes) > 0 else None

def visit_edge_pairs_starting_at_node(G, start):
    # so here the idea is that for drawing we need edge pairs
    # we find edge pairs by starting at a node
    edge_pairs = []
    stack = []
    current_edge=None

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

        u,v = current_edge
        next_edges = list([(w,x) for w,x in G.edges(v) if (w,x) != (v,u)])

        match len(next_edges):
            case 0:
                # there is no next edge
                # we're done for now, let's hope there's something on the stack in the next iteration
                current_edge = None
                continue
            case 1:
                # there is 1 next edge
                next_edge =  next_edges[0]
                edge_pairs.append([current_edge,next_edge])
                current_edge = next_edge
                continue
            case _:
                next_edge =  next_edges[0]
                pairs = product([current_edge], next_edges)
                edge_pairs = edge_pairs + list(pairs)
                stack = stack + next_edges[1:]
                current_edge = next_edge
                continue




    pass
