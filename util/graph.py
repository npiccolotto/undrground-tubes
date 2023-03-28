import numpy as np
import networkx as nx
from util.geometry import dist_euclidean
from itertools import product, pairwise


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


def calculate_path_length(G, path):
    # path is a list of edges
    if len(path) == 0:
        return float("inf")

    length = 0
    for e in path:
        a, b = e
        length += dist_euclidean(G.nodes[a]["pos"], G.nodes[b]["pos"])
    return length


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


def get_shortest_path_between(G, list1, list2):
    shortest = float("inf"), []
    for i in list1:
        for j in list2:
            p = nx.shortest_path(G, i, j)
            p = list(pairwise(p))
            l = calculate_path_length(G, p)
            if l < shortest[0]:
                shortest = (l, p)
    return shortest
