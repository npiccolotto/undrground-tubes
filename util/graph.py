import numpy as np
import networkx as nx
import math
import networkx.algorithms.approximation.traveling_salesman as tsp
import networkx.algorithms.approximation.steinertree as steinertree
from util.geometry import dist_euclidean
from util.perf import timing
from itertools import product, pairwise, combinations
from util.enums import NodeType, EdgePenalty, PortDirs, EdgeType


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
    G.edges[edge]["weight"] += EdgePenalty.IN_SUPPORT

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
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING
        if is_right_tilt:
            e = get_port_edge_between_centers(G, (ux - 1, uy), (ux, uy + 1))
            G.edges[e]["weight"] = G.edges[e]["weight"] + EdgePenalty.CROSSING

    return G


def are_port_edges_crossing(us, ut, vs, vt):
    cw_dirs = PortDirs + PortDirs

    ccw_dirs = (
        [PortDirs[0]]
        + list(reversed(PortDirs[1:]))
        + [PortDirs[0]]
        + list(reversed(PortDirs[1:]))
    )

    # start at us and go clockwise to ut, collect all ports in between
    cw_set = []
    i = cw_dirs.index(us["port"]) + 1
    j = cw_dirs.index(ut["port"])
    if j > i:
        cw_set = cw_dirs[i:j]

    # then start at us and go ccw to ut, collect all ports in between
    ccw_set = []
    i = ccw_dirs.index(ut["port"]) + 1
    j = ccw_dirs.index(us["port"])
    if j > i:
        ccw_set = ccw_dirs[i:j]

    # edges cross iff vs in former and vt in latter set (or vice versa)
    port1 = vs["port"]
    port2 = vt["port"]

    are_crossing = True
    if port1 in cw_set and port2 in cw_set and port1 in ccw_set and port2 in ccw_set:
        are_crossing = False
    if port1 in [ut["port"], us["port"]] or port2 in [ut["port"], us["port"]]:
        are_crossing = False
    return are_crossing


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


def get_ports(G, u, v):
    """for two center nodes u and v in G, returns their respective ports that an edge would use"""
    ux, uy = G.nodes[u]["logpos"]
    vx, vy = G.nodes[v]["logpos"]

    dx = ux - vx
    dy = uy - vy
    delta = (dx, dy)

    # dx = -1 -> ux - vx = -1 -> vx > ux -> v is west of u
    # dx = 1 -> ux - vx = 1 -> ux > vx -> v is east of u
    # dy = -1 -> uy - vy = -1 -> vy > uy -> v is south of u

    match delta:
        # v is west
        case (-1, -1):
            # v is in se of u
            return (get_port(G, u, "sw"), get_port(G, v, "ne"))
        case (-1, 0):
            # v is east of u
            return (get_port(G, u, "w"), get_port(G, v, "e"))
        case (-1, 1):
            # v is ne of u
            return (get_port(G, u, "nw"), get_port(G, v, "se"))

        case (0, -1):
            # v is south of u
            return (get_port(G, u, "s"), get_port(G, v, "n"))
        # (0,0) not possible!
        case (0, 1):
            # v is north of u
            return (get_port(G, u, "n"), get_port(G, v, "s"))

        # v is east
        case (1, -1):
            return (get_port(G, u, "se"), get_port(G, v, "nw"))
        case (1, 0):
            return (get_port(G, u, "e"), get_port(G, v, "w"))
        case (1, 1):
            return (get_port(G, u, "ne"), get_port(G, v, "sw"))

    raise BaseException(f"something went wrong: dx={dx}, dy={dy}, u={u}, v={v}")


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
