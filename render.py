import networkx as nx
import numpy as np
from collections import defaultdict
from itertools import combinations, pairwise, product
import math
from enum import IntEnum
import drawsvg as svg

from util.collections import merge_alternating
from util.geometry import (
    get_side,
    get_closest_point,
    centroid,
    are_faces_adjacent,
    interpolate_biarcs,
    get_angle,
    offset_edge,
    is_point_inside_circle,
    get_segment_circle_intersection,
)
from util.graph import (
    incident_edges,
    get_longest_simple_paths,
    get_closest_pair,
    get_shortest_path_between,
)


class EdgeType(IntEnum):
    ANCHOR = -3  # logical edges that run in the margins between glyphs
    SET = (
        -2
    )  # logical set edges, i.e., these edges connect elements of the same set _completely_
    PHYSICAL = -1  # possible edge that can appear in the drawing
    NEIGHBOR = 0  # logical edge that tells neighborhood relations of the glyphs
    DRAW = 1  # physical edge that is actually drawn


class NodeType(IntEnum):
    CENTER = 0
    CORNER = 1
    SIDE = 2
    ANCHOR = 3


def logical_coords_to_physical(x, y, lattice_type="sqr"):
    if lattice_type == "hex":
        if y % 2 == 1:
            return (x + 0.5, -y)
    return (x, -y)


DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 100,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
}

DEBUG = False


def determine_lattice_size(glyph_positions):
    # TODO implement
    return [10, 10]


def make_hex_graph(m, n):
    # each hex is 17 nodes and...  136 edges :zipface:
    # = 1 center, 6 corners and 6 sides.
    # fully connect corners, sides and center
    # neighbor edges run between adjacent sides of hexes
    side_length = 0.25

    # start with square graph
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.NEIGHBOR
    for logpos, n in G.nodes(data=True):
        x, y = logpos
        n["node"] = NodeType.CENTER
        n["belongs_to"] = logpos
        n["logpos"] = logpos
        n["pos"] = logical_coords_to_physical(x, y, "hex")
        # TODO support glyph spacing in x and y separately

    G_ = nx.MultiGraph(incoming_graph_data=G)

    for logpos, hex_center in G_.nodes(data=True):
        x, y = logpos

        # add diagonal edges
        if y % 2 != 0 and y > 0 and x < m - 1:
            # tilting to right in odd rows
            G.add_edge(
                (x, y), (x + 1, y - 1), EdgeType.NEIGHBOR, edge=EdgeType.NEIGHBOR
            )
        if y % 2 == 0 and y > 0 and x > 0:
            # tilting to left in even rows
            G.add_edge(
                (x, y), (x - 1, y - 1), EdgeType.NEIGHBOR, edge=EdgeType.NEIGHBOR
            )

        hex_corners = [
            (0, side_length),  # N
            (side_length * math.sqrt(3) / 2, side_length / 2),  # NE
            (side_length * math.sqrt(3) / 2, -side_length / 2),  # SE
            (0, -side_length),  # S
            (-side_length * math.sqrt(3) / 2, -side_length / 2),  # SW
            (-side_length * math.sqrt(3) / 2, side_length / 2),  # NW
        ]
        pos = hex_center["pos"]
        hex_corners = [(x + pos[0], y + pos[1]) for x, y in hex_corners]
        for corner in hex_corners:
            G.add_node(corner, node=NodeType.CORNER, belongs_to=logpos, pos=corner)

        sides = pairwise(hex_corners + [hex_corners[0]])
        hex_sides = []
        for c1, c2 in sides:
            n = centroid([c1, c2])
            hex_sides.append(n)
            G.add_node(n, node=NodeType.SIDE, belongs_to=logpos, pos=n)

        # fully connect corners and sides
        corners_and_sides = hex_corners + hex_sides
        for n1, n2 in combinations(corners_and_sides, 2):
            G.add_edge(n1, n2, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL)

        # connect center to corners and sides
        for n in corners_and_sides:
            G.add_edge(logpos, n, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL)

    for logpos, n in G.nodes(data=True):
        x, y = logpos
        if "pos" not in n:
            n["pos"] = logpos

    return G


def make_sqr_graph(m, n):
    # TODO extend similar to hex graph (9 nodes total: 1 center, 4 corners, 4 sides)

    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for pos, u in G.nodes(data=True):
        x, y = pos
        u["logpos"] = (x, y)
        u["node"] = NodeType.CENTER
        u["pos"] = logical_coords_to_physical(x, y)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = EdgeType.NEIGHBOR
    return G


def convert_logical_to_physical_graph(G: nx.Graph, lattice_type, lattice_size):
    """Adds physical paths and routing nodes."""
    # idea: we know where faces are based on grid position and grid type, so we just compute them
    # hex and sqr have the same "base" grid, hex is just shifted and has a few more edges
    # every 2x2 quad [A,B,C,D] cw in the sqr grid is a face
    #   A - B
    #   | . |
    #   D - C
    # hex faces are double - there is an additional edge in the quad to make it 2 triangles
    # the edge is bottom left to top right (B-D) in even rows and bottom right to top left (A-C)
    #   A - B   A - B
    #   |.\.|   |./.|
    #   D - C   D - C
    m, n = lattice_size
    for y in range(0, n - 1):
        for x in range(0, m - 1):

            quad = [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)]
            if lattice_type == "sqr":
                u, v = centroid(quad)
                G.add_node((u, v), node=NodeType.ANCHOR, face=quad, pos=(u, v))

            elif lattice_type == "hex":
                a, b, c, d = quad
                if y % 2 == 0:
                    # ABD and CBD triangles
                    u, v = centroid(
                        list(
                            map(
                                lambda x: logical_coords_to_physical(x[0], x[1], "hex"),
                                [a, b, d],
                            )
                        )
                    )
                    G.add_node((u, v), node=NodeType.ANCHOR, face=[a, b, d], pos=(u, v))
                    u, v = centroid(
                        list(
                            map(
                                lambda x: logical_coords_to_physical(x[0], x[1], "hex"),
                                [c, b, d],
                            )
                        )
                    )
                    G.add_node((u, v), node=NodeType.ANCHOR, face=[c, b, d], pos=(u, v))
                else:
                    # ABC and ACD triangles
                    u, v = centroid(
                        list(
                            map(
                                lambda x: logical_coords_to_physical(x[0], x[1], "hex"),
                                [a, b, c],
                            )
                        )
                    )
                    G.add_node((u, v), node=NodeType.ANCHOR, face=[a, b, c], pos=(u, v))
                    u, v = centroid(
                        list(
                            map(
                                lambda x: logical_coords_to_physical(x[0], x[1], "hex"),
                                [a, c, d],
                            )
                        )
                    )
                    G.add_node((u, v), node=NodeType.ANCHOR, face=[a, c, d], pos=(u, v))

    # this is now a bit inefficient but should be fine in practice i guess
    # to add all anchor edges we have to insert all anchor nodes first, so we can't do it in the previous procedure
    # here we loop quadratically over all anchor nodes (faces), connect it to its glyph nodes and adjacent anchor nodes (faces)
    anchors = [(i, u) for i, u in G.nodes(data=True) if u["node"] == NodeType.ANCHOR]

    # anchors connect to nearest corner nodes of adjacent glyphs
    for i, a in anchors:
        face = a["face"]
        for center in face:
            G.add_edge(i, center, EdgeType.ANCHOR, edge=EdgeType.ANCHOR)
            # find corner nodes of center
            corner_nodes = nx.subgraph_view(
                G,
                filter_node=lambda p: G.nodes[p]["node"] == NodeType.CORNER
                and G.nodes[p]["belongs_to"] == center,
            )
            closest = get_closest_point(i, list(corner_nodes.nodes))
            G.add_edge(closest, i, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL)

    # and to neighboring anchors
    for a, b in combinations(anchors, 2):
        k, u = a
        l, v = b
        face_u = u["face"]
        face_v = v["face"]
        if are_faces_adjacent(face_u, face_v):
            G.add_edge(k, l, EdgeType.ANCHOR, edge=EdgeType.ANCHOR)
            G.add_edge(k, l, EdgeType.PHYSICAL, edge=EdgeType.PHYSICAL)

    # make physical connections for neighboring glyphs
    logical_edges = [
        (u, v, k) for u, v, e in G.edges(data=True) if e["edge"] == EdgeType.NEIGHBOR
    ]
    for u, v, k in logical_edges:
        u_sides = [w for w in G.neighbors(u) if G.nodes[w]["node"] == NodeType.SIDE]
        v_sides = [w for w in G.neighbors(v) if G.nodes[w]["node"] == NodeType.SIDE]

        closest_u, closest_v = get_closest_pair(u_sides, v_sides)
        G.add_edge(
            u_sides[closest_u],
            v_sides[closest_v],
            EdgeType.PHYSICAL,
            edge=EdgeType.PHYSICAL,
        )

    return G


def get_routing_graph(lattice_type, lattice_size):
    """Returns a networkx graph. It contains two types of nodes and two types of edges.
    Nodes can be 'glyph' nodes, which correspond to spots where glyphs will be placed.
    Nodes can also be 'anchor' nodes, which corresponds to anchor points in the margins of the layout along which we trace lines and anchor polygons. We place them in the center of 'glyph' node faces.
    Edges can be 'neighbor' edges, which corresponds to direct neighbor relations of glyphs, e.g., in a sqr grid most nodes have 4 neighbors.
    Edges can also be 'anchor' edges. They connect both 'glyph' and 'anchor' nodes. In a hex lattice, faces (polygons between 'glyph' nodes) are triangles. So each 'anchor' node has 6 'anchor' edges incident: 3 for neighboring anchors and 3 for glyphs on the boundary of the face."""
    m, n = lattice_size
    G = None
    match lattice_type:
        case "hex":
            G = make_hex_graph(m, n)
        case "sqr":
            G = make_sqr_graph(m, n)
        case _:
            raise Exception(f"unknown lattice type {lattice_type}")

    H = convert_logical_to_physical_graph(G, lattice_type, lattice_size)
    return H


def bundle_edges(G):
    # several ideas how we could do this...

    # greedy heuristic
    # the idea would be roughly the following
    # first, compute *all* shortest paths from a to b and store
    # having done that for all set edges, we take the unprocessed edge with the longest shortest path p
    # set p as processed and select any of the shortest paths (TODO selection of that path is actually quite influential for the drawing, in general we want to try it with all options)
    # find all p' that share source or target node with p and sort by length
    # for each p', compute alternatives p'', which are paths that use "one more" anchor of p
    # for each p'', compute ratios of length and shared edges in relation to shortest path p'
    # somehow parametrize this decision and choose an alternative that balances additional length with additional shared edges
    # set p' as processed and start again

    # winding roads
    # basically the same but use a weighted graph where edge weight decreases the more paths use it
    # obv processing order is important then? idk

    # for now only pick one of the shortest graphs without actually bundling anything
    return G


def add_glyphs_to_nodes(instance, G):
    for i, glyph in enumerate(instance["glyph_ids"]):
        logpos = instance["glyph_positions"][i]
        if G.nodes[logpos]["node"] != NodeType.CENTER:
            raise Exception("node to position glyph on is somehow not a glyph center")
        G.nodes[logpos]["occupied"] = True
        G.nodes[logpos]["glyph"] = glyph
    return G


def add_set_edges(instance, G):
    for set_id, elements in instance["set_system"].items():
        # +1 because 0 are neighbor edges. this way we have -1 = anchor, 0 = neighbor, everything after: set in ftb order
        set_ftb_order = instance["set_ftb_order"].index(set_id) + 1
        if set_ftb_order <= 0:
            raise Exception(f"set {set_id} unknown?")
        for a, b in combinations(elements, 2):
            ai = instance["glyph_ids"].index(a)
            bi = instance["glyph_ids"].index(b)
            logpos_a = instance["glyph_positions"][ai]
            logpos_b = instance["glyph_positions"][bi]

            G_ = nx.subgraph_view(G, filter_edge=lambda u, v, k: k == EdgeType.PHYSICAL)

            shortest_paths = [
                list(pairwise(p)) for p in nx.all_shortest_paths(G_, logpos_a, logpos_b)
            ]

            # add the set edge along with physical shortest paths
            G.add_edge(
                logpos_a,
                logpos_b,
                set_ftb_order,
                set_id=set_id,
                edge=EdgeType.SET,
                shortest_paths=shortest_paths,
            )
    return G


def embed_to_routing_graph(instance, G):
    """Routing graph is a graph with positioned nodes. This function embeds glyph nodes and set relations into the host graph.
    Specifically, it adds i) the glyph to their respective nodes as property and ii) additional edges that correspond to the sets."""
    # 1) distribute glyphs onto glyph nodes
    G = add_glyphs_to_nodes(instance, G)
    # 2) add logical set edges between all pairs of elements in a set - this is the reachability graph of kelpfusion
    G = add_set_edges(instance, G)

    return G


def bundle_edges(instance, M):
    # here we re-route the drawn edges in M (hopefully) such that we reduce the number of unique paths
    # what may and may not happen:
    # - edges (directly) connecting neighboring glyphs MUST NOT take a detour through the margins. reason being that these will be the visually most prominent lines.
    # - edges (indirectly) connecting glyphs through anchor nodes MAY take a detour.
    # - this detour MUST NOT include glyph nodes (not even of the same set b/c that may fuck up the general routing connections - e.g. when there are 3 components and the route connecting component 1 and 2 detours over 3, then we could have duplicate connections.)

    # so rough process:
    # - after converting M into a list of simple paths
    # - process them in order by longest first
    # - temporarily increase weight of edges incident to non-source/target nodes to Inf.
    # - find the shortest path
    # - replace current with new shortest path
    # - increase counter for used edges incident to an anchor
    # - set weight to some function of counter. use y= a*x^b+c and find some good parameters a,b,c
    # - repeat until all paths are processed
    return M


def find_edge_ordering(G, p1, p2, u, v, initial_edge, exploring_other_way=False):
    # walk along their common edges starting from u
    # find edges from u to any node but v. G contains two paths, that were up until now
    # coincident, so either we find 0 edges (paths end) or 2 edges (paths continue).
    edges = incident_edges(G, u, avoid_node=v)
    order = 0
    while len(edges) > 0:
        edges = dict([(G.edges[(u, v, k)]["path_id"], v) for u, v, k in edges])
        next_p1 = edges[p1]
        next_p2 = edges[p2]

        if next_p1 != next_p2:
            # the two paths fork at u. path going to left takes precendence over the other.
            # left of what? left of the first edge of left of this edge?
            # let's assume left of this edge.
            side = get_side(u, next_p2, next_p1)
            if not exploring_other_way:
                # we search backwards first
                # so if we find an order when looking back we have to reverse it
                side *= -1
            return side
        else:
            # paths both use same edge
            # check if that edge has an ordering and apply it here, if so, and break. else walk further.
            e = [
                (w, x, k)
                for w, x, k in G.edges(nbunch=u, keys=True)
                if (w, x) == (u, next_p1)
            ]
            e = e[0]
            edata = G.edges[e]
            if "oeb_order" in edata:
                ordering = edata["oeb_order"][(u, next_p1)]
                if ordering.index(p1) < ordering.index(p2):
                    order = -1
                elif ordering.index(p2) < ordering.index(p1):
                    order = 1

                # we search backwards first
                # so if we find an order when looking back we have to reverse it
                if not exploring_other_way:
                    order *= -1
                return order

        edges = incident_edges(G, next_p1, avoid_node=u)
        u = next_p1
        v = u
    # if 0: both paths end at u. repeat process from v and exclude u
    if not exploring_other_way:
        w, x = initial_edge
        return find_edge_ordering(
            G, p1, p2, x, w, initial_edge, exploring_other_way=True
        )

    # in case we did all of the above and still didn't find an ordering, the two paths are coincident (=the same)
    return order


def get_linear_order(O):
    n, _ = O.shape
    # "dependency" graph. edge (i,j) means i must precede j
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    triu = list(zip(*np.triu_indices(n, 1)))
    for u, v in triu:
        if O[u, v] > 0:
            G.add_edge(u, v)
        if O[u, v] < 0:
            G.add_edge(v, u)
    linearized = []
    for comp in nx.connected_components(nx.Graph(incoming_graph_data=G)):
        # start at the node with no in-edges
        # hopefully there is only one
        num_in_edges = [(n, len(G.in_edges(nbunch=n))) for n in comp]
        num_in_edges = list(sorted(num_in_edges, key=lambda x: x[1]))
        src, l_src = num_in_edges[0]
        if l_src > 0:
            raise Exception("cyclellelellele")
        dep_order = nx.dfs_postorder_nodes(G, src)
        for n in dep_order:
            linearized.append(n)

    return linearized


def order_bundles(instance, M):
    """M is a nx MultiGraph with all edges to draw. Returned is an ordering of sets."""

    # we have to preprocess M a little bit so that Pupyrev et al.'s algorithm (ordered edge bundles, OEB) works.
    # 1. drawn edges for a set in M must not contain cycles.
    #   - detect cycles, throw (for now) arbitrary edge away
    #   - since drawn edges are built out of MSTs there can't currently be any cycles so it's fine
    # 2. OEB expects that i) nodes are either terminals or intermediates, never both, and ii) lines are simple paths, i.e., don't fork/merge.
    #   - convert M to separate simple paths. splitting it at nodes with deg > 2 should be sufficient. so each simple path starts and ends at a glyph center (=terminal).
    #   - for each simple path, remove glyph centers that are not the first or last node in the path: replace adjacent edges (u,v),(v,w) with (u,w) if v is glyph center.
    # after these steps, OEB should work.

    paths = defaultdict(list)
    # step 2: split forking paths
    for set_id in instance["set_ftb_order"]:
        # get edges for this set
        G_ = nx.subgraph_view(
            M, filter_edge=lambda i, j, e: M.edges[i, j, e]["set_id"] == set_id
        )
        # then, start at any node with deg=1
        first = list([i for i in G_.nodes() if G_.degree[i] == 1])[0]
        next = list(G_.neighbors(first))[0]
        paths[set_id] = paths[set_id] + get_longest_simple_paths(G_, first, next)

    # step 2: remove intermediate center nodes from paths
    for set_id in instance["set_ftb_order"]:
        paths[set_id] = [
            [
                n
                for i, n in enumerate(path)
                if i in [0, len(path) - 1] or M.nodes[n]["node"] != NodeType.CENTER
            ]
            for path in paths[set_id]
        ]

    # TODO unsure if good to modify this here, but whatevs
    M.remove_edges_from(list(M.edges()))
    for set_id in instance["set_ftb_order"]:
        set_ftb_order = instance["set_ftb_order"].index(set_id) + 1
        for i, path in enumerate(paths[set_id]):
            for u, v in pairwise(path):
                M.add_edge(
                    u,
                    v,
                    set_ftb_order,
                    set_id=set_id,
                    path_id=f"{set_id}-{i}",
                    edge=EdgeType.DRAW,
                )

    # then do the actual OEB algorithm
    unique_edges = list(set(M.edges()))
    for u, v in unique_edges:
        # find paths using this edge,ie. all edges (u,v,k) for any k
        this_edge = nx.subgraph_view(M, filter_edge=lambda w, x, _: (u, v) == (w, x))
        p_uv = [k for u, v, k in this_edge.edges(keys=True)]
        if len(p_uv) < 2:
            # nothing to do here! only one path using this edge.
            k = p_uv[0]
            order = [M.edges[(u, v, k)]["path_id"]]
            M.edges[(u, v, k)]["oeb_order"] = {
                (u, v): order,
                (v, u): order,
            }
            continue

        path_ids = [
            (k, this_edge.edges[(u, v, k)]["path_id"])
            for u, v, k in this_edge.edges(keys=True)
        ]

        # make a square matrix holding the relative orderings: O[a,b] = -1 -> a precedes b. 0 = don't matter, 1 = succeeds, else = don't know
        O = np.empty(shape=(len(p_uv), len(p_uv)))
        np.fill_diagonal(O, 0)

        # direction of (u,v): u -> v
        # (acc. to paper apparently it doesn't matter)

        # for each pair of paths
        for pi1, pi2 in combinations(path_ids, 2):
            k1, p1 = pi1
            k2, p2 = pi2
            # filter here to path ids identified earlier so that we don't deal here with the path itself forking
            M_ = nx.subgraph_view(
                M, filter_edge=lambda w, x, k: M.edges[(w, x, k)]["path_id"] in [p1, p2]
            )
            order = find_edge_ordering(M_, p1, p2, u, v, (u, v))
            # print((u, v), p1, p2, order)
            k1i = p_uv.index(k1)
            k2i = p_uv.index(k2)
            O[k1i, k2i] = order
            O[k2i, k1i] = -order

        # linearize O to something like (b,a,d,...,z)
        path_id_dict = dict(path_ids)
        linear_O = [path_id_dict[p_uv[i]] for i in get_linear_order(O)]
        # save order on all multiedges between u and v
        for k in p_uv:
            M.edges[(u, v, k)]["oeb_order"] = {
                (u, v): linear_O,
                (v, u): list(reversed(linear_O)),
            }

    return M


def add_routes_of_set_systems(instance, G):
    M = nx.MultiGraph(incoming_graph_data=G)
    M.remove_edges_from(list(G.edges()))

    for set_id, elements in instance["set_system"].items():
        set_idx = 1 + instance["set_ftb_order"].index(set_id)
        # 1) filter graph to neighbor edges between elements
        G_ = nx.subgraph_view(
            G,
            filter_edge=lambda u, v, k: k == EdgeType.NEIGHBOR
            and G.nodes[u]["glyph"] in elements
            and G.nodes[v]["glyph"] in elements,
        )

        # 2) compute a spanner or MST, call it S
        # S = nx.spanner(nx.Graph(incoming_graph_data=G_), 100)
        S = nx.minimum_spanning_tree(nx.Graph(incoming_graph_data=G_))

        # components with size > 1
        components = [list(c) for c in nx.connected_components(S) if len(c) > 1]
        # but we also may have components with size = 1, ie., elements as islands
        # by logic these must be the set elements that not in components with size > 1
        element_positions = list(
            map(
                lambda e: instance["glyph_positions"][instance["glyph_ids"].index(e)],
                elements,
            )
        )
        elements_in_large_components = [
            item for sublist in components for item in sublist
        ]
        element_islands = [
            e for e in element_positions if e not in elements_in_large_components
        ]

        for e in element_islands:
            components.append([e])

        # 3) if S has 1 component, cool. replace the neighbor edges the physical path from the appropriate set edge and return
        if len(components) > 1:
            # 4) else, make new graph G' where V = components of S and E = V x V with weight = length of shortest path from v1 to v2
            G_ = nx.Graph()
            G_.add_nodes_from(range(len(components)))
            for i, j in combinations(range(len(components)), 2):
                nodeset_a = components[i]
                nodeset_b = components[j]
                G_path = nx.subgraph_view(
                    G,
                    filter_edge=lambda u, v, k: k == EdgeType.PHYSICAL
                    and (
                        (
                            G.nodes[u]["node"] == NodeType.ANCHOR
                            or G.nodes[v]["node"] == NodeType.ANCHOR
                        )
                        or (
                            G.nodes[u].get("belongs_to") in nodeset_a + nodeset_b
                            or G.nodes[v].get("belongs_to") in nodeset_a + nodeset_b
                        )
                    ),
                )
                weight, shortest_path = get_shortest_path_between(
                    G_path, nodeset_a, nodeset_b
                )
                G_.add_edge(i, j, weight=weight, shortest_path=shortest_path)
            # 5) compute spanner or MST again, call it P
            # P = nx.spanner(G_, 100)
            P = nx.minimum_spanning_tree(G_)
            # 6) merge S and P, replace logical with physical paths
            # add P
            for u, v in P.edges:
                d = G_.edges[(u, v)]
                M.add_edges_from(
                    d["shortest_path"],
                    edge=EdgeType.DRAW,
                    set_id=set_id,
                )

        # add S
        for u, v in S.edges():
            M.add_edges_from(
                G.edges[(u, v, set_idx)]["shortest_paths"][0],
                edge=EdgeType.DRAW,
                set_id=set_id,
            )
    return M


def geometrize(instance, M):
    one_unit_px = DEFAULT_PARAMS["unit_size_in_px"]
    margins = np.array((0.5, 0.5)) * one_unit_px
    factor = one_unit_px
    geometries = []
    mx, my = margins

    # project nodes
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        M.nodes[i]["pos"] = (x * factor + mx, -y * factor + my)

    # glyph nodes
    for i, n in M.nodes(data=True):
        if n["node"] == NodeType.CENTER:
            corners = [
                M.nodes[m]["pos"]
                for m in nx.subgraph_view(
                    M,
                    filter_node=lambda p: G.nodes[p]["node"] == NodeType.CORNER
                    and G.nodes[p]["belongs_to"] == i,
                )
            ]
            xs, ys = zip(*corners)
            # so if we wanted to render to something else than svg, we could
            # replace the drawsvg elements with some dicts and make drawsvg
            # elements in draw_svg function. convert to something else in
            # draw_smth_else.
            # for now it saves time to not do that
            glyph = svg.Lines(*merge_alternating(xs, ys), close=True)
            geometries.append(glyph)

    set_colors = ["red", "blue", "orange", "green", "magenta"]
    uniq_edges = list(set(M.edges()))
    for u, v in uniq_edges:
        ks = [k for w, x, k in M.edges(keys=True) if (w, x) == (u, v)]
        src = M.nodes[u]["pos"]
        tgt = M.nodes[v]["pos"]
        edge_angle = get_angle(src, tgt)

        # print normal of edge vector just for debugging
        # cx1, cy1 = centroid([src, tgt])
        # cx2, cy2 = offset_point((cx1, cy1), edge_angle + math.pi / 2, 10)
        # geometries.append(
        #    svg.Line(cx1, cy1, cx2, cy2, stroke="lightgray", stroke_width=1)
        # )

        paths_at_edge = M.edges[(u, v, ks[0])]["oeb_order"][(u, v)]
        offset = factor / 20  # TODO idk some pixels
        centering_offset = ((len(paths_at_edge) - 1) / 2) * -offset
        for i, path_id in enumerate(paths_at_edge):
            set_id, _ = path_id.split("-")  # TODO meh
            set_idx = instance["set_ftb_order"].index(set_id)
            offset_dir = math.pi / 2
            offset_length = centering_offset + i * offset
            o_u, o_v = offset_edge((src, tgt), edge_angle - offset_dir, offset_length)

            M.edges[(u, v, set_idx + 1)]["edge_pos"] = {
                (u, v): (o_u, o_v),
                (v, u): (o_v, o_u),
            }

    paths = list(set([d["path_id"] for u, v, k, d in M.edges(keys=True, data=True)]))
    for path_id in paths:
        set_id, _ = path_id.split("-")  # TODO meh
        set_idx = instance["set_ftb_order"].index(set_id)
        M_ = nx.subgraph_view(
            M, filter_edge=lambda u, v, k: M.edges[(u, v, k)]["path_id"] == path_id
        )
        endpoints = [n for n in M_.nodes() if M_.degree(n) == 1]
        path = nx.shortest_path(M_, endpoints[0], endpoints[1])
        # idk, like keyframe... find better word
        keypoints = [
            # (type of connection to following point, keypoint coordinates)
        ]

        for i, pair in enumerate(pairwise(path)):
            u, v = pair
            ux, uy = M_.nodes[u]["pos"]
            vx, vy = M_.nodes[v]["pos"]

            # at each node, we find segments from the last key point to the next node's hub circle
            hub_radius = 1 / 16 * factor
            uhub_center = (ux, uy)
            uhub_circle = (uhub_center, hub_radius)

            vhub_center = (vx, vy)
            vhub_circle = (vhub_center, hub_radius)

            edge = M_.edges[(u, v, set_idx + 1)]["edge_pos"][(u, v)]
            a, b = edge

            if is_point_inside_circle(a, uhub_circle):
                # if start point of edge is inside hub circle,
                # we want a straight line from the hub circle circumference to the following point
                keypoints.append(
                    ("L", get_segment_circle_intersection((a, b), uhub_circle))
                )
            else:
                # start point of edge is not inside the hub circle
                # so we need an arc from the hub circle circumference to the edge start
                keypoints.append(
                    (
                        "B",
                        get_segment_circle_intersection((uhub_center, a), uhub_circle),
                    ),
                )
                # and then a straight line
                keypoints.append(("L", a))

            if is_point_inside_circle(b, vhub_circle):
                keypoints.append(
                    ("B", get_segment_circle_intersection((a, b), vhub_circle))
                )
            else:
                # edge to hub
                keypoints.append(("B", b))
                keypoints.append(
                    (
                        "B" if i + 1 < len(path) - 1 else "L",
                        get_segment_circle_intersection((vhub_center, b), vhub_circle),
                    ),
                )

        kwargs = {
            "close": False,
            "stroke_width": 2,
            "fill": "none",
            "stroke": set_colors[set_idx],
            "z": set_idx,
        }
        line = interpolate_biarcs(keypoints, **kwargs)
        geometries.append(line)

    if DEBUG:
        hubs = [n for n in M.nodes() if M.degree[n] > 0]
        for hub in hubs:
            cx, cy = M.nodes[hub]["pos"]
            r = 1 / 16 * factor
            geometries.append(
                svg.Circle(cx, cy, r, fill="none", stroke="gray", stroke_width=1)
            )
    return geometries


def draw_svg(geometries):
    d = svg.Drawing(2000, 2000, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


def render_line(instance, G, p):

    M = add_routes_of_set_systems(instance, G)

    # optional: bundle edges more
    # idea would be to allow bundling of spanner sub-paths between two element nodes
    # so that even after bundling the modified spanner connects all set elements (but is not a t-spanner anymore necessarily)
    M = bundle_edges(instance, M)

    # next step: ordering of bundles
    # for each edge and each pair of paths that use it, find relative order of paths (before, after) acc to paper
    # set_order_in_bundles = get_bundle_order(M, p)
    M = order_bundles(instance, M)

    # next step: defining line segments to draw, this is the annoying part
    # the first thing is to define the hubs, i.e., circles with given radius at anchor, corner, and glyph center nodes. move them a bit closer so that the actual glyph can be drawn over them.
    # along each hub's outline is a segment where lines from a direction can go... no idea how to specify that (very asymmetric for corner and side nodes)
    # for each bundle b connecting hubs u and v we'd like to find |b| straight lines parallel to the line connecting u and v centers, spaced such that they are within the appointed segment on each hub
    # for each anchor hub, connect incident lines of the same path with a biarc (or, simpler but ugly, a straight line)

    # then return the geometries with parameters and hand it off to somewhere for drawing

    return M


def render_envelope(instance, G, p):
    pass


def render_kelpfusion(instance, G, p):
    # sketch:
    # for each set S_i
    # acc to meulemans 2013
    #   compute reachability graph R of S_i, where an edge (u,v) runs in the margins (along 'anchor' edges) if (u,v) is not in host graph. use geometric length of e as weight. R contains an edge for each a,b in S_i a!=b
    #   compute shortest path graph SPG of R
    #   find faces = cycles in SPG
    #   for each face: define if it gets filled according to paper
    #   SOMEHOW render faces and edges
    pass


INSTANCE = {
    "lattice_type": "hex",
    "glyph_ids": [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
    ],  # list of strings (could be file urls)
    # A B C
    #  D E F
    # G H I
    "glyph_positions": [
        (0, 0),
        (1, 0),
        (2, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (0, 2),
        (1, 2),
        (2, 2),
    ],  # a panda df with two columns (x and y) corresponding to logical grid position
    "set_system": {
        "set0": ["A", "B", "H", "D", "E"],
        "set1": ["B", "I"],
        "set2": ["G", "C", "D", "E"],
        "set3": ["A", "F", "D", "E"],
        "set4": ["A", "I", "D", "E"],
    },
    "set_ftb_order": [
        "set4",
        "set3",
        "set2",
        "set1",
        "set0",
    ],  # a list of set ids that defines front to back ordering (index = 0 is most front)
}

if __name__ == "__main__":
    m = 3
    n = 3
    lattice_type = INSTANCE["lattice_type"]
    G = get_routing_graph(lattice_type, (m, n))
    G = embed_to_routing_graph(INSTANCE, G)
    G = render_line(INSTANCE, G, DEFAULT_PARAMS)
    geometries = geometrize(INSTANCE, G)
    img = draw_svg(geometries)

    with open("drawing.svg", "w") as f:
        f.write(img)
        f.flush()
