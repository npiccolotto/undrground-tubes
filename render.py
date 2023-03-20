import networkx as nx
import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from itertools import combinations, pairwise, product
import math
from enum import IntEnum
import drawsvg as svg


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


def are_faces_adjacent(face1, face2):
    # faces are a list of tuples (points)
    # faces are adjacent obv iff they share an edge, ie, two consecutive points
    # however, the way our grid is constructed, any two shared points must be a shared edge
    # thus we count the occurrences of points in both faces and if there are 2 points appearing 2 times, that's an edge
    c = Counter(face1 + face2)
    counter = Counter(list(c.values()))
    return counter[2] == 2


def centroid(points):
    x, y = zip(*points)
    n = len(x)
    return (sum(x) / n, sum(y) / n)


def logical_coords_to_physical(x, y, lattice_type="sqr"):
    if lattice_type == "hex":
        if y % 2 == 1:
            return (x + 0.5, -y)
    return (x, -y)


def dist_euclidean(p1, p2):
    # p1,p2 are tuples
    x1, y1 = p1
    x2, y2 = p2
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def get_closest_point(reference, options):
    dists = [dist_euclidean(reference, option) for option in options]
    closest = np.argmin(np.array(dists))
    return options[closest]


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


def calculate_path_length(G, path):
    # path is a list of edges
    if len(path) == 0:
        return float("inf")

    length = 0
    for e in path:
        a, b = e
        length += dist_euclidean(G.nodes[a]["pos"], G.nodes[b]["pos"])
    return length


DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 50,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
}


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


def make_tri_graph(m, n):
    # TODO implement
    raise Exception("not implemented")


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
        case "tri":
            G = make_tri_graph(m, n)
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


def get_longest_simple_paths(G, src, tgt, visited=[], current_path=[]):
    # follow the path as long as nodes have deg = 2
    # when paths fork (deg>2), add path up until here to list
    # restart procedure along the forks
    # eventually we have visited all nodes

    next_visited = visited + [src]
    if G.degree[tgt] > 2:
        # fork, split here
        current_path = current_path + [tgt]
        next_tgts = [i for i in G.neighbors(tgt) if i not in next_visited]
        next_paths = [
            get_longest_simple_paths(
                G, tgt, next_tgt, visited=next_visited, current_path=[]
            )
            for next_tgt in next_tgts
        ]
        ret = [current_path]
        for p in next_paths:
            ret += p
        return ret

    if G.degree[tgt] == 1:
        # end of the simple path
        current_path = current_path + [tgt]
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
        glsp = get_longest_simple_paths(G_, first, next)
        paths[set_id] = paths[set_id] + glsp

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
        for path in paths[set_id]:
            for i, j in pairwise(path):
                M.add_edge(i, j, set_ftb_order, set_id=set_id, edge=EdgeType.DRAW)

    # then do the actual OEB algorithm

    # bundle order probably will be attribute of a physical edge

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
    margins = (50, 50)
    factor = 100
    geometries = []
    mx, my = margins

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
            xs = [x * factor + mx for x in xs]
            ys = [-y * factor + my for y in ys]
            # https://stackoverflow.com/a/3678938/490524
            flat_corners = [None] * (len(xs) + len(ys))
            flat_corners[::2] = xs
            flat_corners[1::2] = ys
            # so if we wanted to render to something else than svg, we could
            # replace the drawsvg elements with some dicts and make drawsvg
            # elements in draw_svg function.
            # for now it saves time to not do that
            glyph = svg.Lines(*flat_corners, close=True)
            geometries.append(glyph)

    for i, j, e in M.edges(data=True):
        if e["edge"] == EdgeType.DRAW:
            s = M.nodes[i]
            z = M.nodes[j]
            sx, sy = s.get("pos")
            ex, ey = z.get("pos")
            geometries.append(
                svg.Line(
                    sx=sx * factor + mx,
                    sy=-sy * factor + my,
                    ex=ex * factor + mx,
                    ey=-ey * factor + my,
                    stroke="gray",
                    stroke_width=5,
                )
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


def render(instance, p):
    p = DEFAULT_PARAMS | p
    lattice_size = (
        determine_lattice_size(glyph_positions)
        if p["lattice_size"] == "auto"
        else p["lattice_size"]
    )
    G = get_routing_graph(p["lattice_type"], lattice_size)
    G = embed_to_routing_graph(instance, G)

    match p["render_style"]:
        case "line":
            return render_line(instance, G, p)
        case "envelope":
            return render_envelope(instance, G, p)
        case "kelpfusion":
            return render_kelpfusion(instance, G, p)
        case _:
            raise Exception(f'unknown render style {p["render_style"]}')


def draw_rendering(rendering):
    # TODO maybe with PIL for now
    pass


INSTANCE = {
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
        # "set0": ["A", "B", "D", "E"],
        # "set1": ["A", "B", "E"],
        # "set2": ["G", "C"],
        "set3": ["D", "G", "H"],
        "set4": ["A", "H", "G", "F", "E"],
    },
    "set_bundle_order": [
        "set2",
        "set0",
        "set1",
        "set3",
    ],
    "set_ftb_order": [
        "set4",
        "set3",
        # "set2",
        # "set1",
        # "set0"
    ],  # a list of set ids that defines front to back ordering (index = 0 is most front)
}

if __name__ == "__main__":
    m = 3
    n = 3
    lattice_type = "hex"
    G = get_routing_graph(lattice_type, (m, n))
    G = embed_to_routing_graph(INSTANCE, G)
    G = render_line(INSTANCE, G, DEFAULT_PARAMS)
    geometries = geometrize(INSTANCE, G)
    img = draw_svg(geometries)
    print(img)

    if False:
        pos = nx.get_node_attributes(G, "pos")
        node_color_map = {
            NodeType.CENTER: "#882200",
            NodeType.ANCHOR: "#123456",
            NodeType.SIDE: "#123456",
            NodeType.CORNER: "#123456",
        }
        node_size_map = {
            NodeType.CENTER: 300,
            NodeType.ANCHOR: 100,
            NodeType.CORNER: 25,
            NodeType.SIDE: 25,
        }
        nx.draw_networkx_nodes(
            G,
            pos,
            node_shape="h" if lattice_type == "hex" else "s",
            node_size=[node_size_map[node[1]["node"]] for node in G.nodes(data=True)],
            node_color=[node_color_map[node[1]["node"]] for node in G.nodes(data=True)],
        )
        edge_color_map = {
            EdgeType.SET: "#829",
            EdgeType.NEIGHBOR: "#820",
            EdgeType.PHYSICAL: "#000",
            EdgeType.DRAW: "#000",
            EdgeType.ANCHOR: "#128",
        }
        edge_alpha_map = {
            EdgeType.SET: 0,
            EdgeType.NEIGHBOR: 0,
            EdgeType.PHYSICAL: 0,
            EdgeType.DRAW: 1,
            EdgeType.ANCHOR: 0,
        }
        nx.draw_networkx_edges(
            G,
            pos,
            alpha=[edge_alpha_map[e[2]["edge"]] for e in G.edges(data=True)],
            edge_color=[edge_color_map[e[2]["edge"]] for e in G.edges(data=True)],
        )
        nx.draw_networkx_labels(
            G,
            pos,
            dict(zip(INSTANCE["glyph_positions"], INSTANCE["glyph_ids"])),
            font_color="w",
        )
        plt.show()
