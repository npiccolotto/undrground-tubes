import networkx as nx
import matplotlib as mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from itertools import combinations, pairwise, product
import math
from enum import IntEnum


class EdgeType(IntEnum):
    ANCHOR = -3
    SET = -2
    PHYSICAL = -1
    NEIGHBOR = 0
    # everything from here is a 1-based set index


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
    # TODO extend such that each hex is 17 nodes and...  136 edges :zipface:
    # = 1 center, 6 corners and 6 sides.
    # fully connect corners, sides and center
    # neighbor edges run between adjacent sides of hexes
    # maybe useful to add one more edge type that directly connects centers
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

        G.add_nodes_from(hex_corners, node=NodeType.CORNER, belongs_to=logpos)
        sides = pairwise(hex_corners + [hex_corners[0]])
        hex_sides = []
        for c1, c2 in sides:
            n = centroid([c1, c2])
            hex_sides.append(n)
            G.add_node(n, node=NodeType.SIDE, belongs_to=logpos)

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


def embed_to_routing_graph(instance, G):
    """Routing graph is a graph with positioned nodes. This function embeds glyph nodes and set relations into the host graph.
    Specifically, it adds i) the glyph to their respective nodes as property and ii) additional edges that correspond to the sets."""
    # 1) distribute glyphs onto glyph nodes
    for i, glyph in enumerate(instance["glyph_ids"]):
        logpos = instance["glyph_positions"][i]
        if G.nodes[logpos]["node"] != NodeType.CENTER:
            raise Exception("node to position glyph on is somehow not a glyph center")
        G.nodes[logpos]["occupied"] = True
        G.nodes[logpos]["glyph"] = glyph

    # 2) add logical set edges between all pairs of elements in a set - this is the reachability graph
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


def get_bundle_order(M, p):
    """M is a nx MultiGraph with all edges to draw. Returned is an ordering of sets."""

    # M is a spanner. we have to preprocess it a little bit so that Pupyrev et al.'s algorithm (ordered edge bundles, OEB) works.
    # 1. M must not contain cycles.
    #   - detect cycles, throw (for now) arbitrary edge away
    # 2. OEB expects that i) nodes are either terminals or intermediates, never both, and ii) lines are simple paths, i.e., don't fork/merge.
    #   - convert now-cycle-free spanner M to separate simple paths. splitting it at nodes with deg > 2 should be sufficient. so each simple path starts and ends at a glyph center (=terminal).
    #   - for each simple path, remove glyph centers that are not the first or last node in the path: replace adjacent edges (u,v),(v,w) with (u,w) if v is glyph center.
    # after these steps, OEB should work.

    return p[
        "set_bundle_order"
    ]  # btw still wondering in relation to WHAT that order actually is


def render_line(instance, G, p):
    # sketch:
    # for each set:
    # subgraph: keep all anchor edges that incide on elements of set plus neighbor edges that are also set edges
    # find t-spanner in that subgraph
    # add spanner to G
    # return G

    def filter_edge(set_idx, elements):
        def inner(u, v, k):

            is_anchor_edge = k == EdgeType.ANCHOR
            u_is_anchor = G.nodes[u]["node"] == NodeType.ANCHOR
            v_is_anchor = G.nodes[v]["node"] == NodeType.ANCHOR
            both_ends_are_anchors = u_is_anchor and v_is_anchor
            u_is_glyph = G.nodes[u]["node"] == NodeType.CENTER
            v_is_glyph = G.nodes[v]["node"] == NodeType.CENTER
            u_in_set = u_is_glyph and G.nodes[u]["glyph"] in elements
            v_in_set = v_is_glyph and G.nodes[v]["glyph"] in elements
            either_end_in_set = u_in_set or v_in_set

            is_neighbor_edge = k == int(EdgeType.NEIGHBOR)
            has_set_edge = G.has_edge(u, v, set_idx)

            return (
                is_anchor_edge and (both_ends_are_anchors or either_end_in_set)
            ) or (is_neighbor_edge and has_set_edge)

        return inner

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
        S = nx.minimum_spanning_tree(nx.Graph(incoming_graph_data=G_))

        components = [c for c in nx.connected_components(S) if len(c) > 1]
        # 3) if S has 1 component, cool. replace the neighbor edges the physical path from the appropriate set edge and return
        if len(components) == 1:
            for u, v in S.edges():
                M.add_edges_from(
                    G.edges[(u, v, set_idx)]["shortest_paths"][0],
                    edge=EdgeType.SET,
                    set_id=set_id,
                )
        else:
            # 4) else, make new graph G' where V = components of S and E = V x V with weight = length of shortest path from v1 to v2
            # 5) compute spanner or MST again, call it P
            # 6) merge S and P, replace logical with physical paths, retu
            pass

    # optional: bundle edges more
    # idea would be to allow bundling of spanner sub-paths between two element nodes
    # so that even after bundling the modified spanner connects all set elements (but is not a t-spanner anymore necessarily)

    # next step: ordering of bundles
    # for each edge and each pair of paths that use it, find relative order of paths (before, after) acc to paper
    # set_order_in_bundles = get_bundle_order(M, p)

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
        # "set3": ["G", "C"],
        "set4": ["A", "B", "E"]
    },
    "set_bundle_order": [
        "set2",
        "set0",
        "set1",
        "set3",
    ],
    "set_ftb_order": [
        "set4"
        # "set3",
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
        EdgeType.ANCHOR: "#128",
    }
    edge_alpha_map = {
        EdgeType.SET: 1,
        EdgeType.NEIGHBOR: 1,
        EdgeType.PHYSICAL: 1,
        EdgeType.ANCHOR: 1,
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
