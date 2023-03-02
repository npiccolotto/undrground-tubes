import networkx as nx
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from itertools import combinations, pairwise
import math
from enum import Enum


class EdgeType(Enum):
    NEIGHBOR = 0
    ANCHOR = 1


def are_faces_adjacent(face1, face2):
    # faces are a list of tuples (points)
    # faces are adjacent obv iff they share an edge, ie, two consecutive points
    # however, the way our grid is constructed, any two shared points must be a shared edge
    # thus we count the occurrences of points in both faces and if there are 2 points appearing 2 times, that's an edge
    c = Counter(face1 + face2)
    counter = Counter(list(c.values()))
    return counter[2] == 2


def centroid(points, lattice_type="sqr"):
    points = list(
        map(lambda p: logical_coords_to_physical(p[0], p[1], lattice_type), points)
    )
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
    # start with square graph
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for pos, u in G.nodes(data=True):
        x, y = pos

        # add diagonal edges
        if y % 2 != 0 and y > 0 and x < m - 1:
            # tilting to right in odd rows
            G.add_edge((x, y), (x + 1, y - 1), EdgeType.NEIGHBOR)
        if y % 2 == 0 and y > 0 and x > 0:
            # tilting to left in even rows
            G.add_edge((x, y), (x - 1, y - 1), EdgeType.NEIGHBOR)

        u["node"] = "glyph"
        u["logpos"] = (x, y)
        u["pos"] = logical_coords_to_physical(x, y, "hex")

    for _1, _2, e in G.edges(data=True):
        e["edge"] = "neighbor"

    return G


def make_tri_graph(m, n):
    # TODO implement
    raise Exception("not implemented")


def make_sqr_graph(m, n):
    G = nx.grid_2d_graph(m, n)
    G = nx.MultiGraph(incoming_graph_data=G)
    for pos, u in G.nodes(data=True):
        x, y = pos
        u["logpos"] = (x, y)
        u["node"] = "glyph"
        u["pos"] = logical_coords_to_physical(x, y)
    for _1, _2, e in G.edges(data=True):
        e["edge"] = "neighbor"
    return G


def convert_host_to_routing_graph(G: nx.Graph, lattice_type, lattice_size):
    """Given a bare-bones host graph (glyph nodes and neighbor edges), extend it by anchor nodes and edges."""
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
                G.add_node((u, v), node="anchor", face=quad, pos=(u, v))

            elif lattice_type == "hex":
                a, b, c, d = quad
                if y % 2 == 0:
                    # ABD and CBD triangles
                    u, v = centroid([a, b, d], "hex")
                    G.add_node((u, v), node="anchor", face=[a, b, d], pos=(u, v))
                    u, v = centroid([c, b, d], "hex")
                    G.add_node((u, v), node="anchor", face=[c, b, d], pos=(u, v))
                else:
                    # ABC and ACD triangles
                    u, v = centroid([a, b, c], "hex")
                    G.add_node((u, v), node="anchor", face=[a, b, c], pos=(u, v))
                    u, v = centroid([a, c, d], "hex")
                    G.add_node((u, v), node="anchor", face=[a, c, d], pos=(u, v))

    # this is now a bit inefficient but should be fine in practice i guess
    # to add all anchor edges we have to insert all anchor nodes first, so we can't do it in the previous procedure
    # here we loop quadratically over all anchor nodes (faces), connect it to its glyph nodes and adjacent anchor nodes (faces)
    anchors = [(i, u) for i, u in G.nodes(data=True) if u["node"] == "anchor"]

    for i, a in anchors:
        face = a["face"]
        for glyph_node in face:
            G.add_edge(i, glyph_node, EdgeType.ANCHOR, edge="anchor")

    for a, b in combinations(anchors, 2):
        k, u = a
        l, v = b
        face_u = u["face"]
        face_v = v["face"]
        if are_faces_adjacent(face_u, face_v):
            G.add_edge(k, l, EdgeType.ANCHOR, edge="anchor")

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

    H = convert_host_to_routing_graph(G, lattice_type, lattice_size)
    return H


def filter_to_anchor_path(G, a, b):
    def inner(u, v, k):
        is_anchor_edge = k == EdgeType.ANCHOR
        u_is_anchor = G.nodes[u]["node"] == "anchor"
        v_is_anchor = G.nodes[v]["node"] == "anchor"
        both_ends_are_anchors = u_is_anchor and v_is_anchor
        u_is_ab = u in [a, b]
        v_is_ab = v in [a, b]
        either_end_is_ab = u_is_ab or v_is_ab
        return is_anchor_edge and (both_ends_are_anchors or either_end_is_ab)

    return inner


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
    for u, v, k in nx.edges(G):
        e = G.edges[(u, v, k)]
        if k == EdgeType.ANCHOR and e["edge"] == "set":
            e["path"] = e["shortest_paths"][0]

    return G


def embed_to_routing_graph(instance, G):
    """Routing graph is a graph with positioned nodes. This function embeds glyph nodes and set relations into the host graph.
    Specifically, it adds i) the glyph to their respective nodes as property and ii) additional edges that correspond to the sets."""
    # 1) distribute glyphs onto glyph nodes
    for i, glyph in enumerate(instance["glyph_ids"]):
        logpos = instance["glyph_positions"][i]
        if G.nodes[logpos]["node"] != "glyph":
            raise Exception("glyph node is somehow not a glyph node")
        G.nodes[logpos]["occupied"] = True
        G.nodes[logpos]["glyph"] = glyph

    # 2) add logical set edges between all pairs of elements in a set - this is the reachability graph
    for setid, elements in instance["set_system"].items():
        # +1 because 0 are neighbor edges. this way we have -1 = anchor, 0 = neighbor, everything after: set in ftb order
        set_ftb_order = instance["set_ftb_order"].index(setid) + 1
        if set_ftb_order <= 0:
            raise Exception(f"set {setid} unknown?")
        for a, b in combinations(elements, 2):
            ai = instance["glyph_ids"].index(a)
            bi = instance["glyph_ids"].index(b)
            logpos_a = instance["glyph_positions"][ai]
            logpos_b = instance["glyph_positions"][bi]

            phys_paths = []

            if G.has_edge(logpos_a, logpos_b, EdgeType.NEIGHBOR):
                phys_paths = [[(logpos_a, logpos_b)]]
            else:
                # no neighbor edge between these glyphse
                # so make a subgraph of all anchor edges
                # since anchors don't connect two glyph nodes, the shortest path will run just along anchor edges
                G_path = nx.subgraph_view(
                    G, filter_edge=filter_to_anchor_path(G, logpos_a, logpos_b)
                )
                if not nx.has_path(G_path, logpos_a, logpos_b):
                    raise Exception("you dummy")
                phys_paths = [
                    p for p in nx.all_shortest_paths(G_path, logpos_a, logpos_b)
                ]
                # the path has to be longer than two nodes for it to i) run only along anchors ii) have a glyph node at source and target
                # thus safe to use `pairwise`
                phys_paths = [list(pairwise(p)) for p in phys_paths]

            # add the set edge along with physical shortest path and length
            G.add_edge(
                logpos_a,
                logpos_b,
                set_ftb_order,
                edge="set",
                set_id=setid,
                shortest_paths=phys_paths,
            )
    return G


def get_bundle_order(M, p):
    """M is a nx MultiGraph with all edges to draw. Returned is an ordering of sets."""

    # ok so here's the thing.
    # you'd think that ordered bundles by Pupyrev et al. (2016) are a nice an simple option
    # but no. our subgraphs, be it a spanner or TSP tour, 1) connect more than two nodes 2) path terminal property doesn't hold (some paths end at glyph nodes, others don't)
    # one could try to do some splitting into proper paths, treat them separately, put them back together in the most favorable ordering... but requires fast crossing computation unless you obtain additional findings so that you don't have to compute...

    # then okay so maybe let's treat these things as metro lines, i.e., as an MLCM problem, but it's generally NP-hard
    # quite some algorithms were proposed, e.g., by nöllenburg, okamoto...
    # most make some assumptions that are not super useful for me here. e.g., the host graph is a path (nope), lines terminate only at deg=1 nodes (nope), no two lines terminate at same node (nope), 2-sides model (lines enter left and exit right, nope), lines don't fork (could do that)

    # now here's a stupid idea. i expect we generally have few sets, like 6 or so, which amounts to 720 permutations - not THAT many
    # if counting crossings with a given order is fast, we could brute-force this thing and have provably optimal solution

    # OR, and this is kind of funny, maybe don't do it at all? since the sets we're gonna use are ordered categories (e.g., quantiles), there could even be an upside to making this user-steerable. then otoh by definition an element can't be in two categories of the same, like, group (e.g., two quantiles of the same variable). so an ordering of the bundle by, e.g., quantile order of same variable, won't show a nice transition of colors on one node. idk.

    # but if we would know a nice way to get a crossing-minimizing order, we could do it here

    # OKAY I THINK I KNOW
    # so basically, make the host graph more detailed and add the corners of the shape. have centers still connected bei neighbor edges so that we can do the logical routing easily (MST, spanner and such). the corners of the shape are fully connected by anchor shapes, the center connects to all corners but should probs have a different edge type so that we don't accidentally through the center? althouth it would be an additional hop and thus never the shortest path i guess. but the idea is that the glyph CENTER is the terminal node always, period. all other nodes are intermediate. then we can cut up the MST in simple paths and route them separately according to the algorithm.

    return p["set_bundle_order"]


def render_line(instance, G, p):
    # sketch:
    # for each set:
    # subgraph: keep all anchor edges that incide on elements of set plus neighbor edges that are also set edges
    # find t-spanner in that subgraph
    # add spanner to G
    # return G
    def filter_edge(setidx, setid, elements):
        def inner(u, v, k):
            is_anchor_edge = k == EdgeType.ANCHOR
            u_is_anchor = G.nodes[u]["node"] == "anchor"
            v_is_anchor = G.nodes[v]["node"] == "anchor"
            both_ends_are_anchors = u_is_anchor and v_is_anchor
            u_is_glyph = G.nodes[u]["node"] == "glyph"
            v_is_glyph = G.nodes[v]["node"] == "glyph"
            u_in_set = u_is_glyph and G.nodes[u]["glyph"] in elements
            v_in_set = v_is_glyph and G.nodes[v]["glyph"] in elements
            either_end_in_set = u_in_set or v_in_set

            is_neighbor_edge = k == EdgeType.NEIGHBOR
            has_set_edge = G.has_edge(u, v, setidx)

            return (
                is_anchor_edge and (both_ends_are_anchors or either_end_in_set)
            ) or (is_neighbor_edge and has_set_edge)

        return inner

    M = nx.MultiGraph(incoming_graph_data=G)
    M.remove_edges_from(list(M.edges))

    for setid, elements in instance["set_system"].items():
        set_idx = 1 + instance["set_ftb_order"].index(setid)
        G_ = nx.subgraph_view(G, filter_edge=filter_edge(set_idx, setid, elements))
        G_ = nx.Graph(G_)
        S = nx.spanner(G_, 10)  # TODO make configurable
        for u, v in S.edges:
            M.add_edge(u, v, set_idx)

    # optional: bundle edges more
    # idea would be to allow bundling of spanner sub-paths between two element nodes
    # so that even after bundling the modified spanner connects all set elements (but is not a t-spanner anymore necessarily)

    # next step: ordering of bundles
    # for each edge and each pair of paths that use it, find relative order of paths (before, after) acc to paper
    set_order_in_bundles = get_bundle_order(M, p)

    # next step: defining line segments to draw, this is the annoying part
    # we have all paths ordered clockwise from the previous step (i think)
    # the first thing is to define the hubs, i.e., circles with given radius at anchor nodes, regular polygons with given segment length as glyph nodes
    # then for each edge incident on a hub there's a segment along the hub's outline where the bundle can go (8 cones for sqr, 16 cones for hex)
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
    # D E F
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
        "set0": ["A", "B", "D", "E"],
        "set1": ["A", "B", "E"],
        "set2": ["G", "C"],
        "set3": ["G", "C"],
    },
    "set_bundle_order": [
        "set2",
        "set0",
        "set1",
        "set3",
    ],
    "set_ftb_order": [
        "set3",
        "set2",
        "set1",
        "set0",
    ],  # a list of set ids that defines front to back ordering (index = 0 is most front)
}

if __name__ == "__main__":
    m = 3
    n = 3
    lattice_type = "sqr"
    G = get_routing_graph(lattice_type, (m, n))
    G = embed_to_routing_graph(INSTANCE, G)
    G = render_line(INSTANCE, G, DEFAULT_PARAMS)

    pos = nx.get_node_attributes(G, "pos")
    node_color_map = {"glyph": "#882200", "anchor": "#123456"}
    node_size_map = {"glyph": 300, "anchor": 50}
    nx.draw_networkx_nodes(
        G,
        pos,
        node_shape="h" if lattice_type == "hex" else "s",
        node_size=[node_size_map[node[1]["node"]] for node in G.nodes(data=True)],
        node_color=[node_color_map[node[1]["node"]] for node in G.nodes(data=True)],
    )
    edge_color_map = {"set": "#1a6", "neighbor": "#882200", "anchor": "#123456"}
    edge_alpha_map = {"set": 0, "neighbor": 1, "anchor": 1}
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
