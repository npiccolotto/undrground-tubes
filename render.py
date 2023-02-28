import networkx as nx
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
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
        one_end_is_anchor = u_is_anchor or v_is_anchor
        u_is_ab = u in [a, b]
        v_is_ab = v in [a, b]
        either_end_is_ab = u_is_ab or v_is_ab
        return is_anchor_edge and (
            both_ends_are_anchors or (one_end_is_anchor and either_end_is_ab)
        )

    return inner


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
        set_ftb_order = instance["set_ordering"].index(setid) + 1
        if set_ftb_order <= 0:
            raise Exception(f"set {setid} unknown?")
        for a, b in combinations(elements, 2):
            ai = instance["glyph_ids"].index(a)
            bi = instance["glyph_ids"].index(b)
            logpos_a = instance["glyph_positions"][ai]
            logpos_b = instance["glyph_positions"][bi]

            phys_path = None
            phys_length = float("inf")

            if G.has_edge(logpos_a, logpos_b, EdgeType.NEIGHBOR):
                phys_path = [(logpos_a, logpos_b)]
                phys_length = calculate_path_length(G, phys_path)
            else:
                # no neighbor edge between these glyphse
                # so make a subgraph of all anchor edges
                # since anchors don't connect two glyph nodes, the shortest path will run just along anchor edges
                G_path = nx.subgraph_view(
                    G, filter_edge=filter_to_anchor_path(G, logpos_a, logpos_b)
                )
                if not nx.has_path(G_path, logpos_a, logpos_b):
                    raise Exception("you dummy")
                phys_path = nx.shortest_path(G_path, logpos_a, logpos_b)
                # the path has to be longer than two nodes for it to i) run only along anchors ii) have a glyph node at source and target
                phys_path = list(pairwise(phys_path))
                phys_length = calculate_path_length(G, phys_path)

            # add the set edge along with physical shortest path and length
            G.add_edge(
                logpos_a,
                logpos_b,
                set_ftb_order,
                edge="set",
                set_id=setid,
                length=phys_length,
                path=phys_path,
            )
    return G


def render_line(instance, G, p):
    # sketch:
    # for each set S_i:
    #   find nodes of host graph that are in S_i
    #   if these are disconnected, connect them by finding shortest paths between components along 'anchor' edges.
    #   compute MST on that G' -> this is the 'routing graph' G~ in "Edge routing with ordered bundles" (Pupyrev et al., 2016 https://doi.org/10.1016/j.comgeo.2015.10.005)
    #   proceed with path ordering step in aforementioned paper
    #   for each segment:
    #     render 'anchor' edges in MST
    #     render 'neighbor' edges in MST
    pass


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
    "set_ordering": [
        "set3",
        "set2",
        "set1",
        "set0",
    ],  # a list of set ids that defines front to back ordering (index = 0 is most front)
}

if __name__ == "__main__":
    m = 3
    n = 3
    lattice_type = "hex"
    G = get_routing_graph(lattice_type, (m, n))
    G = embed_to_routing_graph(INSTANCE, G)

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
