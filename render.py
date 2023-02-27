import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt

DEFAULT_PARAMS = {
    "render_style": "kelpfusion",  # kelpfusion, line, envelope
    "unit_size_in_px": 50,
    "margin_size": 0.5,  # % of glyph size (which is 1 unit)
    "lane_width": 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    "lattice_type": "sqr",  # hex, tri, sqr
    "lattice_size": "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
    "glyph_positions": None,  # i guess a numpy array?
    "set_hypergraph": None,  # a networkx hypergraph
    "set_ordering": None,  # a list of set ids that defines front to back ordering (index = 0 is most front)
}


def determine_lattice_size(glyph_positions):
    # TODO implement
    return [10, 10]


def make_hex_graph(m, n):
    # start with square graph
    G = nx.grid_2d_graph(m, n)
    for pos, u in G.nodes(data=True):
        x, y = pos

        # add diagonal edges
        if y % 2 != 0 and y > 0 and x < m - 1:
            # tilting to right in odd rows
            G.add_edge((x, y), (x + 1, y - 1))
        if y % 2 == 0 and y > 0 and x > 0:
            # tilting to left in even rows
            G.add_edge((x, y), (x - 1, y - 1))

        u["x"] = x
        u["y"] = y

        # shift every odd row by 0.5 in X
        if y % 2 != 0:
            x += 0.5

        u["pos"] = (x, y)
    return G


def make_tri_graph(m, n):
    # TODO implement
    pass


def make_sqr_graph(m, n):
    G = nx.grid_2d_graph(m, n)
    for pos, u in G.nodes(data=True):
        x, y = pos
        u["pos"] = (x, y)
    return G


def convert_host_to_routing_graph(G):
    """Given a bare-bones host graph (glyph nodes and neighbor edges), extend it by anchor nodes and edges."""
    pass


def get_host_graph(lattice_type, lattice_size):
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
            raise f"unknown lattice type {lattice_type}"

    H = convert_host_to_routing_graph(G)
    return H


def embed_to_host_graph(G, p):
    """Host graph is a graph with positioned nodes. This function embeds glyph nodes and set relations into the host graph.
    Specifically, it adds i) the glyph to their respective nodes as property and ii) additional edges that correspond to the sets."""
    # TODO implement
    pass


def render_line(p):
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


def render_envelope(p):
    pass


def render_kelpfusion(p):
    # sketch:
    # for each set S_i
    # acc to meulemans 2013
    #   compute reachability graph R of S_i, where an edge (u,v) runs in the margins (along 'anchor' edges) if (u,v) is not in host_graph. use geometric length of e as weight. R contains an edge for each a,b in S_i a!=b
    #   compute shortest path graph SPG of R
    #   find faces = cycles in SPG
    #   for each face: define if it gets filled according to paper
    #   SOMEHOW render faces and edges
    pass


def render(p):
    p = DEFAULT_PARAMS | p
    lattice_size = (
        determine_lattice_size(glyph_positions)
        if p["lattice_size"] == "auto"
        else p["lattice_size"]
    )
    p["host_graph"] = get_host_graph(p["lattice_type"], lattice_size)
    match p["render_style"]:
        case "line":
            return render_line(p)
        case "envelope":
            return render_envelope(p)
        case "kelpfusion":
            return render_kelpfusion(p)
        case _:
            raise f'unknown render style {p["render_style"]}'


if __name__ == "__main__":
    G = make_hex_graph(7, 4)
    pos = nx.get_node_attributes(G, "pos")
    nx.draw_networkx_nodes(G, pos)
    nx.draw_networkx_edges(G, pos)
    plt.show()
