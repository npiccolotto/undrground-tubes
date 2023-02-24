import networkx as nx


DEFAULT_PARAMS = {
    render_style: "kelpfusion",  # kelpfusion, line, envelope
    unit_size_in_px: 50,
    margin_size: 0.5,  # % of glyph size (which is 1 unit)
    lane_width: 0.1,  # width of lanes as % of glyph size (which is 1 unit)
    lattice_type: "sqr",  # hex, tri, sqr
    lattice_size: "auto",  # counts glyphs, makes square lattice that fits all. otherwise provide [widht,height] array
    glyph_positions: None,  # i guess a numpy array?
    set_hypergraph: None,  # a networkx hypergraph
    set_ordering: None,  # a list of set ids that defines front to back ordering (index = 0 is most front)
}


def determine_lattice_size(glyph_positions):
    # TODO implement
    return [10, 10]


def make_hex_graph(m, n):
    # TODO implement
    pass


def make_tri_graph(m, n):
    # TODO implement
    pass


def make_sqr_graph(m, n):
    return nx.grid_2d_graph(m, n)


def get_host_graph(lattice_type, lattice_size):
    """Returns a networkx graph. It contains two types of nodes and two types of edges.
    Nodes can be 'glyph' nodes, which correspond to spots where glyphs will be placed.
    Nodes can also be 'anchor' nodes, which corresponds to anchor points in the margins of the layout along which we trace lines and anchor polygons. We place them in the center of 'glyph' node faces.
    Edges can be 'neighbor' edges, which corresponds to direct neighbor relations of glyphs, e.g., in a sqr grid most nodes have 4 neighbors.
    Edges can also be 'anchor' edges. They connect both 'glyph' and 'anchor' nodes. Each 'anchor' glyph has 6 'anchor' edges incident: 3 for neighboring anchors and 3 for glyphs on the boundary of the face."""
    m, n = lattice_size
    match lattice_type:
        case "hex":
            return make_hex_graph(m, n)
        case "sqr":
            return make_sqr_graph(m, n)
        case "tri":
            return make_tri_graph(m, n)
        case _:
            raise f"unknown lattice type {lattice_type}"


def render_line(p):
    # sketch:
    # for each set S_i:
    #   find nodes of lattice graph that are in S_i
    #   if these are disconnected, connect them by finding shortest paths between components along 'anchor' edges.
    #   compute MST
    #   render 'anchor' edges in MST
    #   render 'neighbor' edges in MST
    # TODO: lanes, lane changes, bezier curves / rounded corners?
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
