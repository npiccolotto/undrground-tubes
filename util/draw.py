import drawsvg as svg
import numpy as np
import networkx as nx
from util.enums import NodeType, EdgeType

DARK_MODE = True  # the vs code svg previewer has dark bg, so...
STROKE_WIDTH = 4
LINE_GAP = STROKE_WIDTH * 1.5
CELL_SIZE_PX = 75
GLYPH_SIZE_PX = 0.75 * CELL_SIZE_PX
MARGINS = np.array((0.5, 0.5)) * CELL_SIZE_PX
NODE_CIRCLE_RADIUS = (
    0.3 * CELL_SIZE_PX
)  # must be smaller than 1/2 cell size or we get artefacts in line bends
DRAW_GLYPHS_OVER_LINES = False
DRAW_DEG1_MARKS = True
DRAW_DEG1_MARK_SIZE_PX = STROKE_WIDTH
DRAW_HUBS = False
DRAW_GLYPHS = True
SET_COLORS = [
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
    "#b15928",
    "#666666",
    "#a6cee3",
    "#fb9a99",
    "#b2df8a",
    "#ffff33",
]  # chosen according to set front to back order


def draw_embedding(X, path, **kwargs):
    # X is a
    r = kwargs.get("r", 0.5)
    w = kwargs.get("width", 100)
    h = kwargs.get("height", 100)
    mx = kwargs.get("margin_x", 5)
    my = kwargs.get("margin_y", 5)

    Y = (mx, my) + ((X - np.min(X)) / np.ptp(X)) * (w - mx, h - my)

    geometries = []
    for i in range(len(Y)):
        x, y = Y[i, :]
        geometries.append(
            svg.Circle(cx=x, cy=y, r=r, fill="white" if DARK_MODE else "black")
        )
    with open(path, "w") as f:
        f.write(draw_svg(geometries, w + 2 * mx, h + 2 * my))


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


def draw_svg(geometries, width, height):
    d = svg.Drawing(width, height, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()


def draw_support(instance, M, path="./support.svg"):
    geometries = []
    # project nodes
    mx, my = MARGINS
    for i in M.nodes():
        (x, y) = M.nodes[i]["pos"]
        px, py = (x * CELL_SIZE_PX + mx, -y * CELL_SIZE_PX + my)
        M.nodes[i]["pos"] = (px, py)
        if M.nodes[i]["node"] == NodeType.CENTER and M.nodes[i]["occupied"]:
            c = svg.Circle(
                cx=px,
                cy=py,
                r=NODE_CIRCLE_RADIUS,
                fill="white" if DARK_MODE else "black",
            )
            c.append_title(M.nodes[i]["label"])
            geometries.append(c)

    M_ = nx.subgraph_view(
        M,
        filter_edge=lambda u, v, k: k == EdgeType.SUPPORT
        and edge_filter_ports(M, u, v, same_centers=False),
    )
    for u, v in M_.edges():
        src = M.nodes[M.nodes[u]["belongs_to"]]["pos"]
        tgt = M.nodes[M.nodes[v]["belongs_to"]]["pos"]

        ux, uv = src
        vx, vy = tgt

        geometries.append(
            svg.Line(
                sx=ux,
                sy=uv,
                ex=vx,
                ey=vy,
                stroke_width=STROKE_WIDTH,
                stroke="white" if DARK_MODE else "black",
            )
        )

    w = instance["grid_x"] * CELL_SIZE_PX
    h = instance["grid_y"] * CELL_SIZE_PX
    with open(path, "w") as f:
        f.write(draw_svg(geometries, w, h))
