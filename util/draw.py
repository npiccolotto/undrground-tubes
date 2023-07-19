import drawsvg as svg
import numpy as np


def draw_embedding(X, path, **kwargs):
    # X is a
    r = kwargs.get("r", 0.5)
    w = kwargs.get("width", 100)
    h = kwargs.get("height", 100)
    mx = kwargs.get("margin_x", 5)
    my = kwargs.get("margin_y", 5)

    print(X)
    print((X - X.min()) / (X.max() - X.min()))

    Y = (mx, my) + ((X - np.min(X)) / np.ptp(X)) * (w - mx, h - my)

    geometries = []
    for i in range(len(Y)):
        x, y = Y[i, :]
        geometries.append(svg.Circle(cx=x, cy=y, r=r, fill="white"))
    with open(path, "w") as f:
        f.write(draw_svg(geometries, w + 2 * mx, h + 2 * my))


def draw_svg(geometries, width, height):
    d = svg.Drawing(width, height, origin=(0, 0))

    for e in geometries:
        d.append(e)

    return d.as_svg()
