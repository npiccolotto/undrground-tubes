import math
import json
from collections import Counter, defaultdict
from itertools import pairwise, product, combinations

import os
import gurobipy as gp
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial.distance as scidist
import seaborn as sns
from gurobipy import GRB
from sklearn.manifold import MDS
from sklearn.preprocessing import MinMaxScaler, normalize
from scipy.spatial import procrustes

from util.draw import draw_embedding
from util.DGrid import DGrid
from util.config import config_vars
from util.mip import write_status, write_fake_status
from util.metrics import compute_trustworthyness_EA, compute_trustworthyness_SA, compute_stress, average_local_error


def get_bounds(P):
    xmin = np.min(P[:, 0])
    xmax = np.max(P[:, 0])
    ymin = np.min(P[:, 1])
    ymax = np.max(P[:, 1])

    return {"x0": xmin, "width": xmax - xmin, "y0": ymin, "height": ymax - ymin}


def remap_point(p, bbox1, bbox2, do_round=False):
    x, y = p
    x1, y1, w1, h1 = bbox1
    x2, y2, w2, h2 = bbox2

    fun = lambda x: round(x) if do_round else x

    return [fun((x - x1) / w1 * w2 + x2), fun((y - y1) / h1 * h2 + y2)]


def hilbert_encode(p, size):
    px, py = p
    n = 0
    s = size // 2
    while s > 0:
        rx = (px & s) > 0
        ry = (py & s) > 0
        n += (s**2) * ((3 * rx) ^ ry)
        s = s // 2
    return n


def rotate(size, px, py, rx, ry):
    if ry == 0:
        if rx == 1:
            px = size - 1 - px
            py = size - 1 - py
        return py, px
    return px, py


def hilbert_decode(p, size):
    t = p
    px, py = (0, 0)
    s = 1
    while s < size:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        px, py = rotate(s, px, py, rx, ry)
        px = px + s * rx
        py = py + s * ry
        s = s * 2
        t = t >> 2

    return px, py


def gridify_square(P, level="auto", layer=0):
    """
    Hagrid from Cutura et al. [1] using Hilbert curve and optimal collision resolution.

    R. Cutura, C. Morariu, Z. Cheng, Y. Wang, D. Weiskopf, and M. Sedlmair, “Hagrid — Gridify Scatterplots with Hilbert and Gosper Curves,” in The 14th International Symposium on Visual Information Communication and Interaction, Potsdam Germany: ACM, Sep. 2021, pp. 1–8. doi: 10.1145/3481549.3481569.
    """
    n, m = P.shape
    if m != 2:
        raise Exception("only 2D, sorry")
    if level == "auto":
        level = math.ceil(math.log2(n) / math.log2(4))
    size = 2**level
    bounds = get_bounds(P)
    original_bbox = [
        bounds["x0"],
        bounds["y0"],
        bounds["width"],
        bounds["height"],
    ]
    grid_bbox = [0, 0, size - 1, size - 1]

    Pr = np.array(
        [
            remap_point(P[i, :], original_bbox, grid_bbox, do_round=True)
            for i in range(n)
        ]
    )
    curve_length = 4**level

    Ph = [hilbert_encode(Pr[i, :], size) for i in range(n)]

    cntr = Counter(Ph)
    collision_detected = any([v > 1 for v in cntr.values()])
    if collision_detected:
        Ph = solve_hagrid_optimal(curve_length, Ph, layer=layer)

    result = np.array([hilbert_decode(p, size) for p in Ph])

    return result


def solve_qsap(elements, D, m, n, norm=2, weight=0.5, layer=0):
    """compact formulation of QSAP compared to the linearized version.

    it's still a QSAP though, so the only problems this solves is that
    the model actually builds and a (bad) solution is quickly found.

    takes a `norm` parameter wich tells whether to measure distances as
    manhattan/block (1) or euclidean (2). both variants seem to converge
    equally slowly."""

    if norm not in [1, 2]:
        raise BaseException(f"not an accepted norm: {norm}")

    model = gp.Model("qsap")
    model.params.timeLimit = config_vars["layout.ilptimeoutsecs"].get()
    model.params.MIPGap = config_vars["layout.ilpmipgap"].get()
    if norm == 2:
        model.params.NonConvex = 2

    combs = [(i, j) for i, j in combinations(range(len(elements)), 2)]

    x = model.addVars(range(len(elements)), vtype=GRB.INTEGER, lb=0, ub=m - 1)
    y = model.addVars(range(len(elements)), vtype=GRB.INTEGER, lb=0, ub=n - 1)

    model.update()

    xij_diff = model.addVars(
        combs,
        vtype=GRB.INTEGER,
        lb=-m,
    )
    xij_diff_abs = model.addVars(combs, vtype=GRB.INTEGER, ub=m)
    yij_diff = model.addVars(
        combs,
        vtype=GRB.INTEGER,
        lb=-n,
    )
    yij_diff_abs = model.addVars(combs, vtype=GRB.INTEGER, ub=n)
    if norm == 2:
        dist_ij = model.addVars(combs, vtype=GRB.CONTINUOUS, lb=0)
    discrepancy = model.addVars(combs, vtype=GRB.CONTINUOUS, lb=-m - n, ub=m + n)
    discrepancy_abs = model.addVars(combs, vtype=GRB.CONTINUOUS, ub=m + n)

    for i, j in combs:
        model.addConstr(xij_diff[(i, j)] == x[i] - x[j])
        model.addConstr(xij_diff_abs[(i, j)] == gp.abs_(xij_diff[(i, j)]))

        model.addConstr(yij_diff[(i, j)] == y[i] - y[j])
        model.addConstr(yij_diff_abs[(i, j)] == gp.abs_(yij_diff[(i, j)]))

        # avoid overlaps
        model.addConstr(xij_diff_abs[(i, j)] + yij_diff_abs[(i, j)] >= 1)

        # discrepancy = difference between distance in D and in embedding
        if norm == 2:
            model.addGenConstrNorm(
                dist_ij[(i, j)], [xij_diff[(i, j)], yij_diff[(i, j)]], norm
            )
            model.addConstr(
                discrepancy[(i, j)]
                == dist_ij[(i, j)] - (D[i, j] * math.sqrt((m + n) ** 2))
            )
        else:
            model.addConstr(
                discrepancy[(i, j)]
                == xij_diff_abs[(i, j)] + yij_diff_abs[(i, j)] - (D[i, j] * (m + n))
            )

        model.addConstr(discrepancy_abs[(i, j)] == gp.abs_(discrepancy[(i, j)]))

    model.update()

    model.setObjective(
        gp.quicksum(discrepancy_abs),
        sense=GRB.MINIMIZE,
    )

    model.optimize()

    write_status(f"layout_{layer}", model)

    return np.array([(x[i].X, y[i].X) for i in range(len(elements))])


def solve_qsap_linearized(a, A, b, B):
    # linearize
    # problem: gets very large very quickly
    d = {}
    rla = range(len(a))
    rlb = range(len(b))

    for el1, pos1, el2, pos2 in product(rla, rlb, rla, rlb):
        d[(el1, pos1, el2, pos2)] = (B[pos1, pos2] - A[el1, el2]) ** 2
    mapping, diff = gp.multidict(d)

    model = gp.Model("qsap")
    model.params.timeLimit = config_vars["layout.ilptimeoutsecs"].get()
    model.params.MIPGap = config_vars["layout.ilpmipgap"].get()

    el_to_pos = model.addVars(
        list(product(rla, rlb)),
        vtype=GRB.BINARY,
        name="el_to_pos",
    )

    # at most one assignment to any location
    model.addConstrs((gp.quicksum(el_to_pos[el, p] for el in rla)) <= 1 for p in rlb)
    # assign all elements
    model.addConstrs((gp.quicksum(el_to_pos[el, p] for p in rlb)) == 1 for el in rla)

    model.setObjective(
        (
            gp.quicksum(
                diff[el1, pos1, el2, pos2] * el_to_pos[el1, pos1] * el_to_pos[el2, pos2]
                for el1, pos1, el2, pos2 in mapping
            )
        ),
        sense=GRB.MINIMIZE,
    )

    model.update()
    model.optimize()

    write_status("layout", model)

    result = []
    for i in rla:
        for j in rlb:
            if model.getVarByName(f"el_to_pos[{i},{j}]").X > 0:
                result.append(b[j])
    return np.array(result)


def layout_qsap(elements, D, m=10, n=10, weight=0.5, layer=0):
    """Quadratic assignment onto grid, for now assuming constant space between cells."""

    return solve_qsap(elements, D, m, n, norm=2, weight=weight, layer=layer)
    # 2) make host distances H

    # add dummy elements - qap does a bijective mapping
    # num_cells = len(grid)
    # num_els = len(elements)
    # num_diff = num_cells - num_els
    # if num_diff > 0:
    #    print("adding dummy elements")

    # 0 = "no flow from here to anywhere" = dummies are on the outside, elements in circular shape in the middle
    # 1 = "max flow from here to everywhere" = heuristic gets slowww
    # 0.5 = ?? = heuristic still quick, elements not in single cluster, but not a good layout nonetheless
    # 0.1 = "a little flow everywhere" = circular cluster
    #    dummy_val = 0.01
    #    D = np.c_[D, np.full((num_els, num_diff), dummy_val)]
    #    D = np.r_[D, np.full((num_diff, num_els + num_diff), dummy_val)]

    # 3) put D and H into QAP
    # TODO always puts actual elements in circular shape in center of the grid
    # res = quadratic_assignment(D, H, method="faq")  # '2opt' is too slow
    # grid = list(product(range(m), range(n)))
    # H = normalize(scidist.squareform(scidist.pdist(np.array(grid), "euclidean")))
    # pos = solve_qsap_linearized(elements, D, grid, H)

    # col_ind has the indices in grid to match elements
    # pos = [grid[res.col_ind[i]] for i, el in enumerate(elements)]
    return pos


def solve_hagrid_optimal(max_domain, pos, layer=0):
    model = gp.Model("hagrid")

    rle = list(range(len(pos)))

    # discretized position
    p = model.addVars(
        rle,
        vtype=GRB.INTEGER,
        lb=0,
        ub=max_domain - 1,
        name="p",
    )
    model._p = p

    #  diff p - x
    diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=-10 * max_domain, name="diff")
    abs_diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=0, name="diff_abs")

    for i, x in enumerate(pos):
        model.addConstr(diff[i] == x - p[i])
        model.addConstr(abs_diff[i] == gp.abs_(diff[i]))

    model.Params.LazyConstraints = 1

    def addViolatedConstraints(m, pValues):
        d = defaultdict(list)
        for i, p in pValues.items():
            d[p].append(i)
        for p, indices in d.items():
            if len(indices) > 1:
                for i in indices:
                    for j in indices:
                        if j <= i:
                            continue
                        diffvar = m.getVarByName(f"diff_abs_p{i}_p{j}")
                        m.cbLazy(diffvar >= 1)

    def callback(m, where):
        if where == GRB.Callback.MIPSOL:
            # get solution values for variables x
            pValues = m.cbGetSolution(m._p)
            addViolatedConstraints(m, pValues)
        # check fractional solutions to find violated CECs/DCCs to strengthen the bound
        elif (
            where == GRB.Callback.MIPNODE
            and m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
        ):
            # get solution values for variables x
            pValues = m.cbGetNodeRel(m._p)
            addViolatedConstraints(m, pValues)

    # no double assignment
    # inequality constraints are lazily added, but need the variables already
    for i in rle:
        for j in rle:
            if j <= i:
                continue
            p_between = model.addVar(
                vtype=GRB.INTEGER, lb=-10 * max_domain, name=f"diff_p{i}_p{j}"
            )
            p_between_abs = model.addVar(
                vtype=GRB.INTEGER, lb=0, name=f"diff_abs_p{i}_p{j}"
            )
            model.addConstr(p_between == (p[i] - p[j]))
            model.addConstr(p_between_abs == gp.abs_(p_between))

    model.setObjective(gp.quicksum(abs_diff), sense=GRB.MINIMIZE)

    model.update()
    model.write("hagrid.lp")
    model.optimize(callback)

    write_status(f"overlapremoval_{layer}", model)

    # for i in rle:
    #    print(
    #        pos[i],
    #        "->",
    #        model.getVarByName(f"p[{i}]").X,
    #        "(",
    #        model.getVarByName(f"diff_abs[{i}]").X,
    #        ")",
    #    )

    return [int(model.getVarByName(f"p[{i}]").X) for i in rle]


def compute_error(L1, L2):
    error = 0
    for i in range(len(L1)):
        for j in range(i + 1, len(L2)):
            x_d = L1[i][0] - L2[j][0]
            y_d = L1[i][1] - L2[j][1]

            error += x_d**2 + y_d**2

    return error


def naive_matching(L1, L2):
    rot = 0
    scale = 0

    error = 10000000000000000000000000
    for _ in range(1000):
        rot = np.random.uniform(0.0, 2 * np.pi)
        scale = np.random.uniform(0.1, 2)

        L2_new = []

        for l in L2:
            continue

        new_error = compute_error(L1, L2)

        if new_error < error:
            error = new_error
            best_rot = rot
            best_scale = scale

    return rot, scale


def layout_mds(D, m=10, n=10, weight=0.5, layer=0):
    write_fake_status(f"layout_{layer}")
    return MDS(
        n_components=2,
        metric=True,
        random_state=2,
        dissimilarity="precomputed",
        normalized_stress="auto",
    ).fit_transform(D)


def layout_single(elements, D_EA, D_SR, m=10, n=10, weight=0.5, layer=0):
    layouter = config_vars["layout.layouter"].get()

    if layouter == "auto":
        strat = config_vars["general.strategy"].get()
        layouter = "mds" if strat == "heuristic" else "qsap"

    DE = (D_EA - np.min(D_EA)) / (np.max(D_EA) - np.min(D_EA))
    DS = np.array(D_SR)
    D = (1 - weight) * DE + weight * DS

    match layouter:
        case "mds":
            return layout_mds(D, m=m, n=n, weight=weight, layer=layer)
        case "qsap":
            return layout_qsap(elements, D, m=m, n=n, weight=weight, layer=layer)


def align_layouts(layouts):
    output_pos = []
    output_pos.append(layouts[0])
    mtx_old = layouts[0]

    for i in range(1, len(layouts)):
        mtx1, mtx2, disparity = procrustes(mtx_old, layouts[i])

        mtx_old = mtx2
        output_pos.append(mtx2)

    return output_pos


def remove_overlaps(layout, m=10, n=10, layer=0):
    x_min = layout.min(axis=0)[0]
    y_min = layout.min(axis=0)[1]
    x_max = layout.max(axis=0)[0]
    y_max = layout.max(axis=0)[1]

    w = x_max - x_min
    h = y_max - y_min

    layout[:, 0] = ((layout[:, 0] - x_min) / w) * (m)
    layout[:, 1] = ((layout[:, 1] - y_min) / h) * (n)

    remover = config_vars["layout.overlapremover"].get()

    is_square = m == n
    power2_exp = math.log2(m)
    is_power_of_2_square = is_square and power2_exp == int(power2_exp)
    if remover == "hagrid" and not is_power_of_2_square:
        print(
            "WARN: Hagrid used as overlap remover, but grid size is not a power of 2 square. falling back to DGrid."
        )
        remover = "dgrid"

    match remover:
        case "hagrid":
            return gridify_square(layout, int(power2_exp), layer=layer).astype(int)
        case "dgrid":
            write_fake_status(f"overlapremoval_{layer}")
            return (
                DGrid(glyph_width=1, glyph_height=1, delta=1)
                .fit_transform(layout)
                .astype(int)
            )


def pad_layout(layout, pad):
    return pad + layout


def layout(elements, D_EA, D_SR, m=10, n=10, pad=1):
    num_weights = config_vars["general.numweights"].get()
    weight = config_vars["general.weight"].get()
    layouts = [
        layout_single(elements, D_EA, D_SR, m=m, n=n, weight=w, layer=l)
        for l, w in enumerate(np.linspace(weight, 1, num_weights))
    ]

    layouts = align_layouts(layouts)
    layouts = [
        remove_overlaps(layout, m=m, n=n, layer=l) for l, layout in enumerate(layouts)
    ]

    layouts = [pad_layout(layout, pad) for layout in layouts]

    layouts = list(map(lambda layout: list(map(tuple, layout)), layouts))

    write_metrics(elements, D_EA, D_SR, layouts)

    return layouts


def write_metrics(elements, D_EA, D_SR, layouts):
    if config_vars["general.computemetrics"].get():
        embedding_metrics = []

        for i, l in enumerate(layouts):
            embedding_metrics.append({
                "EA: local": average_local_error(elements, D_EA, l),
                "SR: local": average_local_error(elements, D_SR, l),
                "EA: stress": compute_stress(elements, D_EA, l),
                "SR: stress": compute_stress(elements, D_SR, l),
                "EA: (M1, M2) [k=3]": compute_trustworthyness_EA(elements, D_EA, l),
                "EA: (M1, M2) [k=5]": compute_trustworthyness_EA(elements, D_EA, l, k=5),
                "EA: (M1, M2) [k=7]": compute_trustworthyness_EA(elements, D_EA, l, k=7),
                "EA: (M1, M2) [k=10%]": compute_trustworthyness_EA(elements, D_EA, l, k=(int)(np.ceil(len(elements)/10))),
                "SR: (M1, M2) [k=3]": compute_trustworthyness_SA(elements, D_SR, l, k=3),
                "SR: (M1, M2) [k=5]": compute_trustworthyness_SA(elements, D_SR, l, k=5),
                "SR: (M1, M2) [k=7]": compute_trustworthyness_SA(elements, D_SR, l, k=7),
                "SR: (M1, M2) [k=10%]": compute_trustworthyness_SA(elements, D_SR, l, k=(int)(np.ceil(len(elements)/10)))
            })


        writedir = config_vars["general.writedir"].get()
        with open(f"{writedir}/metrics_embedding.json", "w") as f:
            json.dump(embedding_metrics, f)
        with open(f"{writedir}/embedding.json", "w") as f:
            # https://stackoverflow.com/a/65151218/490524
            def np_encoder(object):
                if isinstance(object, np.generic):
                    return object.item()

            json.dump(layouts, f, default=np_encoder)
