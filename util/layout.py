import math
from collections import Counter, defaultdict
from itertools import pairwise, product

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

# import umap
# import umap.plot
# import umap.utils as utils
# import umap.aligned_umap


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


def gridify_square(P, level="auto"):
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
        Ph = solve_hagrid_optimal(curve_length, Ph)

    result = np.array([hilbert_decode(p, size) for p in Ph])

    return result


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
    # model.params.nonConvex = 2
    model.params.timeLimit = 10
    model.params.MIPGap = 0.5

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

    result = []
    for i in rla:
        for j in rlb:
            if model.getVarByName(f"el_to_pos[{i},{j}]").X > 0:
                result.append(b[j])
    return result


def layout_qsap(elements, D_EA, D_SR, m=10, n=10, weight=0.5):
    """Quadratic assignment onto grid, for now assuming constant space between cells."""

    # 1) make distance matrix D by combining D_EA and D_SR using weight
    # TODO try what ranking instead of minmax norm does
    DE = normalize(np.array(D_EA))
    DS = normalize(np.array(D_SR))
    D = (1 - weight) * DE + weight * DS

    # 2) make host distances H
    grid = list(product(range(m), range(n)))
    H = normalize(scidist.squareform(scidist.pdist(np.array(grid), "euclidean")))

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
    pos = solve_qsap_linearized(elements, D, grid, H)
    # col_ind has the indices in grid to match elements
    # pos = [grid[res.col_ind[i]] for i, el in enumerate(elements)]
    return pos


def solve_hagrid_optimal_comb(max_domain, pos):
    model = gp.Model("hagrid_comb")

    el_count = len(pos)

    # indices of elements
    rle = list(range(el_count))
    # indices of positions
    rlp = list(range(max_domain))
    # discretized positions
    # avoid having 0 factors somewhere, so shift by +1
    domain = range(1, max_domain + 1)

    # discretized position
    p = model.addVars(
        product(rle, rlp),
        vtype=GRB.BINARY,
        name="p",
    )

    #  diff p - x
    diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=-10 * max_domain, name="diff")
    abs_diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=0, name="diff_abs")

    for i, x in enumerate(pos):
        model.addConstr(
            diff[i] == 1 + x - gp.quicksum(p[i, j] * domain[j] for j in rlp)
        )
        model.addConstr(abs_diff[i] == gp.abs_(diff[i]))

    # use each position at most once
    model.addConstrs((gp.quicksum(p[i, j] for i in rle)) <= 1 for j in rlp)
    # assign all elements
    model.addConstrs((gp.quicksum(p[i, j] for j in rlp)) == 1 for i in rle)

    model.setObjective(gp.quicksum(abs_diff), sense=GRB.MINIMIZE)

    model.update()
    model.write("hagrid_comb.lp")
    model.optimize()

    for i in range(el_count):
        for j in range(max_domain):
            if model.getVarByName(f"p[{i},{j}]").X > 0:
                print(
                    pos[i],
                    "->",
                    domain[j] - 1,
                    "(",
                    model.getVarByName(f"diff[{i}]").X,
                    ")",
                )


def solve_hagrid_optimal(max_domain, pos):
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

    for i in rle:
        print(
            pos[i],
            "->",
            model.getVarByName(f"p[{i}]").X,
            "(",
            model.getVarByName(f"diff_abs[{i}]").X,
            ")",
        )

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


def layout_dr_umap(
    elements, D_EA, D_SR, m=10, n=10, weight=0.5, skip_overlap_removal=False
):
    DE = (D_EA - np.min(D_EA)) / (np.max(D_EA) - np.min(D_EA))
    DS = np.array(D_SR)
    D = (1 - weight) * DE + weight * DS

    mds = MDS(n_components=2, metric=True, random_state=2, dissimilarity="precomputed")
    H_mds = mds.fit_transform(D)

    # sns.scatterplot(x=H_mds[:,0], y=H_mds[:,1], palette='Set1')
    # plt.show()

    if not skip_overlap_removal:
        x_min = np.min(H_mds[:, 0])
        y_min = np.min(H_mds[:, 1])
        x_max = np.max(H_mds[:, 0])
        y_max = np.max(H_mds[:, 1])

        w = x_max - x_min
        h = y_max - y_min

        H_mds[:, 0] = ((H_mds[:, 0] - x_min) / w) * (n)
        H_mds[:, 1] = ((H_mds[:, 1] - y_min) / h) * (m)

        h_overlap_removed = DGrid(glyph_width=1, glyph_height=1, delta=1).fit_transform(
            H_mds
        )

        pos = []
        for i in range(len(DE)):
            pos.append((int(h_overlap_removed[i, 0]), int(h_overlap_removed[i, 1])))
    else:
        pos = []
        for i in range(len(DE)):
            pos.append((H_mds[i, 0], H_mds[i, 1]))

    return pos


def layout_dr(elements, D_EA, D_SR, m=10, n=10, weight=0.5, skip_overlap_removal=False):
    """dimensionality reduction onto grid, for now assuming constant space between cells."""

    # 1) make distance matrix D by combining D_EA and D_SR using weight
    # TODO try what ranking instead of minmax norm does
    DE = (D_EA - np.min(D_EA)) / (np.max(D_EA) - np.min(D_EA))
    DS = np.array(D_SR)
    D = (1 - weight) * DE + weight * DS

    mds = MDS(n_components=2, metric=True, random_state=2, dissimilarity="precomputed")
    H_mds = mds.fit_transform(D)

    draw_embedding(H_mds, "./embedding_raw.svg")

    if not skip_overlap_removal:
        x_min = np.min(H_mds[:, 0])
        y_min = np.min(H_mds[:, 1])
        x_max = np.max(H_mds[:, 0])
        y_max = np.max(H_mds[:, 1])

        w = x_max - x_min
        h = y_max - y_min

        H_mds[:, 0] = ((H_mds[:, 0] - x_min) / w) * (m)
        H_mds[:, 1] = ((H_mds[:, 1] - y_min) / h) * (n)

        h_overlap_removed = DGrid(glyph_width=1, glyph_height=1, delta=1).fit_transform(
            H_mds
        )

        draw_embedding(h_overlap_removed, "./embedding_gridded.svg")

        pos = []
        for i in range(len(DE)):
            pos.append((int(h_overlap_removed[i, 0]), int(h_overlap_removed[i, 1])))
    else:
        pos = []
        for i in range(len(DE)):
            pos.append((H_mds[i, 0], H_mds[i, 1]))

    return pos


def layout_dr_multiple(D_EA, D_SR, m=10, n=10, num_samples=10):
    N = len(D_EA)

    pos_mtx = []
    for weight in np.linspace(0, 1, num_samples):
        print(f'weight={weight}')
        DE = (D_EA - np.min(D_EA)) / (np.max(D_EA) - np.min(D_EA))
        DS = np.array(D_SR)
        D = (1 - weight) * DE + weight * DS

        mds = MDS(
            n_components=2, metric=True, random_state=2, dissimilarity="precomputed", normalized_stress='auto'
        )
        H_mds = mds.fit_transform(D)

        # sns.scatterplot(x=H_mds[:,0], y=H_mds[:,1], palette='Set1')
        # plt.show()

        pos = np.zeros((len(DE), 2))
        for i in range(len(DE)):
            pos[i][0] = H_mds[i, 0]
            pos[i][1] = H_mds[i, 1]

        pos_mtx.append(pos)

    transformation_mtx = []

    output_pos = []
    output_pos.append(pos_mtx[0])
    mtx_old = pos_mtx[0]

    for i in range(1, num_samples):
        mtx1, mtx2, disparity = procrustes(mtx_old, pos_mtx[i])

        mtx_old = mtx2
        output_pos.append(mtx2)

    layouts = []

    for i in range(0, num_samples):
        H_mds = output_pos[i]
        x_min = H_mds.min(axis=0)[0]
        y_min = H_mds.min(axis=0)[1]
        x_max = H_mds.max(axis=0)[0]
        y_max = H_mds.max(axis=0)[1]

        w = x_max - x_min
        h = y_max - y_min

        H_mds[:, 0] = ((H_mds[:, 0] - x_min) / w) * (m)
        H_mds[:, 1] = ((H_mds[:, 1] - y_min) / h) * (n)

        h_overlap_removed = DGrid(glyph_width=1, glyph_height=1, delta=1).fit_transform(
            H_mds
        )

        pos = []
        for i in range(len(DE)):
            pos.append((int(h_overlap_removed[i, 0]), int(h_overlap_removed[i, 1])))

        layouts.append(pos)

    for i,l in enumerate(layouts):
        draw_embedding(np.array(l), f"./embedding_gridded_{i}.svg")

    return layouts


if __name__ == "__main__":
    print(gridify_square(np.array([[0.2, 0.1], [0.1, 0.3], [4, 1.2]]), level=1))
