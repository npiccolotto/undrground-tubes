import numpy as np
from scipy.optimize import quadratic_assignment
from sklearn.preprocessing import normalize
import scipy.spatial.distance as scidist
from itertools import permutations, product
import gurobipy as gp
from gurobipy import GRB


def gridify(P, width, height):
    # TODO translate code in these files:
    #   https://github.com/saehm/hagrid/blob/master/src/gridify.js
    #   https://github.com/saehm/hagrid/blob/master/src/hilbert.js
    # using `solve_hagrid_optimal` functions below to avoid ordering issues and non-optimal assignments like in Hagrid paper

    raise Exception("not implemented")


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

    #  diff p - x
    diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=-10 * max_domain, name="diff")
    abs_diff = model.addVars(rle, vtype=GRB.CONTINUOUS, lb=0, name="diff_abs")

    for i, x in enumerate(pos):
        model.addConstr(diff[i] == x - p[i])
        model.addConstr(abs_diff[i] == gp.abs_(diff[i]))

    # no double assignment
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
            model.addConstr(p_between_abs >= 1)

    model.setObjective(gp.quicksum(abs_diff), sense=GRB.MINIMIZE)

    model.update()
    model.write("hagrid.lp")
    model.optimize()

    for i in rle:
        print(
            pos[i],
            "->",
            model.getVarByName(f"p[{i}]").X,
            "(",
            model.getVarByName(f"diff_abs[{i}]").X,
            ")",
        )


if __name__ == "__main__":
    print(solve_hagrid_optimal(4, [2.5, 1.05, 0.9]))
