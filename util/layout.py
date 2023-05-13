import numpy as np
from scipy.optimize import quadratic_assignment
from sklearn.preprocessing import normalize
import scipy.spatial.distance as scidist
from itertools import permutations, product


def solve_gsap_linearized(a, A, b, B):
    import gurobipy as gp
    from gurobipy import GRB

    # linearize
    # problem: gets very large very quickly
    d = {}
    rla = range(len(a))
    rlb = range(len(b))

    for el1, pos1, el2, pos2 in product(rla, rlb, rla, rlb):
        d[(el1, pos1, el2, pos2)] = (B[pos1, pos2] - A[el1, el2]) ** 2
    print(d)
    mapping, diff = gp.multidict(d)

    model = gp.Model("gsap")
    # model.params.nonConvex = 2
    model.params.timeLimit = 60
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


def layout_gsap(elements, D_EA, D_SR, m=10, n=10, weight=0.5):
    """Quadratic assignment onto grid, for now assuming constant space between cells."""

    # 1) make distance matrix D by combining D_EA and D_SR using weight
    # TODO try what ranking instead of minmax norm does
    DE = normalize(np.array(D_EA))
    DS = normalize(np.array(D_SR))
    D = (1 - weight) * DE + weight * DS

    # 2) make host distances H
    grid = list(product(range(m), range(n)))
    H = normalize(scidist.squareform(scidist.pdist(np.array(grid), "cityblock")))

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
    pos = solve_gsap_linearized(elements, D, grid, H)
    # col_ind has the indices in grid to match elements
    # pos = [grid[res.col_ind[i]] for i, el in enumerate(elements)]
    return pos
