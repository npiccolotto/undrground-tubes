from collections import defaultdict
from typing import Dict
from itertools import combinations, chain
import scipy.spatial.distance as scidist
import numpy as np
from util.config import config_vars


def flatten(list):
    return [item for sublist in list for item in sublist]


def merge_alternating(xs, ys=None):
    if ys is None and len(xs) == 2:
        xs, ys = xs
    # https://stackoverflow.com/a/3678938/490524
    result = [None] * (len(xs) + len(ys))
    result[::2] = xs
    result[1::2] = ys
    return result


def invert_dict_of_lists(d: Dict):
    new_d = defaultdict(list)
    for k, v in d.items():
        for e in v:
            new_d[e].append(k)

    return new_d


def invert_dict(d: Dict):
    return dict(zip(d.values(), d.keys()))


def merge_keys_with_same_values(d: Dict):
    """
    {
        "set0": "A",
        "set1": "A",
        "set2": "B",
        "set3": "C",
    }
    ->
    {
        A: set0,set1
        B: set2,
        C: set3
    }
    """
    new_d = defaultdict(list)
    for k, v in d.items():
        new_d[v].append(k)
    return new_d


def get_elements_in_same_lists(d: Dict):
    """
    {
        "set0": ["A", "B", "H", "D", "E"],
        "set1": ["A", "I"],
        "set2": ["G", "C", "D", "E"],
        "set3": ["A", "F", "D", "E"],
        "set4": ["A", "I", "D", "E"],
    }
    ->
    {
        (D,E): (set0,set2,set3,set4),
        (I):(set1,set4),
        (A):(set0,set1,set3,set4),
        ...
    }
    """
    inv = invert_dict_of_lists(d)
    for k in inv.keys():
        inv[k] = tuple(sorted(inv[k]))
    inv = merge_keys_with_same_values(inv)
    for k in inv.keys():
        inv[k] = tuple(sorted(inv[k]))
    inv = invert_dict(inv)
    return inv


def list_of_lists_to_set_system_dict(elements, lol):
    """Returns a dict with sets as keys and their elements as list value."""
    d = defaultdict(list)
    for i, element in enumerate(elements):
        sets_i = lol[i]
        for s in sets_i:
            d[s].append(element)
    return d


# def group_by_intersection_group(set_system_dict):
#    intersection_groups = get_intersection_groups(set_system_dict)
#    sorted_intersection_groups = sorted(
#        zip(intersection_groups.keys(), intersection_groups.values()),
#        key=lambda x: len(x[1]),
#        reverse=True,
#    )
#    return sorted_intersection_groups


def group_by_set(set_system_dict):
    sets, elements = zip(*set_system_dict.items())

    sets = list(
        map(lambda s: s if isinstance(s, list) or isinstance(s, set) else [s], sets)
    )

    return sorted(zip(elements, sets), key=lambda s: len(s[0]), reverse=True)


def invert_list(l):
    """returns a dict that holds the index of each list element"""
    d = {}
    for i, e in enumerate(l):
        d[e] = i
    return d


def set_contains(a, b):
    """returns true if set b is contained in set a"""
    return a.intersection(b) == b


def group_by_intersection_group(d: Dict):
    sets, _ = zip(*d.items())
    igroups = list(
        chain.from_iterable(combinations(sets, l) for l in range(1, len(sets) + 1))
    )

    result = []
    for igroup in igroups:
        igroup_elements = set.intersection(*[set(d[s]) for s in igroup])
        # currently we can't display sets with 1 element in them
        if len(igroup_elements) > 1:
            result.append((tuple(sorted(igroup_elements)), tuple(sorted(igroup))))

    return result


def select_sets(instance, selection):
    sets = list(filter(lambda s: s in selection, instance["sets"]))
    set_system = dict([(s, v) for s, v in instance["set_system"].items() if s in sets])
    SM = np.zeros((len(instance["elements"]), len(sets)))
    for s, els in set_system.items():
        for el in els:
            SM[instance["elements"].index(el), sets.index(s)] = 1

    D_SR = scidist.squareform(scidist.pdist(SM, "jaccard"))
    return instance | {
        "sets": sets,
        "set_system": set_system,
        "D_SR": D_SR,
        "set_ftb_order": sets,
        "set_colors": dict(
            zip(
                sets,
                instance["SC"]
                if "SC" in instance
                else config_vars["draw.setcolors"].get(),
            )
        ),
    }


if __name__ == "__main__":
    print(
        group_by_intersection_group(
            {"horror": ["a", "b"], "comedy": ["a"], "drama": ["a", "b", "c"]}
        )
    )
