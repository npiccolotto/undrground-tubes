from collections import defaultdict
from typing import Dict


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
    d = defaultdict(list)
    for i, element in enumerate(elements):
        sets_i = lol[i]
        for s in sets_i:
            d[s].append(element)
    return d
