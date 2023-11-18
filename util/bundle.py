import networkx as nx
import numpy as np
import os
import json
import subprocess
from itertools import combinations, pairwise
from collections import defaultdict

from util.graph import incident_edges, get_longest_simple_paths, extract_support_layer
from util.geometry import get_side, get_linear_order, get_center_of_mass
from util.enums import NodeType, EdgeType
from util.config import config_vars


def find_edge_ordering_tree(G, k1, k2, u, v, initial_edge, exploring_other_way=False):
    edges = incident_edges(G, u, avoid_node=v)
    order = 0

    while len(edges) > 0:
        edges_ = defaultdict(list)
        for w, x, y in edges:
            edges_[y].append((w, x))

        # first difference: there could be multiple edges per k1/k2
        nexts_k1 = edges_[k1]
        nexts_k2 = edges_[k2]

        edges_common = set(nexts_k1).intersection(set(nexts_k2))
        edges_excl_k1 = set(nexts_k1).difference(set(nexts_k2))
        edges_excl_k2 = set(nexts_k2).difference(set(nexts_k1))

        # now... we don't care about common edges, except if they already have an ordering
        for w, x in edges_common:
            edges_wx = [
                (a, b, c, d)
                for a, b, c, d in G.edges(keys=True, data=True)
                if (a, b) == (w, x)
            ]
            for w, x, k, d in edges_wx:
                if "oeb_order" in d:
                    # found an order, return it
                    ordering = d["oeb_order"][(u, w)]
                    if ordering.index(k1) < ordering.index(k2):
                        order = -1
                    elif ordering.index(k2) < ordering.index(k1):
                        order = 1

                    # we search backwards first
                    # so if we find an order when looking back we have to reverse it
                    if not exploring_other_way:
                        order *= -1
                    return order

        # we get here when there are no common edges or they don't have an ordering
        # so before we move on, check the exclusive edges if they suggest an ordering
        if len(edges_excl_k1) > 0 or len(edges_excl_k2) > 0:
            center = v
            points = [b for (a, b) in edges_excl_k1] + [b for (a, b) in edges_excl_k2]
            categories = [k1] * len(edges_excl_k1) + [k2] * len(edges_excl_k2)
            com = get_center_of_mass(center, points, categories)

            if k1 in com and k2 in com:
                return -1 if com[k1] < com[k2] else (1 if com[k2] < com[k2] else 0)

            elif k1 in com:
                return -1

            elif k2 in com:
                return 1

        further_edges_to_check = (
            list(edges_common) + list(edges_excl_k1) + list(edges_excl_k2)
        )

        for w, x in further_edges_to_check:
            o = find_edge_ordering_tree(
                G, k1, k2, w, x, initial_edge=(u, v), exploring_other_way=False
            )
            if o != 0:
                return o

        return 0


def find_edge_ordering(G, p1, p2, u, v, initial_edge, exploring_other_way=False):
    # walk along their common edges starting from u
    # find edges from u to any node but v. G contains two paths, that were up until now
    # coincident, so either we find 0 edges (paths end) or 2 edges (paths continue).
    edges = incident_edges(G, u, avoid_node=v)
    order = 0
    while len(edges) > 0:
        edges = dict((k, v) for u, v, k in edges)

        # TODO
        # ok so it breaks here in the following two lines. because i removed the
        # segmentation into paths, there are collinear paths that go to a
        # center node, where one of the two paths ends.

        # the interesting question is what to do about it
        # we can't just post-process the terminals and replace P-C-P edges with
        # P-P because it can be ambiguous what to do...

        if p1 not in edges or p2 not in edges:
            # ok, one ends - quick fix just don't decide
            break

        next_p1 = edges[p1]
        next_p2 = edges[p2]

        if next_p1 != next_p2:
            # the two paths fork at u. path going to left takes precendence over the other.
            # left of what? left of the first edge of left of this edge?
            # let's assume left of this edge.
            side = get_side(u, next_p2, next_p1)
            if not exploring_other_way:
                # we search backwards first
                # so if we find an order when looking back we have to reverse it
                side *= -1
            return side
        else:
            # paths both use same edge
            # check if that edge has an ordering and apply it here, if so, and break. else walk further.
            e = [
                (w, x, k)
                for w, x, k in G.edges(nbunch=u, keys=True)
                if (w, x) == (u, next_p1)
            ]
            e = e[0]
            edata = G.edges[e]
            if "oeb_order" in edata:
                ordering = edata["oeb_order"][(u, next_p1)]
                if ordering.index(p1) < ordering.index(p2):
                    order = -1
                elif ordering.index(p2) < ordering.index(p1):
                    order = 1

                # we search backwards first
                # so if we find an order when looking back we have to reverse it
                if not exploring_other_way:
                    order *= -1
                return order

        edges = incident_edges(G, next_p1, avoid_node=u)
        u = next_p1
        v = u
    # if 0: both paths end at u. repeat process from v and exclude u
    if not exploring_other_way:
        w, x = initial_edge
        return find_edge_ordering(
            G, p1, p2, x, w, initial_edge, exploring_other_way=True
        )

    # in case we did all of the above and still didn't find an ordering, the two paths are coincident (=the same)
    return order


def order_bundles(instance, M):
    """M is a nx MultiGraph with all edges to draw. Returned is an ordering of sets."""
    # so... in the general case (trees) we don't really know how bundles should be ordered
    # it's actually enough, i think, if parallel bundles remain parallel in the presence of forks/merge
    # so we do what's simplest first: a generalized pupyrev approach where we basically follow edges and when they fork, we determine the order by where the schwerpunkt of the non-parallel bits is
    # i don't know how if this always leads to good bundles, but let's see

    unique_edges = list(set(M.edges()))
    for u, v in unique_edges:
        # find paths using this edge,ie. all edges (u,v,_)
        this_edge = nx.subgraph_view(M, filter_edge=lambda w, x, _: (u, v) == (w, x))
        p_uv = [k for u, v, k in this_edge.edges(keys=True)]
        if len(p_uv) < 2:
            # nothing to do here! only one path using this edge.
            k = p_uv[0]
            order = [k]
            M.edges[(u, v, k)]["oeb_order"] = {
                (u, v): order,
                (v, u): order,
            }
            continue

        # make a square matrix holding the relative orderings: O[a,b] = -1 -> a precedes b. 0 = don't matter, 1 = succeeds, else = don't know
        O = np.empty(shape=(len(p_uv), len(p_uv)))
        np.fill_diagonal(O, 0)

        # direction of (u,v): u -> v
        # (acc. to paper apparently it doesn't matter)

        # for each pair of paths
        for k1, k2 in combinations(p_uv, 2):
            M_ = nx.subgraph_view(M, filter_edge=lambda w, x, k: k in [k1, k2])
            # print(k1, k2, list(M_.edges(keys=True, data=True)))
            # filter here to path ids identified earlier so that we don't deal here with the path itself forking
            order = find_edge_ordering(M_, k1, k2, u, v, (u, v))
            # print((u, v), p1, p2, order)
            k1i = p_uv.index(k1)
            k2i = p_uv.index(k2)
            O[k1i, k2i] = order
            O[k2i, k1i] = -order

        # linearize O to something like (b,a,d,...,z)
        linear_O = [p_uv[i] for i in get_linear_order(O)]
        # save order on all multiedges between u and v
        for k in p_uv:
            M.edges[(u, v, k)]["oeb_order"] = {
                (u, v): linear_O,
                (v, u): list(reversed(linear_O)),
            }

    return M


def convert_to_line_graph(G):
    """G is a grid graph with support and physical edges. Construct a new graph without port nodes.
    Edges are only between centers."""
    G_ = nx.Graph()

    for u, v in G.edges():
        utype = G.nodes[u]["node"]
        vtype = G.nodes[v]["node"]
        uparent = G.nodes[u]["belongs_to"] if utype == NodeType.PORT else None
        vparent = G.nodes[v]["belongs_to"] if vtype == NodeType.PORT else None

        sets = G.edges[(u, v)]["sets"]

        if uparent is None and vparent is None:
            # both are centers
            # this is a mistake if we get a grid graph but fine for a line graph
            raise BaseException("support edge connecting two centers?")
        else:
            if uparent is not None and vparent is not None:
                # both are ports
                # keep their center if they belong to the same port
                if uparent == vparent:
                    G_.add_node(vparent, **G.nodes[vparent])
                else:
                    G_.add_edge(uparent, vparent, sets=sets, port1=u, port2=v)
                    G_.add_node(uparent, **G.nodes[uparent])
                    G_.add_node(vparent, **G.nodes[vparent])
            else:
                # just one is a port, which must belong to the other node (the center)
                centernode = u if uparent is None else v
                G_.add_node(centernode, **G.nodes[centernode])

    return G_


def convert_to_geojson(G):
    """Takes a line graph, converts to fake GeoJSON"""
    cx = 16.3
    cy = 48.25
    dx = 0.025
    dy = 0.025

    embed_to_fake_geo = lambda p: (cx + p[0] * dx, cy + p[1] * dy)

    stations = []
    connections = []

    for n in G.nodes():
        stations.append(
            {
                "type": "Feature",
                "properties": {
                    "id": str(n),
                    "station_label": str(n),
                    "station_id": str(n),
                },
                "geometry": {
                    "type": "Point",
                    "coordinates": list(embed_to_fake_geo(G.nodes[n]["pos"])),
                },
            }
        )

    for u, v, d in G.edges(data=True):
        connections.append(
            {
                "type": "Feature",
                "properties": {
                    "from": str(u),
                    "to": str(v),
                    "lines": list(map(lambda s: {"id": s, "label": s}, d["sets"])),
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": [
                        list(embed_to_fake_geo(G.nodes[u]["pos"])),
                        list(embed_to_fake_geo(G.nodes[v]["pos"])),
                    ],
                },
            }
        )

    geojson = {"type": "FeatureCollection", "features": stations + connections}

    return json.dumps(geojson)


def str_tuple_to_tuple(s):
    """Convert  "(1,1)" to (1,1)"""
    return tuple(map(lambda x: int(x), s[1:-1].split(",")))


def read_loom_output(output, G):
    geojson_dict = json.loads(output)

    points = [
        (str_tuple_to_tuple(f["properties"]["station_id"]), f)
        for f in geojson_dict["features"]
        if f["geometry"]["type"] == "Point"
    ]
    points_new_to_old_id_dict = dict(
        [(f["properties"]["id"], old) for old, f in points]
    )
    lines = [
        (
            (
                points_new_to_old_id_dict[f["properties"]["from"]],
                points_new_to_old_id_dict[f["properties"]["to"]],
            ),
            f,
        )
        for f in geojson_dict["features"]
        if f["geometry"]["type"] == "LineString"
    ]

    # now go through the old graph, find the new data, set the line order

    for lid, feature in lines:
        u, v = lid
        line_order = feature["properties"]["dbg_lines"].split(",")
        G.edges[(u, v)]["oeb_order"] = {
            (u, v): line_order,
            (v, u): list(reversed(line_order)),
        }
    return G

def bundle_lines_lg(instance, M):
    for layer in range(config_vars["general.numweights"].get()):
        G = extract_support_layer(M, layer)
        G_for_loom = convert_to_geojson(G)

        with open(
            f"{config_vars['general.writedir'].get()}/loom_input_{layer}.json", "w"
        ) as f:
            f.write(G_for_loom)
            loom = subprocess.run(
                [
                    "loom",
                    "-m",
                    "ilp",
                    "--ilp-solver",
                    config_vars["loom.solver"].get(),
                    "--ilp-num-threads",
                    "8",
                    "--ilp-time-limit",
                    str(config_vars["loom.timeoutsecs"].get()),
                ],
                input=G_for_loom.encode(),
                check=True,
                capture_output=True,
            )
        G = read_loom_output(loom.stdout.decode(), G)
        for u, v, d in G.edges(data=True):
            order = {(u, v): d["oeb_order"][(u, v)], (v, u): d["oeb_order"][(v, u)]}
            M.edges[(u,v, (layer, EdgeType.SUPPORT))]["oeb_order"] = order
        with open(
            f"{config_vars['general.writedir'].get()}/loom_output_{layer}.json", "w"
        ) as f:
            f.write(loom.stdout.decode())
    return M

def bundle_lines_gg(instance, M):
    # bundling
    # LOOM is quite fast in bundling and can handle trees/cycles, so...
    # 1) from the grid graph with support, make a line graph again (basically keep centers of used nodes, drop ports)
    # 2) map nodes onto a small geographic area and export as GeoJSON
    # 3) feed it to LOOM
    # 4) read result back in - bam we have an ordering

    for layer in range(config_vars["general.numweights"].get()):
        G = extract_support_layer(M, layer)
        G = convert_to_line_graph(G)
        G_for_loom = convert_to_geojson(G)

        with open(
            os.path.join(
                config_vars["general.writedir"].get(), f"loom_input_{layer}.json"
            ),
            "w",
        ) as f:
            f.write(G_for_loom)
            loom = subprocess.run(
                [
                    "loom",
                    "-m",
                    "ilp",
                    "--ilp-solver",
                    config_vars["loom.solver"].get(),
                    "--ilp-num-threads",
                    "8",
                    "--ilp-time-limit",
                    str(config_vars["loom.timeoutsecs"].get()),
                ],
                input=G_for_loom.encode(),
                check=True,
                capture_output=True,
            )
        G = read_loom_output(loom.stdout.decode(), G)
        with open(
            os.path.join(
                config_vars["general.writedir"].get(), f"loom_output_{layer}.json"
            ),
            "w",
        ) as f:
            f.write(loom.stdout.decode())

        for u, v, d in G.edges(data=True):
            w = d["port1"] if M.nodes[d["port1"]]["belongs_to"] == u else d["port2"]
            x = d["port1"] if M.nodes[d["port1"]]["belongs_to"] == v else d["port2"]
            order = {(w, x): d["oeb_order"][(u, v)], (x, w): d["oeb_order"][(v, u)]}
            # order = {(u, v): d["oeb_order"][(u, v)], (v, u): d["oeb_order"][(v, u)]}

            M.edges[(w, x, (layer, EdgeType.SUPPORT))]["oeb_order"] = order
    return M
