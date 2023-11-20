from flask import Flask, render_template, url_for, request
import os
import json
from collections import Counter
import subprocess
import csv
import numpy as np
import pandas as pd
from util.collections import flatten

app = Flask(__name__)


def determine_type(s):
    if "low" in s:
        return "low"
    if "high" in s:
        return "high"
    if "kernel" in s:
        return "kernel"
    return "mid"


with open("data/sbss/moss.json") as f:
    file = json.load(f)

sets = file["S"]
set_sizes = [
    {
        "name": k,
        "count": v,
        "type": determine_type(k),
        "sort_order": v + (0 if determine_type(k) == "mid" else 1000),
    }
    for k, v in Counter(flatten(file["SR"])).items()
]

with open("data/sbss/kola.csv") as f:
    base_dataset = pd.read_csv(f)
with open("data/sbss/components.csv") as f:
    components = pd.read_csv(f)
with open("data/sbss/loadings.csv") as f:
    loadings = pd.read_csv(f)

coords = base_dataset[["latitude", "longitude"]]
comp_names = [(k) for k in loadings["component"]]
map_center = [coords.mean().latitude, coords.mean().longitude]


@app.route("/variables/<name>")
def variable(name):
    map_data = coords.join(base_dataset[[name]]).to_dict()
    return {
        "map_data": map_data,
    }


@app.route("/components/<name>")
def comp(name):
    map_data = coords.join(components[[name]]).to_dict("records")
    comp_loadings = loadings.loc[comp_names.index(name)].to_dict()
    sets_for_comp = file["SR"][file["E"].index(name)]
    p5 = np.percentile(components[[name]], 5)
    p25 = np.percentile(components[[name]], 25)
    p75 = np.percentile(components[[name]], 75)
    p95 = np.percentile(components[[name]], 95)

    return {
        "p5": p5,
        "p25": p25,
        "p75": p75,
        "p95": p95,
        "map_data": map_data,
        "loadings": comp_loadings,
        "sets": sets_for_comp,
    }


@app.route("/")
def index():
    rows = int(request.args.get("rows", 16))
    cols = int(request.args.get("cols", 16))
    selection = request.args.getlist("selection")
    optRouter = request.args.get("optRoute", "off") == "on"
    optConnect = request.args.get("optConnect", "off") == "on"
    treesupport = request.args.get("support", "supporttree") == "supporttree"
    weight = request.args.get("weight", 0)
    if len(selection) == 0:
        selection = ["kernel: 25"]

    subprocess.run(
        [
            "python",
            "render2.py",
            "--read-dir",
            "data/sbss",
            "--dataset",
            "moss",
            "--weight",
            str(weight),
            "--compute-metrics",
            "0",
            "--connecter",
            "opt" if optConnect else "heuristic",
            "--router",
            "opt" if optRouter else "heuristic",
            "--support-type",
            "steiner-tree" if treesupport else "path",
            "-w",
            str(cols),
            "-h",
            str(rows),
            "--sets",
            ",".join(selection),
        ]
    )
    with open("./drawing_0.svg", "r") as f:
        svg_inline = f.read()
        svg_inline = svg_inline.replace("data/sbss/img/", "/static/data/sbss/img/")
        svg_inline = svg_inline.replace("<svg ", '<svg class="undrground" ')

    return render_template(
        "index.html",
        svg_inline=svg_inline,
        rows=rows,
        cols=cols,
        sets=sets,
        selection=selection,
        weight=weight,
        optRouter=optRouter,
        optConnect=optConnect,
        map_center=map_center,
        treesupport=treesupport,
        set_size_json=json.dumps(set_sizes),
    )
