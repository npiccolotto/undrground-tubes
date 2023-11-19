from flask import Flask, render_template, url_for, request
import os
import json
from collections import Counter
import subprocess
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

    print(selection)
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
        treesupport=treesupport,
        set_size_json=json.dumps(set_sizes),
    )
