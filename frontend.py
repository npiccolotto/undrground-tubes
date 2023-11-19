from flask import Flask, render_template, url_for, request
import os
import subprocess

app = Flask(__name__)


@app.route("/")
def hello_world():
    rows = int(request.args.get("rows", 16))
    cols = int(request.args.get("cols", 16))
    subprocess.run(
        [
            "python",
            "render2.py",
            "--read-dir",
            "data/sbss",
            "--dataset",
            "moss",
            "--compute-metrics",
            "0",
            "-w",
            str(cols),
            "-h",
            str(rows),
            "--sets",
            "kernel: 25",
        ]
    )
    with open("./drawing_0.svg", "r") as f:
        svg_inline = f.read()
        svg_inline = svg_inline.replace("data/sbss/img/", "/static/data/sbss/img/")
    return render_template("index.html", svg_inline=svg_inline, rows=rows, cols=cols)
