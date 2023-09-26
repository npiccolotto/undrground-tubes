from gurobipy import GRB
import os
import json
from util.config import config_vars


def get_mip_status(model):
    is_mip = bool(model.IsMIP)

    if not is_mip:
        return {"mip": False}

    has_mipgap = model.MIPGap != GRB.INFINITY
    return {
        "mip": True,
        "has_mipgap": has_mipgap,
        "mipgap": model.MIPGap if has_mipgap else None,
        "runtime_ms": model.Runtime * 1000,
        "status": model.Status,  # https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html
    }


def write_json(name, data):
    with open(
        os.path.join(config_vars["general.writedir"].get(), f"ilp_status_{name}.json"),
        "w",
    ) as f:
        json.dump(data, f)


def write_status(name, model):
    write_json(name, get_mip_status(model))


def write_fake_status(name):
    write_json(name, {"mip": False})
