import configparser
import numpy as np
import re
from contextvars import ContextVar

config = configparser.ConfigParser()
config.read("config.ini")

config_vars = dict()
for sec in config.sections():
    for k in config[sec]:
        # maybe not ideal but whatevs
        # try guessing the appropriate type of the var
        val = config.get(sec, k)
        if val in ["True", "False"]:
            # a bool
            val = val == "True"
        elif "," in val:
            # a list
            val = val.split(",")
        elif bool(re.match("[0-9]+\.[0-9]+", val)):
            # a float
            val = float(val)
        elif bool(re.match("[0-9]+", val)):
            # int
            val = int(val)

        config_vars[f"{sec}.{k}"] = ContextVar(k, default=val)
