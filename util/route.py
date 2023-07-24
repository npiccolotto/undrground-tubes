import networkx as nx
import gurobipy as gp
from gurobipy import GRB

from itertools import product

def route_multilayer(instance, G, element_set_partition,  support_type='tree'):
    num_layers = instance.get('num_layers',2)

    model = gp.Model('multilayer-route')
    #model.params.timeLimit = 10
    model.params.MIPGap = 0.5

    model.addVars(list(product(range(num_layers), G.edges())))

