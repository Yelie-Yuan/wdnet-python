import wdnet
from copy import deepcopy

from wdnet import hello

hello()

print(dir(wdnet))
from wdnet import WDNet
import numpy as np
from numpy import random
from igraph import Graph


def check_conversion_igraph(g0):
    netwk0 = WDNet.from_igraph(g0)
    g1 = netwk0.to_igraph()
    netwk1 = WDNet.from_igraph(g1)
    g2 = netwk1.to_igraph()
    if not g1.isomorphic(g2):
        raise ValueError("Conversion failed; igraph")
    return True


def check_conversion_wdnet(netwk0):
    g0 = netwk0.to_igraph()
    netwk1 = WDNet.from_igraph(g0)
    g1 = netwk1.to_igraph()
    if not g1.isomorphic(g0):
        raise ValueError("Conversion failed; wdnet")
    return True


def check_conversion_both(g0, netwk0):
    g0 = deepcopy(g0)
    netwk0 = deepcopy(netwk0)
    check_conversion_igraph(g0)
    g0.es["weight"] = random.rand(g0.ecount())
    check_conversion_igraph(g0)
    g0 = Graph.Barabasi(1000, 3, directed=True)
    check_conversion_igraph(g0)
    g0.es["weight"] = random.rand(g0.ecount())
    check_conversion_igraph(g0)

    check_conversion_wdnet(netwk0)
    check_conversion_wdnet(netwk0.to_undirected())
    check_conversion_wdnet(netwk0.to_unweighted())


g0 = Graph.Barabasi(1000, 3, directed=True)
g0.es["weight"] = random.rand(g0.ecount())
netwk0 = WDNet.from_igraph(g0)
check_conversion_both(g0, netwk0)

edgelist0, edgeweight0 = netwk0.to_edgelist()
adj0 = netwk0.to_adjacency()

netwk1 = WDNet(edgelist=edgelist0, edgeweight=edgeweight0, directed=netwk0.directed)
g1 = netwk1.to_igraph()
check_conversion_both(g1, netwk1)
netwk2 = WDNet.from_adjacency(adj0, weighted=True, directed=netwk0.directed)
g2 = netwk2.to_igraph()
check_conversion_both(g2, netwk2)

sorted_rows_netwk0 = np.sort(netwk0.edgelist, axis=1)
sorted_rows_netwk1 = np.sort(netwk1.edgelist, axis=1)

# Sort the rows themselves
sorted_netwk0 = sorted_rows_netwk0[np.lexsort(sorted_rows_netwk0.T[::-1])]
sorted_netwk1 = sorted_rows_netwk1[np.lexsort(sorted_rows_netwk1.T[::-1])]

# Check if the sorted arrays are equal
np.array_equal(sorted_netwk0, sorted_netwk1)

np.allclose(
    sorted(
        [
            (*edge, netwk0.edge_attr["weight"][i])
            for i, edge in enumerate(netwk0.edgelist)
        ]
    ),
    sorted(
        [
            (*edge, netwk1.edge_attr["weight"][i])
            for i, edge in enumerate(netwk1.edgelist)
        ]
    ),
)
print(netwk0.directed, netwk1.directed, netwk2.directed)
print(
    len(netwk0.edge_attr["weight"]),
    len(netwk1.edge_attr["weight"]),
    len(netwk2.edge_attr["weight"]),
)
np.allclose(
    sorted(
        [
            (*edge, netwk0.edge_attr["weight"][i])
            for i, edge in enumerate(netwk0.edgelist)
        ]
    ),
    sorted(
        [
            (*edge, netwk2.edge_attr["weight"][i])
            for i, edge in enumerate(netwk2.edgelist)
        ]
    ),
)

np.allclose(netwk0.node_attr, netwk1.node_attr)
np.allclose(netwk0.node_attr, netwk2.node_attr)
print(netwk0.weighted, netwk1.weighted, netwk2.weighted)
netwk0.save_edgelist(file="test.csv")
netwk0.to_unweighted().save_adjacency(file="test1.csv")

print(netwk0.copy().to_undirected().directed, netwk0.directed)

print(netwk0)

print(netwk1)


import numpy as np
from numpy import random

# set random seed
random.seed(123)
random.uniform(0, 1, 10)
random.uniform(0, 1, 10)

# set random seed using python default random number generator
import random as pyrandom

pyrandom.seed(123)
pyrandom.uniform(0, 1)
pyrandom.uniform(0, 1)
pyrandom.uniform(0, 1)
