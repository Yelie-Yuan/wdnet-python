import wdnet
wdnet.hello()
print(dir(wdnet))
print(wdnet.fib(10))
from wdnet._utils import node_strength_py as s
print(s([(1, 2), (2, 3)], [0.5, 10]))
from wdnet import WDNet
a = WDNet(edgelist=[(0, 1), (1, 2)], edgeweight=[2, 10])
print(a.edgelist)
print(a.edge_attr)
print(a.node_attr)
print(a.to_undirected().node_attr)
print(a.to_unweighted().node_attr)
