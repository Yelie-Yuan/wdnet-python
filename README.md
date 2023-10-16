# wdnet

This is a Python version of the R package
[wdnet](https://cran.r-project.org/web/packages/wdnet/index.html).


Implementations of network analysis including:
+ assortativity coefficient of weighted and directed networks, Yuan,
  Yan and Zhang (2021) <doi:10.1093/comnet/cnab017>;
+ centrality measures for weighted and directed networks, Opsahl,
  Agneessens and Skvoretz (2010) <doi:10.1016/j.socnet.2010.03.006>,
  Zhang, Wang and Yan (2022) <doi:10.1016/j.physa.2021.126438>;
+ clustering coefficient of weighted and directed networks, Fagiolo
  (2007) <doi:10.1103/PhysRevE.76.026107> and Clemente and Grassi
  (2018) <doi:10.1016/j.chaos.2017.12.007>;
+ rewiring networks with given assortativity coefficients, Wang, Yan,
  Yuan and Zhang (2022) <doi:10.1007/s11222-022-10161-8>;
+ preferential attachment network generation, Yuan, Wang, Yan and
  Zhang (2023) <doi:10.6339/23-JDS1110>.


# Development

This package is under development.

- [x] Assortativity coefficient; `assortcoef()`
- [ ] Centrality measures
- [ ] Clustering coefficient
- [x] Network rewiring
- [x] Assortativity coefficient range `dprewire_range()`
- [ ] PA network generation; `rpanet()`
  - [ ] Binary method
  - [ ] Linear method
  - [ ] Bag method
  - [ ] Bagx method
- [x] `WDNet` class
  - [x] `edgelist`, `numpy` 2-D array with data type `numpy.int_`
  - [x] `edgeweight`, `numpy` 1-D array with data type `numpy.float64`
    if `weighted` and `numpy.int_` if `unweighted`.
- [ ] `rpacontrol` class

# Structure (tentative)

## Classes

+ WDNet
  - assortcoef()
  - centrality()
  - to_undirected()
  - to_unweighted()
  - to_igraph()
  - from_igraph()
  - from_edgelist()
  - from_adj()
  - to_edgelist()
  - save_edgelist()
  - to_adjacency()
  - save_adjacency()
+ RPACtrl
  - scenario()
  - edgeweight()
  - newedge()
  - reciprocal()
  - preference()

## Modules

+ rewire
  - dprewire()
  - dprewire_range()
  - sprewire()
  - sprewire_range()

## Functions

+ rpanet


## Example Usage

```python
from numpy import random as nprandom
from igraph import Graph
from wdnet import WDNet, rewire
import matplotlib.pyplot as plt
import random
random.seed(123) # for igraph

# ER network
g = Graph.Erdos_Renyi(n=200, m=5000, directed=True)
netwk = WDNet.from_igraph(g)
netwk.assortcoef()
netwk.to_undirected()
netwk.to_igraph()

# Use dprewire and dprewire_range functions
# rewire.dprewire(net) # under development
rewire.dprewire_range(netwk, which_range="inin")

nprandom.seed(123) # for wdnet
_, assort_trace, _ = rewire.dprewire(
    netwk, target_assortcoef={"outout": 0.2, "outin": 0.2, "inout": -0.2, "inin": -0.2}
)

_, axs = plt.subplots(2, 2, figsize=(10, 10))
axs[0, 0].plot(assort_trace["outout"])
axs[0, 0].set_title("outout")
axs[0, 1].plot(assort_trace["outin"])
axs[0, 1].set_title("outin")
axs[1, 0].plot(assort_trace["inout"])
axs[1, 0].set_title("inout")
axs[1, 1].plot(assort_trace["inin"])
axs[1, 1].set_title("inin")
plt.show()
# Generate a PA network
# rpanet(100, control=RPACtrl(), initial_network=WDNet()) # under development
```
