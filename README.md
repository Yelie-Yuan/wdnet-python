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
- [ ] Network rewiring
- [x] Assortativity coefficient range `dprewire_range()`
- [ ] PA network generation; `rpanet()`
  - [ ] Binary method
  - [ ] Linear method
  - [ ] Bag method
  - [ ] Bagx method
- [x] `WDNet` class
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
  - to_adj()
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
from igraph import Graph
from wdnet import WDNet, rewire

# ER network
g = Graph.Erdos_Renyi(n=100, m=1000, directed=True)
netwk = WDNet.from_igraph(g)
netwk.assortcoef()
netwk.to_undirected()
netwk.to_igraph()

# Use dprewire and dprewire_range functions
# rewire.dprewire(net) # under development
rewire.dprewire_range(netwk, which_range="inin")
rewire.dprewire_range(netwk, which_range="inin", target_assortcoef={"outout": 0.5, "outin": 0.5})

# Generate a PA network
# rpanet(100, control=RPACtrl(), initial_network=WDNet()) # under development
```