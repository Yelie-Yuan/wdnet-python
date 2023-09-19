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

- [ ] Assortativity coefficient; `assortcoef()`
- [ ] Centrality measures
- [ ] Clustering coefficient
- [ ] Network rewiring
- [ ] PA network generation; `rpanet()`
  - [ ] Binary method
  - [ ] Linear method
  - [ ] Bag method
  - [ ] Bagx method
- [ ] `wdnet` and `rpacontrol` classes


# Structure (tentative)

## Directory Structure

```
wdnet/
├── __init__.py
├── wdnet/
│   ├── __init__.py
│   ├── wdnet.py
│   ├── rpacontrol.py
│   ├── rpanet.py
│   ├── dprewire.py
│   └── _internal/
│       ├── __init__.py
│       └── utils.py  # Internal Python utility functions
└── _wdnet/
    ├── some_cpp_file.cpp
    ├── another_cpp_file.cpp
    └── ...

```

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


# Files

## `rpacontrol.pyx`

```python
class RPACtrl:
    def __init__(self, some_setting, another_setting):
        # Initialization code here
    
    def scenario(self):
        # Method implementation

    def edgeweight(self):
        # Method implementation
    
    def newedge(self):
        # Method implementation
    
    def reciprocal(self):
        # Method implementation

    def preference(self):
        # Method implementation

    def __str__(self):
        return "Summary of RPA control settings..."
```

## `rpanet.pyx`

```python
def rpanet(arg1, arg2, control=RPACtrl()):
    # Generate WDNet object
```

## `rewire.pyx`

```python
def dprewire(wdnet_obj):
    # Function logic here

def dprewire_range(wdnet_obj):
    # Function logic here

def sprewire(wdnet_obj):
    # Function logic here

def sprewire_range(wdnet_obj):
    # Function logic here
```

## `utils.pyx` and `utils.cpp` are omitted for now

## Example Usage

```python
from wdnet import WDNet, RPACtrl, rpanet, rewire

# Create and use a WDNet object
net = WDNet(edgelist=[(0, 1), (1, 2), [2, 0]], edgeweight=[1, 2, 5])
net.assortcoef()
net.centrality()
net.to_undirected()
net.to_igraph()
net.to_edgelist()

# Use dprewire and dprewire_range functions
rewire.dprewire(net)
rewire.dprewire_range(net)

# Generate a PA network
rpanet(100, control=RPACtrl(), initial_network=WDNet())
```