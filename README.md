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

## `__init__.py`

```python
from .wdnet import WDNet
from .rpacontrol import RPACtrl
from .rpanet import rpanet
from .dprewire import dprewire, dprewire_range
```

## `wdnet.py`

```python
class WDNet:
    def __init__(self, directed=False, weighted=False, edgeweight=None, edgelist=None, edge_attr=None, node_attr=None):
        self.directed = directed
        self.weighted = weighted
        self.edgeweight = edgeweight if edgeweight is not None else []
        self.edgelist = edgelist if edgelist is not None else []
        self.edge_attr = edge_attr if edge_attr is not None else {}
        self.node_attr = node_attr if node_attr is not None else {}

    def assortcoef(self):
        # Method implementation

    def centrality(self):
        # Method implementation

    def save_edgelist(self, file_path):
        # Save the edgelist to a file or return as an array

    def to_undirected(self):
        # Convert to an undirected network
    
    def to_unweighted(self):
        # Convert to an unweighted network
    
    def to_igraph(self):
        # Convert to an igraph.Graph object

    def from_igraph(self, igraph_obj):
        # Convert from an igraph.Graph object
      
    def to_edgelist(self):
        # Convert to an edgelist
    
    def to_adj(self, file=None):
        adj = self.convert_to_adj()

        if file is None:
            print(adj)
        else:
            # Save to file
            np.save(file, adj, delimiter=",")

    def __str__(self):
        return f"RPACtrl Object: {self.some_setting}, {self.another_setting}"
```

## `rpacontrol.py`

```python
class RPACtrl:
    def __init__(self, some_setting, another_setting):
        # Initialization code here
    def scenario(self):
        # Method implementation

    def preference(self):
        # Method implementation

    def __str__(self):
        return "Summary of RPA control settings..."
```

## `rpanet.py`

```python
def rpanet(arg1, arg2, control=RPACtrl()):
    # Generate WDNet object
```

## `dprewire.py`

```python
def dprewire(wdnet_obj):
    # Function logic here

def dprewire_range(wdnet_obj):
    # Function logic here
```

## Example Usage

```python
from wdnet import WDNet, RPACtrl, rpanet, dprewire, dprewire_range

# Create and use a WDNet object
net = WDNet()

# Use dprewire and dprewire_range functions
dprewire(net)
dprewire_range(net)
```