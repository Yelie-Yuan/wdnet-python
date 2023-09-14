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
