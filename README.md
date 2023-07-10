# draw-ensemble-set-vis

BUGS

- lines are not really going the way i expect? weird turns, long detours... probs something wrong with the weights
- line bundling doesn't work 100% the way the routing is implemented (paths can end at centers while others don't)
- and the biarcs also look like shit

TODO

- think about how line graph would look like if we were to apply bast et al 2019
  - connect port edges to centers, save port data at edge (direction we'll need for ILP later)
  - delete port-center edges
  - throw out all nodes and ports with deg = 0 (unused)
  - that should be it: this graph can be pruned, untangled, cut
  - then ILP
  - then rendering
- probably need some small test dataset again that has actual collinear lines (maybe a subset of WL tramways?)
- make pipeline and functions configurable
  - support for single set: steiner tree, path
  - DR method: MDS, UMAP, t-SNE, optimal (= QSAP)
  - gridification: Dgrid, hagrid, none (if `optimal` chosen in prior step)
  - grid size: (depends on above, hagrid can only do square grids)
  - edge penalties
- implement bundle ordering for trees
- add other test datasets
- finally get rid of pygame dependency, which screws the upstart time of our scripts
