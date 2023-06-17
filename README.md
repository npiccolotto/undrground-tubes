# draw-ensemble-set-vis

TODO

- make pipeline and functions configurable
  - support for single set: steiner tree, path
  - DR method: MDS, UMAP, t-SNE, optimal (= QSAP)
  - gridification: Dgrid, hagrid, none (if `optimal` chosen in prior step)
  - grid size: (depends on above, hagrid can only do square grids)
  - edge penalties
- implement path support
- implement bundle ordering for trees
- add other test datasets
- finally get rid of pygame dependency, which screws the upstart time of our scripts
