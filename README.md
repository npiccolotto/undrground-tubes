# draw-ensemble-set-vis

BUGS

- lines are not really going the way i expect? weird turns, long detours...
- not all nodes are connected? eg in wienerlinien_sm, stubentor should be connected to other U3 stations but it's not? but it seems also not to be a problem of the intersection groups? since neubaugasse and ottakring are in the same group... FUN
- line bundling doesn't work the way the routing is implemented (paths can end at centers while others don't)
- and probably line drawing is also broken given that we have tons of tiny ports all over the place? idk

TODO

- probably need some small test dataset again
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
