# draw-ensemble-set-vis

## Installation

* Needs everything in `requirements.txt`
* Needs [GLPK](https://www.gnu.org/software/glpk/) or alternatively don't specify solver in `bundle.py`
* Needs [LOOM](https://github.com/ad-freiburg/loom)

## TODO

- make pipeline and functions configurable
  - support for single set: steiner tree, path
  - DR method: MDS, UMAP, t-SNE, optimal (= QSAP)
  - gridification: Dgrid, hagrid, none (if `optimal` chosen in prior step)
  - grid size: (depends on above, hagrid can only do square grids)
  - edge penalties
- add other test datasets
- finally get rid of pygame dependency, which screws the upstart time of our scripts
