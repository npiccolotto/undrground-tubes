# UnDRground Tubes

## Installation

Requires Python 3.10.

* Needs `pyenv`
* Needs everything in `requirements.txt`. Install with `python -m pip install -r requirements.txt`
* Needs [Gurobi](https://www.gurobi.com/)
* Needs [LOOM](https://github.com/ad-freiburg/loom) and LOOM needs to be in your PATH variable.

## Usage

`python render2.py`

The result will be written as `drawing_{layer}.svg`, usually `drawing_0.svg`. For available options, see `python render2.py --help`. The defaults and additional configuration options can be changed in the `config.ini` file.

## Datasets

The required data is expected to be provided in a JSON file. The schema of the JSON file is available in `schema.json` as [JSON Schema](https://json-schema.org/). The **glyphs** are expected to be images and the path to them provided in that JSON file. See `data/imdb/imdb10` for an example.

## Folders in this Repository

* `data` holds several example datasets, e.g., IMDb movies, the Viennese tram and metro network, a synthetic example, R's `mtcars` dataset, and the SBSS example used in the paper. Note that some example datasets require additional steps before they work. Most notably regarding their glyphs, which we did not want to add to the git repository.
* `designspace-result` holds the results from our computational experiments, both raw data and a Jupyter notebook (`merge.ipynb`) interpreting them.
* `designspace-test` holds the example datasets generated for the computational experiments.
* `eval` has some initial and currently irrelevant experiments regarding a controlled user study.
* `figures` has some figures used in the paper.
* `static`, `templates` are required by the Web frontend described in Section 5 of the paper. Use `frontend.py` to start it, but it won't have glyphs unless you generate them with the R file.
* `tests`, `util` and `render2.py` are the implementation of UnDRground Tubes.

# License

GPL v3 https://choosealicense.com/licenses/gpl-3.0/
