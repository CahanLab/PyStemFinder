# Changelog

All notable changes to PyStemFinder are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/). While the version is below 1.0,
minor releases may change the API.

## [0.2.0] - 2026-09-28

### Changed
- **Breaking:** `run_stemFinder(adata, markers, thresh=0.0, neighbors_key=None)`
  now writes the same columns as the R package: `stemFinder_raw` (higher = less
  differentiated) and `stemFinder` (`1 - raw / max(raw)`, lower = less
  differentiated, like pseudotime). Previously `stemFinder` held the raw score,
  and `stemFinder_invert` / `stemFinder_comp` were also written; both are gone.
- **Breaking:** `run_stemFinder` no longer takes `k`. The neighborhood size is
  read from the kNN graph, so it can no longer disagree with the graph (a
  mismatched `k` used to produce negative scores).
- `run_stemFinder` and `diffOmeter` are vectorized (sparse matrix products
  instead of per-cell Python loops).
- `run_stemFinder` and `diffOmeter` accept `neighbors_key` to pick a graph made
  with `sc.pp.neighbors(key_added=...)`, and ignore absent genes with a warning.
- `generate_scRNAseq_test_data` takes `random_state` (default 42, same data as
  before) and no longer reseeds numpy's global random state.
- `from PyStemFinder import *` exports only the package's functions.

### Fixed
- `run_stemFinder` crashed on sparse `.X` with scipy >= 1.14.
- `generate_scRNAseq_test_data` crashed with anndata >= 0.12 (`dtype` argument).
- `diffOmeter` crashed with `weight_by='expression'` or `'presence'` on dense
  `.X`, and on datasets with fewer than 10 cells.
- Neighbors at distance 0 (duplicate cells) were dropped from neighborhoods.
- `sf_norm_hvg_scale_pca` raised `KeyError` when the blacklist contained genes
  absent from the data.

## [0.1.0] - 2023-08-24

Initial Python port: `run_stemFinder` (Gini method), `diffOmeter`, toy data
generator, and preprocessing helpers.

[0.2.0]: https://github.com/pcahan1/PyStemFinder/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/pcahan1/PyStemFinder/releases/tag/v0.1.0
