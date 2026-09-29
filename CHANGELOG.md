# Changelog

All notable changes to PyStemFinder are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/). While the version is below 1.0,
minor releases may change the API.

## [0.4.0] - 2026-09-29

### Changed
- **Breaking:** E2f8 (mouse) and E2F8 (human) are S phase genes only. They
  are removed from the bundled G2M lists; the R package lists them under both
  phases, while Seurat's `cc.genes` lists them under S.
  `cell_cycle_genes(phase='both')` lists each gene once, which also covers the
  11 genes shared by the C. elegans S and G2M lists.
- **Breaking:** `run_stemFinder`, `diffOmeter`, and `gene_set_score` count a
  repeated gene once, with a warning. R counts it once per occurrence.
- Scores on the R vignette's bone marrow data therefore differ slightly from
  the published ones, which count E2f8 twice (rank correlation 0.9999, maximum
  difference 0.0065 in `stemFinder`). The R parity reference
  (`tests/data/bmmc_r_reference.h5ad`) was regenerated in R with E2f8 counted
  once, and Python still matches it to 1e-12.

## [0.3.2] - 2026-09-29

### Added
- MIT license (`LICENSE`).
- GitHub Actions workflow running the tests (Python 3.10 and 3.12) and a
  strict docs build.

### Changed
- The repository moved to https://github.com/CahanLab/PyStemFinder and is
  public, so PyStemFinder can be installed with
  `pip install git+https://github.com/CahanLab/PyStemFinder.git`.

## [0.3.1] - 2026-09-28

### Added
- Documentation: the table of contents works again, plus an API reference, a
  changelog page, and a notebook that reproduces the R vignette and its
  benchmarks (`docs/notebooks/benchmarking.ipynb`).
- `environment.yml` for a portable conda environment (replaces the exported
  personal environment files).
- README with installation, usage, differences from R, and citation.

### Changed
- `setup.py`: dependencies trimmed to what the package imports (matplotlib
  and seaborn dropped), `python_requires >= 3.9`, README as long
  description, fixed author email.
- Read the Docs installs the package so the API reference can be generated;
  Sphinx toolchain updated to current releases.

## [0.3.0] - 2026-09-28

### Added
- `run_stemFinder(method='stdev' | 'variance')`, ported from R: the summed
  sample sd / variance of marker expression across each neighborhood
  (including the cell). R uses log-normalized data for these, so
  `run_stemFinder` takes a `layer` argument.
- `compute_performance_single` and `pct_recover`, ported from R, for
  benchmarking scores against ground truth (single-cell Spearman, phenotype
  correlation, AUC; optional competitor method).
- `gene_set_score`, ported from R.
- `cell_cycle_genes(species, phase)` and `transcription_factors(species)`
  with the S/G2M and TF lists bundled with the R package (human, mouse,
  C. elegans).
- Parity tests: all three methods match R to 1e-12 on the bone marrow data,
  and so do the metrics wherever R computes them correctly (see below).

### Changed
- The gini score is summed from integer counts, so mathematically equal scores
  are bit-identical and tie in rank-based metrics. R's own single-cell
  Spearman on the vignette data (0.742814) differs from the tie-exact value
  (0.742816) because round-off splits some ties.

### Differences from R
- AUC uses every (most, least differentiated) pair. R's `auc_probability`
  receives 0/1 numeric labels, so `scores[labels]` repeats the first positive
  cell's score and the AUC compares only that one cell with the negatives
  (vignette: R 0.9725, exact 0.9662; CCAT: R 0.9299, exact 0.9530).
- The phenotype correlation defaults to Spearman, as its name says. R's
  `cor.test` call uses the default Pearson (vignette: R 0.888, Spearman
  0.865); pass `pheno_method='pearson'` to reproduce R.

## [0.2.1] - 2026-09-28

### Added
- R parity test (`tests/test_r_parity.py`) on the R vignette's Tabula Muris
  bone marrow data (3,427 cells, R's own kNN graph). `run_stemFinder` matches
  R's published scores (`bmmc_sF_results.csv`) to 1e-12. The fixture
  `tests/data/bmmc_r_reference.h5ad` is rebuilt with
  `tests/data/export_bmmc_reference.R` and `tests/data/make_bmmc_fixture.py`.

### Fixed
- A marker listed twice counts twice again, as in R (0.2.0 removed repeats).
  The R vignette's mouse S + G2M list contains E2f8 twice, so the published
  scores depend on this.

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

[0.4.0]: https://github.com/CahanLab/PyStemFinder/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/CahanLab/PyStemFinder/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/CahanLab/PyStemFinder/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/CahanLab/PyStemFinder/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/CahanLab/PyStemFinder/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/CahanLab/PyStemFinder/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/CahanLab/PyStemFinder/releases/tag/v0.1.0
