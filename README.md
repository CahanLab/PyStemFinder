# PyStemFinder

[![tests](https://github.com/CahanLab/PyStemFinder/actions/workflows/tests.yml/badge.svg)](https://github.com/CahanLab/PyStemFinder/actions/workflows/tests.yml)
[![docs](https://readthedocs.org/projects/pystemfinder/badge/?version=latest)](https://pystemfinder.readthedocs.io/en/latest/)

PyStemFinder estimates how far single cells have progressed through differentiation from
scRNA-seq data. Less differentiated cells vary more in their expression of cell cycle genes than
their neighbors do, and stemFinder scores that heterogeneity within each cell's k-nearest-neighbor
neighborhood.

PyStemFinder is a Python port of the R package [stemFinder](https://github.com/CahanLab/stemfinder)
that works on [AnnData](https://anndata.readthedocs.io) objects and
[scanpy](https://scanpy.readthedocs.io) neighbor graphs. Given the same inputs, its scores match R's
to within floating point precision (checked in `tests/test_r_parity.py` against the R vignette's
bone marrow data).

## Installation

```bash
pip install git+https://github.com/CahanLab/PyStemFinder.git
```

To run the tutorials or the tests, clone the repository and install into the conda environment
defined in `environment.yml`:

```bash
git clone https://github.com/CahanLab/PyStemFinder.git
cd PyStemFinder
conda env create -f environment.yml
conda activate pystemfinder
pip install -e .
```

## Usage

```python
import numpy as np
import scanpy as sc
import pystemfinder as psf

# adata: raw counts (cells x genes)
markers = psf.cell_cycle_genes("mouse")  # S and G2M phase genes; also "human", "celegans"

# normalize, log-transform, scale, and run PCA, leaving the markers out of the highly variable genes
psf.recipe_stemfinder(adata, exclude=markers, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=int(round(np.sqrt(adata.n_obs))), n_pcs=32)

psf.stemfinder(adata, markers)
# adata.obs["stemfinder_raw"]: higher = less differentiated
# adata.obs["stemfinder"]:     1 - raw / max(raw), oriented like pseudotime
```

With ground truth differentiation stages in `adata.obs["Ground_truth"]` and cell types in
`adata.obs["Phenotype"]`, `psf.compute_performance(adata)` and `psf.pct_recover(adata)` benchmark the
scores as in the R package.

The [documentation](https://pystemfinder.readthedocs.io) includes notebooks (in `docs/notebooks`) that
walk through a toy dataset (`quickstart.ipynb`) and reproduce the R vignette (`benchmarking.ipynb`).

## Coming from the R package

| R stemFinder | pystemfinder |
| --- | --- |
| `run_stemFinder(adata, nn, k, thresh, markers, method)` | `stemfinder(adata, markers, threshold, method)` |
| `compute_performance_single(adata, competitor, comp.inverted)` | `compute_performance(adata, competitor_key=..., competitor_inverted=...)` |
| `pct_recover(adata)` | `pct_recover(adata)` |
| `gene_set_score(gene.set, adata)` | `gene_set_score(adata, genes)` |
| `c(s_genes_mouse, g2m_genes_mouse)`, `mmTFs` | `cell_cycle_genes("mouse")`, `transcription_factors("mouse")` |

Differences:

- `compute_performance` computes the AUC over every pair of most and least differentiated cells. R's
  `auc_probability` compares only the first most differentiated cell with the least differentiated ones.
- The phenotype-level correlation is Spearman, as its name says; R's uses Pearson. Pass
  `pheno_method="pearson"` to reproduce R's value.
- `stemfinder` reads the neighborhood size from the kNN graph, so it takes no `k` argument.
- E2f8 (E2F8 in human) is an S phase gene only. The R package also lists it under G2M, so the R
  vignette counts it twice, and R counts any repeated marker once per occurrence; pystemfinder counts
  each marker once. Scores on the vignette's data therefore differ slightly from the published ones
  (rank correlation 0.9999).

See [CHANGELOG.md](CHANGELOG.md) for details.

## Documentation

https://pystemfinder.readthedocs.io

To build the docs locally:

```bash
pip install -r docs/requirements.txt
sphinx-build docs docs/_build/html
```

## Tests

```bash
pytest
```

## Citation

Noller K, Cahan P. Cell cycle expression heterogeneity predicts degree of differentiation.
*Briefings in Bioinformatics* 25(6):bbae536 (2024). https://doi.org/10.1093/bib/bbae536

## License

MIT; see [LICENSE](LICENSE).
