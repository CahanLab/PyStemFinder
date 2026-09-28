"""Parity with the R stemFinder package on the R vignette's Tabula Muris bone marrow data.

tests/data/bmmc_r_reference.h5ad holds R's inputs (scaled marker expression and Seurat's kNN graph) and
R's outputs; see tests/data/export_bmmc_reference.R. R's gini scores there are identical to the published
https://cnobjects.s3.amazonaws.com/stemFinder/bmmc_sF_results.csv.
"""
from pathlib import Path

import anndata
import numpy as np
import pytest

import PyStemFinder as psf

FIXTURE = Path(__file__).parent / "data" / "bmmc_r_reference.h5ad"


@pytest.fixture
def bmmc():
    return anndata.read_h5ad(FIXTURE)


def test_gini_scores_match_r(bmmc):
    psf.run_stemFinder(bmmc, markers=bmmc.uns["R_markers"])

    np.testing.assert_allclose(bmmc.obs["stemFinder_raw"], bmmc.obs["R_gini_stemFinder_raw"], rtol=1e-12)
    np.testing.assert_allclose(bmmc.obs["stemFinder"], bmmc.obs["R_gini_stemFinder"], rtol=1e-12, atol=1e-12)
