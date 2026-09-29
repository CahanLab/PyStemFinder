"""Parity with the R stemFinder package on the R vignette's Tabula Muris bone marrow data.

tests/data/bmmc_r_reference.h5ad holds R's inputs (scaled marker expression and Seurat's kNN graph) and
R's outputs; see tests/data/export_bmmc_reference.R. Both sides use the 91 mouse cell cycle genes with E2f8
counted once, so R's scores differ slightly from the vignette's published ones, which count it twice.
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


@pytest.mark.parametrize("method", ["stdev", "variance"])
def test_dispersion_scores_match_r(bmmc, method):
    # R computes these from log-normalized data (Seurat @data) rather than scale.data
    psf.run_stemFinder(bmmc, markers=bmmc.uns["R_markers"], method=method, layer="data")

    np.testing.assert_allclose(bmmc.obs["stemFinder_raw"], bmmc.obs[f"R_{method}_stemFinder_raw"], rtol=1e-12)
    np.testing.assert_allclose(bmmc.obs["stemFinder"], bmmc.obs[f"R_{method}_stemFinder"], rtol=1e-12, atol=1e-12)


def test_bundled_mouse_cell_cycle_genes_are_the_r_markers(bmmc):
    # R: unique(c(s_genes_mouse, g2m_genes_mouse)), all of which are in the data
    assert psf.cell_cycle_genes("mouse") == list(bmmc.uns["R_markers"])

    psf.run_stemFinder(bmmc, markers=psf.cell_cycle_genes("mouse"))

    np.testing.assert_allclose(bmmc.obs["stemFinder_raw"], bmmc.obs["R_gini_stemFinder_raw"], rtol=1e-12)


def test_performance_metrics_match_r(bmmc):
    bmmc.obs["stemFinder"] = bmmc.obs["R_gini_stemFinder"]
    r = bmmc.uns["R_performance"]

    # R's "phenotypic Spearman" is cor.test's default method, Pearson
    result = psf.compute_performance_single(bmmc, pheno_method="pearson")

    # R: cor.test(<these scores>, Ground_truth, method = "spearman"). R's in-session value, 0.741172286, differs by
    # 1e-6 because round-off in R's sums splits some tied scores, which the 15-digit export then re-ties.
    assert result.loc["stemFinder", "Spearman_SingleCell"] == pytest.approx(0.741173472392783, rel=1e-12)
    assert result.loc["stemFinder", "Spearman_Pheno"] == pytest.approx(r["Spearman_Pheno"], rel=1e-12)
    assert psf.pct_recover(bmmc) == pytest.approx(r["pct_recover"], rel=1e-12)


def test_equal_scores_tie_exactly(bmmc):
    # Mathematically equal gini scores must be bit-identical for rank-based metrics to treat them as ties.
    # Expected: R's cor.test on round(stemFinder, 12), where distinct scores differ by > 1e-5.
    psf.run_stemFinder(bmmc, markers=bmmc.uns["R_markers"])

    result = psf.compute_performance_single(bmmc)

    assert result.loc["stemFinder", "Spearman_SingleCell"] == pytest.approx(0.741173311223115, rel=1e-12)


def test_auc_and_spearman_pheno_are_exact(bmmc):
    # R reports AUC 0.9724572 because auc_probability compares only the first most-differentiated cell with the
    # least differentiated ones. Expected values from R: wilcox.test(pos, neg)$statistic / (364 * 399) and
    # cor.test(..., method = "spearman") on phenotype means.
    bmmc.obs["stemFinder"] = bmmc.obs["R_gini_stemFinder"]

    result = psf.compute_performance_single(bmmc, competitor_key="ccat_invert", competitor_inverted=True)

    assert result.loc["stemFinder", "AUC"] == pytest.approx(0.966258365694456, rel=1e-12)
    assert result.loc["stemFinder", "Spearman_Pheno"] == pytest.approx(0.865048975764109, rel=1e-12)
    np.testing.assert_allclose(
        result.loc["ccat_invert"], [0.675567206638321, 0.84115259521814, 0.953014404142224], rtol=1e-12
    )
