import anndata
import numpy as np
import pandas as pd
import pytest

import PyStemFinder as psf

# Ground truth: higher = more differentiated. Phenotype means: A (0.3, 1), B (0.2, 2), C (0.9, 3), D (0.625, 4).
GROUND_TRUTH = [1, 1, 2, 3, 4, 4]
SCORES = [0.1, 0.5, 0.2, 0.9, 0.3, 0.95]
PHENOTYPE = ["A", "A", "B", "C", "D", "D"]

# Spearman = Pearson on average ranks: ground truth ranks [1.5, 1.5, 3, 4, 5.5, 5.5], score ranks [1, 4, 2, 5, 3, 6]
SPEARMAN_SINGLE_CELL = 9.5 / np.sqrt(16.5 * 17.5)
# phenotype mean score ranks [2, 1, 4, 3] vs [1, 2, 3, 4]: 1 - 6 * 4 / (4 * 15)
SPEARMAN_PHENO = 0.6
PEARSON_PHENO = 0.8375 / np.sqrt(0.30546875 * 5)
# most differentiated (4) [0.3, 0.95] vs least (1) [0.1, 0.5]: 3 of 4 pairs ranked correctly
AUC = 0.75


def _adata(scores=SCORES, ground_truth=GROUND_TRUTH, **obs):
    n = len(scores)
    df = pd.DataFrame(
        {"stemFinder": scores, "Ground_truth": ground_truth, "Phenotype": PHENOTYPE[:n], **obs},
        index=[f"c{i}" for i in range(n)],
    )
    return anndata.AnnData(np.zeros((n, 1)), obs=df)


def test_performance_metrics():
    result = psf.compute_performance_single(_adata())

    assert list(result.columns) == ["Spearman_SingleCell", "Spearman_Pheno", "AUC"]
    np.testing.assert_allclose(result.loc["stemFinder"], [SPEARMAN_SINGLE_CELL, SPEARMAN_PHENO, AUC])


def test_auc_uses_every_positive_cell():
    # R's auc_probability compares only the first positive cell (0.3), which gives 0.5 here
    result = psf.compute_performance_single(_adata())

    assert result.loc["stemFinder", "AUC"] == pytest.approx(0.75)


def test_auc_gives_half_credit_for_ties():
    # positives [0.5, 0.9] vs negatives [0.2, 0.5]: 3 wins + 1 tie out of 4 pairs
    adata = _adata(scores=[0.2, 0.5, 0.5, 0.9], ground_truth=[1, 1, 2, 2])

    result = psf.compute_performance_single(adata)

    assert result.loc["stemFinder", "AUC"] == pytest.approx(0.875)


def test_pearson_phenotype_correlation_reproduces_r():
    result = psf.compute_performance_single(_adata(), pheno_method="pearson")

    assert result.loc["stemFinder", "Spearman_Pheno"] == pytest.approx(PEARSON_PHENO)


def test_competitor_scores_are_inverted_unless_already_inverted():
    competitor = 1 - np.array(SCORES)  # opposite orientation to stemFinder
    adata = _adata(ccat=competitor)

    inverted = psf.compute_performance_single(adata, competitor_key="ccat")
    as_is = psf.compute_performance_single(adata, competitor_key="ccat", competitor_inverted=True)

    # 1 - x / max(x) of the competitor ranks cells exactly like stemFinder
    np.testing.assert_allclose(inverted.loc["ccat"], [SPEARMAN_SINGLE_CELL, SPEARMAN_PHENO, AUC])
    np.testing.assert_allclose(as_is.loc["ccat"], [-SPEARMAN_SINGLE_CELL, -SPEARMAN_PHENO, 1 - AUC])


def test_custom_keys():
    adata = _adata()
    adata.obs = adata.obs.rename(columns={"stemFinder": "score", "Ground_truth": "gt", "Phenotype": "celltype"})

    result = psf.compute_performance_single(adata, score_key="score", ground_truth_key="gt", phenotype_key="celltype")

    np.testing.assert_allclose(result.loc["score"], [SPEARMAN_SINGLE_CELL, SPEARMAN_PHENO, AUC])


def test_pct_recover():
    # 4 ground truth levels -> 25th percentile of scores = 0.225; least differentiated cells score [0.1, 0.5]
    assert psf.pct_recover(_adata()) == pytest.approx(50.0)
