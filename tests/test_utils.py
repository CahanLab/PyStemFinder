import warnings

import numpy as np
import pytest

import PyStemFinder as psf


def test_generator_default_output_is_unchanged():
    # reference values recorded from PyStemFinder 0.1 (np.random.seed(42))
    adata = psf.generate_scRNAseq_test_data()

    assert adata.shape == (300, 50)
    assert int(adata.X.sum()) == 30881
    assert adata.X[0, :6].tolist() == [10, 8, 5, 4, 6, 8]


def test_generator_quickstart_parameters_output_is_unchanged():
    adata = psf.generate_scRNAseq_test_data(
        dropout_rate=0.25, stochastic_dropout_by_population=(0.25, 0.4, 0.6), stochastic_expression_diff=(0, 2, 4)
    )

    assert int(adata.X.sum()) == 40076
    assert int(adata.X[:, 40:45].sum()) == 3286
    assert adata.X[299, 40:45].tolist() == [0, 4, 0, 3, 4]


def test_generator_leaves_global_random_state_alone():
    np.random.seed(0)
    expected = np.random.rand()
    np.random.seed(0)

    psf.generate_scRNAseq_test_data()

    assert np.random.rand() == expected


def test_generator_random_state_changes_output():
    a = psf.generate_scRNAseq_test_data(random_state=1)
    b = psf.generate_scRNAseq_test_data(random_state=2)

    assert not np.array_equal(a.X, b.X)


def test_generator_emits_no_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        psf.generate_scRNAseq_test_data()


def test_blacklisted_genes_leave_hvgs_even_if_some_are_absent_from_data():
    adata = psf.generate_scRNAseq_test_data()
    hvgs = psf.sf_norm_hvg_scale_pca(adata, [], n_comps=10).var.query("highly_variable").index[:3].tolist()
    assert len(hvgs) == 3  # guard: the blacklist must hit genes that would otherwise be HVGs

    out = psf.sf_norm_hvg_scale_pca(adata, hvgs + ["not_a_gene"], n_comps=10)

    assert not out.var.loc[hvgs, "highly_variable"].any()
    assert out.obsm["X_pca"].shape == (300, 10)


def test_count_high_expr_genes(toy):
    psf.count_high_expr_genes(toy, ["g0", "g1", "not_a_gene"], 0, "n_on")

    assert toy.obs["n_on"].tolist() == [1, 2, 1, 1]


def test_gene_set_score_is_mean_expression(toy):
    with pytest.warns(UserWarning, match="not_a_gene"):
        psf.gene_set_score(toy, ["g0", "g1", "not_a_gene"])

    np.testing.assert_allclose(toy.obs["gene_set_score"], [1.0, 2.0, 1.0, 0.5])


def test_gene_set_score_layer_and_key(toy_dense):
    toy_dense.layers["lognorm"] = toy_dense.X.copy()
    toy_dense.X = np.zeros_like(toy_dense.X)

    psf.gene_set_score(toy_dense, ["g0", "g1"], layer="lognorm", key_added="cc")

    np.testing.assert_allclose(toy_dense.obs["cc"], [1.0, 2.0, 1.0, 0.5])
