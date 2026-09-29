import warnings

import numpy as np
import pytest

import pystemfinder as psf


def test_generator_default_output_is_unchanged():
    # reference values recorded from PyStemFinder 0.1 (np.random.seed(42))
    adata = psf.simulate_data()

    assert adata.shape == (300, 50)
    assert int(adata.X.sum()) == 30881
    assert adata.X[0, :6].tolist() == [10, 8, 5, 4, 6, 8]


def test_generator_quickstart_parameters_output_is_unchanged():
    adata = psf.simulate_data(
        dropout_rate=0.25, stochastic_dropout_by_population=(0.25, 0.4, 0.6), stochastic_expression_diff=(0, 2, 4)
    )

    assert int(adata.X.sum()) == 40076
    assert int(adata.X[:, 40:45].sum()) == 3286
    assert adata.X[299, 40:45].tolist() == [0, 4, 0, 3, 4]


def test_generator_leaves_global_random_state_alone():
    np.random.seed(0)
    expected = np.random.rand()
    np.random.seed(0)

    psf.simulate_data()

    assert np.random.rand() == expected


def test_generator_random_state_changes_output():
    a = psf.simulate_data(random_state=1)
    b = psf.simulate_data(random_state=2)

    assert not np.array_equal(a.X, b.X)


def test_generator_emits_no_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        psf.simulate_data()


def test_recipe_excluded_genes_leave_hvgs_even_if_some_are_absent_from_data():
    adata = psf.simulate_data()
    hvgs = psf.recipe_stemfinder(adata, n_comps=10, copy=True).var.query("highly_variable").index[:3].tolist()
    assert len(hvgs) == 3  # guard: the excluded genes must be ones that would otherwise be HVGs

    out = psf.recipe_stemfinder(adata, exclude=hvgs + ["not_a_gene"], n_comps=10, copy=True)

    assert not out.var.loc[hvgs, "highly_variable"].any()
    assert out.obsm["X_pca"].shape == (300, 10)


def test_recipe_scales_every_gene_in_place():
    adata = psf.simulate_data()

    result = psf.recipe_stemfinder(adata, n_comps=10)

    assert result is None
    # the stochastic genes 40-44 are not all HVGs, but stemfinder needs them scaled too
    np.testing.assert_allclose(adata.X.mean(axis=0), 0, atol=1e-6)
    assert "X_pca" in adata.obsm


def test_recipe_copy_leaves_input_untouched():
    adata = psf.simulate_data()

    out = psf.recipe_stemfinder(adata, n_comps=10, copy=True)

    assert adata.X.sum() == 30881  # raw counts
    assert "X_pca" not in adata.obsm and "X_pca" in out.obsm


def test_recipe_without_scaling_keeps_log_normalized_data():
    adata = psf.simulate_data()

    psf.recipe_stemfinder(adata, scale=False, n_comps=10)

    assert adata.X.min() == 0  # log1p of normalized counts is never negative


def test_count_expressed_genes(toy):
    with pytest.warns(UserWarning, match="not_a_gene"):
        result = psf.count_expressed_genes(toy, ["g0", "g1", "not_a_gene"], key_added="n_on")

    assert result is None
    assert toy.obs["n_on"].tolist() == [1, 2, 1, 1]


def test_count_expressed_genes_threshold(toy_dense):
    # g0 [2, 3, 0, 1] and g1 [0, 1, 2, 0] above 1.5
    psf.count_expressed_genes(toy_dense, ["g0", "g1"], threshold=1.5)

    assert toy_dense.obs["n_expressed_genes"].tolist() == [1, 1, 1, 0]


def test_binarize():
    np.testing.assert_array_equal(psf.binarize(np.array([[0.5, -1.0], [2.0, 0.0]])), [[1, 0], [1, 0]])
    np.testing.assert_array_equal(psf.binarize(np.array([[0.5, -1.0], [2.0, 0.0]]), threshold=1.0), [[0, 0], [1, 0]])


def test_gini_impurity():
    assert psf.gini_impurity([1, 1, 0, 0]) == 0.5
    assert psf.gini_impurity([True, True, True]) == 0.0
    with pytest.raises(ValueError, match="binary"):
        psf.gini_impurity([0, 1, 2])


def test_gene_set_score_is_mean_expression(toy):
    with pytest.warns(UserWarning, match="not_a_gene"):
        psf.gene_set_score(toy, ["g0", "g1", "not_a_gene"])

    np.testing.assert_allclose(toy.obs["gene_set_score"], [1.0, 2.0, 1.0, 0.5])


def test_gene_set_score_layer_and_key(toy_dense):
    toy_dense.layers["lognorm"] = toy_dense.X.copy()
    toy_dense.X = np.zeros_like(toy_dense.X)

    psf.gene_set_score(toy_dense, ["g0", "g1"], layer="lognorm", key_added="cc")

    np.testing.assert_allclose(toy_dense.obs["cc"], [1.0, 2.0, 1.0, 0.5])
