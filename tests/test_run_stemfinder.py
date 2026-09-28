import anndata
import numpy as np
import pytest
import scanpy as sc

import PyStemFinder as psf
from conftest import knn_graph

MARKERS = ["g0", "g1"]
# Hand-derived from the toy data in conftest.py (see its docstring):
# per gene, p = fraction of the 2 neighbors whose binarized state matches the cell; score = sum p * (1 - p)
RAW = [0.25, 0.0, 0.25, 0.5]
INVERTED = [0.5, 1.0, 0.5, 0.0]  # 1 - RAW / max(RAW)


def test_scores_match_hand_computed_values(toy):
    psf.run_stemFinder(toy, markers=MARKERS)

    np.testing.assert_allclose(toy.obs["stemFinder_raw"], RAW)
    np.testing.assert_allclose(toy.obs["stemFinder"], INVERTED)


def test_threshold_controls_binarization(toy_dense):
    # at 1.5, g0 -> [1, 1, 0, 0] and g1 -> [0, 0, 1, 0]
    psf.run_stemFinder(toy_dense, markers=MARKERS, thresh=1.5)

    np.testing.assert_allclose(toy_dense.obs["stemFinder_raw"], [0.5, 0.25, 0.25, 0.5])
    np.testing.assert_allclose(toy_dense.obs["stemFinder"], [0.0, 0.5, 0.5, 0.0])


def test_neighbors_at_distance_zero_still_count(toy_dense):
    # duplicate cells sit at distance 0; scanpy stores those edges as explicit zeros
    toy_dense.obsp["distances"] = knn_graph(
        {0: [1, 2], 1: [0, 3], 2: [1, 3], 3: [0, 2]},
        n_obs=4,
        values={0: [1.0, 1.0], 1: [0.0, 1.0], 2: [1.0, 1.0], 3: [1.0, 1.0]},
    )

    psf.run_stemFinder(toy_dense, markers=MARKERS)

    np.testing.assert_allclose(toy_dense.obs["stemFinder_raw"], RAW)


def test_self_edges_in_graph_are_ignored(toy_dense):
    # Seurat-style kNN graphs list each cell as its own neighbor
    toy_dense.obsp["distances"] = knn_graph({0: [0, 1, 2], 1: [0, 1, 3], 2: [1, 2, 3], 3: [0, 2, 3]}, n_obs=4)

    psf.run_stemFinder(toy_dense, markers=MARKERS)

    np.testing.assert_allclose(toy_dense.obs["stemFinder_raw"], RAW)


def test_uses_graph_named_by_neighbors_key(toy_dense):
    graph = toy_dense.obsp["distances"]
    del toy_dense.obsp["distances"]
    toy_dense.obsp["alt_distances"] = graph
    toy_dense.uns["alt"] = {"distances_key": "alt_distances"}

    psf.run_stemFinder(toy_dense, markers=MARKERS, neighbors_key="alt")

    np.testing.assert_allclose(toy_dense.obs["stemFinder_raw"], RAW)


def test_missing_neighbor_graph_raises(toy_dense):
    del toy_dense.obsp["distances"]

    with pytest.raises(ValueError, match="sc.pp.neighbors"):
        psf.run_stemFinder(toy_dense, markers=MARKERS)


def test_absent_markers_are_dropped_with_warning(toy_dense):
    with pytest.warns(UserWarning, match="not_a_gene"):
        psf.run_stemFinder(toy_dense, markers=MARKERS + ["not_a_gene"])

    np.testing.assert_allclose(toy_dense.obs["stemFinder_raw"], RAW)


def test_no_markers_present_raises(toy_dense):
    with pytest.raises(ValueError, match="markers"):
        psf.run_stemFinder(toy_dense, markers=["x", "y"])


def test_matches_literal_translation_of_r_code():
    rng = np.random.default_rng(0)
    adata = anndata.AnnData(rng.normal(size=(150, 12)))
    k = 10
    sc.pp.neighbors(adata, n_neighbors=k, use_rep="X")
    markers = list(adata.var_names[:8])

    psf.run_stemFinder(adata, markers=markers)

    # R/run_stemFinder.R (method = 'gini'); Seurat's kNN graph includes the cell itself
    exp_dat = adata[:, markers].X
    graph = adata.obsp["distances"]
    expected = []
    for i in range(adata.n_obs):
        neigh = np.concatenate([[i], graph.indices[graph.indptr[i] : graph.indptr[i + 1]]])
        exp = exp_dat[neigh] > 0
        n_match = (exp == (exp_dat[i] > 0)).sum(axis=0) - 1
        p_g = n_match / (k - 1)
        expected.append(np.sum(p_g * (1 - p_g)))
    expected = np.array(expected)
    np.testing.assert_allclose(adata.obs["stemFinder_raw"], expected)
    np.testing.assert_allclose(adata.obs["stemFinder"], 1 - expected / expected.max())
