import numpy as np
import pytest

import PyStemFinder as psf
from conftest import knn_graph

MARKERS = ["g0", "g1"]
# Hand-derived from the toy data in conftest.py. With self included each neighborhood has 3 cells,
# so per-gene impurity 2p(1-p) is 0 (all agree) or 4/9 (2 vs 1). Only c1 has a pure gene (g0).
EQUAL = [4 / 9, 2 / 9, 4 / 9, 4 / 9]


def test_equal_weights(toy):
    psf.diffOmeter(toy, MARKERS, 0.0)

    np.testing.assert_allclose(toy.obs["diffOmeter"], EQUAL)


@pytest.mark.parametrize(
    "weight_by, c1_score",
    [
        ("expression", 4 / 27),  # gene weights = mean expression = [1.5, 0.75]
        ("presence", 8 / 45),  # gene weights = fraction of cells above threshold = [0.75, 0.5]
    ],
)
def test_weighted(toy, weight_by, c1_score):
    psf.diffOmeter(toy, MARKERS, 0.0, weight_by=weight_by)

    np.testing.assert_allclose(toy.obs["diffOmeter"], [4 / 9, c1_score, 4 / 9, 4 / 9])


def test_excluding_self(toy_dense):
    psf.diffOmeter(toy_dense, MARKERS, 0.0, include_self=False)

    np.testing.assert_allclose(toy_dense.obs["diffOmeter"], [0.25, 0.0, 0.25, 0.5])


def test_neighbors_at_distance_zero_still_count(toy_dense):
    toy_dense.obsp["distances"] = knn_graph(
        {0: [1, 2], 1: [0, 3], 2: [1, 3], 3: [0, 2]},
        n_obs=4,
        values={0: [1.0, 1.0], 1: [0.0, 1.0], 2: [1.0, 1.0], 3: [1.0, 1.0]},
    )

    psf.diffOmeter(toy_dense, MARKERS, 0.0)

    np.testing.assert_allclose(toy_dense.obs["diffOmeter"], EQUAL)


def test_invalid_weight_by_raises(toy_dense):
    with pytest.raises(ValueError, match="weight_by"):
        psf.diffOmeter(toy_dense, MARKERS, 0.0, weight_by="nope")
