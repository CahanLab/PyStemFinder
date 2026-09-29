"""stemfinder(method='stdev' | 'variance'): R sums, over markers, the sample sd / variance of expression
across each cell's kNN neighborhood (Seurat's graph, which includes the cell itself)."""
import numpy as np
import pytest

import pystemfinder as psf

MARKERS = ["g0", "g1"]
# Hand-derived from the toy data in conftest.py. Neighborhoods including self:
#   c0 {c0, c1, c2}: g0 [2, 3, 0] var 7/3, g1 [0, 1, 2] var 1
#   c1 {c1, c0, c3}: g0 [3, 2, 1] var 1,   g1 [1, 0, 0] var 1/3
#   c2 {c2, c1, c3}: g0 [0, 3, 1] var 7/3, g1 [2, 1, 0] var 1
#   c3 {c3, c0, c2}: g0 [1, 2, 0] var 1,   g1 [0, 0, 2] var 4/3
VARIANCE_RAW = [10 / 3, 4 / 3, 10 / 3, 7 / 3]
STDEV_RAW = [np.sqrt(7 / 3) + 1, 1 + np.sqrt(1 / 3), np.sqrt(7 / 3) + 1, 1 + np.sqrt(4 / 3)]


def test_variance(toy):
    psf.stemfinder(toy, markers=MARKERS, method="variance")

    np.testing.assert_allclose(toy.obs["stemfinder_raw"], VARIANCE_RAW)
    np.testing.assert_allclose(toy.obs["stemfinder"], [0.0, 0.6, 0.0, 0.3], atol=1e-12)


def test_stdev(toy):
    psf.stemfinder(toy, markers=MARKERS, method="stdev")

    np.testing.assert_allclose(toy.obs["stemfinder_raw"], STDEV_RAW)
    np.testing.assert_allclose(toy.obs["stemfinder"], 1 - np.array(STDEV_RAW) / max(STDEV_RAW), atol=1e-12)


def test_layer_is_used_instead_of_x(toy_dense):
    toy_dense.layers["lognorm"] = toy_dense.X.copy()
    toy_dense.X = np.zeros_like(toy_dense.X)

    psf.stemfinder(toy_dense, markers=MARKERS, method="variance", layer="lognorm")

    np.testing.assert_allclose(toy_dense.obs["stemfinder_raw"], VARIANCE_RAW)


def test_invalid_method_raises(toy_dense):
    with pytest.raises(ValueError, match="method"):
        psf.stemfinder(toy_dense, markers=MARKERS, method="entropy")
