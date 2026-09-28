import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse


def _toy_adata(sparse_x=False):
    """Four cells, three genes, and a hand-built kNN graph (2 neighbors per cell, self excluded).

    Binarized at 0, the marker genes g0/g1 are:
        c0 [1, 0]   neighbors c1, c2
        c1 [1, 1]   neighbors c0, c3
        c2 [0, 1]   neighbors c1, c3
        c3 [1, 0]   neighbors c0, c2
    g2 is heterogeneous so that scoring it by mistake changes the result.
    """
    X = np.array(
        [
            [2.0, 0.0, 5.0],
            [3.0, 1.0, 0.0],
            [0.0, 2.0, 5.0],
            [1.0, 0.0, 0.0],
        ]
    )
    adata = anndata.AnnData(
        X=scipy.sparse.csr_matrix(X) if sparse_x else X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    adata.obsp["distances"] = knn_graph({0: [1, 2], 1: [0, 3], 2: [1, 3], 3: [0, 2]}, n_obs=4)
    return adata


def knn_graph(neighbors, n_obs, values=None):
    """CSR distance matrix from {cell: [neighbor, ...]}; explicit zeros in `values` are kept."""
    indptr, indices, data = [0], [], []
    for i in range(n_obs):
        js = neighbors.get(i, [])
        indices.extend(js)
        data.extend(values[i] if values is not None else [1.0] * len(js))
        indptr.append(len(indices))
    return scipy.sparse.csr_matrix((np.array(data), np.array(indices), np.array(indptr)), shape=(n_obs, n_obs))


@pytest.fixture(params=[False, True], ids=["dense", "sparse"])
def toy(request):
    return _toy_adata(sparse_x=request.param)


@pytest.fixture
def toy_dense():
    return _toy_adata(sparse_x=False)
