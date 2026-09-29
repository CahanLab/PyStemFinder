"""stemFinder scores and related per-cell measures of expression heterogeneity."""
import warnings

import numpy as np
import scipy.sparse


def _genes_present(adata, genes, label):
    """Return the unique `genes` that are in adata.var_names, warning about any that are absent or repeated.

    Unlike R, where a repeated marker is counted once per occurrence, each gene is counted once.
    """
    genes = list(genes)
    repeated = list(dict.fromkeys(g for i, g in enumerate(genes) if g in genes[:i]))
    if repeated:
        warnings.warn(f"Repeated {label} are counted once: {', '.join(repeated)}", UserWarning, stacklevel=3)
    genes = list(dict.fromkeys(genes))
    var_names = set(adata.var_names)
    present = [g for g in genes if g in var_names]
    missing = [g for g in genes if g not in var_names]
    if not present:
        raise ValueError(f"None of the {label} are in adata.var_names.")
    if missing:
        shown = ", ".join(missing[:10]) + (", ..." if len(missing) > 10 else "")
        warnings.warn(
            f"{len(missing)} of {len(genes)} {label} are not in adata.var_names and were ignored: {shown}",
            UserWarning,
            stacklevel=3,
        )
    return present


def _expression(adata, genes, layer=None):
    """Dense (cells x genes) expression from adata.X or adata.layers[layer]."""
    view = adata[:, genes]
    X = view.X if layer is None else view.layers[layer]
    return X.toarray() if scipy.sparse.issparse(X) else np.asarray(X)


def _binarized(adata, genes, threshold, layer=None):
    """Dense boolean (cells x genes) matrix of expression > threshold."""
    return _expression(adata, genes, layer) > threshold


def _neighbor_graph(adata, neighbors_key=None):
    """Binary (cells x cells) kNN adjacency, self excluded, from a graph computed by sc.pp.neighbors."""
    key = "neighbors" if neighbors_key is None else neighbors_key
    default_key = "distances" if neighbors_key is None else f"{neighbors_key}_distances"
    distances_key = adata.uns.get(key, {}).get("distances_key", default_key)
    if distances_key not in adata.obsp:
        raise ValueError(f"No neighbor graph in adata.obsp['{distances_key}']. Run sc.pp.neighbors first.")

    graph = scipy.sparse.csr_matrix(adata.obsp[distances_key])
    # Use the sparsity structure rather than the values: neighbors at distance 0 (duplicate cells)
    # are stored as explicit zeros. Self edges (as in Seurat's kNN graph) are dropped.
    rows = np.repeat(np.arange(graph.shape[0]), np.diff(graph.indptr))
    keep = rows != graph.indices
    return scipy.sparse.csr_matrix((np.ones(keep.sum()), (rows[keep], graph.indices[keep])), shape=graph.shape)


def stemfinder(adata, markers, threshold=0.0, method='gini', layer=None, neighbors_key=None, key_added='stemfinder'):
    """Score how differentiated each cell is (port of ``run_stemFinder`` from the R package).

    Each method measures how heterogeneous the expression of marker genes (typically S and G2M phase cell cycle
    genes, see :func:`cell_cycle_genes`) is within each cell's kNN neighborhood, which is high in less
    differentiated cells:

    - ``'gini'`` (default): expression is binarized at ``threshold`` and p_g is the fraction of the cell's
      neighbors (excluding the cell itself) whose binarized state matches the cell's. The raw score is the sum over
      markers of p_g * (1 - p_g). As in R, use scaled expression (see :func:`recipe_stemfinder`), so that the
      default threshold of 0 splits each gene at its mean.
    - ``'stdev'`` / ``'variance'``: the sum over markers of the sample standard deviation / variance of expression
      across the neighborhood including the cell itself. R uses log-normalized (not scaled) expression for these
      methods; pass the layer that holds it. ``threshold`` is not used.

    The neighborhood size is read from the kNN graph built by ``sc.pp.neighbors``. The R package uses
    k = sqrt(number of cells).

    Args:
        adata (anndata.AnnData): Annotated data matrix (cells x genes) with a kNN graph from ``sc.pp.neighbors``.
        markers (list of str): Marker genes. Genes absent from ``adata.var_names`` are ignored with a warning, and
            a gene listed twice counts once (in R it counts twice).
        threshold (float, optional): Expression threshold for binarizing, used by method ``'gini'``. Defaults to 0.
        method (str, optional): ``'gini'``, ``'stdev'``, or ``'variance'``. Defaults to ``'gini'``.
        layer (str, optional): Layer holding the expression to use instead of ``adata.X``. Defaults to None.
        neighbors_key (str, optional): Key of the neighbors graph, as passed to
            ``sc.pp.neighbors(key_added=...)``. Defaults to the graph in ``adata.obsp['distances']``.
        key_added (str, optional): Name of the score columns in ``adata.obs``. Defaults to ``'stemfinder'``.

    Returns:
        None: Adds ``adata.obs[key_added + '_raw']``, the raw score (higher means less differentiated), and
        ``adata.obs[key_added]``, ``1 - raw / max(raw)``, which is oriented like pseudotime (lower means less
        differentiated).
    """
    if method not in ['gini', 'stdev', 'variance']:
        raise ValueError("Invalid value for 'method'. Expected one of 'gini', 'stdev', 'variance'.")

    markers = _genes_present(adata, markers, "markers")
    nn = _neighbor_graph(adata, neighbors_key)

    if method == 'gini':
        n = np.asarray(nn.sum(axis=1)).ravel()
        # number of each cell's neighbors above threshold, per marker
        m = nn @ _binarized(adata, markers, threshold, layer).astype(float)
        # sum over markers of p * (1 - p) with p = m / n. This is symmetric in p and 1 - p, so it does not matter
        # whether p counts the neighbors that match the cell's own state (as in R) or those above threshold.
        # Summing the integer numerators keeps equal scores bit-identical, so rank-based metrics see them as ties.
        raw = np.sum(m * (n[:, None] - m), axis=1) / n**2
    else:
        nn = nn + scipy.sparse.identity(adata.n_obs, format='csr')
        n = np.asarray(nn.sum(axis=1)).ravel()[:, None]
        X = _expression(adata, markers, layer)
        mean = (nn @ X) / n
        variance = np.clip(((nn @ X**2) / n - mean**2) * n / (n - 1), 0, None)
        raw = np.sum(variance if method == 'variance' else np.sqrt(variance), axis=1)

    adata.obs[f'{key_added}_raw'] = raw
    adata.obs[key_added] = 1 - raw / raw.max()


def diffometer(adata, genes, threshold=0.0, weight_by='equal', include_self=True, layer=None, neighbors_key=None,
               key_added='diffometer'):
    """Score each cell by the Gini impurity of binarized gene expression across its kNN neighborhood.

    Expression of ``genes`` is binarized at ``threshold``, the Gini impurity 2p(1 - p) of each gene is computed
    across each cell's neighborhood (p = fraction of the neighborhood above threshold), and the weighted mean over
    genes is stored. Higher means more heterogeneous, i.e. less differentiated. With ``include_self=False`` and
    equal weights this equals ``2 * stemfinder_raw / len(genes)``.

    Args:
        adata (anndata.AnnData): Annotated data matrix (cells x genes) with a kNN graph from ``sc.pp.neighbors``.
        genes (list of str): Genes to consider. Genes absent from ``adata.var_names`` are ignored with a warning.
        threshold (float, optional): Expression threshold for binarizing. Defaults to 0.
        weight_by (str, optional): How to weight genes: ``'equal'``, ``'expression'`` (mean expression across
            cells), or ``'presence'`` (fraction of cells above threshold). Defaults to ``'equal'``.
        include_self (bool, optional): Whether the neighborhood includes the cell itself. Defaults to True.
        layer (str, optional): Layer holding the expression to use instead of ``adata.X``. Defaults to None.
        neighbors_key (str, optional): Key of the neighbors graph, as passed to
            ``sc.pp.neighbors(key_added=...)``. Defaults to the graph in ``adata.obsp['distances']``.
        key_added (str, optional): Name of the score column in ``adata.obs``. Defaults to ``'diffometer'``.

    Returns:
        None: The score is stored in ``adata.obs[key_added]``.
    """
    if weight_by not in ['equal', 'expression', 'presence']:
        raise ValueError("Invalid value for 'weight_by'. Expected one of 'equal', 'expression', 'presence'.")

    genes = _genes_present(adata, genes, "genes")
    expressed = _binarized(adata, genes, threshold, layer)

    if weight_by == 'expression':
        gene_weights = _expression(adata, genes, layer).mean(axis=0)
    elif weight_by == 'presence':
        gene_weights = expressed.mean(axis=0)
    else:
        gene_weights = np.ones(len(genes))

    nn = _neighbor_graph(adata, neighbors_key)
    if include_self:
        nn = nn + scipy.sparse.identity(adata.n_obs, format='csr')
    neighborhood_size = np.asarray(nn.sum(axis=1)).ravel()

    p = (nn @ expressed.astype(float)) / neighborhood_size[:, None]
    impurities = 2 * p * (1 - p)
    adata.obs[key_added] = impurities @ gene_weights / np.sum(gene_weights)


def gene_set_score(adata, genes, layer=None, key_added='gene_set_score'):
    """Compute the mean expression of a gene set in each cell (port of ``gene_set_score`` from the R package).

    Args:
        adata (anndata.AnnData): Annotated data matrix (cells x genes). R uses log-normalized expression.
        genes (list of str): Genes in the set. Genes absent from ``adata.var_names`` are ignored with a warning.
        layer (str, optional): Layer holding the expression to use instead of ``adata.X``. Defaults to None.
        key_added (str, optional): Name of the score column in ``adata.obs``. Defaults to ``'gene_set_score'``.

    Returns:
        None: The score is stored in ``adata.obs[key_added]``.
    """
    genes = _genes_present(adata, genes, "genes")
    adata.obs[key_added] = _expression(adata, genes, layer).mean(axis=1)


def count_expressed_genes(adata, genes, threshold=0.0, layer=None, key_added='n_expressed_genes'):
    """Count how many of the given genes each cell expresses above a threshold.

    Args:
        adata (anndata.AnnData): Annotated data matrix (cells x genes).
        genes (list of str): Genes to count. Genes absent from ``adata.var_names`` are ignored with a warning.
        threshold (float, optional): A gene counts if its expression is above this value. Defaults to 0.
        layer (str, optional): Layer holding the expression to use instead of ``adata.X``. Defaults to None.
        key_added (str, optional): Name of the count column in ``adata.obs``. Defaults to ``'n_expressed_genes'``.

    Returns:
        None: The counts are stored in ``adata.obs[key_added]``.
    """
    genes = _genes_present(adata, genes, "genes")
    adata.obs[key_added] = _binarized(adata, genes, threshold, layer).sum(axis=1)


def binarize(data, threshold=0.0):
    """Binarize expression values: 1 where above ``threshold``, 0 elsewhere.

    Args:
        data (numpy.ndarray or scipy.sparse.csr_matrix): Expression values.
        threshold (float, optional): Values above this are set to 1. Defaults to 0.

    Returns:
        numpy.ndarray or scipy.sparse.csr_matrix: Integer array of the same shape with values 0 and 1.
    """
    return (data > threshold).astype(int)


def gini_impurity(labels):
    """Compute the Gini impurity 2p(1 - p) of binary labels, where p is the fraction of ones.

    Args:
        labels (list of bool or int): Binary labels (0 or 1, True or False).

    Returns:
        float: The Gini impurity, from 0 (all labels equal) to 0.5 (half ones).

    Raises:
        ValueError: If the labels contain non-binary values.
    """
    if not all(label in [0, 1, True, False] for label in labels):
        raise ValueError("Labels list contains non-binary values.")

    p = np.mean(labels)
    return 2 * p * (1 - p)
