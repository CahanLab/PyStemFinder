import warnings

import numpy as np
import pandas as pd
import scipy.sparse
import anndata
import scanpy as sc

__all__ = [
    'binary_gini_impurity',
    'binarize_data',
    'diffOmeter',
    'run_stemFinder',
    'generate_scRNAseq_test_data',
    'sf_norm_hvg_scale_pca',
    'count_high_expr_genes',
]

def binary_gini_impurity(labels):
    """
    Compute the Gini impurity of a binary list of labels.
    
    Args:
        labels (List[bool or int]): A list of binary labels (0 or 1, True or False).
        
    Returns:
        float: The computed Gini impurity.
        
    Raises:
        ValueError: If the labels list contains non-binary values.
    """
    # Check for non-binary values
    if not all(label in [0, 1, True, False] for label in labels):
        raise ValueError("Labels list contains non-binary values.")
    
    p = np.mean(labels)
    return 2 * p * (1 - p)


def binarize_data(data, threshold):
    """
    Binarize data based on a given threshold.
    
    Args:
        data (numpy.ndarray): A numpy array of data values to be binarized.
        threshold (float): The threshold value. Data values greater than this threshold will be set to 1, 
                           otherwise they will be set to 0.
                           
    Returns:
        numpy.ndarray: A binarized numpy array where values greater than the threshold are set to 1, 
                       and all others are set to 0.
    """
    return (data > threshold).astype(int)

def _genes_present(adata, genes, label):
    """Return `genes` that are in adata.var_names, warning about any that are absent.

    Repeated genes are kept, as in R, so they are counted once per occurrence.
    """
    genes = list(genes)
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


def _binarized(adata, genes, threshold):
    """Dense boolean (cells x genes) matrix of expression > threshold."""
    X = adata[:, genes].X
    X = X.toarray() if scipy.sparse.issparse(X) else np.asarray(X)
    return X > threshold


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


def diffOmeter(adata, genes, threshold, weight_by='equal', include_self=True, neighbors_key=None):
    """
    Compute a 'diffOmeter' score for each cell in an anndata object based on the Gini impurity of binary gene expression.
    
    This function binarizes the expression of `genes` at `threshold`, computes the Gini impurity 2p(1 - p) of each 
    gene across the kNN neighborhood of each cell (p = fraction of the neighborhood above threshold), and stores the 
    weighted mean of these impurities per cell.
    
    Args:
        adata (anndata.AnnData): The annotated data matrix of shape (n_obs, n_vars), with a kNN graph from 
                                 sc.pp.neighbors.
        genes (list of str): Genes to consider. Genes absent from adata.var_names are ignored with a warning.
        threshold (float): The threshold value used to binarize gene expression data.
        weight_by (str, optional): Method to weight genes. 'equal', 'expression' (mean expression across cells), or 
                                   'presence' (fraction of cells above threshold). Defaults to 'equal'.
        include_self (bool, optional): Whether to include the cell itself when considering its neighborhood. 
                                       Defaults to True.
        neighbors_key (str, optional): Key of the neighbors graph, as passed to sc.pp.neighbors(key_added=...). 
                                       Defaults to the graph in adata.obsp['distances'].
                                       
    Returns:
        None: The results are stored in the 'diffOmeter' column of the .obs attribute of the input anndata object.
        
    Raises:
        ValueError: If the 'weight_by' parameter is not one of the expected values ('equal', 'expression', 'presence').
    """
    if weight_by not in ['equal', 'expression', 'presence']:
        raise ValueError("Invalid value for 'weight_by'. Expected one of 'equal', 'expression', 'presence'.")

    genes = _genes_present(adata, genes, "genes")
    expressed = _binarized(adata, genes, threshold)

    if weight_by == 'expression':
        gene_weights = np.asarray(adata[:, genes].X.mean(axis=0)).ravel()
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
    adata.obs['diffOmeter'] = impurities @ gene_weights / np.sum(gene_weights)


def run_stemFinder(adata, markers, thresh=0.0, neighbors_key=None):
    """
    Compute the 'stemFinder' scores for each cell in an anndata object.
    
    Port of run_stemFinder (method = 'gini') from the R package (https://github.com/CahanLab/stemfinder). For each 
    cell and marker gene, expression is binarized at `thresh` and p_g is the fraction of the cell's kNN neighbors 
    (excluding the cell itself) whose binarized state matches the cell's. The raw score is the sum over markers of 
    p_g * (1 - p_g): heterogeneous marker expression within a neighborhood, which is high in less differentiated cells.
    
    As in R, the input should be scaled expression (e.g. sc.pp.scale, or sf_norm_hvg_scale_pca(gene_scale=True)), so 
    that the default threshold of 0 splits each gene at its mean. The markers are typically S and G2M phase cell 
    cycle genes. The neighborhood size (k - 1) is read from the kNN graph built by sc.pp.neighbors.
    
    Args:
        adata (anndata.AnnData): The annotated data matrix of shape (n_obs, n_vars), with a kNN graph from 
                                 sc.pp.neighbors.
        markers (list of str): Marker genes. Markers absent from adata.var_names are ignored with a warning. As in 
                               R, a marker listed twice counts twice.
        thresh (float, optional): The threshold value used to binarize gene expression data. Defaults to 0.
        neighbors_key (str, optional): Key of the neighbors graph, as passed to sc.pp.neighbors(key_added=...). 
                                       Defaults to the graph in adata.obsp['distances'].
                                       
    Returns:
        None: Adds two columns to adata.obs, matching the R package:
              'stemFinder_raw': raw score; higher = less differentiated.
              'stemFinder': 1 - stemFinder_raw / max(stemFinder_raw); lower = less differentiated (like pseudotime).
    """
    markers = _genes_present(adata, markers, "markers")
    nn = _neighbor_graph(adata, neighbors_key)
    n_neighbors = np.asarray(nn.sum(axis=1)).ravel()

    # fraction of each cell's neighbors above threshold, per marker
    q = (nn @ _binarized(adata, markers, thresh).astype(float)) / n_neighbors[:, None]
    # p * (1 - p) is symmetric in p and 1 - p, so it does not matter whether p counts the neighbors that
    # match the cell's own state (as in R) or those above threshold
    raw = np.sum(q * (1 - q), axis=1)

    adata.obs['stemFinder_raw'] = raw
    adata.obs['stemFinder'] = 1 - raw / raw.max()


def generate_scRNAseq_test_data(n_cells=300, n_genes=50, lambda_val=2, dropout_rate=0.6, 
                                n_populations=3, stochastic_gene_range=(40, 45),
                                marker_expression_diff=5, stochastic_expression_diff=(0, 0, 0),
                                stochastic_mean=None, stochastic_dropout_by_population=None, random_state=42):
    """
    Generate synthetic single-cell RNA sequencing (scRNAseq) data.
    
    This function creates synthetic scRNAseq data by simulating expression values for a given number of cells and genes.
    Base expression values are sampled from a Poisson distribution, and general dropout is introduced to simulate the
    zero-inflation observed in real scRNAseq data.
    
    Args:
        n_cells (int, optional): Total number of cells in the generated dataset. Defaults to 300.
        n_genes (int, optional): Number of genes in the dataset. Defaults to 50.
        lambda_val (float, optional): Mean of the Poisson distribution used to generate base expression values. Defaults to 2.
        dropout_rate (float, optional): General dropout rate for the data. Defaults to 0.6.
        n_populations (int, optional): Number of cell populations in the dataset. Defaults to 3.
        stochastic_gene_range (tuple of int, optional): Range (inclusive) of gene indices that have stochastic expression. Defaults to (40, 45).
        marker_expression_diff (float, optional): Difference in expression for marker genes. Defaults to 5.
        stochastic_expression_diff (tuple of float, optional): Difference in expression for stochastic genes. Defaults to (0, 0, 0).
        stochastic_mean (float, optional): Mean expression for stochastic genes. If None, uses the global mean. Defaults to None.
        stochastic_dropout_by_population (list or array-like, optional): Dropout rates for stochastic genes for each population. 
                                                                        If None, uses the global dropout rate. Defaults to None.
        random_state (int, optional): Seed for the random number generator. The global numpy random state is not 
                                      touched. Defaults to 42.

    Returns:
        anndata.AnnData: An AnnData object containing the generated synthetic scRNAseq data.
    """
    rng = np.random.RandomState(random_state)  # same stream as np.random.seed(random_state)
    cells_per_population = n_cells // n_populations

    # Base expression matrix
    data = rng.poisson(lambda_val, (n_cells, n_genes))

    # Introducing general dropout
    dropout_mask = (rng.rand(n_cells, n_genes) < dropout_rate)
    data[dropout_mask] = 0

    # If stochastic dropout by population is not provided, use the general dropout rate
    if stochastic_dropout_by_population is None:
        stochastic_dropout_by_population = [dropout_rate] * n_populations

    if stochastic_mean is None:
        stochastic_mean = lambda_val

    # Add unique expression patterns for marker genes of each population
    for i in range(n_populations):
        start_idx = i * cells_per_population
        end_idx = (i + 1) * cells_per_population
        
        marker_gene_counts = rng.poisson(lambda_val + marker_expression_diff, (cells_per_population, 10))
        data[start_idx:end_idx, i*10:i*10+10] = marker_gene_counts

        # Stochastic expression for genes in the specified range
        stochastic_gene_indices = range(*stochastic_gene_range)
        
        for gene_idx in stochastic_gene_indices:
            stochastic_expr = rng.poisson(stochastic_mean + stochastic_expression_diff[i], cells_per_population)
            dropout_mask = (rng.rand(cells_per_population) < stochastic_dropout_by_population[i])
            stochastic_expr[dropout_mask] = 0
            data[start_idx:end_idx, gene_idx] = stochastic_expr


    genes = [f'gene{i}' for i in range(n_genes)]
    cells = [f'cell{i}' for i in range(n_cells)]

    return anndata.AnnData(X=data, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))



def sf_norm_hvg_scale_pca(
    adQ: anndata.AnnData,
    blacklist,
    tsum: float = 1e4,
    min_mean: float = 0.0125,
    max_mean: float = 6,
    min_disp: float = 0.25,
    scale_max: float = 10,
    n_comps: int = 100,
    gene_scale: bool = False
) -> anndata.AnnData:
    """
    Normalize, detect highly variable genes, optionally scale, and perform PCA on an AnnData object.
    
    This function takes a copy of an AnnData object, normalizes and log-transforms it, identifies highly variable genes 
    (excluding blacklist genes), optionally scales all genes, and finally performs PCA on the highly variable genes.
    
    Args:
        adQ (anndata.AnnData): Annotated data matrix with observations (cells) and variables (features).
        blacklist (list of str): Genes to exclude from the highly variable genes used for PCA (e.g. cell cycle genes). 
                                 Genes absent from the data are ignored.
        tsum (float, optional): The total count to which data is normalized. Defaults to 1e4.
        min_mean (float, optional): Minimum mean expression value of genes to be considered as highly variable. Defaults to 0.0125.
        max_mean (float, optional): Maximum mean expression value of genes to be considered as highly variable. Defaults to 6.
        min_disp (float, optional): Minimum dispersion value of genes to be considered as highly variable. Defaults to 0.25.
        scale_max (float, optional): Maximum value of scaled expression data. Defaults to 10.
        n_comps (int, optional): Number of principal components to compute. Defaults to 100.
        gene_scale (bool, optional): Whether to scale every gene to zero mean and unit variance (clipped at scale_max), 
                                     as run_stemFinder expects. Defaults to False.
    
    Returns:
        anndata.AnnData: Annotated data matrix with normalized, highly variable, optionally scaled, and PCA-transformed data.
    """

    # Create a copy of the input data
    adata = adQ.copy()

    # Normalize the data to a target sum
    sc.pp.normalize_total(adata, target_sum=tsum)

    # Log-transform the data
    sc.pp.log1p(adata)

    # Detect highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp
    )

    # Optionally scale all genes (the stemFinder markers are usually not HVGs)
    if gene_scale:
        sc.pp.scale(adata, max_value=scale_max)

    adata.var.loc[adata.var_names.isin(blacklist), 'highly_variable'] = False

    # Perform PCA on the data
    sc.tl.pca(adata, n_comps=n_comps)

    return adata


def count_high_expr_genes(adata, genes, threshold, column_name):
    """
    Add a new .obs column to an AnnData object counting genes with expression above a given threshold.
    
    This function counts the number of genes from a provided list that have expression values greater than the specified 
    threshold for each cell. The counts are stored in a new column in the .obs attribute of the AnnData object.
    
    Args:
        adata (anndata.AnnData): Annotated data matrix with observations (cells) and variables (features).
        genes (list of str): List of genes to consider for counting.
        threshold (float): Expression threshold above which genes are counted.
        column_name (str): Name for the new column in the .obs dataframe of the AnnData object.
    
    Returns:
        anndata.AnnData: Modified AnnData object with the new column added to its .obs dataframe.
    """
    
    # Ensure genes in list are in the AnnData object
    valid_genes = [gene for gene in genes if gene in adata.var.index]

    # Extract the subset of the data matrix (cells x genes)
    submatrix = adata[:, valid_genes].X

    # Count genes with expression > threshold for each cell
    high_expr_count = (submatrix > threshold).sum(axis=1)

    # Convert to 1D array or Series if it's not (depends on the backend of AnnData's X)
    if isinstance(high_expr_count, pd.DataFrame):
        high_expr_count = high_expr_count.squeeze()

    # Add to the .obs dataframe
    adata.obs[column_name] = high_expr_count

    return adata





