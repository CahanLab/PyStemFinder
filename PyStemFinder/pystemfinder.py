import numpy as np
import pandas as pd
import scipy.sparse
import anndata
import scanpy as sc

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

def diffOmeter(adata, genes, threshold, weight_by='equal', include_self=True):
    """
    Compute a 'diffOmeter' score for each cell in an anndata object based on the Gini impurity of binary gene expression.
    
    This function binarizes the gene expression data of each cell based on a provided threshold. Then, it calculates 
    the Gini impurity for each gene across the neighborhood of each cell. Finally, a weighted mean of these impurities 
    is computed for each cell, based on the specified weighting scheme, and stored in the anndata object.
    
    Args:
        adata (anndata.AnnData): The annotated data matrix of shape (n_obs, n_vars). Rows correspond to cells 
                                 and columns to genes.
        genes (list of str): List of gene names to consider for the computation.
        threshold (float): The threshold value used to binarize gene expression data.
        weight_by (str, optional): Method to weight genes. Can be 'equal', 'expression', or 'presence'. 
                                   Defaults to 'equal'.
        include_self (bool, optional): Whether to include the cell itself when considering its neighborhood. 
                                       Defaults to True.
                                       
    Returns:
        None: The results are stored in the 'diffOmeter' column of the .obs attribute of the input anndata object.
        
    Raises:
        ValueError: If the 'weight_by' parameter is not one of the expected values ('equal', 'expression', 'presence').
    """
    # Check for valid weight_by value
    if weight_by not in ['equal', 'expression', 'presence']:
        raise ValueError("Invalid value for 'weight_by'. Expected one of 'equal', 'expression', 'presence'.")
    
    # Convert gene expression data to binary based on threshold
    binarized_data = (adata[:, genes].X > threshold).astype(int)
    
    # Precompute gene weights if necessary
    if weight_by == 'expression':
        gene_weights = np.mean(adata[:, genes].X, axis=0).A1  # .A1 to convert to 1D array
    elif weight_by == 'presence':
        gene_weights = np.mean(binarized_data.toarray(), axis=0)
    else:
        gene_weights = np.ones(len(genes))
    
    n_cells = adata.shape[0]
    impurities = np.zeros(n_cells)
    
    nn = adata.obsp['distances']

    # Progress tracking
    ten_percent = n_cells // 10

    for i in range(n_cells):
        neighbors = nn[i].nonzero()[1]
        if include_self:
            neighbors = np.append(i, neighbors)
        
        if isinstance(binarized_data, np.ndarray):
            all_data = binarized_data[neighbors]
        else:
            all_data = binarized_data[neighbors].toarray()  # Convert only the slice to dense
        
        # Calculate the weighted mean of impurities
        impurities_per_gene = np.apply_along_axis(binary_gini_impurity, 0, all_data)
        weighted_impurities = np.sum(impurities_per_gene * gene_weights) / np.sum(gene_weights)
        impurities[i] = weighted_impurities

        # Progress tracking
        if i % ten_percent == 0:
            print(".", end="", flush=True)

    print()  # To ensure newline after dots

    # Store the results in the anndata object
    adata.obs['diffOmeter'] = impurities
    


def run_stemFinder(adata, k, thresh, markers):
    """
    Compute the 'stemFinder' scores for each cell in an anndata object.
    
    This function calculates a Gini-based metric to determine the stemness of each cell based on the expression 
    of marker genes. The marker gene expression is binarized based on a provided threshold, and the Gini impurity 
    is computed for each marker gene across the neighborhood of each cell. The results are stored in the input 
    anndata object.
    
    Args:
        adata (anndata.AnnData): The annotated data matrix of shape (n_obs, n_vars). Rows correspond to cells 
                                 and columns to genes.
        k (int): Number of nearest neighbors considered.
        thresh (float): The threshold value used to binarize gene expression data.
        markers (list of str): List of gene names (markers) to consider for the computation.
                                       
    Returns:
        None: The results are stored in the 'stemFinder', 'stemFinder_invert', and 'stemFinder_comp' columns of 
              the .obs attribute of the input anndata object.
    """
    # Assuming adata.X is the equivalent of adata@assays$RNA@scale.data
    expDat = adata[:, markers].X.copy()

    gini_agg = pd.DataFrame(index=adata.obs.index, 
                            columns=['gini_index_agg'], 
                            data=np.nan)

    # get nearest neighbors for each cell from precomputed neighbors graph
    nn = adata.obsp['distances']

    for i, cell in enumerate(adata.obs.index):
        neigh = nn[i].nonzero()[1]
        exp = expDat[neigh, :] > thresh
        exp_cell = expDat[i,:] > thresh
        if isinstance(exp_cell, scipy.sparse.csr_matrix):
            exp_cell = scipy.sparse.csr_matrix(np.repeat(exp_cell.A, exp.shape[0], axis=0))
        exp_match = exp == exp_cell
        n_match = np.sum(exp_match, axis=0) # sum across cells for each gene
        if isinstance(n_match, np.matrix):
            n_match = n_match.A.flatten()
        p_g = n_match / (k - 1) # gene-specific match rate for neighboring cells {0 - 1}
        gini_g = p_g * (1 - p_g) # if p_g == 1 (all matches) ||  0 (no matches) --> gini_g == 0
        gini_agg.loc[cell, 'gini_index_agg'] = np.sum(gini_g)

    adata.obs['stemFinder'] = gini_agg['gini_index_agg']
    adata.obs['stemFinder_invert'] = 1 - adata.obs['stemFinder'] / np.max(adata.obs['stemFinder'])
    adata.obs['stemFinder_comp'] = adata.obs['stemFinder'] / len(markers)


def generate_scRNAseq_test_data(n_cells=300, n_genes=50, lambda_val=2, dropout_rate=0.6, 
                                n_populations=3, stochastic_gene_range=(40, 45),
                                marker_expression_diff=5, stochastic_expression_diff=(0, 0, 0),
                                stochastic_mean=None, stochastic_dropout_by_population=None):
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

    Returns:
        anndata.AnnData: An AnnData object containing the generated synthetic scRNAseq data.
    """
    np.random.seed(42)  # for reproducibility
    cells_per_population = n_cells // n_populations

    # Base expression matrix
    data = np.random.poisson(lambda_val, (n_cells, n_genes))

    # Introducing general dropout
    dropout_mask = (np.random.rand(n_cells, n_genes) < dropout_rate)
    data[dropout_mask] = 0

    marker_genes = []

    # If stochastic dropout by population is not provided, use the general dropout rate
    if stochastic_dropout_by_population is None:
        stochastic_dropout_by_population = [dropout_rate] * n_populations

    if stochastic_mean is None:
        stochastic_mean = lambda_val

    # Add unique expression patterns for marker genes of each population
    for i in range(n_populations):
        start_idx = i * cells_per_population
        end_idx = (i + 1) * cells_per_population
        
        marker_gene_counts = np.random.poisson(lambda_val + marker_expression_diff, 
                                               (cells_per_population, 10))
        data[start_idx:end_idx, i*10:i*10+10] = marker_gene_counts
        marker_genes.extend(range(i*10, i*10+10))

        # Stochastic expression for genes in the specified range
        stochastic_gene_indices = range(*stochastic_gene_range)
        
        for gene_idx in stochastic_gene_indices:
            stochastic_expr = np.random.poisson(stochastic_mean + stochastic_expression_diff[i], cells_per_population)
            dropout_mask = (np.random.rand(cells_per_population) < stochastic_dropout_by_population[i])
            stochastic_expr[dropout_mask] = 0
            data[start_idx:end_idx, gene_idx] = stochastic_expr


    genes = [f'gene{i}' for i in range(n_genes)]
    cells = [f'cell{i}' for i in range(n_cells)]

    return anndata.AnnData(X=data, dtype=data.dtype, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))



def sf_norm_hvg_scale_pca(
    adQ: anndata,
    blacklist, 
    tsum: float = 1e4,
    min_mean: float = 0.0125,
    max_mean: float = 6,
    min_disp: float = 0.25,
    scale_max: float = 10,
    n_comps: int = 100,
    gene_scale: bool = False
) -> anndata:
    """
    Normalize, detect highly variable genes, optionally scale, and perform PCA on an AnnData object.
    
    This function takes an AnnData object, normalizes it, identifies highly variable genes (excluding blacklist genes), 
    optionally scales the expression values of these genes, and finally performs PCA on the data.
    
    Args:
        adQ (anndata.AnnData): Annotated data matrix with observations (cells) and variables (features).
        blacklist (list of str): List of genes to exclude from HVG detection and PCA.
        tsum (float, optional): The total count to which data is normalized. Defaults to 1e4.
        min_mean (float, optional): Minimum mean expression value of genes to be considered as highly variable. Defaults to 0.0125.
        max_mean (float, optional): Maximum mean expression value of genes to be considered as highly variable. Defaults to 6.
        min_disp (float, optional): Minimum dispersion value of genes to be considered as highly variable. Defaults to 0.25.
        scale_max (float, optional): Maximum value of scaled expression data. Defaults to 10.
        n_comps (int, optional): Number of principal components to compute. Defaults to 100.
        gene_scale (bool, optional): Whether to scale the expression values of highly variable genes. Defaults to False.
    
    Returns:
        anndata.AnnData: Annotated data matrix with normalized, highly variable, optionally scaled, and PCA-transformed data.
    """

    # genes to keep for HVG and PCA
    whitelist = [x for x in adQ.var_names if x not in blacklist]

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

    # Optionally scale the expression values of highly variable genes
    if gene_scale:
        sc.pp.scale(adata, max_value=scale_max)

    adata.var.loc[blacklist, 'highly_variable'] = False

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





