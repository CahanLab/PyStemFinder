"""Synthetic scRNA-seq data for exploring stemFinder scores."""
import anndata
import numpy as np
import pandas as pd


def simulate_data(n_cells=300, n_genes=50, lambda_val=2, dropout_rate=0.6, n_populations=3,
                  stochastic_gene_range=(40, 45), marker_expression_diff=5, stochastic_expression_diff=(0, 0, 0),
                  stochastic_mean=None, stochastic_dropout_by_population=None, random_state=42):
    """Generate a toy scRNA-seq count matrix with cell populations and stochastically expressed genes.

    Counts are drawn from a Poisson distribution with dropout. Each population (a consecutive block of cells)
    highly expresses its own 10 marker genes (genes 0-9 for the first population, 10-19 for the second, ...).
    The genes in ``stochastic_gene_range`` stand in for cell cycle genes: their mean and dropout rate can differ
    between populations, which makes their expression more or less heterogeneous.

    Args:
        n_cells (int, optional): Number of cells. Defaults to 300.
        n_genes (int, optional): Number of genes. Defaults to 50.
        lambda_val (float, optional): Mean of the Poisson distribution for background expression. Defaults to 2.
        dropout_rate (float, optional): Fraction of counts set to zero. Defaults to 0.6.
        n_populations (int, optional): Number of cell populations. Defaults to 3.
        stochastic_gene_range (tuple of int, optional): Start (inclusive) and end (exclusive) index of the
            stochastic genes. Defaults to (40, 45).
        marker_expression_diff (float, optional): How much higher the marker genes' mean is than
            ``lambda_val``. Defaults to 5.
        stochastic_expression_diff (tuple of float, optional): Per population, how much higher the stochastic
            genes' mean is than ``stochastic_mean``. Defaults to (0, 0, 0).
        stochastic_mean (float, optional): Base mean of the stochastic genes. Defaults to ``lambda_val``.
        stochastic_dropout_by_population (list of float, optional): Per population, the dropout rate of the
            stochastic genes. Defaults to ``dropout_rate`` for every population.
        random_state (int, optional): Seed for the random number generator. numpy's global random state is not
            touched. Defaults to 42.

    Returns:
        anndata.AnnData: Integer counts, cells named ``cell0...`` and genes ``gene0...``.
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
