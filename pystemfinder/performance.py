"""Benchmark stemFinder scores against ground truth (port of the R package's evaluation functions)."""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, spearmanr


def _auc(positive, negative):
    """Probability that a positive outranks a negative, with half credit for ties (Mann-Whitney U / n_pos n_neg)."""
    ranks = rankdata(np.concatenate([positive, negative]))
    n_pos, n_neg = len(positive), len(negative)
    return (ranks[:n_pos].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def _metrics(scores, ground_truth, phenotype, pheno_method):
    single_cell = spearmanr(scores, ground_truth)[0]

    means = pd.DataFrame({'score': scores, 'ground_truth': ground_truth, 'phenotype': phenotype})
    means = means.groupby('phenotype', observed=True).mean()
    correlate = spearmanr if pheno_method == 'spearman' else pearsonr
    pheno = correlate(means['score'], means['ground_truth'])[0]

    most = ground_truth == ground_truth.max()
    least = ground_truth == ground_truth.min()
    return [single_cell, pheno, _auc(scores[most], scores[least])]


def compute_performance(adata, score_key='stemfinder', ground_truth_key='Ground_truth', phenotype_key='Phenotype',
                        competitor_key=None, competitor_inverted=False, pheno_method='spearman'):
    """Quantify how well a score recovers ground truth differentiation on one dataset.

    Port of ``compute_performance_single`` from the R package. Scores should be oriented like pseudotime (lower
    means less differentiated), as in the ``'stemfinder'`` column written by :func:`stemfinder`, and ground truth
    should increase with differentiation. Three metrics are computed:

    - ``Spearman_SingleCell``: Spearman correlation between scores and ground truth across cells.
    - ``Spearman_Pheno``: correlation between the mean score and the mean ground truth of each phenotype.
    - ``AUC``: probability that a most differentiated cell (maximum ground truth) scores higher than a least
      differentiated cell (minimum ground truth), with half credit for ties.

    Two results differ from the R function by design. R's ``auc_probability`` compares only the first most
    differentiated cell with the least differentiated ones, while this AUC uses every pair. R's "phenotypic
    Spearman" correlation is computed with ``cor.test``'s default, Pearson; pass ``pheno_method='pearson'`` to
    reproduce it.

    Args:
        adata (anndata.AnnData): Annotated data with the score, ground truth, and phenotype columns in ``.obs``.
        score_key (str, optional): ``.obs`` column with the scores. Defaults to ``'stemfinder'``.
        ground_truth_key (str, optional): ``.obs`` column with numeric ground truth. Defaults to
            ``'Ground_truth'``.
        phenotype_key (str, optional): ``.obs`` column with cell type annotations. Defaults to ``'Phenotype'``.
        competitor_key (str, optional): ``.obs`` column with another method's scores to evaluate the same way.
            Defaults to None.
        competitor_inverted (bool, optional): Whether the competitor scores are already oriented like pseudotime.
            If False, they are inverted as ``1 - x / max(x)``. Defaults to False.
        pheno_method (str, optional): ``'spearman'`` or ``'pearson'`` for the phenotype-level correlation.
            Defaults to ``'spearman'``.

    Returns:
        pandas.DataFrame: One row per method (``score_key``, and ``competitor_key`` if given) with columns
        ``Spearman_SingleCell``, ``Spearman_Pheno``, and ``AUC``.
    """
    if pheno_method not in ['spearman', 'pearson']:
        raise ValueError("Invalid value for 'pheno_method'. Expected 'spearman' or 'pearson'.")

    obs = adata.obs
    ground_truth = obs[ground_truth_key].to_numpy(dtype=float)
    phenotype = obs[phenotype_key].to_numpy()

    rows = {score_key: _metrics(obs[score_key].to_numpy(dtype=float), ground_truth, phenotype, pheno_method)}
    if competitor_key is not None:
        competitor = obs[competitor_key].to_numpy(dtype=float)
        if not competitor_inverted:
            competitor = 1 - competitor / competitor.max()
        rows[competitor_key] = _metrics(competitor, ground_truth, phenotype, pheno_method)

    return pd.DataFrame.from_dict(rows, orient='index', columns=['Spearman_SingleCell', 'Spearman_Pheno', 'AUC'])


def pct_recover(adata, score_key='stemfinder', ground_truth_key='Ground_truth'):
    """Percentage of the least differentiated cells that a score ranks among the least differentiated.

    Port of ``pct_recover`` from the R package: the percentage of cells with the minimum ground truth whose score
    is below the 1 / (number of ground truth levels) quantile of all scores.

    Args:
        adata (anndata.AnnData): Annotated data with the score and ground truth columns in ``.obs``.
        score_key (str, optional): ``.obs`` column with pseudotime-oriented scores. Defaults to ``'stemfinder'``.
        ground_truth_key (str, optional): ``.obs`` column with numeric ground truth. Defaults to
            ``'Ground_truth'``.

    Returns:
        float: Percentage (0-100) of least differentiated cells recovered.
    """
    scores = adata.obs[score_key].to_numpy(dtype=float)
    ground_truth = adata.obs[ground_truth_key].to_numpy(dtype=float)
    threshold = np.quantile(scores, 1 / len(np.unique(ground_truth)))
    least = ground_truth == ground_truth.min()
    return 100 * np.mean(scores[least] < threshold)
