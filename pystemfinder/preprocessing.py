"""Preprocessing for stemFinder."""
import scanpy as sc


def recipe_stemfinder(adata, exclude=(), target_sum=1e4, min_mean=0.0125, max_mean=6, min_disp=0.25, scale=True,
                      max_value=10, n_comps=100, copy=False):
    """Normalize, log-transform, scale, and run PCA on raw counts, ready for a kNN graph and :func:`stemfinder`.

    Runs ``sc.pp.normalize_total``, ``sc.pp.log1p``, ``sc.pp.highly_variable_genes``, and ``sc.pp.scale`` (on every
    gene, because the stemFinder markers are usually not highly variable), then ``sc.tl.pca`` on the highly variable
    genes. Genes in ``exclude`` are removed from the highly variable genes, so that, as in the R vignette, the
    markers do not shape the PCA and the kNN graph they are scored on.

    Args:
        adata (anndata.AnnData): Raw counts (cells x genes).
        exclude (list of str, optional): Genes to leave out of the highly variable genes, typically the markers
            passed to :func:`stemfinder`. Genes absent from the data are ignored. Defaults to none.
        target_sum (float, optional): Total count each cell is normalized to. Defaults to 1e4.
        min_mean (float, optional): Minimum mean expression of highly variable genes. Defaults to 0.0125.
        max_mean (float, optional): Maximum mean expression of highly variable genes. Defaults to 6.
        min_disp (float, optional): Minimum dispersion of highly variable genes. Defaults to 0.25.
        scale (bool, optional): Whether to scale every gene to zero mean and unit variance. Defaults to True.
        max_value (float, optional): Scaled values are clipped at this value. Defaults to 10.
        n_comps (int, optional): Number of principal components. Defaults to 100.
        copy (bool, optional): Return a processed copy instead of modifying ``adata``. Defaults to False.

    Returns:
        anndata.AnnData or None: The processed copy if ``copy`` is True, otherwise None.
    """
    adata = adata.copy() if copy else adata

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    adata.var.loc[adata.var_names.isin(list(exclude)), 'highly_variable'] = False
    if scale:
        sc.pp.scale(adata, max_value=max_value)
    sc.tl.pca(adata, n_comps=n_comps)

    return adata if copy else None
