API reference
=============

.. currentmodule:: PyStemFinder

Scoring
-------

.. autofunction:: run_stemFinder
.. autofunction:: diffOmeter
.. autofunction:: gene_set_score
.. autofunction:: count_high_expr_genes

Benchmarking
------------

.. autofunction:: compute_performance_single
.. autofunction:: pct_recover

Gene lists
----------

.. autofunction:: cell_cycle_genes
.. autofunction:: transcription_factors

Preprocessing and toy data
--------------------------

.. autofunction:: sf_norm_hvg_scale_pca
.. autofunction:: generate_scRNAseq_test_data
.. autofunction:: binarize_data
.. autofunction:: binary_gini_impurity
