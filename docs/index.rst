PyStemFinder
============

PyStemFinder estimates how far single cells have progressed through differentiation from scRNA-seq data.
Less differentiated cells vary more in their expression of cell cycle genes than their neighbors do, and
stemFinder scores that heterogeneity within each cell's k-nearest-neighbor neighborhood.

PyStemFinder is a Python port of the R package `stemFinder <https://github.com/CahanLab/stemfinder>`_ that works
on `AnnData <https://anndata.readthedocs.io>`_ objects and `scanpy <https://scanpy.readthedocs.io>`_ neighbor
graphs. Its scores match R's to within floating point precision; :doc:`notebooks/benchmarking` shows the comparison.

If you use stemFinder, please cite:

   Noller K, Cahan P. Cell cycle expression heterogeneity predicts degree of differentiation.
   *Briefings in Bioinformatics* 25(6):bbae536 (2024). https://doi.org/10.1093/bib/bbae536

.. toctree::
   :maxdepth: 2
   :caption: Contents

   notebooks/quickstart
   notebooks/benchmarking
   api
   changelog
