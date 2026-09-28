import os
import re

from setuptools import setup

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'PyStemFinder', '_version.py')) as f:
    version = re.search(r'__version__ = "(.+)"', f.read()).group(1)

setup(name='PyStemFinder',
      version=version,
      description='Tools to infer extent of differentiation and cell fate potential from single cell omics data',
      url='http://github.com/pcahan1/PyStemFinder/',
      author='Kathleen Noller, Patrick Cahan',
      author_email='patrick.cahan@jhmi.ed',
      license='MIT',
      packages=['PyStemFinder'],
      package_data={'PyStemFinder': ['data/*.txt']},
      install_requires=[
          'pandas',
          'numpy',
          'scipy',
          'matplotlib',
          'scanpy',
          'anndata',
          'seaborn',
      ],
)
