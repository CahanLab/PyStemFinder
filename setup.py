import os
import re

from setuptools import setup

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'PyStemFinder', '_version.py')) as f:
    version = re.search(r'__version__ = "(.+)"', f.read()).group(1)
with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(name='PyStemFinder',
      version=version,
      description='Tools to infer extent of differentiation and cell fate potential from single cell omics data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/pcahan1/PyStemFinder',
      author='Kathleen Noller, Patrick Cahan',
      author_email='patrick.cahan@jhmi.edu',
      license='MIT',
      packages=['PyStemFinder'],
      package_data={'PyStemFinder': ['data/*.txt']},
      python_requires='>=3.9',
      install_requires=[
          'pandas',
          'numpy',
          'scipy',
          'scanpy',
          'anndata',
      ],
)
