from setuptools import setup

setup(name='PyStemFinder',
      version='0.1',
      description='Tools to infer extent of differentiation and cell fate potential from single cell omics data',
      url='http://github.com/pcahan1/PyStemFinder/',
      author='Kathleen Noller, Patrick Cahan',
      author_email='patrick.cahan@jhmi.ed',
      license='MIT',
      packages=['PyStemFinder'],
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
