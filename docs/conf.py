# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyStemFinder'
copyright = '2023-2026, Patrick Cahan, Kathleen Noller'
author = 'Patrick Cahan, Kathleen Noller'
_version_file = Path(__file__).parent.parent / 'pystemfinder' / '_version.py'
release = re.search(r'__version__ = "(.+)"', _version_file.read_text()).group(1)
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx_copybutton",
    "sphinx_inline_tabs",
    "sphinx_design",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon"
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable', None),
}
napoleon_preprocess_types = True
napoleon_type_aliases = {'array-like': ':term:`array-like <array_like>`'}

# Notebooks are stored with their outputs; do not re-execute them when building the docs
nb_execution_mode = 'off'

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']

html_theme_options = {
    "light_logo": "stemFinder_logo_light.png",
    "dark_logo": "stemFinder_logo_dark.png",
    "sidebar_hide_name": True
}
