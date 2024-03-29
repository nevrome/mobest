# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'mobest'
copyright = '2023, Clemens Schmid'
author = 'Clemens Schmid'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
  'myst_parser',
  'sphinx_rtd_theme',
  'sphinx_multiversion',
  'sphinxcontrib.bibtex',
  'sphinx.ext.autosectionlabel',
  'sphinx_toolbox.sidebar_links',
  'sphinx_toolbox.github'
]

myst_enable_extensions = [
  "dollarmath"
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

exclude_patterns = ['README.md']

# https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html
bibtex_bibfiles = ['references.bib']
bibtex_reference_style = 'author_year'
bibtex_default_style = 'plain'

# make sure the autosectionlabel targets are unique
autosectionlabel_prefix_document = True

# https://sphinx-toolbox.readthedocs.io/en/latest/extensions/sidebar_links.html
github_username = 'nevrome'
github_repository = 'mobest'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
