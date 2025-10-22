# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
import datetime
# sys.path.insert(0, os.path.abspath('.'))


# Configuration file for the Sphinx documentation builder.

# -- Project information --------------------------------

project = 'chemdiskpy'
copyright = '{0}, {1}'.format(datetime.datetime.now().year, 'Sacha Gavino')
authors = 'Sacha Gavino'

# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- General configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.graphviz",
    #"autoapi.extension",
    # custom extentions
    "_extension.gallery_directive",
    "_extension.component_directive",
    # For extension examples and demos
    #"myst_parser",
    #"ablog",
    "jupyter_sphinx",
    #"sphinxcontrib.youtube",
    #"nbsphinx",
    "numpydoc",
    #"sphinx_togglebutton",
    "jupyterlite_sphinx",
    "sphinx_favicon",
]


jupyterlite_config = "jupyterlite_config.json"



intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']


exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# specifying the natural language populates some key tags
language = "en"

# -- Options for LaTeX output 

latex_elements = {
    'sphinxsetup': 'VerbatimColor={rgb}{0.95,0.95,0.95}',
}


# -- Options for HTML output
html_theme = "pydata_sphinx_theme"
#html_theme = "furo"
#html_theme_path = ["_themes", ]


html_sidebars = {
    "community/index": [
        "sidebar-nav-bs",
        "custom-template",
    ],  # This ensures we test for custom sidebars
    "examples/no-sidebar": [],  # Test what page looks like with no sidebar items
    "examples/persistent-search-field": ["search-field"],

}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for EPUB output
epub_show_urls = 'footnote'