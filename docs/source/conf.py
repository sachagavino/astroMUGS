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
import os
# import sys
import datetime
# sys.path.insert(0, os.path.abspath('.'))


# Configuration file for the Sphinx documentation builder.

# -- Project information --------------------------------

project = 'chemdiskpy'
copyright = '{0}, {1}'.format(datetime.datetime.now().year, 'Sacha Gavino')
authors = 'Sacha Gavino'

# The full version, including alpha/beta/rc tags
#release = '1.0.0'

# -- General configuration

extensions = [
    #"sphinx.ext.napoleon",
    #"sphinx.ext.autodoc",
    #"sphinx.ext.autosummary",
    #"sphinx.ext.todo",
    #"sphinx.ext.viewcode",
    #"sphinx.ext.intersphinx",
    #"sphinx.ext.graphviz",
    #"autoapi.extension",
    # custom extentions
    #"_extension.gallery_directive",
    #"_extension.component_directive",
    # For extension examples and demos
    "myst_parser",
    #"ablog",
    #"jupyter_sphinx",
    #"sphinxcontrib.youtube",
    #"nbsphinx",
    #"numpydoc",
    #"sphinx_togglebutton",
    #"jupyterlite_sphinx",
    #"sphinx_favicon",
]


source_suffix = ['.rst', '.md']

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

html_logo = "_static/logo-dark.png"


# Define the version we use for matching in the version switcher.
version_match = os.environ.get("READTHEDOCS_VERSION")
release = '1.0.5'
# If READTHEDOCS_VERSION doesn't exist, we're not on RTD
# If it is an integer, we're in a PR build and the version isn't correct.
# If it's "latest" → change to "dev" (that's what we want the switcher to call it)
if not version_match or version_match.isdigit() or version_match == "latest":
    # For local development, infer the version to match from the package.
    if "dev" in release or "rc" in release:
        version_match = "dev"
        # We want to keep the relative reference if we are in dev mode
        # but we want the whole url if we are effectively in a released version
        json_url = "_static/switcher.json"
    else:
        version_match = f"v{release}"
elif version_match == "stable":
    version_match = f"v{release}"


html_theme_options = {
    "header_links_before_dropdown": 4,

    "logo": {
        "text": "astroMUGS",
        "image_dark": "_static/logo-dark.png",
    },
    "show_toc_level": 1,
    # [left, content, right] For testing that the navbar items align properly
    "navbar_align": "left",
    # "show_nav_level": 2,
    "show_version_warning_banner": True,
    "navbar_center": ["version-switcher", "navbar-nav"],
    # "navbar_start": ["navbar-logo"],
    # "navbar_end": ["theme-switcher", "navbar-icon-links"],
    # "navbar_persistent": ["search-button"],
    # "primary_sidebar_end": ["custom-template", "sidebar-ethical-ads"],
    # "article_footer_items": ["test", "test"],
    # "content_footer_items": ["test", "test"],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
    "secondary_sidebar_items": {
        "**/*": ["page-toc", "edit-this-page", "sourcelink"],
        "examples/no-sidebar": [],
    },
    # "back_to_top_button": False,
    "search_as_you_type": True,
}



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
#html_static_path = ['_static']

# -- Options for EPUB output
epub_show_urls = 'footnote'