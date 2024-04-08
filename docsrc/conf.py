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
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'autoDEER'
copyright = '2021-2024, Hugo Karas, Gunnar Jeschke, Stefan Stoll'
author = 'Hugo Karas'

# The full version, including alpha/beta/rc tags
release = '0.7'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
            #   'sphinx.ext.autodoc',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx',
              #'sphinx.ext.autosummary',
              'autoapi.extension',
              'sphinx_toolbox.collapse',
              'sphinx_toolbox.code',
              'sphinx_copybutton',
              'numpydoc',
              'sphinx.ext.graphviz',
              'myst_parser']

# Add any paths that contain templates here, relative to this directory.
# Configuration of Sphinx-Autosymmary
# ----------------------------------------------------------------------
add_module_names = True
# Turn on sphinx.ext.autosummary
autoapi_dirs = ['../autodeer']
autodoc_typehints = "description"
autoapi_template_dir = "_templates/autoapi"
# autoapi_options = [
#     "members",
#     "undoc-members",
#     "show-inheritance",
#     # "show-module-summary",
#     "imported-members",
# ]
autoapi_keep_files = True
autoapi_add_toctree_entry = False
autoapi_python_class_content= "both"
# suppress_warnings = ["autoapi"]
autoapi_python_use_implicit_namespaces = True
autoapi_own_page_level = 'class'


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**_old**']
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
intersphinx_mapping = {
    "deerlab": ("https://jeschkelab.github.io/DeerLab/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}
intersphinx_disabled_reftypes = ["*"]
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
html_title = " "
html_static_path = ["_static"]
html_theme_options = {
    # "announcement": "Version 0.7 out now!",
    "light_logo": "autoDEER_EPR_light.svg",
    "dark_logo": "autoDEER_EPR_dark.svg",

}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
# html_logo = "autoDEER_name_purple.svg"

