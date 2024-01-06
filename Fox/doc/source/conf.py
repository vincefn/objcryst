# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Fox'
copyright = '2024, Vincent Favre-Nicolin'
author = 'Vincent Favre-Nicolin'
release = '2024.x'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_logo = "_static/Fox.png"
html_css_files = ['css/custom.css']

# Not sure what all these actually do
html_theme_options = {
    "show_nav_level": 2,
    "navigation_depth": 1,
    "navbar_align": "left",
    # "primary_sidebar_end": ["indices.html", "sidebar-ethical-ads.html"]
}

# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-html_sidebars
html_sidebars = {
    "**": ["globaltoc.html", "sidebar-nav-bs"],
    # "**": ["localtoc.html"],
    # "**": ["sidebar-nav-bs"],
    # "<page_pattern>": ["index", "manual-intro", "tutorials", "manual"]
}
