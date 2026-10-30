import pathlib

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'dgcode_projects'
copyright = 'ESS ERIC and the dgcode_projects developers'
author = 'The dgcode_projects developers'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

nitpicky = True
extensions = [
    'sphinx_toolbox.sidebar_links',
    'sphinx_toolbox.github',
    ]

def _reporoot():
    p=pathlib.Path(__file__).parent.parent.parent
    assert ( p / 'simplebuild.cfg' ).is_file()
    return p.absolute()

# for :sbpkg:`MyPkg` links
extensions += '_simple_build_system._sphinxext',
def _sbbundles():
    import simplebuild_dgcode
    sbdgdata = pathlib.Path(simplebuild_dgcode.__file__).parent / 'data'
    import _simple_build_system as _
    sbdata = pathlib.Path(_.__file__).parent / 'data'
    sbdgdata_online = 'https://github.com/mctools/simplebuild-dgcode/[blob|tree]/main/src/simplebuild_dgcode/data'
    sbdata_online = 'https://github.com/mctools/simplebuild/[blob|tree]/main/src/_simple_build_system/data'

    bundles = { 'dgcode_projects' : ( _reporoot(),
                                      'https://github.com/ess-dg/dgcode_projects/[blob|tree]/main' ),
                'dgcode' : ( sbdgdata / 'pkgs',
                             sbdgdata_online + '/pkgs' ),
                'core' :  ( sbdata / 'pkgs-core',
                            sbdata_online + '/pkgs-core' ),
               }
    return bundles
simplebuild_pkgbundles = _sbbundles()

def setup(app):
    #Let table cells wrap (the theme keeps them on one line):
    app.add_css_file('custom.css')

#For the toolbox.sidebar_links:
github_username = 'ess-dg'
github_repository = 'dgcode_projects'

templates_path = []
exclude_patterns = []

#Links that can not be checked from outside ESS (internal hosts, private
#repositories, shares needing a login):
linkcheck_ignore = [ r'https://confluence\.esss\.lu\.se/.*',
                     r'https://project\.esss\.dk/owncloud/.*',
                     r'https://github\.com/ess-loki-sim/.*' ]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

github_url = 'https://github.com/ess-dg/dgcode_projects'

html_theme_options = {
    'logo_only': False,
    'collapse_navigation' : False,
    'sticky_navigation': True, #sidebar stays in place while contents scrolls
    'navigation_with_keys': True, #navigate with left/right arrow keys
    'navigation_depth': 4,
}
