# Configuration file for the Sphinx documentation of dgcode_projects.
#
# The site has a short main page (index.md) and one section per project. A project contributes its documentation
# as a doc/ directory in its own project directory (e.g. LOKI/doc/, starting page LOKI/doc/index.md), see
# contributing.md. At each build, the doc/ directories of the projects are copied here as projects/<Project> (not
# committed; edit the originals), and the table of the projects on the main page is generated.
#
# Build: make html (in the doc/ directory of the repository), see doc/README.md.
# Sphinx configuration reference: https://www.sphinx-doc.org/en/master/usage/configuration.html

import pathlib
import shutil

_SOURCE = pathlib.Path(__file__).parent.absolute()
_REPO = _SOURCE.parent.parent
_GITHUB = 'https://github.com/ess-dg/dgcode_projects'
_GITHUB_BRANCH = 'main'
# top level directories of the repository that are not projects
_NOT_PROJECTS = {'doc', 'install', 'simplebuild_cache', '.git', '.github', '.githooks'}

# -- Project information -----------------------------------------------------

project = 'dgcode_projects'
copyright = 'ESS ERIC and the dgcode_projects developers'
author = 'The dgcode_projects developers'

# -- The projects ------------------------------------------------------------


def _project_dirs():
    """The project directories of the repository: top level directories with packages (pkg.info files), in git
    (local directories that are not in the repository are left out), or all of them without git."""
    import subprocess
    try:
        tracked = subprocess.run(['git', 'ls-files', '*/pkg.info'], cwd=_REPO, capture_output=True, text=True,
                                 check=True).stdout.split()
        names = {pathlib.PurePosixPath(p).parts[0] for p in tracked}
    except (OSError, subprocess.CalledProcessError):
        names = {p.relative_to(_REPO).parts[0] for p in _REPO.glob('*/**/pkg.info')}
    return [_REPO / n for n in sorted(names, key=str.lower) if n not in _NOT_PROJECTS and not n.startswith('.')]


def _doc_index(project_dir):
    for name in ('index.md', 'index.rst'):
        if (project_dir / 'doc' / name).is_file():
            return project_dir / 'doc' / name
    return None


def _copy_project_docs():
    """Copies <Project>/doc to projects/<Project> for each project with documentation (the file times are kept, so
    Sphinx rebuilds only the changed pages), and writes the table of the projects (_generated/projects_table.md)."""
    copies = _SOURCE / 'projects'
    shutil.rmtree(copies, ignore_errors=True)
    copies.mkdir()
    rows = []
    for project_dir in _project_dirs():
        index = _doc_index(project_dir)
        github = f'{_GITHUB}/tree/{_GITHUB_BRANCH}/{project_dir.name}'
        if index is None:
            rows.append(f'| {project_dir.name} | (no documentation yet) | [source]({github}) |')
            continue
        shutil.copytree(project_dir / 'doc', copies / project_dir.name)
        rows.append(f'| {project_dir.name} | [documentation](projects/{project_dir.name}/{index.stem}) | [source]({github}) |')
    generated = _SOURCE / '_generated'
    generated.mkdir(exist_ok=True)
    table = ['| Project | Documentation | Code |', '| --- | --- | --- |'] + rows
    (generated / 'projects_table.md').write_text('\n'.join(table) + '\n')


_copy_project_docs()

# -- General configuration ---------------------------------------------------

extensions = ['myst_parser']
try:
    # :sbpkg:`G4GeoLoki` links to the code of a package (from simplebuild, as in the dgcode documentation)
    import _simple_build_system._sphinxext  # noqa: F401
    extensions.append('_simple_build_system._sphinxext')
    simplebuild_pkgbundles = {'dgcode_projects': (_REPO, f'{_GITHUB}/[blob|tree]/{_GITHUB_BRANCH}')}
except ImportError:
    print('WARNING: simple-build-system not found: the :sbpkg: links will not work (use the dgcode environment)')

source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
myst_enable_extensions = ['colon_fence', 'deflist', 'dollarmath', 'attrs_inline']
myst_heading_anchors = 3
templates_path = ['_templates']
exclude_patterns = ['build', '_generated', 'projects/*/figures']
nitpicky = True

# -- HTML output -------------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
html_title = 'dgcode_projects documentation'
html_context = {
    'display_github': True,
    'github_user': 'ess-dg',
    'github_repo': 'dgcode_projects',
    'github_version': _GITHUB_BRANCH,
    'conf_py_path': '/doc/source/',
}
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_with_keys': True,
    'navigation_depth': 4,
}


def _github_url_of_project_pages(app, pagename, templatename, context, doctree):
    """'Edit on GitHub' of the pages of a project: the file in <Project>/doc/, not its copy in projects/."""
    parts = pagename.split('/')
    if len(parts) >= 3 and parts[0] == 'projects':
        suffix = context.get('page_source_suffix', '.md')
        path = '/'.join([parts[1], 'doc'] + parts[2:])
        context['meta'] = dict(context.get('meta') or {},
                               github_url=f'{_GITHUB}/blob/{_GITHUB_BRANCH}/{path}{suffix}')


def setup(app):
    app.connect('html-page-context', _github_url_of_project_pages)
