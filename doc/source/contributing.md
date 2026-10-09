# Documenting a project

Every project documents itself in a `doc/` directory of its own project directory, next to its packages:

```text
dgcode_projects/
├── doc/                  the site: main page, this page, configuration (doc/source/conf.py)
├── LOKI/                 a project
│   ├── doc/              its documentation
│   │   ├── index.md      its starting page (required)
│   │   ├── geometry.md   further pages, listed in a toctree of index.md
│   │   └── images/       its figures
│   ├── G4GeoLoki/        its packages
│   └── ...
└── <YourProject>/
    └── doc/index.md
```

The site finds the documentation by itself: a project with a `doc/index.md` (or `doc/index.rst`) appears on the main
page and in the menu; the other projects are listed without a link. Nothing outside the project directory needs to be
changed.

## Writing the pages

- Pages are Markdown ([MyST](https://myst-parser.readthedocs.io/en/latest/), `.md`) or reStructuredText (`.rst`).
- The first heading of `index.md` is the title of the project in the menu. Further pages are listed in a `toctree`
  in `index.md`:

  ````md
  ```{toctree}
  :maxdepth: 2

  installation
  geometry
  ```
  ````

- Link to the code of a package with `` {sbpkg}`G4GeoLoki` `` or to a file in it with
  `` {sbpkg}`G4GeoLoki/libinc/BcsBanks.hh` `` (links to GitHub, checked at the build).
- Figures go in the project's `doc/images/`. The repository rejects binary files above 20 kB (text files: 200 kB):
  prefer SVG for drawings, and WebP or JPEG for screenshots. A figure made by a script (e.g.
  `LOKI/doc/figures/make_numbering_figure.py`) can be regenerated when the code changes.
- Write for users of the project: what it does, how to run it, what the options mean. Keep commands copyable and
  check them against the current code.

## Building the documentation

See `doc/README.md`: in the dgcode environment, `python3 -m pip install -r doc/requirements.txt` once, then
`make html` in `doc/`, and open `doc/build/html/index.html`. The build must finish without warnings
(`make html SPHINXOPTS=-W`).
