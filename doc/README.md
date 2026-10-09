# dgcode_projects documentation

The documentation website of dgcode_projects (Sphinx, like the
[simplebuild](https://mctools.github.io/simplebuild/) and [dgcode](https://mctools.github.io/simplebuild-dgcode/)
documentation). This directory holds the site: the main page and the configuration (`source/`). The projects keep
their documentation in their own directories (`<Project>/doc/`), see `source/contributing.md`.

Build (in the dgcode conda environment, or any Python 3.9+ environment):

```sh
conda activate dgcode
python3 -m pip install -r doc/requirements.txt   # once
cd doc
make html                  # or: make html SPHINXOPTS=-W  (warnings are errors)
open build/html/index.html
```

`make clean` removes the build. `make linkcheck` checks the external links.
