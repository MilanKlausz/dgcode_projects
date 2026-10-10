# dgcode_projects

[dgcode_projects](https://github.com/ess-dg/dgcode_projects) is a collection of simulation projects of the Detector
Group at the European Spallation Source ERIC (ESS) and its collaborators. The projects are
[simplebuild](https://mctools.github.io/simplebuild/) packages built on
[dgcode](https://mctools.github.io/simplebuild-dgcode/), the Geant4-based simulation framework of the ESS Detector
Group. Install dgcode as described in its [documentation](https://mctools.github.io/simplebuild-dgcode/install.html),
then build all projects by running `sb` in a clone of the repository:

```sh
git clone https://github.com/ess-dg/dgcode_projects.git
cd dgcode_projects
conda activate dgcode
sb
eval "$(sb --env-setup)"    # in every new terminal: makes the sb_<package>_<script> commands available
```

## Projects

Each project is a directory of the repository and documents itself (see [](contributing.md)).

```{include} _generated/projects_table.md
```

```{toctree}
:hidden:
:glob:
:caption: Projects

projects/*/index
```

```{toctree}
:hidden:
:caption: About

contributing
```
