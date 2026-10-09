# Installation

## Geant4 simulation (dgcode)

1. Install dgcode as described in its
   [documentation](https://mctools.github.io/simplebuild-dgcode/install.html): create the `dgcode` conda environment
   and activate it in every new terminal:

   ```sh
   conda env create -f conda_dgcode.yml
   conda activate dgcode
   ```

2. Clone [dgcode_projects](https://github.com/ess-dg/dgcode_projects) and build it:

   ```sh
   git clone https://github.com/ess-dg/dgcode_projects.git
   cd dgcode_projects
   sb
   ```

   Run `sb` again after changing C++ code; changes of Python scripts need no rebuild. The commands of the packages
   are then available as `sb_<package>_<script>`, e.g. `sb_loki_sim` (the `sim` script of the `LOKI` package).

Useful `sb` options ([sb documentation](https://mctools.github.io/simplebuild/cmdline.html#the-sb-command)):

| Command | Effect |
| --- | --- |
| `sb --help` | All options. |
| `sb -t` | Build and run the tests (see [](testing.md)). |
| `sb -i` (`--insist`) | Reconfigure and rebuild everything. |
| `sb -c` (`--clean`) | Remove the build and install directories. |
| `sb -s` (`--show`) | Show the configuration. |

Every simulation and analysis command has a `-h` option listing its options.

## McStas

The McStas models of the instruments are in separate repositories (see [](mcstas.md)):

- LoKI: [loki-mcstas-master](https://bitbucket.org/MKlausz/loki-mcstas-master/src/master/)
- Larmor (only for the Larmor tests): [larmor-mcstas](https://github.com/ess-loki-sim/larmor-mcstas)

Install McStas from [McCode](https://github.com/mccode-dev/McCode) (versions 2.7 and 3.x were used). The model
repositories have branches for McStas 2.x and 3.x, with and without the components for Geant4 (`git branch -a`):

| Branch | McStas | Geant4 (MCPL) output |
| --- | --- | --- |
| `master` | 2.x | no |
| `mcstas2.x-g4` | 2.x | yes |
| `mcstas3.x` | 3.x | no |
| `mcstas3.x-g4` | 3.x | yes |

## Mantid

The {sbpkg}`LokiMantid` package converts the detection events to NeXus files with Mantid (see [](mantid.md)). It is
built by `sb` only if the command `mantidpython` is available. Install Mantid from
[mantidproject.org](https://www.mantidproject.org/installation/index) (version 5.0 or newer), then make `mantidpython`
available:

- **Older Mantid versions** (e.g. 5.1) ship an IPython executable named `mantidpython`. Find it and add its
  directory to `PATH`:

  ```sh
  find /Applications/Mantid* -name "mantidpython"
  export PATH="/Applications/MantidWorkbench.app/Contents/MacOS/:$PATH"
  ```

- **Newer Mantid versions** (e.g. 6.7) have no `mantidpython`: link the Python interpreter of Mantid under that
  name, and add the Mantid and dgcode_projects Python directories to `PYTHONPATH`:

  ```sh
  mkdir -p $HOME/mantidpythonsymlink
  ln -sf /Applications/MantidWorkbench.app/Contents/Resources/bin/python $HOME/mantidpythonsymlink/mantidpython
  export PATH="$HOME/mantidpythonsymlink:$PATH"
  export PYTHONPATH="$PYTHONPATH:/Applications/MantidWorkbench.app/Contents/Resources/bin:<path to dgcode_projects>/install/python"
  ```

Put the commands in `~/.bashrc` (or `~/.bash_profile`) to have them in every terminal, and run `sb` again to build
LokiMantid.

The conversion is also available as a Mantid algorithm, `LoadLokiDetectionEvents`: in Mantid, add the
`dgcode_projects/install/python/` directory to the Python script directories (File → Manage User Directories →
Python Script Directories; without it: `ModuleNotFoundError: No module named 'LokiMantid'`), then open and run
`dgcode_projects/LOKI/LokiMantid/python/lokiAlgorithm.py` to register the algorithm.

## DMSC cluster

:::{note}
This section is from the earlier LOKI simulation guide and has not been checked against the current cluster setup
(node names, modules, operating systems).
:::

Connect as described in the [ESS SSH guide](https://confluence.esss.lu.se/display/DCC/SSH+Best+Practices). For the
quark nodes, add to `$HOME/.ssh/config`:

```text
Host quarkcompile
    User <user name>
    ProxyCommand ssh -W %h:%p login
    ForwardAgent yes
    ForwardX11 yes
```

and connect with `ssh quarkcompile`.

**McStas**: no installation needed, load the modules:

| Nodes | Modules |
| --- | --- |
| quark (`ssh quarkcompile`) | `module load mcstas/2.7.2 gcc/10.2.0 openmpi/3.0_gcc1020` or `module load mcstas/3.3 gcc/10.2.0 openmpi/4.0_gcc1020` |
| newlong, short (`ssh compile`) | `module load mcstas/2.7.1 intel/17.0 openmpi/3_0_mellanox_intel170` |

**Geant4 (dgcode)**: on a quark node, once:

```sh
module load miniconda3/4.8.3
conda init          # then disconnect and connect again
conda env create -f conda_dgcode.yml
```

then clone dgcode_projects and build it as above (`conda activate dgcode`, `sb`). The environment also works on the
compile node. If both the compile and the quarkcompile nodes are used, keep one clone of dgcode_projects for each
(different operating systems). If git does not work on these nodes, run the git commands on the login node.

**Mantid**: available only on the quark nodes, as a module. Run `conda init` after loading the miniconda module (as
above), then in each session load Mantid first and activate the environment after it (the Mantid module changes
`PATH`, and would hide the CMake of the environment otherwise):

```sh
module load mantid/6.1.2
conda activate dgcode
cd dgcode_projects
sb                  # builds LokiMantid as well
```

A function in `~/.bashrc` can do this on both node types:

```sh
setup_dgcode() {
    if module avail 2>&1 | grep -q "mantid/6.1.2"; then  # quark node
        module load mantid/6.1.2
    fi
    conda activate dgcode
}
```
