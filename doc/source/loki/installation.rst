Installation
============

Geant4 simulation (dgcode)
--------------------------

1. Install dgcode as described in its
   `documentation <https://mctools.github.io/simplebuild-dgcode/install.html>`__: download ``conda_dgcode.yml`` from
   there and create the ``dgcode`` conda environment:

   .. code:: sh

      conda env create -f conda_dgcode.yml

2. Clone `dgcode_projects <https://github.com/ess-dg/dgcode_projects>`__ and build it:

   .. code:: sh

      git clone https://github.com/ess-dg/dgcode_projects.git
      cd dgcode_projects
      conda activate dgcode
      sb

3. Set up the shell, so that the commands of the packages are found, as ``sb_<package>_<script>`` (e.g.
   ``sb_loki_sim``, the ``sim`` script of the ``LOKI`` package):

   .. code:: sh

      eval "$(sb --env-setup)"       # in the dgcode_projects directory
      sb_loki_sim -h                 # check

   Do this in every new terminal, after ``conda activate dgcode``. Without it, ``sb_loki_sim`` gives
   ``command not found``. (Alternatively, run single commands as ``sbenv sb_loki_sim ...``.)

Run ``sb`` again after changing C++ code; changes of Python scripts need no rebuild. A first simulation is in
:ref:`loki-simulation-a-first-simulation`.

Useful ``sb`` options (`sb documentation <https://mctools.github.io/simplebuild/cmdline.html#the-sb-command>`__):

.. list-table::
   :header-rows: 1

   * - Command
     - Effect
   * - ``sb --help``
     - All options.
   * - ``sb -t``
     - Build and run the tests (see :doc:`testing`).
   * - ``sb -i`` (``--insist``)
     - Reconfigure and rebuild everything.
   * - ``sb -c`` (``--clean``)
     - Remove the build and install directories.
   * - ``sb -s`` (``--show``)
     - Show the configuration.

The simulation commands and the Python tools have a ``-h`` option listing their options (the analysis programs
``sb_loki_ana`` and ``sb_lokimasking_ana`` take only GRIFF files; ``sb_loki_ana`` also ``--nodetfile``).

.. _loki-installation-loki-geometry-optional:

loki-geometry (optional)
------------------------

`loki-geometry <https://github.com/MilanKlausz/loki-geometry>`__ is a separate Python package (Python ≥ 3.9 and numpy)
that computes the surveyed bank calibrations (the ``G4GeoLoki/data`` files) and the pixel geometry of the instrument's
NeXus-JSON file (``loki-pixel-json``, used in :ref:`loki-geometry-verification`). It is needed only for those tasks, and can
be installed in its own environment (the Geant4 tools only read its output files):

.. code:: sh

   git clone https://github.com/MilanKlausz/loki-geometry.git
   python3 -m venv lg-venv && lg-venv/bin/pip install ./loki-geometry

McStas
------

The McStas models of the instruments are in separate repositories (see :doc:`mcstas`):

-  LoKI: `LOKI <https://git.esss.dk/dmsc-instrumentmodels/LOKI>`__ (DMSC instrument models), branch ``geant4`` for the
   MCPL output for Geant4. It needs McStas 3.9.1 or later; the repository has its conda environment:

   .. code:: sh

      git clone -b geant4 https://git.esss.dk/dmsc-instrumentmodels/LOKI.git loki_mcstas
      conda env create -n loki -f loki_mcstas/conda.yml

-  Larmor (only for the Larmor tests): `larmor-mcstas <https://github.com/ess-loki-sim/larmor-mcstas>`__

The Larmor model and the earlier LoKI model
(`loki-mcstas-master <https://bitbucket.org/MKlausz/loki-mcstas-master/src/master/>`__, no longer developed) have
branches for McStas 2.x and 3.x (from `McCode <https://github.com/mccode-dev/McCode>`__; versions 2.7 and 3.x were
used), with and without the components for Geant4 (``git branch -a``):

.. list-table::
   :header-rows: 1

   * - Branch
     - McStas
     - Geant4 (MCPL) output
   * - ``master``
     - 2.x
     - no
   * - ``mcstas2.x-g4``
     - 2.x
     - yes
   * - ``mcstas3.x``
     - 3.x
     - no
   * - ``mcstas3.x-g4``
     - 3.x
     - yes

.. _loki-installation-mantid:

Mantid
------

The :sbpkg:`LokiMantid` package converts the detection events to NeXus files with Mantid (see :doc:`mantid`). It is
built by ``sb`` only if the command ``mantidpython`` is available. Install Mantid from
`mantidproject.org <https://www.mantidproject.org/installation/index>`__ (version 5.0 or newer), then make ``mantidpython``
available:

-  **Older Mantid versions** (e.g. 5.1) ship an IPython executable named ``mantidpython``. Find it and add its
   directory to ``PATH``:

   .. code:: sh

      find /Applications/Mantid* -name "mantidpython"
      export PATH="/Applications/MantidWorkbench.app/Contents/MacOS/:$PATH"

-  **Newer Mantid versions** (e.g. 6.7) have no ``mantidpython``: link the Python interpreter of Mantid under that
   name, and add the Mantid and dgcode_projects Python directories to ``PYTHONPATH``:

   .. code:: sh

      mkdir -p $HOME/mantidpythonsymlink
      ln -sf /Applications/MantidWorkbench.app/Contents/Resources/bin/python $HOME/mantidpythonsymlink/mantidpython
      export PATH="$HOME/mantidpythonsymlink:$PATH"
      export PYTHONPATH="$PYTHONPATH:/Applications/MantidWorkbench.app/Contents/Resources/bin:<path to dgcode_projects>/simplebuild_cache/install/python"

Put the commands in ``~/.bashrc`` (or ``~/.bash_profile``) to have them in every terminal, and run ``sb`` again to build
LokiMantid.

The conversion is also available as a Mantid algorithm, ``LoadLokiDetectionEvents``: in Mantid, add the
``dgcode_projects/simplebuild_cache/install/python/`` directory to the Python script directories (File → Manage User Directories →
Python Script Directories; without it: ``ModuleNotFoundError: No module named 'LokiMantid'``), then open and run
``dgcode_projects/LOKI/LokiMantid/python/lokiAlgorithm.py`` to register the algorithm.

.. _loki-installation-dmsc-cluster:

DMSC cluster
------------

.. note::

   This section is from the earlier LOKI simulation guide and has not been checked against the current cluster setup
   (node names, modules, operating systems).

Connect as described in the `ESS SSH guide <https://confluence.esss.lu.se/display/DCC/SSH+Best+Practices>`__. For the
quark nodes, add to ``$HOME/.ssh/config``:

.. code:: text

   Host quarkcompile
       User <user name>
       ProxyCommand ssh -W %h:%p login
       ForwardAgent yes
       ForwardX11 yes

and connect with ``ssh quarkcompile``.

**McStas**: no installation needed, load the modules:

.. list-table::
   :header-rows: 1

   * - Nodes
     - Modules
   * - quark (``ssh quarkcompile``)
     - ``module load mcstas/2.7.2 gcc/10.2.0 openmpi/3.0_gcc1020`` or ``module load mcstas/3.3 gcc/10.2.0 openmpi/4.0_gcc1020``
   * - newlong, short (``ssh compile``)
     - ``module load mcstas/2.7.1 intel/17.0 openmpi/3_0_mellanox_intel170``

**Geant4 (dgcode)**: on a quark node, once:

.. code:: sh

   module load miniconda3/4.8.3
   conda init          # then disconnect and connect again
   conda env create -f conda_dgcode.yml

then clone dgcode_projects and build it as above (``conda activate dgcode``, ``sb``). The environment also works on the
compile node. If both the compile and the quarkcompile nodes are used, keep one clone of dgcode_projects for each
(different operating systems). If git does not work on these nodes, run the git commands on the login node.

**Mantid**: available only on the quark nodes, as a module. Run ``conda init`` after loading the miniconda module (as
above), then in each session load Mantid first and activate the environment after it (the Mantid module changes
``PATH``, and would hide the CMake of the environment otherwise):

.. code:: sh

   module load mantid/6.1.2
   conda activate dgcode
   cd dgcode_projects
   sb                  # builds LokiMantid as well

A function in ``~/.bashrc`` can do this on both node types:

.. code:: sh

   setup_dgcode() {
       if module avail 2>&1 | grep -q "mantid/6.1.2"; then  # quark node
           module load mantid/6.1.2
       fi
       conda activate dgcode
   }
