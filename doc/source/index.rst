*********************************************
Welcome to the dgcode_projects documentation!
*********************************************

`dgcode_projects <https://github.com/ess-dg/dgcode_projects>`__ is a collection of simulation projects of the
Detector Group at the European Spallation Source ERIC (ESS) and its collaborators. The projects are `simplebuild
<https://mctools.github.io/simplebuild>`__ packages built on `dgcode <https://mctools.github.io/simplebuild-dgcode>`__,
the Geant4-based simulation framework of the ESS Detector Group. Install dgcode as described in its `documentation
<https://mctools.github.io/simplebuild-dgcode/install.html>`__, then build all projects by running ``sb`` in a clone of
the repository:

.. code-block:: sh

   git clone https://github.com/ess-dg/dgcode_projects.git
   cd dgcode_projects
   conda activate dgcode
   sb
   eval "$(sb --env-setup)"    # in every new terminal: makes the sb_<package>_<script> commands available

The documented projects are listed in the menu. A project adds its pages in a directory of its own under ``doc/source/``
(e.g. ``doc/source/loki/``), and its starting page to the table of contents below.

.. toctree::
   :hidden:
   :caption: Projects

   loki/index

.. sidebar-links::
   :github:

   Issue tracker <https://github.com/ess-dg/dgcode_projects/issues>
   dgcode <https://mctools.github.io/simplebuild-dgcode>
   simplebuild <https://mctools.github.io/simplebuild>
