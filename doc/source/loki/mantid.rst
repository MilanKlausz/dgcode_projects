Mantid processing
=================

.. note::

   This page is from the earlier LOKI simulation guide. LokiMantid is not built without Mantid, and its tests do not run
   in the usual test setup; the options below follow the current code.

The :sbpkg:`LokiMantid` package converts the detection file of the analysis (``detectionEvents.mcpl``) to NeXus files
that Mantid loads for the data reduction: a detector file with the time-of-flight spectra of the pixels, and a
monitor file. See :ref:`loki-installation-mantid` for the setup.

.. important::

   LokiMantid does not use the weights of the detection events: every event counts as one neutron. Process only
   detection files of unit-weight events, e.g. of a McStas MCPL file written with ``mcpl_min_weight=1
   mcpl_max_weight=1`` (see :ref:`loki-mcstas-loki-model`) or of a flood source simulation; with weighted events (the
   default MCPL files of the LoKI McStas model) the result is wrong.

.. important::

   The instrument definition files of LokiMantid are nominal: process only simulations made with
   ``bank_calibration=nominal-geant4-geometry``. For a detection file with another bank calibration (recorded in its
   metadata) the processing prints a warning: its pixel positions in Mantid would be wrong (by up to ~33 mm for
   2026-September-SAM-606). Processing calibrated geometries would need instrument definition files made from the
   calibration (see ``checkBankCalibration`` in :sbpkg:`LokiMantid/python/workspaceCreator.py`).

   This concerns only simulations: for measured data Mantid takes the pixel geometry from the instrument's NeXus file,
   which carries the surveyed geometry once its JSON template is updated with loki-geometry.

Contents
--------

.. list-table::
   :header-rows: 1

   * - File
     - Contents
   * - :sbpkg:`LokiMantid/scripts/loki`
     - The command line script (``sb_lokimantid_loki``).
   * - :sbpkg:`LokiMantid/python/lokiAlgorithm.py`
     - ``LoadLokiDetectionEvents``: the processing, also usable as a Mantid algorithm.
   * - :sbpkg:`LokiMantid/python/paramHandler.py`
     - The parameters: defaults < metadata of the detection file < command line.
   * - :sbpkg:`LokiMantid/python/workspaceCreator.py`
     - Instrument definition files from templates, and the workspaces.
   * - :sbpkg:`LokiMantid/python/monitorSpectrum.py`
     - McStas monitor files, and the calculated pre-sample monitor of flood source simulations.
   * - :sbpkg:`LokiMantid/python/mantidApi.py`
     - The calls of the Mantid Python API.
   * - :sbpkg:`LokiMantid/data`
     - The instrument definition templates (``LOKI_Definition_template_detectors.xml``, ``..._monitors.xml``).

The scripts run with the Python of Mantid (``mantidpython``), as the Mantid Python package is built for it. Therefore
the C++ modules of G4GeoLoki (built for the dgcode Python) cannot be used in them. With the older-Mantid setup,
``mantidpython`` is IPython and the script arguments must follow ``--``, e.g. ``sb_lokimantid_loki -- -h``; the examples
below use the newer setup, without ``--`` (``sb_lokimantid_loki -h`` works with it, ``-- -h`` does not, which shows the
setup in use).

Processing
----------

.. code:: sh

   sb_lokimantid_loki -h
   sb_lokimantid_loki detectionEvents.mcpl --showMetadata
   sb_lokimantid_loki detectionEvents.mcpl -s testLokiMcStas --mcstasDir DIR    # DIR: the McStas output directory
   sb_lokimantid_loki detectionEvents.mcpl -s testFlood -n100000                # flood source: the simulated neutrons

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
   * - ``filename``
     - The detection file (first argument), relative to ``$G4PROC_MCPL_BASEDIR_LOKI`` unless a full path.
   * - ``-s``, ``--savename``
     - Output name base: ``<name>_detector.nxs`` and ``<name>_monitor.nxs``, relative to ``$G4PROC_SAVEDIR_LOKI`` unless a full path.
   * - ``--mcstasDir``
     - McStas output directory with the monitor files (McStas + Geant4), relative to ``$G4PROC_MCSTAS_BASEDIR_LOKI`` unless a full path.
   * - ``--mcstas_monitors``
     - The monitor files to load (replaces the default list).
   * - ``-n``, ``--neutronNumber``
     - The total number of simulated neutrons (flood source; all jobs of a merged file).
   * - ``-show``, ``--showMetadata``
     - Print the metadata of the detection file and exit.
   * - ``-det``, ``--detectorOnly``; ``-mon``, ``--monitorOnly``
     - Write only the detector or only the monitor file (the monitor file alone is quick).
   * - ``-i``, ``--idfCreation``
     - Write the instrument definition files instead of the NeXus files.
   * - ``--singleNexus``
     - One file for all banks (default: one per bank).
   * - ``--detectorWorkspaceTofStart``, ``...BinWidth``, ``...End``
     - Time binning of the detectors in µs (default 11000, 1000, 113000).
   * - ``--monitorWorkspaceTofStart``, ``...BinWidth``, ``...End``
     - Time binning of the monitors in µs (default 9000, 1000, 100000).
   * - ``--rear_detector_distance_m``, ``--analysis_straw_pixel_number``, ...
     - Simulation parameters for detection files without them in their metadata (or to override them).
   * - ``-no``, ``--noOutput``
     - Process without writing (for development).
   * - ``-v``, ``--verbose``
     - More output.

The base directories can be set once, e.g. in a file sourced in each session:

.. code:: sh

   export G4PROC_MCSTAS_BASEDIR_LOKI="$HOME/loki-mcstas"
   export G4PROC_MCPL_BASEDIR_LOKI="$HOME/lokiSimulations"
   export G4PROC_SAVEDIR_LOKI="$HOME/nexusFiles/loki"

**How it works.** The parameters are collected (defaults, the metadata of the detection file, the command line). The
instrument definition files are made from the templates for the rear bank distance, the pixels per straw and the
banks (placeholders and per-bank comment blocks in the templates). The detection events are read in chunks and added
to the event lists of the detector workspace, with the detector number = pixel id + 1; the workspace is binned and
saved. Only the monitors differ between the McStas + Geant4 and the flood source simulations.

Monitors of McStas + Geant4 simulations
---------------------------------------

The monitor file holds the spectra of the McStas monitors (see :ref:`loki-mcstas-loki-model`) from ``--mcstasDir``:
``Mon10_PostFOC``, ``PreSampleMonitor``, ``PostSampleMonitor`` and ``beamstopMonitor_<rear distance>m``, the monitors of
the earlier LoKI model. The current model has no ``PreSampleMonitor`` and ``beamstopMonitor``: use ``--detectorOnly``
for its detection files. The earlier model wrote separate MCPL files for the front/middle and the rear banks, which give
two detection files: process them separately (the monitors are needed only once).

Monitors of flood source simulations
------------------------------------

The flood source has no McStas monitors: the script calculates the pre-sample monitor spectrum of an isotropic
scatterer from the generator parameters in the detection file. The simulation sampled only a cone (opening angles α
maximum, β minimum) instead of the full sphere, so the number of neutrons (``-n``) is scaled by the ratio of the sphere
to the sampled spherical segment:

.. math::


   F = \frac{4\pi}{2\pi h} = \frac{2}{\cos\beta - \cos\alpha}

(F = 2/(1 − cos α) without a minimum angle). The scaled number is spread uniformly over the times of flight between
those of the minimum and maximum wavelength over the source to monitor distance (computed with 1 µs bins from 0 to
100 000 µs, then rebinned). Defaults for detection files without the metadata: wavelengths 2-13 Å, cone 0-49.6°,
monitor at 23.599 m from the source.

Running on the DMSC cluster
---------------------------

The detection files are large, so the processing is best done on the cluster (quark nodes, with the Mantid module,
see :ref:`loki-installation-dmsc-cluster`). Submit with ``sbatch submit.batch``:

.. code:: sh

   #!/bin/bash
   #SBATCH --mail-user=<email>
   #SBATCH --mail-type=ALL
   #SBATCH --job-name=lokimantid
   #SBATCH --output=loadMantid.out
   #SBATCH --error=loadMantid.err
   #SBATCH --partition=quark
   #SBATCH --nodes 1-1
   #SBATCH --time=20:00:00
   #SBATCH --exclusive

   sb_lokimantid_loki /full/path/to/detectionEvents.mcpl -s /full/path/to/testOutput -n100000
