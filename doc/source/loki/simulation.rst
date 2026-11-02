Simulation
==========

.. _loki-simulation-a-first-simulation:

A first simulation
------------------

After the :doc:`installation <installation>` (including ``eval "$(sb --env-setup)"``), in an empty directory: a flood
source (uniform scattering from the sample) aimed at bank 3, its analysis, and a look at the result.

.. code:: sh

   mkdir firstrun && cd firstrun
   sb_loki_sim event_gen=flood aiming_bank_id=3 -n2000 -j4 > sim.log 2>&1
   sb_loki_ana *.griff > ana.log 2>&1
   python3 -m LOKI.OutputSummary          # numbers of detections and histogram summaries

.. list-table::
   :header-rows: 1

   * - Step
     - Writes (in the current directory)
   * - ``sb_loki_sim ... -j4``
     - ``lokisim.0.griff`` ... ``lokisim.3.griff``: one GRIFF file per process (``lokisim.griff`` without ``-j``); about 10 kB per event.
   * - ``sb_loki_ana *.griff``
     - ``detectionEvents.mcpl`` (the detection events) and ``bcsloki_sans.shist`` (histograms), see :doc:`analysis`.

The simulation log is long (an overlap check line per volume): keep it in a file and look for ``NOTE``, ``WARNING`` and
``ERROR`` lines. The surveyed bank positions are used by default (see :ref:`loki-geometry-bank-calibrations`).
:ref:`loki-analysis-inspecting-the-results` shows how to read the detection events.

Running a simulation
--------------------

A simulation is started by a script of a package: ``sb_loki_sim`` (:sbpkg:`LOKI/scripts/sim`) for neutrons and
``sb_lokimasking_sim`` (:sbpkg:`LokiMasking/scripts/sim`) for the masking. These scripts take the geometry
(:doc:`geometry`) and a launcher (:sbpkg:`LOKI/python/Launcher.py`, :sbpkg:`LokiMasking/python/Launcher.py`), which
chooses the particle generator. The parameters given on the command line (``name=value``) belong to:

.. list-table::
   :header-rows: 1

   * - Parameters of
     - Defined in
     - Examples
     - List with
   * - the geometry
     - ``GeoBCSBanks``
     - ``rear_detector_distance_m``, ``bank_calibration``, ``beamstop_id``, ``with_calibration_slits``
     - ``-g``
   * - the launcher
     - ``Launcher.py``
     - ``event_gen``, ``analysis_straw_pixel_number``, ``aiming_bank_id``
     - ``-d``
   * - the generator
     - the generator
     - ``neutron_wavelength_min_aangstrom`` (flood), ``input_file`` (MCPL)
     - ``-p``

``-d`` lists them all with their values (also with other parameters set, e.g. ``sb_loki_sim event_gen=flood -d``). In
its output the launcher parameters come after the generator parameters, followed by an empty ``LauncherPars:`` line
(a quirk of the dgcode output). Options of all simulation commands:

.. list-table::
   :header-rows: 1

   * - Option
     - Meaning
   * - ``-n N``
     - Simulate N events in total, shared by the ``-j`` processes (for an MCPL source: all particles in the file if not given).
   * - ``-j N``
     - Run N processes in parallel (one GRIFF file each). Leave a core free for the system.
   * - ``-s S``
     - Random seed (default 123456789).
   * - ``-l PL``
     - Physics list (default ``QGSP_BIC_HP_EMZ``).
   * - ``-o NAME``
     - Name of the GRIFF output (default ``lokisim``; ``loki_masking`` for the masking).
   * - ``--view``, ``--data``, ``--aim``
     - Visualisation instead of a simulation (see :ref:`loki-geometry-visualisation`).
   * - ``-t``
     - Check the geometry for overlaps and exit.

The output goes to the current directory. The pixels are not part of the simulation: ``analysis_straw_pixel_number``
(default 256) is stored in the GRIFF file and used by the analysis; the instrument's NeXus geometry has 512.

Neutron generators
------------------

``sb_loki_sim`` selects the generator with ``event_gen``:

.. list-table::
   :header-rows: 1

   * - ``event_gen``
     - Generator
     - Use
   * - ``mcpl``
     - ``MCPLGen`` (dgcode)
     - Neutrons of a McStas simulation (:ref:`loki-simulation-mcpl-based-simulation`).
   * - ``flood``
     - :sbpkg:`LOKI/python/FloodSourceGen.py`
     - Uniform (isotropic) scattering from the sample (:ref:`loki-simulation-flood-source-simulation`).
   * - ``spheremodel``
     - :sbpkg:`LOKI/python/SansSphereGen.py`
     - Ideal SANS sphere-model distribution.
   * - ``isotheta``
     - :sbpkg:`LOKI/python/IsoThetaGen.py`
     - Uniform polar angle distribution (tests).
   * - (not set)
     - ``FlexGen`` (dgcode)
     - A narrow cone of 3 Å neutrons along the beam (change with ``random_max_polarangle_deg``, ``neutron_wavelength_aangstrom``, ...).

``geantino=1`` replaces the neutrons by geantinos (non-interacting particles, for visualisation), and
``primary_only=True`` keeps only the primary particles in the output.

.. _loki-simulation-mcpl-based-simulation:

MCPL-based simulation
~~~~~~~~~~~~~~~~~~~~~

The neutrons of a McStas simulation (see :doc:`mcstas`) are read from an MCPL file (``input_file``, required; ``-d``
needs it too):

.. code:: sh

   sb_loki_sim event_gen=mcpl input_file=/path/to/file.mcpl.gz -d
   sb_loki_sim -j10 event_gen=mcpl rear_detector_distance_m=10 \
       analysis_straw_pixel_number=512 input_file=/path/to/loki-collen_5.0m-ap1_30mm-ap2_10mm-lmin_3.0A-lmax_11.5A_rear.mcpl.gz

-  ``input_file``: the full path of the MCPL file.
-  Parameters of ``MCPLGen``: ``skip_events`` (skip the first particles, see :doc:`dmsc`), ``allow_zero_weight``,
   ``input_filter``, ``dz_meter`` (set automatically, see the next point).
-  The particles are placed relative to the sample position: the distance of the MCPL output component from the sample
   in the McStas model is read from the file (``sample_mcpl_distance_m``); files without it use 0.2 m (the LoKI McStas
   model).
-  All metadata of the MCPL file (e.g. the McStas command line, with the parameter values of the run) is copied to the
   GRIFF file and to the detection file, and also ``mcpl_universal_weight``, the comments (``mcpl_comments``) and the
   statistics sums (``mcpl_stat_sum_<name>``, e.g. ``mcpl_stat_sum_initial_ray_count`` of McStas).
-  The particle weights of the file are used (particles of weight 0 are allowed), and the detection events keep them
   (for a file without a universal weight, e.g. of the LoKI McStas model): the sum of the McStas weights divided by
   ``initial_ray_count`` is the intensity in neutrons/s.
-  ``gen_x_offset_meters`` shifts the beam horizontally.
-  ``nominal_source_sample_distance_meters`` (default 23.6 m) is used only for approximate histograms in the analysis.

**The two MCPL files of a LoKI McStas simulation** (see :ref:`loki-mcstas-loki-model`): one for the front and middle
banks (independent of the rear bank distance) and one for the rear bank (at the McStas ``reardet``). Together they
hold all neutrons towards the detectors; each can come from a simulation with its own statistics. Simulate and analyse
each with the ``rear_detector_distance_m`` and ``beamstop_id`` of the configuration: its result is the sum of the two
detection files, each normalised with its own ``initial_ray_count``.

.. _loki-simulation-flood-source-simulation:

Flood source simulation
~~~~~~~~~~~~~~~~~~~~~~~

``FloodSourceGen`` emits neutrons from the sample in random directions within a cone, with a uniform wavelength
distribution: a flat (isotropic) scatterer without a McStas simulation. The start time of each neutron is set from
its wavelength and the source to sample distance (a single, instantaneous pulse).

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Meaning
   * - ``neutron_wavelength_min_aangstrom``, ``neutron_wavelength_max_aangstrom``
     - 2, 13
     - Wavelength range.
   * - ``cone_opening_deg``, ``cone_opening_min_deg``
     - 49.6, 0
     - Range of the angle from the reference direction (replaced when ``aiming_bank_id`` is set).
   * - ``ref_dir_x``, ``ref_dir_y``, ``ref_dir_z``
     - 0, 0, 0 (= along the beam)
     - Reference direction of the cone (set by ``aiming_bank_id``).
   * - ``source_monitor_distance_meters``
     - 23.599
     - Position of a calculated pre-sample monitor (stored in the output).
   * - ``gen_x_offset_meters``, ``gen_y_offset_meters``, ``gen_z_offset_meters``
     - 0
     - Source position offset.
   * - ``gen_x_width_meters``, ``gen_y_width_meters``
     - 0.006, 0.008
     - Size of the source (the sample).
   * - ``particle``
     - ``neutron``
     - The particle (``geantino=1`` of the launcher sets it).

The launcher parameter ``nominal_source_sample_distance_meters`` (default 23.6 m) sets the source to sample distance
(time of flight). The generator parameters are stored in the output (e.g. to calculate the incident spectrum).

**``aiming_bank_id``** limits the directions to a bank or a group of banks (more efficient):

.. list-table::
   :header-rows: 1

   * - Value
     - Direction
     - Cone (replaces ``cone_opening_deg``, ``cone_opening_min_deg``)
   * - a bank id, 0-8
     - the centre of that bank (in the bank placement of the geometry in use)
     - 0 to 5.07, 9.9, 4.9, 9.9, 4.9, 31.9, 20.9, 29.5, 22.1° for banks 0-8 (fixed values for the nominal geometry)
   * - ``1234``
     - along the beam
     - 3.5-15.8° (the middle banks)
   * - ``5678``
     - along the beam
     - 10.5-49.6° (the front banks)

The cone also reaches other banks. The analysis writes only the detections in the aimed banks to the detection file;
its histograms count all banks (see :ref:`loki-analysis-neutron-analysis`). ``cone_view=min`` or ``max`` emits only at the
minimum or maximum angle, to check the cone visually:

.. code:: sh

   sb_loki_sim event_gen=flood -d
   sb_loki_sim event_gen=flood primary_only=True --data -n1000 aiming_bank_id=1234 cone_view=min
   sb_loki_sim event_gen=flood --data -n1000 aiming_bank_id=3 cone_view=max geantino=1
   sb_loki_sim event_gen=flood aiming_bank_id=3 -n10000 -j10   # simulation

.. _loki-simulation-masking-simulation:

Masking simulation
------------------

``sb_lokimasking_sim`` (:sbpkg:`LokiMasking/python/MaskingSourceGen.py`) shoots one geantino from the sample at the
centre of each pixel, in the order of the pixel ids. Geantinos travel in straight lines without interacting; the
analysis (:ref:`loki-analysis-masking-analysis`) finds the pixels a geantino reaches without passing an absorber (boron
masks, absorber panels, beamstop, calibration slits). Uses:

1. **Masking**: the pixels not directly visible from the sample are masked in the data reduction.
2. **Calibration slits**: with ``with_calibration_slits=true``, the pixels illuminated through the slits. The
   overlapping layers see a slit at different positions along the straws (parallax), which can be used for position
   corrections. The mask file then lists the pixels that are *not* illuminated; the illuminated ones are the rest.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Meaning
   * - ``aiming_straw_pixel_number``
     - 512
     - Pixels per straw for the aiming only: more than ``analysis_straw_pixel_number`` aims several geantinos at each pixel (better for partly shadowed pixels).
   * - ``analysis_straw_pixel_number``
     - 256
     - Pixels per straw of the analysis (mask files).
   * - ``aiming_bank_id``
     - -1
     - Aim only at this bank (from its first pixel); -1: all banks from pixel 0.
   * - ``aiming_pixel_id_min``
     - -1
     - First pixel to aim at (overrides the start of ``aiming_bank_id``).
   * - ``masking_only``
     - no
     - Store only the volumes needed for the masking (smaller GRIFF files).
   * - ``gen_x_offset_meters``, ``gen_x_width_meters``, ``gen_y_width_meters``
     - 0
     - Position and size of the source.

``-d`` also lists ``geo_rear_detector_distance_m`` and ``geo_bank_calibration``: copies of the geometry parameters used by
the aiming (set the geometry parameters, not these).

**Which settings.** The mask files hold detector numbers, which depend on the pixels per straw and on the bank
placement:

.. list-table::
   :header-rows: 1

   * - Mask for
     - ``analysis_straw_pixel_number``
     - ``bank_calibration``
   * - measured data
     - 512, as the instrument's NeXus geometry (``loki.json``)
     - the surveyed calibration (the default, ``2026-September-SAM-606``)
   * - simulated data
     - the value of the simulation
     - the calibration of the simulation

Use at least as many aiming pixels as analysis pixels (e.g. ``aiming_straw_pixel_number=1024``).

**Number of events.** ``-n`` must cover the pixels to aim at: the pixels of a bank are its tubes × 7 × the aiming pixels
per straw (the tubes are in :ref:`loki-geometry-numbering-and-pixel-ids`), and ``sb_loki_info -s 1024 -d`` lists them (6 422 528
for the full geometry with 1024). The simulation also prints the number, at the end of its log. The generator
supports ``-j`` (the processes share the pixels).

.. code:: sh

   sb_lokimasking_sim -d
   sb_lokimasking_sim aiming_bank_id=5 aiming_straw_pixel_number=1024 -n1024 --data
   sb_lokimasking_sim aiming_straw_pixel_number=1024 analysis_straw_pixel_number=512 -n6422528 -j10
   sb_lokimasking_sim with_calibration_slits=true aiming_straw_pixel_number=1024 analysis_straw_pixel_number=512 -n6422528 -j10

.. _loki-simulation-simulation-and-analysis-in-one-step:

Simulation and analysis in one step
-----------------------------------

The chain commands run a simulation in ``DIR/runsim`` and then the analysis of its GRIFF files in ``DIR/runana``; they
also run a 1-event simulation in ``DIR/runsample`` (to keep an example GRIFF file when ``--cleanup`` is used). Nothing is
printed to the terminal; the logs are in the directories.

.. list-table::
   :header-rows: 1

   * - Command
     - Runs
   * - ``sb_loki_simanachain DIR ...``
     - ``sb_loki_sim ...`` → ``sb_loki_ana``
   * - ``sb_loki_simanachain_det DIR ...``
     - ``sb_loki_sim_det ...`` → ``sb_loki_ana_det``
   * - ``sb_lokimasking_simanachain DIR ...``
     - ``sb_lokimasking_sim ...`` → ``sb_lokimasking_ana``

The arguments after the directory go to the simulation. Options: ``--cleanup`` removes the GRIFF files after the
analysis (they are large: about 10 kB per event), ``--chainmult=M`` runs the chain M times with different seeds and
merges the histograms. ``sb_loki_sim_det`` / ``sb_loki_ana_det`` (also ``det_only=True`` of ``sb_loki_sim``) store only the
converter and counting gas volumes, for detection studies with small files.

Other scripts
-------------

``sb_loki_scan`` and ``sb_loki_scanana`` are an old example of a parameter scan (with dgcode's ScanUtils; ``scan`` calls a
chain command that no longer exists), and ``sb_loki_binmcplneutronwavelength`` writes the wavelength distribution of an
MCPL file to a CSV file.
