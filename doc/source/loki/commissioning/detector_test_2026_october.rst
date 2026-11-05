.. _loki-detector-test-2026-october:

Detector test with an AmBe source (October 2026)
================================================

In October 2026 the installed LoKI detector banks were tested with an AmBe neutron source in an HDPE moderator
(proposal 336349, runs 345-365), at atmospheric pressure, before the detector vessel was evacuated. The source was
placed at eight positions (setups A-H) in front of and between the banks, and B4C sheets were taped on some banks to
cast shadows. The simulation of each setup gives the expected counts per pixel, to be compared with the measured
data pixel by pixel (readout mapping, straw order, bank positions, rear-bank calibration mask).

The setups are stored in the package ``DetectorTest2026October`` (in ``LOKI/Commissioning/LokiCommissioning/``) as
simulation parameters; the source is the package ``LokiAmBe`` (in ``LOKI/Commissioning/``). The source positions,
the B4C sheets and the positions of the middle banks were found by comparing simulations with the measured data
(the shadows of the sheets and the count rates of the banks), setup by setup; the stored setups are the final ones.

The setups
----------

.. list-table::
   :header-rows: 1

   * - Setup
     - Runs
     - Source (centre of the moderator, m)
     - Where
   * - A
     - 345-348
     - (-0.153, 0.064, 1.70)
     - first source position: between the front banks, next to bank 8
   * - B
     - 349
     - (-0.012, 0.067, 1.70)
     - between banks 6 and 8
   * - C
     - 350-358
     - (-0.03, 0.065, 1.30)
     - under bank 5; new sheets on banks 6 and 8
   * - D
     - 359
     - (0.0, -0.07, 3.20)
     - in front of banks 2 and 4, low; bank 2 sheet replaced
   * - E
     - 360
     - (0.0, -0.365, 9.7654)
     - in front of the rear bank (bank 0), for its calibration mask
   * - F
     - 361, 362
     - (-0.02, -0.233, 4.07)
     - between banks 2 and 4, behind them
   * - G
     - 363
     - (-0.008, 0.065, 1.74)
     - between banks 6 and 8
   * - H
     - 365
     - (0.0, -0.365, 9.7654)
     - in front of the rear bank, as in setup E

Common to all setups (the ``[common]`` part of ``data/setups.toml``):

- the surveyed bank positions (``bank_calibration=2026-September-SAM-606``), with the rear bank at the depth of the
  measurement (``rear_detector_distance_m=9.985615``: the surveyed front face with the detector carriage at
  4999.96 mm) and the middle banks 2 and 4 570 mm further downstream than surveyed (``bank_shifts_mm``: they had been
  moved on their rails);
- the design B4C masks of the banks, without the ones that were not on the detector (``boron_mask_exclude``,
  ``with_triangular_masks=no``);
- the rear-bank calibration mask (SI-7615-993, ``with_rear_calibration_mask=yes``) with a 0.30 mm Gd2O3 paint layer;
- air at atmospheric pressure (``world_material=G4_AIR``, ``pack_fill_material=G4_AIR``) and the aluminium detector
  vessel (``vacuum_tank_wall_mm=25``);
- the AmBe source of ``LokiAmBe`` with the table of its 26 cm x 26 cm moderator, all neutron energies;
- 512 pixels per straw in the analysis (as in the measured files).

The B4C sheets on the banks are ``absorber_strips`` and ``absorber_plates`` items of the geometry; the catalogue
lists each sheet with the setups it was in.

Running a setup
---------------

The list of setups, and the parameters of one:

.. code:: sh

   sb_detectortest2026october_params
   sb_detectortest2026october_params C

``setup=<name>`` gives a simulation all parameters of the setup; parameters given on the command line replace those
of the setup:

.. code:: sh

   sb_detectortest2026october_sim setup=C -n 10000
   sb_detectortest2026october_sim setup=C --viewer
   sb_detectortest2026october_sim setup=C -n 10000 z_meters=1.35      # the source 5 cm further downstream

The simulation and its analysis in one step (as with ``sb_loki_simanachain``, see
:ref:`loki-simulation-simulation-and-analysis-in-one-step`):

.. code:: sh

   sb_detectortest2026october_simanachain results_C setup=C -n 100000 --cleanup

The setup name and a checksum of the setups file are written to the GRIFF metadata (user data
``detector_test_setup``, ``detector_test_setups_checksum``).

The statistics of a comparison with the data need 10⁸-10⁹ source neutrons, i.e. a cluster run (see :doc:`../dmsc`): many
jobs of ``sb_detectortest2026october_simanachain results setup=<name> -n<events> -s<seed> --cleanup`` with their
own seeds. On the DMSC cluster, keep a job at about 2.5·10⁵ events: the GRIFF file of a larger job can fill the
node's local disk (``--cleanup`` removes it after the analysis).

.. _loki-detector-test-2026-october-output:

Simulated counts per pixel (HDF5)
---------------------------------

``sb_lokioutput_writehdf5`` (package ``LokiOutput``) adds up the analysis outputs of the jobs and writes the counts
per pixel into an HDF5 file:

.. code:: sh

   sb_lokioutput_writehdf5 jobs/job*/results -o setupC.h5 description="setup C, 4e9 source neutrons"

Per bank (``/entry/loki_detector_<bank>``): ``detector_number`` (as in the measured files), ``counts`` and
``counts_variance``, the pixel positions of the simulated geometry (``x/y/z_position``, m) and the ``tube``,
``straw``, ``pixel`` indices. For an AmBe simulation the file also has ``rate_1Ci``, the counts as a rate (1/s) for a
1 Ci source (2.2·10⁶ neutrons/s): each simulated neutron leaves the moderator, and the yield printed by the source
(neutrons leaving the moderator per source neutron, ``AmBeSourceGen: ... neutrons per source neutron``, read from the
simulation logs or given with ``--ambe-yield``) converts the number of simulated neutrons into a measuring time. The
file loads into scipp with ``LOKI/LokiOutput/data/hdf5_scipp_snippet.py``. ``sb_lokioutput_writehdf5 --help`` lists the
options.

The AmBe source
---------------

``sb_lokiambe_sim`` is ``sb_loki_sim`` with the generator ``event_gen=ambe`` (``LokiAmBe/python/AmBeSourceGen.py``).
The neutrons start on the surface of the moderator cylinder (vertical axis): its curved side, top and bottom, in
proportion to the neutrons leaving through each of them, with the energy spectrum and the angular distribution (to the
surface normal) of the neutrons leaving through that surface. These come from a table in ``LokiAmBe/data``; the
moderator itself is not in the geometry.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``x_meters``, ``y_meters``, ``z_meters``
     - the centre of the moderator (geometry coordinates: z along the beam, y up, origin at the sample)
   * - ``spectrum_file``
     - the table (default ``ambe_hdpe_d26cm_h26cm_spectrum.txt``: 26 cm diameter, 26 cm height), or the path of a table
   * - ``moderator_radius_cm``, ``moderator_height_cm``
     - the moderator size (default: from the table header)
   * - ``min_energy_eV``, ``max_energy_eV``
     - the energy range generated (the yield is that of this range)
   * - ``surfaces``
     - the surfaces to start from (default ``side,top,bottom``)
   * - ``angular_distribution``
     - ``moderator`` (the table) or ``isotropic`` (in the outward hemisphere)
   * - ``aim_half_angle_deg``, ``aim_x/y/z_meters``
     - direction biasing: directions only within a cone around the direction to the aim point, with weights (e.g. for
       the rear bank)

The table (spectra through the side, top and bottom in 10 bins per decade, and angular distributions in five energy
groups, per source neutron) is made with the moderator simulation of ``LokiAmBe``: the bare AmBe source (isotropic,
approximate ISO 8529-1 spectrum) at the centre of the HDPE cylinder (NCrystal ``Polyethylene_CH2.ncmat``, 0.95 g/cm³,
with a steel capsule), and the neutrons leaving the cylinder written to an MCPL file:

.. code:: sh

   sb_lokiambe_moderatorsim -n 200000 -s 5000 -m MINIMAL moderator_radius_cm=13 moderator_height_cm=26 \
       --mcpl="RecordShell where trk.pdgcode==2112 to leak"
   sb_lokiambe_makespectrum leak.mcpl.gz 200000 ambe_hdpe_d26cm_h26cm_spectrum.txt 13 26

The table in ``LokiAmBe/data`` was made from 100 such jobs on the DMSC cluster (2·10⁵ source neutrons each, seeds
5000-5099, each in its own directory ``jobs/job<i>``), with ``sb_lokiambe_makespectrum jobs 20000000 <table> 13 26``
(a directory argument reads ``jobs/job*/leak.mcpl.gz``). A table for another moderator size is made the same way and
given to the source with ``spectrum_file`` (a table name in ``LokiAmBe/data`` or a path).

Geometry parameters used by the setups
--------------------------------------

These geometry parameters (``G4GeoLoki``, ``GeoBCSBanks``) describe the state of the detector during a measurement;
the setups set them, and they can be used in any simulation.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``bank_shifts_mm``
     - banks moved from their calibrated positions: ``"bank:dx,dy,dz;..."`` (mm, geometry frame), e.g.
       ``"2:0,0,570;4:0,0,570"``. The absorbers given for a shifted bank move with it, and the analysis uses the
       shifted positions.
   * - ``boron_mask_banks``, ``with_triangular_masks``, ``boron_mask_exclude``
     - the design B4C masks: the banks that have them, the triangular masks of banks 5 and 7, and individual masks not
       present (``"bank:mask,mask;..."``, the mask ids of ``BoronMasks.cc``)
   * - ``absent_packs``
     - detector packs not on the detector: ``"bank:pack,pack;..."``
   * - ``absorber_strips``
     - absorber strips (``absorber_strip_material``, default B4C) in front of a bank: ``"bank,x1,y1,z1,x2,y2,z2,width,thickness;..."``
       (mm): the strip from point 1 to point 2 (geometry coordinates), projected onto the plane just in front of the
       bank. A 10th value > 0 places the strip inside the bank volume from point 1 to point 2 as given.
   * - ``absorber_plates``
     - absorber plates on a side face of a bank: ``"bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance[,on_pack]];..."``
       (mm): the centre line from point 1 to point 2, projected onto the plane just outside the face nearest to the
       points. A distance > 0 places the plate inside the bank volume at this distance from the points; with
       ``on_pack`` 1 it lies on the outer face of the outermost pack, parallel to the packs.
   * - ``with_rear_calibration_mask``, ``calibration_paint_thickness_mm``
     - the rear-bank calibration mask (SI-7615-993: an aluminium plate painted with Gd2O3 except ten slits) and the
       thickness of its paint
   * - ``world_material``, ``pack_fill_material``
     - the gas around the banks and inside them (``G4_AIR`` for a measurement at atmospheric pressure)
   * - ``vacuum_tank_wall_mm`` (and ``vacuum_tank_*``)
     - the detector vessel: an aluminium cylinder around the banks (0: none)
