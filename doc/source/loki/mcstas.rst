McStas simulation
=================

.. note::

   The McStas models are not part of dgcode_projects: the LoKI model is ``instr/LOKI_main.instr`` of the
   `LOKI <https://git.esss.dk/dmsc-instrumentmodels/LOKI>`__ repository of the DMSC instrument models, with the MCPL
   output for Geant4 on its ``geant4`` branch; the Larmor model is
   `larmor-mcstas <https://github.com/ess-loki-sim/larmor-mcstas>`__. The Larmor part of this page is from the earlier
   LOKI simulation guide and has not been checked against the current version of that repository.

The McStas simulation transports the neutrons from the moderator to the sample and saves the scattered neutrons in an
MCPL file, the source of the Geant4 simulation (``sb_loki_sim event_gen=mcpl``, see
:ref:`loki-simulation-mcpl-based-simulation`).

Particle weights
----------------

McStas uses particle weights; the Geant4 simulation and the data reduction need whole neutrons (weight 1). The
MCPL output component (``MCPL_output_noacc_deweight`` in the LoKI model, ``MCPL_output_deweight`` in the Larmor model)
converts the weights when writing the MCPL file (``deweight=1``):

-  **Splitting**: while the weight p ≥ 1, a neutron of weight 1 is written and p is reduced by 1.
-  **Russian roulette**: for 0 < p < 1, a neutron of weight 1 is written with probability p (``if (p > rand) p = 1;   else p = 0;``).

This increases the variance; it is best if the weights are already below 1 when they are written, so that only the
Russian roulette is applied. The weights depend on the source intensity and the number of simulated neutrons: with
1000 n/s and 100 simulated neutrons the mean weight is 10. To lower the weights without simulating more neutrons than
the intensity, ``SPLIT n`` before the sample makes n copies of each neutron with 1/n of its weight. This improves the
sampling of the sample scattering at the cost of a coarser sampling of the instrument before it, so it should not be
overused. In practice: increase the number of neutrons and the SPLIT without deweighting until the MCPL file has more
particles than the intensity, then use the same numbers with deweighting.

.. _loki-mcstas-mcpl-file-size:

MCPL file size
--------------

Keep the number of particles in an MCPL file below 2³¹ (≈ 2.15·10⁹; the ``skip_events`` limit of the Geant4
submission, see :doc:`dmsc`); that is already 78 GB compressed. To keep only useful neutrons:

-  A slit after the sample (``mcpl_slit`` in the LoKI model) passes only neutrons that can reach the detectors. Check its size
   with the Geant4 visualisation of the MCPL file (a few thousand particles, ``--data primary_only=True``), and compare
   the detection counts of a full simulation with a larger slit.
-  The beamstop removes the transmitted beam (earlier versions also stopped the unscattered neutrons in McStas).

Split larger files with ``mcpltool`` (also in the dgcode environment): ``-e`` extracts, ``-sN`` skips the first N particles,
``-lN`` limits the number:

.. code:: sh

   mcpltool -e larmor_postsample.mcpl.gz larmor_postsample_after2e9 -s2000000000

.. _loki-mcstas-geant4-input:

Geant4 input
------------

The MCPL tools ``mcpltool`` and ``pymcpltool`` come with McStas and with the dgcode environment. The
number of particles of a file is needed for the Geant4 submission:

.. code:: sh

   pymcpltool --stats larmor_postsample.mcpl.gz > stats.txt &   # slow for ~1e9 particles

The MCPL output component also stores the McStas command line and the instrument file in the MCPL file:

.. code:: sh

   pymcpltool file.mcpl.gz -bmccode_cmd_line
   pymcpltool file.mcpl.gz -bmccode_instr_file
   pymcpltool file.mcpl.gz -bsample_mcpl_distance_m

``sample_mcpl_distance_m`` is the distance of the MCPL output component from the sample; ``sb_loki_sim`` reads it and
places the particles accordingly. The particles start at the post-sample slit, not at the sample (seen in the
visualisation).

.. _loki-mcstas-loki-model:

LoKI model
----------

``instr/LOKI_main.instr`` of the `LOKI <https://git.esss.dk/dmsc-instrumentmodels/LOKI>`__ repository, branch
``geant4`` (McStas 3.9.1 or later; the repository has the conda environment, ``conda.yml``). The parameters of the
instrument are described in its README; those for the Geant4 simulation:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``mcplout``
     - 1: write the MCPL file for Geant4 (default 0).
   * - ``acc_power``
     - Accelerator power in MW (default 5), the intensity of the ``ESS_butterfly`` source.
   * - ``deweight``
     - Apply the deweighting (default 1).
   * - ``beamstop``
     - Beamstop in front of the rear bank: 0 = 70 mm × 70 mm (default), 1-5 = the beamstops of ESS-1178830 (the
       ``beamstop_id`` of the Geant4 geometry, see :ref:`loki-geometry-beamstops`).
   * - ``sampletype``
     - The sample: negative = empty beam (no sample), 0 = water, 1-18 = the models of ``SANS_benchmark2`` (e.g. 10 =
       logarithmically spaced sharp peaks in Q; 5, 6, 15 = spheres), 19 = flat scatterer towards all detector banks.
   * - ``collen``, ``reardet``
     - Collimation length (3, 5 or 8 m) and rear bank distance from the sample [m].

**MCPL output.** ``mcpl_out`` writes one file, 0.2 m after the sample: the neutrons that pass the beamstop and the slit
``mcpl_slit`` (towards the detector banks). Its name is made of the instrument settings, e.g.
``loki-collen_5.0m-ap1_30mm-ap2_10mm-lmin_3.0A-lmax_11.5A-reardet_10.0m.mcpl.gz``. Simulate it with the same rear bank
distance (``rear_detector_distance_m``) and beamstop (``beamstop_id``). The SPLIT before the sample is fixed (10), so
the weights are lowered with more neutrons (``-n``) or a lower ``acc_power``; with ``deweight=1`` the file has the
neutrons of one second of beam at ``acc_power``.

**From McStas to Geant4** (the parameter values are examples):

.. code:: sh

   # 1. McStas: writes the MCPL file (and the monitor files) to DIR
   mcrun -y LOKI_main.instr -d DIR -n1e9 collen=5 reardet=10 beamstop=5 mcplout=1 acc_power=0.5
   # 2. Geant4: simulation and analysis of the MCPL file
   sb_loki_simanachain g4 event_gen=mcpl input_file=DIR/<name>.mcpl.gz rear_detector_distance_m=10 beamstop_id=5 \
       analysis_straw_pixel_number=512 -j10

For a quick test of the chain use ``-n1e7 acc_power=5e-3`` (about 2·10⁴ neutrons in the MCPL file). For large files
see :doc:`dmsc`.

**Mantid.** The Mantid processing (see :doc:`mantid`) was written for the earlier LoKI model
(`loki-mcstas-master <https://bitbucket.org/MKlausz/loki-mcstas-master/src/master/>`__, Bitbucket), which wrote three
MCPL files (front and middle banks, rear bank at 5 m and at 10 m) and had the monitors ``PreSampleMonitor`` and
``beamstopMonitor_5m``/``beamstopMonitor_10m``. These two monitors are not in the current model: process its detection
files with ``--detectorOnly``.

**Numbers of the earlier model** on the DMSC cluster (apertures ``sourceapx=0.030 sourceapy=0.025 sampleapx=0.010
sampleapy=0.015``; front/middle banks and rear bank files; the flat scatterer was ``sampletype=0``, the Gaussian coil
``sampletype=18`` is not in the current model, and the SPLIT was a parameter). Local runs with 1/100 of the neutrons
and of the power give 1/100 of these numbers.

.. list-table::
   :header-rows: 1

   * - Sample
     - collen [m]
     - Power [MW]
     - SPLIT
     - Neutrons
     - bs5 / bs10
     - Sum of weights in MCPL
     - Detection events
   * - Empty beam
     - 3
     - 5
     - 1
     - 1e11
     - -
     - 1.4e9
     - 7.7e8
   * - 
     - 5
     - 5
     - 1
     - 1e11
     - -
     - 5.3e8
     - 2.9e8
   * - 
     - 8
     - 5
     - 1
     - 1e11
     - -
     - 2.05e8
     - 1.14e8
   * - Flat scatterer
     - 3
     - 50
     - 1000
     - 1e11
     - 5 / -
     - 1.2e9 / 1.2e7
     - 4.2e8 / 4.8e6
   * - 
     - 5
     - 200
     - 10000
     - 1e11
     - 3 / 5
     - 1.8e9 / 1.8e7
     - 6.4e8 / 7.2e6
   * - 
     - 8
     - 500
     - 10000
     - 1e11
     - 3 / 4
     - 1.8e9 / 1.7e7
     - 6e8 / 7e6
   * - Peaks in Q
     - 3
     - 50k
     - 1000
     - 1e11
     - 5 / -
     - 1.2e9 / 7.4e8
     - 5e9 / 4e8
   * - 
     - 5
     - 200k
     - 10000
     - 1e11
     - 3 / 5
     - 1.8e9 / 1.5e9
     - 7.6e8 / 8.4e8
   * - 
     - 8
     - 500k
     - 10000
     - 1e11
     - 3 / 4
     - 1.7e9 / 1.4-1.7e9
     - 7.2e8 / 7.8-9.5e8
   * - Gaussian coil
     - 3
     - 500
     - 1000
     - 1e11
     - 5 / -
     - 2.5e9 / 1.4e9
     - 1.1e9 / 7.0e8
   * - 
     - 5
     - 1000
     - 1000
     - 1e11
     - 3 / 5
     - 1.8e9 / 1.1e9
     - 8e8 / 5.8e8
   * - 
     - 8
     - 500
     - 1000
     - 1e11
     - 3 / 4
     - 3.5e8 / 2.2e8
     - 1.5e8 / 1.1-1.2e8
   * - Hard spheres
     - 3
     - 200k
     - 4000
     - 1e11
     - 5 / -
     - 3.11e8 / 1.94e9
     - 1.3e8 / 1.1e9
   * - 
     - 5
     - 500k
     - 10000
     - 1e11
     - 3 / 5
     - 2.85e8 / 2.01e9
     - 1.2e8 / 1.1e9
   * - 
     - 8
     - 1250k
     - 10000
     - 2e11
     - 3 / 4
     - 2.70e8 / 1.93-1.98e9
     - 1.16e8 / 1.1e9

.. _loki-mcstas-larmor-model:

Larmor model
------------

``Larmor.instr`` (with ``ViewModISISver1.comp``, which sets the intensity), main parameters:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``sampletype``
     - Negative = empty beam, 0 = flat scatterer, 10 = logarithmically spaced sharp peaks in Q, 18 = Gaussian coil,
       20 = monodisperse hard spheres, R = 50 Å (the numbering of the earlier LoKI model). Samples: 10 × 10 × 1 mm
       (the empty sample slit 50 × 50 mm).
   * - ``deweight``
     - Apply the deweighting.
   * - ``current``
     - Beam current in A (default 2); scales the intensity (``beamCurrent`` of ``ViewModISISver1.comp``).
   * - ``sampdet``
     - Sample to detector distance (default 4.099 m, the 2022 test).
   * - ``s2``
     - A1 aperture size (default 0.02 m); the A2 aperture is fixed at 6 × 8 mm.

The SPLIT is after the ``Samplepoint`` arm. The beamstop is at ``sampdet - 0.09`` (not in empty beam runs), the MCPL output
(``mcpl_out_postsample``, file ``larmor_postsample.mcpl.gz``) at ``sampdet - 0.05`` = 4.049 m from the sample. Monitors:
``TOFmon1`` (after the chopper, 9.8195 m), ``TOFmon2`` (before the A1 aperture, 20.313 m), ``TOFpresamp`` (before the
sample, 25.57 m), ``TOFmon4`` (after the sample, 25.76 m), files ``<name>.dat``.

.. list-table::
   :header-rows: 1

   * - Sample
     - Current [A]
     - SPLIT
     - Neutrons
     - Sum of weights in MCPL
     - Detection events
   * - Empty beam
     - 0.020
     - 1000
     - 1e11
     - 4.8e8
     - 2.2e8
   * - Flat scatterer
     - 20
     - 1000
     - 1e11
     - 1.5e9 (expected)
     - 6e8 (expected)
   * - Peaks in Q
     - 2000
     - 1000
     - 1e11
     - 2.1e9
     - 1.2e9
   * - Gaussian coil
     - 20
     - 1000
     - 1e11
     - 3.9e9
     - 2.1e9

Running McStas on the DMSC cluster
----------------------------------

Load the modules (see :ref:`loki-installation-dmsc-cluster`) and compile for MPI on the right node (``compile`` for the
short and newlong partitions, ``quarkcompile`` for quark): start a run with ``-c --mpi=4`` and stop it (Ctrl+C) when it
asks for the parameters:

.. code:: sh

   mcrun -c --mpi=4 Larmor.instr -n1e6 -dtest

Then submit a batch script with ``sbatch mpirun.batch``:

.. code:: sh

   #!/bin/bash
   #SBATCH --mail-user=<email>
   #SBATCH --job-name=larmor
   #SBATCH --output=job.out
   #SBATCH --error=job.err
   #SBATCH --partition=quark
   #SBATCH --nodes=1-1
   #SBATCH --time=24:00:00

   mpirun Larmor.out -dsamp18_split1000_1e11 sampletype=18 deweight=1 sampdet=4.099 -n1e11

-  Use the compiled ``.out`` file in ``mpirun``. Keep earlier commands in the file as comments.
-  Do not use dots in the output directory name (the merging of the MCPL files fails).
-  The merging and compression of the MCPL files can take longer than the simulation and does not use more nodes:
   use one node (28 cores on newlong, 32 on quark).
-  If the merging fails, repeat the simulation.
-  Partition time limits (``sinfo``): short 4 hours, newlong 7 days, quark 1 day.
-  Larmor on 1 quark node: peaks in Q 5 h + compression (43 GB), Gaussian coil 2.5 h (78 GB), empty beam 1 h (6.6 GB).
