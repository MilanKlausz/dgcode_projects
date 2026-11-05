Commissioning
=============

Simulations of the measurements made while the LoKI detectors are commissioned, with the setups of the measurements
stored in the repository, so that each can be simulated again (with more statistics, or after a code change) and
compared with the measured data pixel by pixel. The packages are in ``LOKI/Commissioning/``:

.. list-table::
   :header-rows: 1

   * - Package
     - Contents
   * - ``LokiAmBe``
     - The AmBe neutron source in an HDPE moderator (``sb_lokiambe_sim event_gen=ambe``) and the moderator simulation
       that makes its source table.
   * - ``LokiCommissioning/DetectorTest2026October``
     - The setups of the detector test with an AmBe source in October 2026 (``sb_detectortest2026october_sim setup=<name>``).

The simulated counts per pixel are written to an HDF5 file readable with scipp by ``sb_lokioutput_writehdf5`` of the
package ``LOKI/LokiOutput``, see :ref:`loki-detector-test-2026-october-output`.

.. toctree::
   :maxdepth: 1

   commissioning/detector_test_2026_october
