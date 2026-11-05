Tests
=====

``sb -t`` builds everything and runs the tests of all packages; run it before pushing changes (GitHub Actions runs the
tests as well). A test is a script ``scripts/test_*`` of a package: it passes if it exits with 0 and its output equals
the reference log next to it (``scripts/test_*.log``, compared ignoring whitespace). Only the LOKI tests:

.. code:: sh

   sb -t --testfilter='sb_g4geolokitests_*,sb_g4geoloki_*,sb_loki_test*,sb_lokimasking_test*,sb_larmor_*,sb_larmortests_*'

``sb -t`` runs the tests in parallel (``-j1``: one at a time, e.g. to time them: ``time sb -t -j1 ...``). The results,
with the output of each test, are in ``simplebuild_cache/bld/testresults/<test>/``. A run with a filter clears the
results of the other tests.

The LOKI tests
--------------

.. list-table::
   :header-rows: 1

   * - Test
     - Checks
   * - **G4GeoLokiTests**
     - 
   * - ``test_geo_construct_overlaps``, ``test_geo_dense_overlaps``
     - No overlaps in the geometry (Geant4 check; dense random sampling), for the nominal geometry and the calibrations, with beamstop and slits; the notch of banks 6 and 8.
   * - ``test_geo_placements_pinned``, ``test_geo_navigation``
     - The placements of all volumes; Geant4 navigation through the banks.
   * - ``test_bank_transform``
     - The bank transforms are rotations; pixel ids at and near the pixel centres (nominal and SAM-606).
   * - ``test_bank_calibration``
     - Loading of calibration files, errors for bad files, the default calibration.
   * - ``test_aimhelper_pixels``, ``test_pixelid_roundtrip``, ``test_pixels_per_straw``
     - Pixel centres and ids, for different numbers of pixels per straw.
   * - ``test_geantino_aiming``
     - Geantinos aimed at pixel centres reach their pixel in the Geant4 geometry (every straw, nominal and SAM-606).
   * - ``test_simana_pixelid``
     - Simulation and analysis give the pixel ids of the aimed pixels (also with a calibration).
   * - ``test_masking_chain``, ``test_masking_full_instrument``, ``test_mask_files``
     - The masking simulation and analysis, the mask files.
   * - ``test_recorded_calibration``
     - The analysis uses the calibration recorded by the simulation.
   * - ``test_pre_calibration_griff``
     - GRIFF files from before the calibrations are analysed with the nominal geometry.
   * - ``test_verify_json``
     - ``sb_g4geoloki_verifyjson`` with JSON files of loki-geometry.
   * - **G4GeoLoki**
     - 
   * - ``test_pixel_position``, ``test_pixel_position_default``, ``test_write_pixel_file``
     - The pixel positions (nominal, default calibration) and ``sb_g4geoloki_writepixelfile``.
   * - **LOKI**, **LokiMasking**
     - 
   * - ``test_info``, ``test_info_pixel``
     - ``sb_loki_info``.
   * - ``test_det``, ``test_flood``, ``test_mcpl``, ``test``, ``test_calib``
     - Short simulations and analyses (summarised with ``python3 -m LOKI.OutputSummary``).
   * - **Larmor**, **LarmorTests**
     - The Larmor setups: geometries, pixel positions and conversions, simulations, masking.

The whole LOKI test set takes about 4.5 minutes when the tests run one at a time; the longest are
``test_geantino_aiming`` and ``test_masking_full_instrument`` (about 30 s each). The LokiMantid tests run only with
Mantid.

Updating a reference log
------------------------

When a change is meant to change the output of a test, check the difference, then copy the new output over the log:

.. code:: sh

   sb -t --testfilter=sb_g4geolokitests_test_bank_transform
   diff -w LOKI/G4GeoLokiTests/scripts/test_bank_transform.log \
       simplebuild_cache/bld/testresults/sb_g4geolokitests_test_bank_transform/output.log
   cp simplebuild_cache/bld/testresults/sb_g4geolokitests_test_bank_transform/output.log \
       LOKI/G4GeoLokiTests/scripts/test_bank_transform.log

The repository rejects text files larger than 200 kB and binary files larger than 20 kB (also changed ones), so keep
logs and test data small.
