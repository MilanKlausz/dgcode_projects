# Larmor detector tests

The LoKI rear bank detectors were tested at the Larmor instrument (ISIS) in December 2020 and in March 2022. The
models of these tests are in `LOKI/Larmor/`:

| Package | Contents |
| --- | --- |
| {sbpkg}`Larmor` | The geometries (`GeoLarmor2020`, `GeoLarmor2022`), simulation, analysis and masking commands, pixel id and mask file conversion. |
| {sbpkg}`LarmorTests` | Tests of the Larmor geometries, aiming and pixel ids. |
| {sbpkg}`LarmorMantid` | Mantid processing of the Larmor simulations (`larmor2020`, `larmor2022`), built only with Mantid. |

The Larmor setups are frozen: they have their own copy of the bank placement and of the tube and pixel numbering
(`Larmor2022Bank`), so changes of the LOKI geometry (e.g. the bank calibrations) do not change them. The detector
hardware (tubes, packs, the rear bank box and its masks) comes from G4GeoLoki.

## 2022: the LoKI rear bank

The test used the rear bank with 16 of its 28 packs. The geometry `GeoLarmor2022` has the full bank at 4.099 m (at the
Larmor beam height); the analysis keeps only the detections in the 16 packs of the test (backscattering from the other
packs is expected to be negligible). The detection area is also limited by the end window of the vacuum tank.

| Command | Does |
| --- | --- |
| `sb_larmor_sim_2022` | Neutron simulation (generators as for LOKI: `event_gen=mcpl`, `flood`, ...). Defaults of the test: `rear_detector_distance_m=4.099`, `gen_x_offset_meters=0.005` (the beam 5 mm to the left), `nominal_source_sample_distance_meters=25.61`, `analysis_straw_pixel_number=512`. |
| `sb_larmor_ana_2022` | Analysis: detection file with the pixel ids of the test. |
| `sb_larmor_simanachain_2022` | Simulation and analysis. |
| `sb_larmor_masking_sim_2022`, `sb_larmor_masking_ana_2022`, `sb_larmor_masking_simanachain_2022` | Masking (as `sb_lokimasking_*`). |
| `sb_larmor_maskconverter_2022` | Merge mask files and convert them to the pixel ids of the test. |
| `sb_larmor_writepixelfile_2022` | Pixel centres of the test geometry (`--fullGeomOutput`: all 28 packs). |
| `sb_larmor_verifyjson_2022` | Compare a NeXus-JSON geometry of the test with the Geant4 pixel centres. |

`GeoLarmor2022` has a cadmium calibration mask (`with_calibration_slits=true`) and the parameter `old_tube_numbering`
(the tube numbering of the data of the test). The pixel id conversions between the LoKI rear bank and the test are in
{sbpkg}`Larmor/python/PixelIdConverters.py`.

```sh
sb_larmor_sim_2022 -d
sb_larmor_sim_2022 with_calibration_slits=true --aim
sb_larmor_sim_2022 event_gen=mcpl input_file=/path/to/larmor_postsample.mcpl.gz -j10
```

The McStas source is the Larmor model (see [](mcstas.md#larmor-model)); its MCPL output is 4.049 m from the sample.

## 2020: BCS detectors

The 2020 test used BCS detectors of different lengths and converter thicknesses. `GeoLarmor2020` is a separate model,
built from the same tube and pack classes. Its rear detector distance is 4.42 m; the MCPL source is 4.356 m from the
sample and the beam is offset by −40 mm.

| Command | Does |
| --- | --- |
| `sb_larmor_sim_2020`, `sb_larmor_ana_2020`, `sb_larmor_simanachain_2020` | Simulation and analysis. |
| `sb_larmor_ana_2020_det`, `sb_larmor_simanachain_2020_det` | The detector-only variant. |

## Mantid

`sb_larmormantid_larmor2022` works like `sb_lokimantid_loki` (see [](mantid.md)) with the instrument definition of the
2022 test. `sb_larmormantid_larmor2020` writes a single NeXus file with the monitors and the detectors (the first 10
detector ids are monitors, so the pixel ids are offset by 11).
