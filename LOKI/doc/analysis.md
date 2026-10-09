# Analysis

The simulations write GRIFF files with the tracks of all particles. An analysis program (C++, in an `app_<name>`
directory of a package, run as `sb_<package>_<name>`) reads them and derives the quantities of interest; it can be
changed and rerun without repeating the simulation. The LOKI analyses use `PixelatedBanks` for the pixel ids, with the
bank placement recorded by the simulation (see [](geometry.md#bank-calibrations)).

## Neutron analysis

`sb_loki_ana` ({sbpkg}`LOKI/app_ana/analysis_program.cc`) analyses the GRIFF files of the MCPL-based and the flood
source simulations:

```sh
sb_loki_ana lokisim.griff
sb_loki_ana *.griff               # the files of a simulation with -j
sb_loki_ana *.griff --nodetfile   # only the histograms
```

**Detection.** The model has no anode wires: neither the charge collection nor the readout is simulated. A neutron
is detected if its conversion products deposit more than 120 keV in the counting gas (a threshold against the gamma
background; 30-200 keV are reported for BCS detectors). Each neutron is simulated independently, so there is no
pile-up. The detection position is the energy-weighted mean of the energy deposits in the gas; its straw, tube and
bank (copy numbers) and its position along the straw give the pixel id (`PixelatedBanks::getPixelId`). With
`aiming_bank_id` (flood source), only the detections in the aimed banks are kept.

**Outputs.**

| File | Contents |
| --- | --- |
| `detectionEvents.mcpl` | The detection events: one particle per detection, with the time of flight in the time field and the pixel id in the kinetic energy field (`DetectionFileCreator`, {sbpkg}`LOKI/libinc/DetectionFileCreator.hh`). Its metadata holds the simulation parameters (pixels per straw, rear bank distance, bank calibration, flood source parameters, McStas metadata). |
| `bcsloki_sans.shist` | Histograms (SimpleHists); view with `sb_simplehists_browse bcsloki_sans.shist`, plot with [](plotting.md). |

The detection file uses the MCPL format for its convenient C++ writer and Python reader. It uses a universal PDG code
and weight (28 bytes per event); older files with the pixel id in the user flags can still be processed.

`sb_loki_ana_det` ({sbpkg}`LOKI/app_ana_det/analysis_program.cc`) is the variant for the detector-only simulations
(`sb_loki_sim_det`).

## Masking analysis

`sb_lokimasking_ana` ({sbpkg}`LokiMasking/app_ana/analysis_program.cc`) follows each geantino of a masking simulation.
A pixel is *visible* if a geantino enters its converter before passing an absorber: a boron mask, an absorber panel
(B4C or aluminium), the beamstop or a calibration slit mask. The other pixels are masked.

```sh
sb_lokimasking_ana loki_masking.griff
sb_lokimasking_ana *.griff        # the files of a simulation with -j
```

| File | Contents |
| --- | --- |
| `bank<N>_maskFile.xml` | One Mantid mask file per bank (only the aimed bank with `aiming_bank_id`): the detector numbers (pixel id + 1) of the masked pixels. |
| `bcsloki_masking.shist` | Histograms, e.g. the entered pixels (`h_geantino_pixel_enter`) and the visible ones (`h_geantino_pixel_enter_masked`). |

The analysis warns if not every pixel was aimed at (`MaskFileCreator`, {sbpkg}`LokiMasking/libinc/MaskFileCreator.hh`).

Separate simulations (e.g. per bank) can be combined with `sb_lokimasking_maskmerger`: it keeps the pixels masked in
all input files (a pixel visible in any simulation is visible).

```sh
sb_lokimasking_maskmerger -s merged_maskFile.xml bank0_maskFile.xml other/bank0_maskFile.xml
```
