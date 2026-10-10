# Analysis

The simulations write GRIFF files with the tracks of all particles. An analysis program (C++, in an `app_<name>`
directory of a package, run as `sb_<package>_<name>`) reads them and derives the quantities of interest; it can be
changed and rerun without repeating the simulation. The LOKI analyses use `PixelatedBanks` for the pixel ids, with the
bank placement recorded by the simulation (see [](geometry.md#bank-calibrations)). They take the GRIFF files as
arguments and write their outputs to the current directory.

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
pile-up. The pixel id of a detection is made of:

- the bank, tube and straw where the neutron **converted**: the copy numbers of the converter volume and of its
  mother volumes (a copy number is the index Geant4 gives a placed volume; here it follows the numbering of
  [](geometry.md#numbering-and-pixel-ids));
- the pixel along that straw from the **detection position**, the energy-weighted mean of the energy deposits in the
  gas: pixel = ⌊d / (L / N)⌋, with d the distance along the straw from its pixel-0 end, L the straw length of the bank
  and N the pixels per straw (clamped to 0 … N−1); `PixelatedBanks::getPixelId`.

**`aiming_bank_id`** (flood source): only the detections in the aimed banks are written to the detection file. The
histograms count the detections in all banks (an aimed cone also reaches other banks), so their totals can be larger
than the number of detection events.

**Outputs.**

| File | Contents |
| --- | --- |
| `detectionEvents.mcpl` | The detection events: one particle per detection, with the time of flight in the time field (ms) and the pixel id in the kinetic energy field (`DetectionFileCreator`, {sbpkg}`LOKI/libinc/DetectionFileCreator.hh`). Its metadata holds the simulation parameters (pixels per straw, rear bank distance, bank calibration, flood source parameters, McStas metadata). |
| `bcsloki_sans.shist` | Histograms (SimpleHists, dgcode's histogram format) and counters, e.g. `neutron_counters` (simulated, converted, detected neutrons). |

The detection file uses the MCPL format for its convenient C++ writer and Python reader. It uses a universal PDG code
and weight (28 bytes per event); older files with the pixel id in the user flags can still be processed.

`sb_loki_ana_det` ({sbpkg}`LOKI/app_ana_det/analysis_program.cc`) is the variant for the detector-only simulations
(`sb_loki_sim_det`).

## Inspecting the results

In the directory of the analysis:

```sh
python3 -m LOKI.OutputSummary          # all histograms and counters, and the number of detection events
mcpltool -j detectionEvents.mcpl       # the header of the detection file: number of events, metadata
sb_simplehists_browse bcsloki_sans.shist   # browse the histograms (see also Plotting)
sb_loki_info -p 622200                 # bank, tube, straw, pixel and position of a pixel id
```

`mcpltool` prints the pixel ids rounded (as energies); read them with Python:

```python
import mcpl
import numpy as np

f = mcpl.MCPLFile("detectionEvents.mcpl")
ids = np.concatenate([block.ekin for block in f.particle_blocks]).astype(int)      # pixel ids (0-based)
tof = np.concatenate([block.time for block in f.particle_blocks])                  # time of flight [ms]
print(f.nparticles, "detections, pixel ids", ids[:5], "times of flight [ms]", tof[:5])
```

`OutputSummary` ends with `detectionEvents.mcpl: N detection events`. In the histograms, `neutron_counters` →
`count_neutrons_hit` counts the detections in all banks.

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
| `bank<N>_maskFile.xml` | One Mantid mask file per bank (only the aimed bank with `aiming_bank_id`): the detector numbers (pixel id + 1) of the masked pixels. With calibration slits these are the pixels *not* illuminated through the slits. |
| `bcsloki_masking.shist` | Histograms: the entered pixels (`h_geantino_pixel_enter`) and the visible ones, i.e. entered without passing an absorber (`h_geantino_pixel_enter_masked`). |

The analysis warns if not every pixel was aimed at (`MaskFileCreator`, {sbpkg}`LokiMasking/libinc/MaskFileCreator.hh`).
The mask files do not record the settings they were made with: keep the pixels per straw and the bank calibration
with them (or in their file names). See [](simulation.md#masking-simulation) for the settings of masks for measured
and for simulated data. In Mantid, a mask file is applied to a workspace of the same instrument geometry and pixel
numbering (e.g. `LoadMask` with the instrument, then `MaskDetectors`).

Separate simulations (e.g. per bank) can be combined with `sb_lokimasking_maskmerger`: it keeps the pixels masked in
all input files (a pixel visible in any simulation is visible).

```sh
sb_lokimasking_maskmerger -s merged_maskFile.xml bank0_maskFile.xml other/bank0_maskFile.xml
```
