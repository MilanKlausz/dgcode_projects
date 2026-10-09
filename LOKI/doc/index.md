# LOKI

Geant4 simulation of the boron-coated straw (BCS) detectors of
[LoKI](https://ess.eu/instruments/loki), the small-angle neutron scattering instrument of the
European Spallation Source (ESS): the detector geometry, the simulation and analysis of neutron detection, masking
of the pixels that are not directly visible from the sample, and the processing of the results with Mantid. It also
contains the models of the detector tests at the Larmor instrument (ISIS) in 2020 and 2022.

```{toctree}
:maxdepth: 1

installation
geometry
simulation
analysis
mcstas
dmsc
mantid
larmor
testing
plotting
```

## Packages

| Package | Contents |
| --- | --- |
| {sbpkg}`G4GeoLoki` | The detector geometry (`GeoBCSBanks`), the geometry classes (banks, packs, tubes, masks, pixels), bank calibrations, `AimHelper` (pixel positions), pixel file export and JSON verification. |
| {sbpkg}`LOKI` | Neutron simulation (`sb_loki_sim`) with the McStas (MCPL), flood source and other generators, the analysis (`sb_loki_ana`), the geometry info tool (`sb_loki_info`), plotting, DMSC submission templates. |
| {sbpkg}`LokiMasking` | Masking simulation with geantinos (`sb_lokimasking_sim`), the analysis writing Mantid mask files (`sb_lokimasking_ana`), merging of mask files. |
| {sbpkg}`LokiMantid` | Conversion of detection events to Mantid NeXus files (`sb_lokimantid_loki`); built only if `mantidpython` is available. |
| {sbpkg}`G4GeoLokiTests` | Tests of the geometry, the pixel ids, the aiming, the masking chain and the bank calibrations. |
| {sbpkg}`Larmor`, {sbpkg}`LarmorTests`, {sbpkg}`LarmorMantid` | The Larmor 2020 and 2022 detector tests (in `LOKI/Larmor/`), see [](larmor.md). |

## Workflows

```{figure} images/workflow.webp
:alt: The simulation and data workflow

The simulation and data workflow
```

The full workflow starts with a McStas simulation of the instrument from the moderator to the sample. The neutrons
scattered by the sample are saved in an MCPL file, which is the source of the Geant4 simulation of the detectors. The
Geant4 simulation writes a GRIFF file; its analysis writes the detection events (pixel id and time of flight) to an
MCPL file, which a Mantid script converts to NeXus files, together with the monitor spectra of the McStas simulation.

| Workflow | Source of the neutrons | Steps | Result |
| --- | --- | --- | --- |
| McStas + Geant4 | MCPL file of a McStas simulation | [McStas](mcstas.md) → `sb_loki_sim event_gen=mcpl` → `sb_loki_ana` → `sb_lokimantid_loki` | NeXus files (detectors, monitors) |
| Flood source | `FloodSourceGen`: uniform scattering from the sample, no McStas | `sb_loki_sim event_gen=flood` → `sb_loki_ana` → `sb_lokimantid_loki` | NeXus files (detectors, a calculated pre-sample monitor) |
| Masking | `MaskingSourceGen`: one geantino per pixel centre | `sb_lokimasking_sim` → `sb_lokimasking_ana` | Mantid mask files (one per bank) |

## Abbreviations

BCS
: Boron-coated straws

DMSC
: Data Management and Software Centre (its computing cluster)

GRIFF
: Geant4 Results In Friendly Format (the dgcode simulation output)

ICD
: Interface Control Document (of the LoKI detectors)

IDF
: Instrument Definition File (Mantid)

MCPL
: Monte Carlo Particle List

TOF
: Time of flight
