# Simulation

## Running a simulation

A simulation is started by a script of a package: `sb_loki_sim` ({sbpkg}`LOKI/scripts/sim`) for neutrons and
`sb_lokimasking_sim` ({sbpkg}`LokiMasking/scripts/sim`) for the masking. These scripts take the geometry
({doc}`geometry`) and a launcher ({sbpkg}`LOKI/python/Launcher.py`, {sbpkg}`LokiMasking/python/Launcher.py`), which
chooses the particle generator. The parameters given on the command line (`name=value`) belong to:

| Parameters of | Defined in | Examples | List with |
| --- | --- | --- | --- |
| the geometry | `GeoBCSBanks` | `rear_detector_distance_m`, `bank_calibration`, `beamstop_id`, `with_calibration_slits` | `-g` |
| the launcher | `Launcher.py` | `event_gen`, `analysis_straw_pixel_number`, `aiming_bank_id` | `-d` (`LauncherPars`) |
| the generator | the generator | `neutron_wavelength_min_aangstrom` (flood), `input_file` (MCPL) | `-p` |

`-d` lists them all with their values (also with other parameters set, e.g. `sb_loki_sim event_gen=flood -d`).
Options of all simulation commands:

| Option | Meaning |
| --- | --- |
| `-n N` | Simulate N events (for an MCPL source: all particles in the file if not given). |
| `-j N` | Run N processes in parallel (one GRIFF file each). Leave a core free for the system. |
| `-s S` | Random seed (default 123456789). |
| `-l PL` | Physics list (default `QGSP_BIC_HP_EMZ`). |
| `-o NAME` | Name of the GRIFF output (default `lokisim`; `loki_masking` for the masking). |
| `--view`, `--data`, `--aim` | Visualisation instead of a simulation (see [](geometry.md#visualisation)). |
| `-t` | Check the geometry for overlaps and exit. |

The pixels are not part of the simulation: `analysis_straw_pixel_number` (default 256) is stored in the GRIFF file
and used by the analysis. Mantid processing needs 256 or 512 (its instrument definition files).

## Neutron generators

`sb_loki_sim` selects the generator with `event_gen`:

| `event_gen` | Generator | Use |
| --- | --- | --- |
| `mcpl` | `MCPLGen` (dgcode) | Neutrons of a McStas simulation ([](#mcpl-based-simulation)). |
| `flood` | {sbpkg}`LOKI/python/FloodSourceGen.py` | Uniform (isotropic) scattering from the sample ([](#flood-source-simulation)). |
| `spheremodel` | {sbpkg}`LOKI/python/SansSphereGen.py` | Ideal SANS sphere-model distribution. |
| `isotheta` | {sbpkg}`LOKI/python/IsoThetaGen.py` | Uniform polar angle distribution (tests). |
| (not set) | `FlexGen` (dgcode) | A narrow cone of 3 Å neutrons along the beam (change with `random_max_polarangle_deg`, `neutron_wavelength_aangstrom`, ...). |

`geantino=1` replaces the neutrons by geantinos (non-interacting particles, for visualisation), and
`primary_only=True` keeps only the primary particles in the output.

### MCPL-based simulation

The neutrons of a McStas simulation (see [](mcstas.md)) are read from an MCPL file:

```sh
sb_loki_sim event_gen=mcpl -d
sb_loki_sim -j10 event_gen=mcpl rear_detector_distance_m=10 input_file=/path/to/loki-collen_5p0m-ap1_30mm-ap2_10mm-lmin_3p0A-lmax_11p5A_rear_10p0m.mcpl.gz
```

- `input_file`: the full path of the MCPL file.
- The particles are placed relative to the sample position: the distance of the MCPL output component from the sample
  in the McStas model is read from the file (`sample_mcpl_distance_m`); files without it use 0.2 m (the LoKI McStas
  model).
- All metadata of the MCPL file (e.g. the McStas command line) is copied to the GRIFF file and to the detection file.
- `gen_x_offset_meters` shifts the beam horizontally.
- `nominal_source_sample_distance_meters` (default 23.6 m) is used only for approximate histograms in the analysis.
- A LoKI McStas simulation writes three MCPL files (front and middle banks, rear bank at 5 m, rear bank at 10 m):
  simulate each separately, with the matching `rear_detector_distance_m`.
- Use `skip_events` and `-n` to simulate a part of a file (see [](dmsc.md)).

### Flood source simulation

`FloodSourceGen` emits neutrons from the sample in random directions within a cone, with a uniform wavelength
distribution: a flat (isotropic) scatterer without a McStas simulation. The start time of each neutron is set from
its wavelength and the source to sample distance (a single, instantaneous pulse).

| Parameter | Default | Meaning |
| --- | --- | --- |
| `neutron_wavelength_min_aangstrom`, `neutron_wavelength_max_aangstrom` | 2, 13 | Wavelength range. |
| `cone_opening_deg`, `cone_opening_min_deg` | 49.6, 0 | Range of the angle from the reference direction (default: the beam). |
| `nominal_source_sample_distance_meters` | 23.6 | Source to sample distance (time of flight). |
| `source_monitor_distance_meters` | 23.599 | Position of the calculated pre-sample monitor (Mantid processing). |
| `gen_x_offset_meters`, `gen_y_offset_meters`, `gen_z_offset_meters` | 0 | Source position offset. |
| `gen_x_width_meters`, `gen_y_width_meters` | 0.006, 0.008 | Size of the source (the sample). |

The generator parameters are stored in the output, so that the Mantid processing can calculate the pre-sample
monitor spectrum (see [](mantid.md)).

`aiming_bank_id` limits the directions to a bank or a group of banks (more efficient): a bank id (0-8) aims a cone at
the centre of that bank (in the bank placement of the geometry in use), `1234` covers the middle banks and `5678` the
front banks with a minimum and maximum angle. The analysis then keeps only the detections in those banks.
`cone_view=min` or `max` emits only at the minimum or maximum angle, to check the cone visually:

```sh
sb_loki_sim event_gen=flood -d
sb_loki_sim event_gen=flood primary_only=True --data -n1000 aiming_bank_id=1234 cone_view=min
sb_loki_sim event_gen=flood --data -n1000 aiming_bank_id=3 cone_view=max geantino=1
sb_loki_sim event_gen=flood aiming_bank_id=3 -n10000 -j10        # simulation
```

## Masking simulation

`sb_lokimasking_sim` ({sbpkg}`LokiMasking/python/MaskingSourceGen.py`) shoots one geantino from the sample at the
centre of each pixel, in the order of the pixel ids. Geantinos travel in straight lines without interacting; the
analysis ([](analysis.md#masking-analysis)) finds the pixels a geantino reaches without passing an absorber (boron
masks, absorber panels, beamstop, calibration slits). Uses:

1. **Masking**: the pixels not directly visible from the sample are masked in the data reduction.
2. **Calibration slits**: with `with_calibration_slits=true`, the pixels illuminated through the slits. The
   overlapping layers see a slit at different positions along the straws (parallax), which can be used for position
   corrections.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `aiming_straw_pixel_number` | 512 | Pixels per straw for the aiming only: more than `analysis_straw_pixel_number` aims several geantinos at each pixel (better for partly shadowed pixels). |
| `analysis_straw_pixel_number` | 256 | Pixels per straw of the analysis (mask files). |
| `aiming_bank_id` | -1 | Aim only at this bank (from its first pixel); -1: all banks from pixel 0. |
| `aiming_pixel_id_min` | -1 | First pixel to aim at (overrides the start of `aiming_bank_id`). |
| `masking_only` | no | Store only the volumes needed for the masking (smaller GRIFF files). |
| `gen_x_offset_meters`, `gen_x_width_meters`, `gen_y_width_meters` | 0 | Position and size of the source. |

The aiming uses the bank placement of the geometry (`bank_calibration`). The number of events `-n` must cover the
pixels to aim at (the simulation prints it): with 1024 aiming pixels per straw, the full geometry has 6 422 528 pixels
(`sb_loki_info -s 1024 -d`). The generator supports `-j` (the processes share the pixels).

```sh
sb_lokimasking_sim -d
sb_lokimasking_sim aiming_bank_id=5 aiming_straw_pixel_number=1024 -n1024 --data
sb_lokimasking_sim aiming_straw_pixel_number=1024 -n6422528 -j10
sb_lokimasking_sim with_calibration_slits=true aiming_straw_pixel_number=1024 -n6422528 -j10
```

## Simulation and analysis in one step

The chain commands run a simulation and then the analysis of its GRIFF files (after a short test run of both), in a
new output directory:

| Command | Runs |
| --- | --- |
| `sb_loki_simanachain DIR ...` | `sb_loki_sim ...` → `sb_loki_ana` |
| `sb_loki_simanachain_det DIR ...` | `sb_loki_sim_det ...` → `sb_loki_ana_det` |
| `sb_lokimasking_simanachain DIR ...` | `sb_lokimasking_sim ...` → `sb_lokimasking_ana` |

The arguments after the directory go to the simulation. `sb_loki_sim_det` / `sb_loki_ana_det` (also `det_only=True`
of `sb_loki_sim`) store only the converter and counting gas volumes, for detection studies with small files.
