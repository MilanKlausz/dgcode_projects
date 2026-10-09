# McStas simulation

:::{note}
The McStas models are not part of dgcode_projects: [loki-mcstas-master](https://bitbucket.org/MKlausz/loki-mcstas-master/src/master/)
(LoKI) and [larmor-mcstas](https://github.com/ess-loki-sim/larmor-mcstas) (Larmor). This page is from the earlier
LOKI simulation guide and has not been checked against the current versions of these repositories.
:::

The McStas simulation transports the neutrons from the moderator to the sample and saves the scattered neutrons in
MCPL files, the source of the Geant4 simulation (`sb_loki_sim event_gen=mcpl`, see
[](simulation.md#mcpl-based-simulation)). The McStas monitors give the incident spectra for the Mantid processing.

## Particle weights

McStas uses particle weights; the Geant4 simulation and the data reduction need whole neutrons (weight 1). The
`MCPL_output_deweight` component converts the weights when writing the MCPL file:

- **Splitting**: while the weight p ≥ 1, a neutron of weight 1 is written and p is reduced by 1.
- **Russian roulette**: for 0 < p < 1, a neutron of weight 1 is written with probability p (`if (p > rand) p = 1;
  else p = 0;`).

This increases the variance; it is best if the weights are already below 1 when they are written, so that only the
Russian roulette is applied. The weights depend on the source intensity and the number of simulated neutrons: with
1000 n/s and 100 simulated neutrons the mean weight is 10. To lower the weights without simulating more neutrons than
the intensity, `SPLIT n` before the sample makes n copies of each neutron with 1/n of its weight. This improves the
sampling of the sample scattering at the cost of a coarser sampling of the instrument before it, so it should not be
overused. In practice: increase the number of neutrons and the SPLIT without deweighting until the MCPL file has more
particles than the intensity, then use the same numbers with deweighting.

## MCPL file size

Keep the number of particles in an MCPL file below 2³¹ (≈ 2.15·10⁹; the `skip_events` limit of the Geant4
submission, see [](dmsc.md)); that is already 78 GB compressed. To keep only useful neutrons:

- A slit after the sample (`focusAroundDetector`) passes only neutrons that can reach the detectors. Check its size
  with the Geant4 visualisation of the MCPL file (a few thousand particles, `--data primary_only=True`), and compare
  the detection counts of a full simulation with a larger slit.
- The beamstop removes the transmitted beam (earlier versions also stopped the unscattered neutrons in McStas).

Split larger files with `mcpltool` (`sb_mcpl_tool` in dgcode): `-e` extracts, `-sN` skips the first N particles,
`-lN` limits the number:

```sh
mcpltool -e larmor_postsample.mcpl.gz larmor_postsample_after2e9 -s2000000000
```

## Geant4 input

The MCPL tools come with McStas (`pymcpltool`, `mcpltool`) and with dgcode (`sb_mcpl_pytool`, `sb_mcpl_tool`). The
number of particles of a file is needed for the Geant4 submission:

```sh
pymcpltool --stats larmor_postsample.mcpl.gz > stats.txt &   # slow for ~1e9 particles
```

`MCPL_output_deweight` also stores the McStas command line and the instrument file in the MCPL file:

```sh
pymcpltool file.mcpl.gz -bmccode_cmd_line
pymcpltool file.mcpl.gz -bmccode_instr_file
pymcpltool file.mcpl.gz -bsample_mcpl_distance_m
```

`sample_mcpl_distance_m` is the distance of the MCPL output component from the sample; `sb_loki_sim` reads it and
places the particles accordingly. The particles start at the post-sample slit, not at the sample (seen in the
visualisation).

## LoKI model

`loki_master_model.instr`, main parameters:

| Parameter | Meaning |
| --- | --- |
| `acc_power` | Accelerator power in MW (default 5), the intensity of the `ESS_butterfly` source. |
| `deweight` | Apply the deweighting. |
| `sampletype` | The sample: negative = empty beam (no sample, no beamstop), 0 = flat scatterer (`Isotropic_Sqw`), 10 = logarithmically spaced sharp peaks in Q (`SANS_benchmark2_RKH`), 18 = Gaussian coil, 20 = monodisperse hard spheres, R = 50 (`Sans_spheres`). |
| `collen` | Collimation length: 3, 5 or 8 m. |
| `bs5`, `bs10` | Beamstop (1-5, see [](geometry.md#beamstops)) for the rear bank at 5 and 10 m. |

**MCPL output.** One simulation writes three files with the components `mcpl_out_frontMidBanks`,
`mcpl_out_frontRearBank_5m` and `mcpl_out_frontRearBank_10m` (0.2 m after the sample): the neutrons towards the
front and middle banks (independent of the rear bank distance), and those reaching the rear bank at 5 m and at 10 m
past their beamstops. An empty beam simulation writes one file.

**Monitors** (loaded by the Mantid processing):

| Monitor | Position (default) | File |
| --- | --- | --- |
| `Mon10_PostFOC` | after the frame overlap chopper (15.4905 m) | `Mon10_PostFOC_*.t` |
| `PreSampleMonitor` | before the sample (23.5921 m) | `PreSampleMonitor_*.t` |
| `PostSampleMonitor` | after the sample (23.8322 m) | `PostSampleMonitor_*.t` |
| `beamstopMonitor_5m`, `beamstopMonitor_10m` | at the beamstops | `beamstopMonitor_5m_*.t`, `beamstopMonitor_10m_*.t` |

**Parameters and resulting numbers** on the DMSC cluster (apertures `sourceapx=0.030 sourceapy=0.025 sampleapx=0.010
sampleapy=0.015`; front/middle banks and rear bank files). Local runs with 1/100 of the neutrons and of the power give
1/100 of these numbers.

| Sample | collen [m] | Power [MW] | SPLIT | Neutrons | bs5 / bs10 | Sum of weights in MCPL | Detection events |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Empty beam | 3 | 5 | 1 | 1e11 | - | 1.4e9 | 7.7e8 |
| | 5 | 5 | 1 | 1e11 | - | 5.3e8 | 2.9e8 |
| | 8 | 5 | 1 | 1e11 | - | 2.05e8 | 1.14e8 |
| Flat scatterer | 3 | 50 | 1000 | 1e11 | 5 / - | 1.2e9 / 1.2e7 | 4.2e8 / 4.8e6 |
| | 5 | 200 | 10000 | 1e11 | 3 / 5 | 1.8e9 / 1.8e7 | 6.4e8 / 7.2e6 |
| | 8 | 500 | 10000 | 1e11 | 3 / 4 | 1.8e9 / 1.7e7 | 6e8 / 7e6 |
| Peaks in Q | 3 | 50k | 1000 | 1e11 | 5 / - | 1.2e9 / 7.4e8 | 5e9 / 4e8 |
| | 5 | 200k | 10000 | 1e11 | 3 / 5 | 1.8e9 / 1.5e9 | 7.6e8 / 8.4e8 |
| | 8 | 500k | 10000 | 1e11 | 3 / 4 | 1.7e9 / 1.4-1.7e9 | 7.2e8 / 7.8-9.5e8 |
| Gaussian coil | 3 | 500 | 1000 | 1e11 | 5 / - | 2.5e9 / 1.4e9 | 1.1e9 / 7.0e8 |
| | 5 | 1000 | 1000 | 1e11 | 3 / 5 | 1.8e9 / 1.1e9 | 8e8 / 5.8e8 |
| | 8 | 500 | 1000 | 1e11 | 3 / 4 | 3.5e8 / 2.2e8 | 1.5e8 / 1.1-1.2e8 |
| Hard spheres | 3 | 200k | 4000 | 1e11 | 5 / - | 3.11e8 / 1.94e9 | 1.3e8 / 1.1e9 |
| | 5 | 500k | 10000 | 1e11 | 3 / 5 | 2.85e8 / 2.01e9 | 1.2e8 / 1.1e9 |
| | 8 | 1250k | 10000 | 2e11 | 3 / 4 | 2.70e8 / 1.93-1.98e9 | 1.16e8 / 1.1e9 |

## Larmor model

`Larmor.instr` (with `ViewModISISver1.comp`, which sets the intensity), main parameters:

| Parameter | Meaning |
| --- | --- |
| `sampletype` | As for LoKI (negative, 0, 10, 18, 20). Samples: 10 × 10 × 1 mm (the empty sample slit 50 × 50 mm). |
| `deweight` | Apply the deweighting. |
| `current` | Beam current in A (default 2); scales the intensity (`beamCurrent` of `ViewModISISver1.comp`). |
| `sampdet` | Sample to detector distance (default 4.099 m, the 2022 test). |
| `s2` | A1 aperture size (default 0.02 m); the A2 aperture is fixed at 6 × 8 mm. |

The SPLIT is after the `Samplepoint` arm. The beamstop is at `sampdet - 0.09` (not in empty beam runs), the MCPL output
(`mcpl_out_postsample`, file `larmor_postsample.mcpl.gz`) at `sampdet - 0.05` = 4.049 m from the sample. Monitors:
`TOFmon1` (after the chopper, 9.8195 m), `TOFmon2` (before the A1 aperture, 20.313 m), `TOFpresamp` (before the
sample, 25.57 m), `TOFmon4` (after the sample, 25.76 m), files `<name>.dat`.

| Sample | Current [A] | SPLIT | Neutrons | Sum of weights in MCPL | Detection events |
| --- | --- | --- | --- | --- | --- |
| Empty beam | 0.020 | 1000 | 1e11 | 4.8e8 | 2.2e8 |
| Flat scatterer | 20 | 1000 | 1e11 | 1.5e9 (expected) | 6e8 (expected) |
| Peaks in Q | 2000 | 1000 | 1e11 | 2.1e9 | 1.2e9 |
| Gaussian coil | 20 | 1000 | 1e11 | 3.9e9 | 2.1e9 |

## Running McStas on the DMSC cluster

Load the modules (see [](installation.md#dmsc-cluster)) and compile for MPI on the right node (`compile` for the
short and newlong partitions, `quarkcompile` for quark): start a run with `-c --mpi=4` and stop it (Ctrl+C) when it
asks for the parameters:

```sh
mcrun -c --mpi=4 Larmor.instr -n1e6 -dtest
```

Then submit a batch script with `sbatch mpirun.batch`:

```sh
#!/bin/bash
#SBATCH --mail-user=<email>
#SBATCH --job-name=larmor
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=quark
#SBATCH --nodes=1-1
#SBATCH --time=24:00:00

mpirun Larmor.out -dsamp18_split1000_1e11 sampletype=18 deweight=1 sampdet=4.099 -n1e11
```

- Use the compiled `.out` file in `mpirun`. Keep earlier commands in the file as comments.
- Do not use dots in the output directory name (the merging of the MCPL files fails).
- The merging and compression of the MCPL files can take longer than the simulation and does not use more nodes:
  use one node (28 cores on newlong, 32 on quark).
- If the merging fails, repeat the simulation.
- Partition time limits (`sinfo`): short 4 hours, newlong 7 days, quark 1 day.
- Larmor on 1 quark node: peaks in Q 5 h + compression (43 GB), Gaussian coil 2.5 h (78 GB), empty beam 1 h (6.6 GB).
