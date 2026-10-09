# Running on the DMSC cluster

:::{note}
This page is from the earlier LOKI simulation guide; the queue names, paths and submission details have not been
checked against the current cluster.
:::

Build dgcode_projects on the cluster (see [](installation.md#dmsc-cluster)). Large simulations are split into many
jobs, each running a [simulation and its analysis](simulation.md#simulation-and-analysis-in-one-step) in its own
directory, submitted with the dgcode utility `sb_dmscutils_submit` from a Python script. Templates are in
{sbpkg}`LOKI/scripts`: `submit_template_flood` and `submit_template_mcpl`.

## Flood source

```python
#!/usr/bin/env python
import os
import subprocess

baseDir = '/mnt/groupdata/detector/workspaces/<user>/lokiSimulations/'
outputdir = baseDir + '/floodSource_bank3_1e9'

nJobs, nEventsPerJob = (10000, 100000)  # total: 1e4 x 1e5 = 1e9
submitJobRange = [0, 500]  # submit only these jobs now (do not flood the cluster)

queue = 'newlong'
cmdbase = ['sb_loki_simanachain', 'results', 'event_gen=flood', 'aiming_bank_id=3']

for i in range(nJobs):
    if i < submitJobRange[0] or i >= submitJobRange[1]:
        continue
    cmd = cmdbase + ['-n%i' % nEventsPerJob, '-s%i' % (7999954 + 483 * i)]
    cmd = ['sb_dmscutils_submit', queue, os.path.join(outputdir, 'job%04i' % i)] + cmd
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    print(stdout, stderr)
```

- Batches of about 10⁵ events per job keep the memory use low.
- Each job has its own seed (`-s`).
- Raise `submitJobRange` step by step instead of submitting all jobs at once.
- A job uses one core unless the queue name ends with `/exclusive` (e.g. `newlong/exclusive`); a node has 28
  (newlong) or 32 (quark) cores.

## MCPL file

The same, with each job taking its part of the MCPL file (`skip_events`); the last job also takes the remainder. The
number of particles comes from `pymcpltool --stats` (see [](mcstas.md#geant4-input)):

```python
mcplFile = baseDir + '/mcstas/samp18/loki-collen_3p0m-ap1_30mm-ap2_10mm-lmin_3p0A-lmax_11p5A_rear_10p0m.mcpl.gz'
outputdir = baseDir + '/geant4/samp18/rear10'
rearDetDistance = 10  # [m]
nJobs, nEventsPerJob, additionalNeutrons = (21126, 100000, 70927)  # total: 2112670927 (from the MCPL file)
submitJobRange = [1000, 2000]

queue = 'newlong'
cmdbase = ['sb_loki_simanachain', 'results', 'event_gen=mcpl', 'allow_zero_weight=y',
           'analysis_straw_pixel_number=512', 'bank_calibration=nominal-geant4-geometry',
           f'rear_detector_distance_m={rearDetDistance}', f'input_file={mcplFile}', '--cleanup']

for i in range(nJobs):
    if i < submitJobRange[0] or i >= submitJobRange[1]:
        continue
    nskip = i * nEventsPerJob
    n = nEventsPerJob + (additionalNeutrons if i == nJobs - 1 else 0)
    cmd = cmdbase + ['skip_events=%i' % nskip, '-n%i' % n, '-s%i' % (7999954 + 483 * i)]
    cmd = ['sb_dmscutils_submit', queue, os.path.join(outputdir, 'job%05i' % i)] + cmd
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    print(stdout, stderr)
```

- `skip_events` is limited to 2³¹: split larger MCPL files first (see [](mcstas.md#mcpl-file-size)).
- A LoKI McStas simulation gives three MCPL files: run one submission for each.
- `bank_calibration=nominal-geant4-geometry` for results processed with Mantid (see [](mantid.md)).

## Checking and merging the results

Each job writes its detection file to `job*/results/runana/`. Count them (in the background, it can be slow) and
compare with the number of jobs:

```sh
ls -l job*/results/runana/*.mcpl | wc -l > countMCPL.txt &
```

Failed jobs have a non-zero `exitcode.txt`:

```sh
for d in job*; do [ "$(cat $d/exitcode.txt)" != "0" ] && echo $d; done
```

Resubmit them with `sb_dmscutils_resubmit <job directories>`. If that fails (e.g. `Ignoring job for which queue
status could not be determined`), delete those job directories carefully and submit them again with the submission
script (existing job directories are not overwritten, only reported). A detection file that was not closed properly
gives `MCPL WARNING: Input file appears to not have been closed properly` when merging.

Merge the detection files and the histograms (more than about 10 000 files at once can be a problem):

```sh
sb_mcpl_tool -m floodSource_10000to19999 job1*/results/runana/detectionEvents.mcpl
sb_simplehists_merge -o floodSource_10000to19999 job1*/results/runana/*.shist
```

Geant4 jobs had problems on the quark partition; the newlong and short partitions were used.
