#!/bin/bash
# Example: the HDF5 files of the final simulations of the detector test (setups A-G, 4,000M source neutrons each; setup H
# has the geometry and source position of setup E). The simulations ran on LUMI on 2026-10-29/30: their jobs are in
# <case>_jobs.tar.gz (40000 *.shist files and logs, no GRIFF files) in /scratch/project_465003341/klauszmi/loki/prod.
#
# The analysis outputs of these runs are older than the recorded geometry: the geometry below (the same for all setups)
# must be given by hand, and sb_lokioutput_writehdf5 warns about it. Without --bank-shifts banks 2 and 4 would be 570 mm
# off. All generated neutrons have weight 1, so counts_variance = counts is right. --ambe-yield: from the simulation logs.
#
# Reading the jobs takes about 1.5 hours per setup: run this on a cluster, in the dgcode environment with the project
# built (sb), in a directory with room for the extracted files. The long list of inputs needs Linux (macOS allows less).

PROD=${1:-/scratch/project_465003341/klauszmi/loki/prod}

for case in fin_A345_4G fin_B349_4G fin_C358_4G fin_D359_4G fin_E360_4G fin_F361_v3 fin_G363_4G; do
    tar xzf "$PROD/${case}_jobs.tar.gz" --wildcards "$case/job*.shist" || exit 1
    sb_lokioutput_writehdf5 "$case"/job*.shist -o "$case.h5" \
        --rear-detector-distance 9.985615 \
        --bank-calibration 2026-September-SAM-606 \
        --bank-shifts '2:0,0,570;4:0,0,570' \
        --ambe-yield 0.4034 \
        "description=LoKI detector test 2026 October, final simulation $case" || exit 1
    rm -r "$case"
done
