"""The counts of an AmBe source simulation as a rate at 1 Ci, for the HDF5 output.

One simulated neutron stands for 1 / yield source neutrons (yield: neutrons leaving the moderator per source neutron,
printed by AmBeSourceGen of LokiAmBe), so N simulated neutrons are N / (yield * 2.2e6) seconds of a 1 Ci source.
"""
import glob
import gzip
import os
import re

NEUTRONS_PER_S_1CI = 2.2e6
RATE_DATASET = 'rate_1Ci'
YIELD_LINE = re.compile(r'AmBeSourceGen: \S+: ([0-9.eE+-]+) neutrons per source neutron')


def yield_from_logs(paths):
    """The yield in the simulation logs of simanachain directories; None if there is no such directory or line."""
    values = set()
    for path in paths:
        for log in glob.glob(os.path.join(path, 'runsim*', 'output.log*')):
            with (gzip.open(log, 'rt') if log.endswith('.gz') else open(log)) as f:
                values.update(float(v) for v in YIELD_LINE.findall(f.read()))
    if len(values) > 1:
        raise ValueError(f'different AmBe yields in the simulation logs: {sorted(values)}; give --ambe-yield')
    return values.pop() if values else None


def rate_per_count(n_events, ambe_yield):
    """The rate [1/s] at 1 Ci of one count."""
    return ambe_yield * NEUTRONS_PER_S_1CI / n_events
