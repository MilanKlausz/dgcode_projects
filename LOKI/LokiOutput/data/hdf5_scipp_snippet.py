# Loads a LOKI simulation HDF5 file (sb_lokioutput_writehdf5, see LokiOutput/python/Hdf5Output.py) into scipp:
# sim[f'loki_detector_{b}'] with dims layer, tube, straw, pixel (the folding of the ess.loki workflow);
# values = counts (sum of weights) with their variances, coords detector_number and position (the pixel
# centres of the simulated geometry).
import h5py
import numpy as np
import scipp as sc

sim_file = 'simulation.h5'                                   # the file of sb_lokioutput_writehdf5
sim = sc.DataGroup()
with h5py.File(sim_file, 'r') as f:
    for key, g in f['entry'].items():
        sizes = {'layer': int(g.attrs['layers']), 'tube': int(g.attrs['tubes_per_layer']),
                 'straw': int(g.attrs['straws_per_tube']), 'pixel': int(g.attrs['pixels_per_straw'])}
        sim[key] = sc.DataArray(
            sc.array(dims=['detector_number'], values=g['counts'][:], variances=g['counts_variance'][:], unit='counts'),
            coords={'detector_number': sc.array(dims=['detector_number'], values=g['detector_number'][:], unit=None),
                    'position': sc.vectors(dims=['detector_number'], unit='m',
                                           values=np.stack([g[f'{a}_position'][:] for a in 'xyz'], axis=1))},
        ).fold(dim='detector_number', sizes=sizes)
