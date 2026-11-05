"""HDF5 output of a LOKI simulation: the simulated counts per pixel, in the format of the AmBe source study files.

  sb_lokioutput_writehdf5 INPUT [INPUT ...] -o OUT.h5 [options] [KEY=VALUE ...]

INPUT: the SimpleHists files of the analysis (sb_loki_ana: *.shist), or the output directories of sb_loki_simanachain
(the *.shist file of their runana directory). Several inputs (e.g. the jobs of a cluster run) are added up; they must
have the same geometry. KEY=VALUE pairs are written as attributes of /entry (e.g. description=...).

The file (h5py and scipp friendly; NeXus-like names, but not a NeXus file):
  /entry                              attributes: description, n_source_events, the geometry, KEY=VALUE, ...
  /entry/loki_detector_<b>/           b = 0..8 (all banks); attributes: layers, tubes_per_layer, straws_per_tube,
                                      pixels_per_straw (the folding of the ess.loki workflow: layer, tube, straw, pixel)
      detector_number   (n,) int64    the pixel IDs (the 0-based pixel id of the Geant4 geometry + 1, as in the ICD and
                                      the measured NeXus files)
      counts            (n,) float64  the simulated counts: the sum of the weights of the hits
      counts_variance   (n,) float64  their variance: the sum of the squared weights (= counts for unweighted events)
      x_position, y_position, z_position (n,) float64 [m]
                                      the pixel centres of the simulated geometry (its bank calibration and rear detector
                                      distance): x horizontal (positive to the left looking along the beam), y up, z
                                      along the beam, origin at the sample position
      tube, straw, pixel (n,) int32   the indices in the bank (detector_number - first = (tube * 7 + straw) * pixels + pixel;
                                      tube = layer * tubes_per_layer + row)
  with --ambe-yield (or found in the simulation log of an AmBe run directory), the AmBe source normalisation:
      rate_1Ci          (n,) float64  the counts as a rate [1/s] for a 1 Ci AmBe source (2.2e6 neutrons/s)
  and the /entry attributes neutrons_per_source_event and rate_1Ci_per_count_per_s.

The files of the AmBe source study (2026) had the same layout, but with the pixel positions of the measured files and
x/y/z_pixel_offset; this output has the positions of the simulated geometry and no pixel offsets.

Loading with scipp (h5py + scipp, outside dgcode): LokiOutput/data/hdf5_scipp_snippet.py.
"""

import datetime
import glob
import gzip
import os
import re
import tempfile

import numpy as np

PIXEL_HIST = 'h_neutron_pixel_hit'
PIXEL_W2_HIST = 'h_neutron_pixel_hit_w2'  # sum of squared weights, with the geometry in its comment
AMBE_NEUTRONS_PER_S_1CI = 2.2e6
NUMBER_OF_BANKS = 9
LAYERS = 4
STRAWS_PER_TUBE = 7
AMBE_YIELD_LINE = re.compile(r'AmBeSourceGen: \S+: ([0-9.eE+-]+) neutrons per source neutron')


class SimulatedCounts:
    """The summed counts per pixel of one or more analysis outputs, with the geometry of their pixel ids."""

    def __init__(self, counts, variances, n_events, geometry, sources):
        self.counts = counts          # (pixels,) sum of weights, pixel-id order
        self.variances = variances    # (pixels,) sum of squared weights, or None (older outputs)
        self.n_events = n_events      # generated events (count_initial_neutrons)
        self.geometry = geometry      # rear_detector_distance_m, straw_pixel_number, bank_calibration[_text]
        self.sources = sources        # the input files


def parse_geometry_comment(comment):
    """The geometry recorded in the comment of h_neutron_pixel_hit_w2 by the analysis."""
    head, _, text = comment.partition('bank_calibration_text:\n')
    geometry = dict(line.split('=', 1) for line in head.splitlines() if '=' in line)
    geometry['rear_detector_distance_m'] = float(geometry['rear_detector_distance_m'])
    geometry['straw_pixel_number'] = int(geometry['straw_pixel_number'])
    geometry['bank_calibration_text'] = text
    return geometry


def _pixel_order(hist):
    """The contents of a (pixel along the straw, straw) histogram in pixel-id order (straw * pixels + pixel)."""
    return np.array(hist.contents(), dtype=np.float64).T.reshape(-1)


def find_shist(path):
    """The analysis output of an input: the file itself, or the *.shist of a sb_loki_simanachain directory."""
    if os.path.isfile(path):
        return path
    files = sorted(glob.glob(os.path.join(path, 'runana', '*.shist'))) or sorted(glob.glob(os.path.join(path, '*.shist')))
    if len(files) != 1:
        raise ValueError(f'{path}: expected one *.shist file (in {path}/runana or {path}), found {len(files)}')
    return files[0]


def read_counts(paths, geometry=None):
    """Add up the pixel counts of the analysis outputs. geometry: for outputs without the recorded geometry (older
    analyses): rear_detector_distance_m and bank_calibration (the straw pixel number is the histogram binning)."""
    import SimpleHists as sh
    counts = variances = None
    n_events = 0.0
    recorded = None
    sources = []
    for path in paths:
        shist = find_shist(path)
        sources.append(shist)
        collection = sh.HistCollection(shist)
        if PIXEL_HIST not in collection.keys:
            raise ValueError(f'{shist}: no {PIXEL_HIST} histogram (not an output of sb_loki_ana)')
        hist = collection.hist(PIXEL_HIST)
        c = _pixel_order(hist)
        if PIXEL_W2_HIST in collection.keys:
            w2hist = collection.hist(PIXEL_W2_HIST)
            this = parse_geometry_comment(w2hist.comment)
            if recorded is not None and this != recorded:
                raise ValueError(f'{shist}: another geometry than {sources[0]} (bank calibration, rear detector distance '
                                 'or pixels per straw); only outputs of the same geometry can be added up')
            recorded = this
            w2 = _pixel_order(w2hist)
        else:
            if recorded is not None or (variances is not None):
                raise ValueError(f'{shist}: older analysis output (no {PIXEL_W2_HIST}) mixed with newer ones')
            w2 = None
        if counts is None:
            counts, variances = c, w2
        else:
            if c.shape != counts.shape:
                raise ValueError(f'{shist}: {c.size} pixels, but {sources[0]} has {counts.size}')
            counts = counts + c
            if w2 is not None:
                variances = variances + w2
        counters = {counter.label: counter.value for counter in collection.hist('neutron_counters').counters}
        n_events += counters['count_initial_neutrons']
    pixels_per_straw = int(sh.HistCollection(sources[0]).hist(PIXEL_HIST).nbinsx)
    if recorded is None:
        if not geometry or 'rear_detector_distance_m' not in geometry or 'bank_calibration' not in geometry:
            raise ValueError(f'{sources[0]}: no recorded geometry (an analysis from before {PIXEL_W2_HIST}): give the '
                             'rear detector distance and the bank calibration of the simulation')
        recorded = dict(rear_detector_distance_m=float(geometry['rear_detector_distance_m']),
                        straw_pixel_number=pixels_per_straw, bank_calibration=geometry['bank_calibration'],
                        bank_calibration_text=None)
    elif geometry:
        raise ValueError(f'{sources[0]}: the geometry is recorded in the analysis output; do not give it')
    return SimulatedCounts(counts, variances, n_events, recorded, sources)


def ambe_yield_from_logs(paths):
    """The AmBe yield (neutrons per source neutron) printed by AmBeSourceGen in the simulation logs of sb_loki_simanachain
    directories; None if there is no such directory or line."""
    values = set()
    for path in paths:
        for log in glob.glob(os.path.join(path, 'runsim', 'output.log*')):
            text = gzip.open(log, 'rt').read() if log.endswith('.gz') else open(log).read()
            values.update(float(v) for v in AMBE_YIELD_LINE.findall(text))
    if len(values) > 1:
        raise ValueError(f'different AmBe yields in the simulation logs: {sorted(values)}; give --ambe-yield')
    return values.pop() if values else None


def pixel_positions(geometry, exact=False):
    """The pixel centres [m] (pixel-id order, (pixels, 3)) and the AimHelper of the simulated geometry."""
    import G4GeoLoki.LokiAimHelper as LokiAim
    rear_mm = geometry['rear_detector_distance_m'] * 1000.0
    pixels = geometry['straw_pixel_number']
    text = geometry['bank_calibration_text']
    if text:  # the placements recorded by the simulation (the same even if the calibration file has changed)
        with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as f:
            f.write(text)
        try:
            helper = LokiAim.AimHelper(rear_mm, pixels, NUMBER_OF_BANKS, f.name)
        finally:
            os.unlink(f.name)
    else:  # the nominal geometry, or an older output (by name)
        helper = LokiAim.AimHelper(rear_mm, pixels, NUMBER_OF_BANKS, geometry['bank_calibration'])
    n = helper.getTotalNumberOfPixels()
    if not exact:  # the pixels of a straw are equally spaced: its first and last pixel centre, linearly interpolated
        first = np.arange(0, n, pixels)
        start = np.array([helper.getPixelCentreCoordinates(int(i)) for i in first], dtype=np.float64)
        end = np.array([helper.getPixelCentreCoordinates(int(i) + pixels - 1) for i in first], dtype=np.float64)
        t = (np.arange(pixels) / max(pixels - 1, 1))[None, :, None]
        positions = (start[:, None, :] + t * (end - start)[:, None, :]).reshape(n, 3) / 1000.0
    else:  # every pixel centre from the geometry (for the tests)
        positions = np.array([helper.getPixelCentreCoordinates(i) for i in range(n)], dtype=np.float64) / 1000.0
    return positions, helper


def write(counts, out, ambe_yield=None, attributes=None):
    """Write the HDF5 file (see the module documentation)."""
    try:
        import h5py
    except ImportError:
        raise SystemExit('ERROR: the HDF5 output needs the Python module h5py (e.g. conda install -c conda-forge h5py)')
    geometry = counts.geometry
    positions, helper = pixel_positions(geometry)
    if positions.shape[0] != counts.counts.size:
        raise ValueError(f'the geometry has {positions.shape[0]} pixels, the analysis output {counts.counts.size}')
    pixels = geometry['straw_pixel_number']
    rate_per_count = (ambe_yield * AMBE_NEUTRONS_PER_S_1CI / counts.n_events
                      if ambe_yield is not None and counts.n_events > 0 else None)
    kw = dict(compression='gzip', compression_opts=4)
    with h5py.File(out, 'w') as f:
        e = f.create_group('entry')
        e.attrs['description'] = ('Geant4 (dgcode LOKI) simulation: simulated counts per pixel ID (the sum of the hit '
                                  'weights), without any readout or calibration effect; pixel positions of the '
                                  'simulated geometry')
        e.attrs['created'] = datetime.datetime.now().isoformat(timespec='seconds')
        e.attrs['analysis_files'] = '\n'.join(os.path.abspath(s) for s in counts.sources)
        e.attrs['n_source_events'] = counts.n_events
        e.attrs['weighted'] = bool(counts.variances is not None and not np.array_equal(counts.variances, counts.counts))
        e.attrs['bank_calibration'] = geometry['bank_calibration']
        if geometry['bank_calibration_text']:
            e.attrs['bank_calibration_text'] = geometry['bank_calibration_text']
        e.attrs['rear_detector_distance_m'] = geometry['rear_detector_distance_m']
        e.attrs['pixels_per_straw'] = pixels
        if counts.variances is None:
            e.attrs['counts_variance_note'] = ('analysis output without the squared weights: counts_variance = counts '
                                               '(right for unweighted events only)')
        if rate_per_count is not None:
            e.attrs['neutrons_per_source_event'] = ambe_yield
            e.attrs['rate_1Ci_per_count_per_s'] = rate_per_count
        for key, value in (attributes or {}).items():
            e.attrs[key] = value
        variances = counts.variances if counts.variances is not None else counts.counts
        for b in range(NUMBER_OF_BANKS):
            first, n = helper.getBankPixelOffset(b), helper.getNumberOfPixels(b)
            tubes = n // (STRAWS_PER_TUBE * pixels)
            local = np.arange(n, dtype=np.int64)
            g = e.create_group(f'loki_detector_{b}')
            g.attrs['layers'] = LAYERS
            g.attrs['tubes_per_layer'] = tubes // LAYERS
            g.attrs['straws_per_tube'] = STRAWS_PER_TUBE
            g.attrs['pixels_per_straw'] = pixels
            g.create_dataset('detector_number', data=first + 1 + local, **kw)
            g.create_dataset('counts', data=counts.counts[first:first + n], **kw)
            g.create_dataset('counts_variance', data=variances[first:first + n], **kw)
            if rate_per_count is not None:
                g.create_dataset('rate_1Ci', data=counts.counts[first:first + n] * rate_per_count, **kw)
            for axis, name in enumerate('xyz'):
                g.create_dataset(f'{name}_position', data=positions[first:first + n, axis], **kw)
            g.create_dataset('tube', data=(local // (STRAWS_PER_TUBE * pixels)).astype(np.int32), **kw)
            g.create_dataset('straw', data=(local // pixels % STRAWS_PER_TUBE).astype(np.int32), **kw)
            g.create_dataset('pixel', data=(local % pixels).astype(np.int32), **kw)
    return rate_per_count


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(prog='sb_lokioutput_writehdf5', description=__doc__.split('\n')[0],
                                     epilog='See the documentation of LokiOutput/python/Hdf5Output.py for the file layout.')
    parser.add_argument('inputs', nargs='+', metavar='INPUT|KEY=VALUE',
                        help='analysis outputs (*.shist or sb_loki_simanachain directories); KEY=VALUE: /entry attributes')
    parser.add_argument('-o', '--output', required=True, help='the HDF5 file to write')
    parser.add_argument('--ambe-yield', type=float, default=None,
                        help='AmBe normalisation: neutrons per source neutron (printed by AmBeSourceGen); by default '
                             'taken from the simulation logs of sb_loki_simanachain directories')
    parser.add_argument('--rear-detector-distance', type=float, default=None, metavar='M',
                        help='only for analysis outputs without the recorded geometry: rear_detector_distance_m')
    parser.add_argument('--bank-calibration', default=None, metavar='NAME',
                        help='only for analysis outputs without the recorded geometry: bank_calibration')
    args = parser.parse_intermixed_args(argv)
    inputs = [i for i in args.inputs if '=' not in os.path.basename(i)]
    attributes = dict(i.split('=', 1) for i in args.inputs if '=' in os.path.basename(i))
    geometry = {}
    if args.rear_detector_distance is not None:
        geometry['rear_detector_distance_m'] = args.rear_detector_distance
    if args.bank_calibration is not None:
        geometry['bank_calibration'] = args.bank_calibration
    try:
        counts = read_counts(inputs, geometry)
        ambe_yield = args.ambe_yield if args.ambe_yield is not None else ambe_yield_from_logs(inputs)
        rate_per_count = write(counts, args.output, ambe_yield, attributes)
    except (ValueError, KeyError) as error:
        parser.exit(1, f'ERROR: {error}\n')
    print(f'Created {args.output}: {len(inputs)} analysis output(s), {counts.n_events:g} source events, '
          f'{counts.counts.sum():g} counts (sum of weights), bank calibration {counts.geometry["bank_calibration"]}, '
          f'rear detector distance {counts.geometry["rear_detector_distance_m"]:g} m, '
          f'{counts.geometry["straw_pixel_number"]} pixels per straw'
          + (f'; AmBe: 1 count = {rate_per_count:.6g} /s at 1 Ci' if rate_per_count is not None else ''))


if __name__ == '__main__':
    main()
