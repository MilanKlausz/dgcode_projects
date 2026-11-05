"""HDF5 output of a LOKI simulation: the simulated counts per pixel.

  sb_lokioutput_writehdf5 INPUT [INPUT ...] -o OUT.h5 [options] [KEY=VALUE ...]

INPUT: the SimpleHists files of the analysis (sb_loki_ana: *.shist), or the output directories of a simanachain command
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
  for an AmBe source (--ambe-yield, or the yield found in the simulation logs of a simanachain directory, see AmBeRate.py):
      rate_1Ci          (n,) float64  the counts as a rate [1/s] for a 1 Ci source
  and the /entry attributes neutrons_per_source_event and rate_1Ci_per_count_per_s.

Loading with scipp (h5py + scipp, outside dgcode): LokiOutput/data/hdf5_scipp_snippet.py.
"""
import datetime
import glob
import os
import re
import tempfile

import h5py
import numpy as np
import SimpleHists as sh

import G4GeoLoki.LokiAimHelper as LokiAim
import LokiOutput.AmBeRate as AmBeRate

PIXEL_HIST = 'h_neutron_pixel_hit'
PIXEL_W2_HIST = 'h_neutron_pixel_hit_w2'  # sum of squared weights, with the geometry in its comment
NUMBER_OF_BANKS = 9
LAYERS = 4
STRAWS_PER_TUBE = 7


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
    """The analysis output of an input: the file itself, or the *.shist of a simanachain directory."""
    if os.path.isfile(path):
        return path
    files = sorted(glob.glob(os.path.join(path, 'runana', '*.shist'))) or sorted(glob.glob(os.path.join(path, '*.shist')))
    if len(files) != 1:
        raise ValueError(f'{path}: expected one *.shist file (in {path}/runana or {path}), found {len(files)}')
    return files[0]


def _read_output(path):
    """(shist file, counts, squared weights or None, events, recorded geometry or None, pixels per straw) of an input."""
    shist = find_shist(path)
    collection = sh.HistCollection(shist)
    if PIXEL_HIST not in collection.keys:
        raise ValueError(f'{shist}: no {PIXEL_HIST} histogram (not an output of sb_loki_ana)')
    hist = collection.hist(PIXEL_HIST)
    w2, geometry = None, None
    if PIXEL_W2_HIST in collection.keys:
        w2hist = collection.hist(PIXEL_W2_HIST)
        w2, geometry = _pixel_order(w2hist), parse_geometry_comment(w2hist.comment)
    counters = {counter.label: counter.value for counter in collection.hist('neutron_counters').counters}
    return shist, _pixel_order(hist), w2, counters['count_initial_neutrons'], geometry, int(hist.nbinsx)


def read_counts(paths, given_geometry=None):
    """Add up the pixel counts of the analysis outputs. given_geometry: for outputs without the recorded geometry (older
    analyses): rear_detector_distance_m and bank_calibration (the straw pixel number is the histogram binning)."""
    counts = variances = recorded = None
    n_events = 0.0
    sources = []
    for path in paths:
        shist, c, w2, events, geometry, pixels_per_straw = _read_output(path)
        n_events += events
        sources.append(shist)
        if counts is None:
            counts, variances, recorded = c, w2, geometry
            continue
        if (geometry is None) != (recorded is None):
            raise ValueError(f'{shist}: older analysis output (no {PIXEL_W2_HIST}) mixed with newer ones')
        if geometry != recorded:
            raise ValueError(f'{shist}: another geometry than {sources[0]} (bank calibration, rear detector distance '
                             'or pixels per straw); only outputs of the same geometry can be added up')
        if c.shape != counts.shape:
            raise ValueError(f'{shist}: {c.size} pixels, but {sources[0]} has {counts.size}')
        counts = counts + c
        if w2 is not None:
            variances = variances + w2
    if recorded is not None:
        if given_geometry:
            raise ValueError(f'{sources[0]}: the geometry is recorded in the analysis output; do not give it')
        return SimulatedCounts(counts, variances, n_events, recorded, sources)
    if not given_geometry or 'rear_detector_distance_m' not in given_geometry or 'bank_calibration' not in given_geometry:
        raise ValueError(f'{sources[0]}: no recorded geometry (an analysis from before {PIXEL_W2_HIST}): give the '
                         'rear detector distance and the bank calibration of the simulation')
    geometry = dict(rear_detector_distance_m=float(given_geometry['rear_detector_distance_m']),
                    straw_pixel_number=pixels_per_straw, bank_calibration=given_geometry['bank_calibration'],
                    bank_calibration_text=None)
    return SimulatedCounts(counts, variances, n_events, geometry, sources)


def aim_helper(geometry):
    """The AimHelper of the simulated geometry."""
    rear_mm = geometry['rear_detector_distance_m'] * 1000.0
    pixels = geometry['straw_pixel_number']
    text = geometry['bank_calibration_text']
    if not text:  # the nominal geometry, or an older output (by name)
        return LokiAim.AimHelper(rear_mm, pixels, NUMBER_OF_BANKS, geometry['bank_calibration'])
    with tempfile.NamedTemporaryFile('w', suffix='.txt') as f:  # the placements recorded by the simulation
        f.write(text)
        f.flush()
        return LokiAim.AimHelper(rear_mm, pixels, NUMBER_OF_BANKS, f.name)


def pixel_positions(helper, pixels):
    """The pixel centres [m], pixel-id order, (pixels in the geometry, 3). The pixels of a straw are equally spaced: the
    first and last pixel centre of each straw from the geometry, linearly interpolated (much faster than every pixel)."""
    n = helper.getTotalNumberOfPixels()
    first = np.arange(0, n, pixels)
    start = np.array([helper.getPixelCentreCoordinates(int(i)) for i in first], dtype=np.float64)
    end = np.array([helper.getPixelCentreCoordinates(int(i) + pixels - 1) for i in first], dtype=np.float64)
    t = (np.arange(pixels) / max(pixels - 1, 1))[None, :, None]
    return (start[:, None, :] + t * (end - start)[:, None, :]).reshape(n, 3) / 1000.0


def _write_bank(entry, bank, helper, counts, variances, positions, pixels, scaled_counts):
    first, n = helper.getBankPixelOffset(bank), helper.getNumberOfPixels(bank)
    tubes = n // (STRAWS_PER_TUBE * pixels)
    local = np.arange(n, dtype=np.int64)
    kw = dict(compression='gzip', compression_opts=4)
    g = entry.create_group(f'loki_detector_{bank}')
    g.attrs['layers'] = LAYERS
    g.attrs['tubes_per_layer'] = tubes // LAYERS
    g.attrs['straws_per_tube'] = STRAWS_PER_TUBE
    g.attrs['pixels_per_straw'] = pixels
    g.create_dataset('detector_number', data=first + 1 + local, **kw)
    g.create_dataset('counts', data=counts[first:first + n], **kw)
    g.create_dataset('counts_variance', data=variances[first:first + n], **kw)
    for name, factor in scaled_counts.items():
        g.create_dataset(name, data=counts[first:first + n] * factor, **kw)
    for axis, name in enumerate('xyz'):
        g.create_dataset(f'{name}_position', data=positions[first:first + n, axis], **kw)
    g.create_dataset('tube', data=(local // (STRAWS_PER_TUBE * pixels)).astype(np.int32), **kw)
    g.create_dataset('straw', data=(local // pixels % STRAWS_PER_TUBE).astype(np.int32), **kw)
    g.create_dataset('pixel', data=(local % pixels).astype(np.int32), **kw)


def write(counts, out, attributes=None, scaled_counts=None):
    """Write the HDF5 file (see the module documentation). scaled_counts: {dataset name: factor}: further datasets of the
    banks, the counts times the factor."""
    geometry = counts.geometry
    pixels = geometry['straw_pixel_number']
    helper = aim_helper(geometry)
    positions = pixel_positions(helper, pixels)
    if positions.shape[0] != counts.counts.size:
        raise ValueError(f'the geometry has {positions.shape[0]} pixels, the analysis output {counts.counts.size}')
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
        for key, value in (attributes or {}).items():
            e.attrs[key] = value
        variances = counts.variances if counts.variances is not None else counts.counts
        for bank in range(NUMBER_OF_BANKS):
            _write_bank(e, bank, helper, counts.counts, variances, positions, pixels, scaled_counts or {})


def _is_attribute(argument):
    """KEY=VALUE, and not an input file with '=' in its name."""
    return re.match(r'\w+=', argument) is not None and not os.path.exists(argument)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(prog='sb_lokioutput_writehdf5', description=__doc__.split('\n')[0],
                                     epilog='See the documentation of LokiOutput/python/Hdf5Output.py for the file layout.')
    parser.add_argument('inputs', nargs='+', metavar='INPUT|KEY=VALUE',
                        help='analysis outputs (*.shist or simanachain directories); KEY=VALUE: /entry attributes')
    parser.add_argument('-o', '--output', required=True, help='the HDF5 file to write')
    parser.add_argument('--ambe-yield', type=float, default=None,
                        help='AmBe source: neutrons per source neutron (printed by AmBeSourceGen), for the rate at 1 Ci; '
                             'by default taken from the simulation logs of simanachain directories')
    parser.add_argument('--rear-detector-distance', type=float, default=None, metavar='M',
                        help='only for analysis outputs without the recorded geometry: rear_detector_distance_m')
    parser.add_argument('--bank-calibration', default=None, metavar='NAME',
                        help='only for analysis outputs without the recorded geometry: bank_calibration')
    args = parser.parse_intermixed_args(argv)
    inputs = [i for i in args.inputs if not _is_attribute(i)]
    attributes = dict(i.split('=', 1) for i in args.inputs if _is_attribute(i))
    if not inputs:
        parser.error('no input')
    given_geometry = {}
    if args.rear_detector_distance is not None:
        given_geometry['rear_detector_distance_m'] = args.rear_detector_distance
    if args.bank_calibration is not None:
        given_geometry['bank_calibration'] = args.bank_calibration
    rate_note = ''
    scaled_counts = {}
    try:
        counts = read_counts(inputs, given_geometry)
        ambe_yield = args.ambe_yield if args.ambe_yield is not None else AmBeRate.yield_from_logs(inputs)
        if ambe_yield is not None and counts.n_events > 0:
            rate = AmBeRate.rate_per_count(counts.n_events, ambe_yield)
            scaled_counts[AmBeRate.RATE_DATASET] = rate
            attributes = dict(neutrons_per_source_event=ambe_yield, rate_1Ci_per_count_per_s=rate, **attributes)
            rate_note = f'; AmBe: 1 count = {rate:.6g} /s at 1 Ci'
        write(counts, args.output, attributes, scaled_counts)
    except (ValueError, KeyError) as error:
        parser.exit(1, f'ERROR: {error}\n')
    print(f'Created {args.output}: {len(inputs)} analysis output(s), {counts.n_events:g} source events, '
          f'{counts.counts.sum():g} counts (sum of weights), bank calibration {counts.geometry["bank_calibration"]}, '
          f'rear detector distance {counts.geometry["rear_detector_distance_m"]:g} m, '
          f'{counts.geometry["straw_pixel_number"]} pixels per straw{rate_note}')


if __name__ == '__main__':
    main()
