"""Metadata of the MCPL files of the LoKI McStas -> Geant4 chain, from the file header only (no particle is read).

  python3 -m LOKI.McplMetadata FILE [--get NAME ... | --json | --intensity]

Reads the McStas output files (LOKI_main, MCPL_output: <name>_frontmid / <name>_rear) and the detection files of the
analysis (detectionEvents.mcpl, which keep the McStas header as blobs next to the Geant4 user data):

- the McStas version and instrument, the output component and which file it is (frontmid / rear),
- the instrument parameters of the McStas run (the mccode_cmd_line blob: the values used),
- the distance of the output plane from the sample (the AT of the output component in the instrument file blob),
- the number of particles, the normalisation of the weights (initial_ray_count: sum of the weights / it = n/s),
- for detection files, the recorded Geant4 settings (the GRIFF user data).

--get prints the values of fields or McStas parameters, one per line (for the submission scripts), e.g.
  python3 -m LOKI.McplMetadata run_frontmid.mcpl.gz --get nparticles reardet beamstop
--intensity reads all particles for the sum of the weights / initial_ray_count (n/s).

Needs only the mcpl Python module, so it can be used outside dgcode as well.
"""

import dataclasses
import json
import re
import shlex
import sys

MCSTAS_SOURCE = re.compile(r'^mcstas\s+(?P<version>\S+)\s+(?P<instrument>\S+)$', re.IGNORECASE)
OUTPUT_COMPONENT = re.compile(r'^Output by COMPONENT:\s*(?P<name>\S+)$')
USERFLAGS_COMMENT = re.compile(r'^userflags:\s*(?P<text>.*)$')
DETECTION_SOURCE = 'LOKI/DetectionFileCreator Class'
OLD_PLANE_DISTANCE_BLOB = 'sample_mcpl_distance_m'  # written by the earlier LoKI models
SAMPLE_COMPONENT = 'Sample'  # the origin of the Geant4 geometry (the nominal sample position)


@dataclasses.dataclass
class McplMetadata:
    path: str
    kind: str                           # 'mcstas', 'detection' or 'other'
    source: str                         # the MCPL source name
    nparticles: int                     # particles (detection events in a detection file)
    universal_weight: float             # 0: each particle has its own weight
    comments: list
    mcstas_version: str = None
    instrument: str = None
    component: str = None               # the MCPL_output component, e.g. mcpl_out_rear
    output: str = None                  # 'frontmid' or 'rear' (the end of the component name)
    initial_ray_count: float = None     # the McStas ncount: sum of the weights / it = n/s
    stat_sum: dict = dataclasses.field(default_factory=dict)
    userflags: str = None
    cmd_line: str = None
    parameters: dict = dataclasses.field(default_factory=dict)  # McStas instrument parameters (as given, str)
    plane_distance_m: float = None      # output plane from the sample (z of the component AT, RELATIVE Sample)
    geant4: dict = dataclasses.field(default_factory=dict)      # detection files: the Geant4 user data

    def parameter(self, name):
        """A McStas parameter as a number (int or float) if it is one, else as given; None if it is not there."""
        return _number(self.parameters[name]) if name in self.parameters else None

    def get(self, name):
        """A field, else a McStas parameter, else a Geant4 setting (for --get)."""
        if name in (f.name for f in dataclasses.fields(self)):
            return getattr(self, name)
        if name in self.parameters:
            return self.parameter(name)
        if name in self.geant4:
            return _number(self.geant4[name])
        raise KeyError(f'{self.path}: no field, McStas parameter or Geant4 setting named "{name}"')

    def summary(self):
        lines = [f'{self.path}: {self.kind} file, {self.nparticles} {"events" if self.kind == "detection" else "particles"}',
                 f'  source: {self.source}']
        if self.instrument:
            lines.append(f'  McStas {self.mcstas_version}, instrument {self.instrument}, output {self.component}'
                         f' ({self.output or "?"}), plane at {self.plane_distance_m} m from the sample')
        if self.initial_ray_count is not None:
            lines.append(f'  initial_ray_count (ncount): {self.initial_ray_count:g}')
        if self.parameters:
            lines.append('  McStas parameters: ' + ' '.join(f'{k}={v}' for k, v in self.parameters.items()))
        if self.geant4:
            lines.append('  Geant4 settings: ' + ' '.join(f'{k}={v}' for k, v in self.geant4.items()
                                                          if '\n' not in v and len(v) < 80))
        return '\n'.join(lines)


def _number(text):
    for kind in (int, float):
        try:
            return kind(text)
        except (TypeError, ValueError):
            pass
    return text


def _text_blobs(mcplfile):
    blobs = {}
    for key, value in mcplfile.blobs.items():
        try:
            blobs[key] = value.decode('utf-8')
        except UnicodeDecodeError:
            pass  # a binary blob (e.g. the configuration stored by the MCPL writer of dgcode)
    return blobs


def user_data(path):
    """The header of an MCPL file as {key: text} for the GRIFF user data; read() reads it back from detection files."""
    import mcpl
    mcplfile = mcpl.MCPLFile(str(path))
    data = _text_blobs(mcplfile)
    # the weights of the file: the universal weight (0: each particle has its own, e.g. the McStas weights), and the
    # normalisation of the McStas weights (stat:sum initial_ray_count: sum of the weights / it = neutrons/s)
    data['mcpl_universal_weight'] = str(mcplfile.opt_universalweight)
    for key, value in mcplfile.stat_sum.items():
        if value is not None:  # None: not set (the run did not finish)
            data[f'mcpl_stat_sum_{key}'] = repr(value)
    data['mcpl_comments'] = '\n'.join(mcplfile.comments)
    data['mcpl_source'] = mcplfile.sourcename
    return data


def parse_cmd_line(cmd_line):
    """The name=value tokens of a McStas command line (the instrument parameters, as given)."""
    return dict(token.split('=', 1) for token in shlex.split(cmd_line) if '=' in token and not token.startswith('-'))


def plane_distance(instr_text, component):
    """The z of the AT of a component placed RELATIVE Sample with a literal position, else None."""
    if not instr_text or not component:
        return None
    code = re.sub(r'//[^\n]*|/\*.*?\*/', '', instr_text, flags=re.DOTALL)
    match = re.search(rf'\bCOMPONENT\s+{re.escape(component)}\s*=.*?\bAT\s*\(([^)]*)\)\s*RELATIVE\s+(\w+)', code,
                      flags=re.DOTALL)
    if not match or match.group(2) != SAMPLE_COMPONENT:
        return None
    position = [_number(c.strip()) for c in match.group(1).split(',')]
    if len(position) != 3 or not all(isinstance(c, (int, float)) for c in position) or position[:2] != [0, 0]:
        return None
    return float(position[2])


def read(path):
    """The metadata of an MCPL file (header only)."""
    import mcpl
    mcplfile = mcpl.MCPLFile(str(path))
    blobs = _text_blobs(mcplfile)
    comments = list(mcplfile.comments)
    kind = 'detection' if mcplfile.sourcename == DETECTION_SOURCE else 'other'
    stat_sum = dict(mcplfile.stat_sum)
    if kind == 'detection':  # the McStas header, as the Launcher stored it in the user data
        comments = blobs.get('mcpl_comments', '').splitlines()
        stat_sum = {key[len('mcpl_stat_sum_'):]: float(value) for key, value in blobs.items()
                    if key.startswith('mcpl_stat_sum_')}
    meta = McplMetadata(path=str(path), kind=kind, source=mcplfile.sourcename, nparticles=mcplfile.nparticles,
                        universal_weight=mcplfile.opt_universalweight, comments=comments, stat_sum=stat_sum)
    for comment in comments:
        if match := OUTPUT_COMPONENT.match(comment):
            meta.component = match['name']
        elif match := USERFLAGS_COMMENT.match(comment):
            meta.userflags = match['text']
    if match := MCSTAS_SOURCE.match(mcplfile.sourcename):
        meta.kind = 'mcstas'
        meta.mcstas_version, meta.instrument = match['version'], match['instrument']
    elif match := MCSTAS_SOURCE.match(blobs.get('mcpl_source', '')):  # detection file: the source of the McStas file
        meta.mcstas_version, meta.instrument = match['version'], match['instrument']
    if meta.component:
        meta.output = next((name for name in ('frontmid', 'rear') if meta.component.endswith(name)), None)
    if stat_sum.get('initial_ray_count'):  # None: not set (the run did not finish)
        meta.initial_ray_count = stat_sum['initial_ray_count']
    meta.cmd_line = blobs.get('mccode_cmd_line')
    if meta.cmd_line:
        meta.parameters = parse_cmd_line(meta.cmd_line)
    meta.plane_distance_m = plane_distance(blobs.get('mccode_instr_file'), meta.component)
    if meta.plane_distance_m is None and OLD_PLANE_DISTANCE_BLOB in blobs:
        meta.plane_distance_m = float(blobs[OLD_PLANE_DISTANCE_BLOB])
    if kind == 'detection':
        meta.geant4 = {key: value for key, value in blobs.items() if not key.startswith(('mccode_', 'mcpl_'))}
    return meta


def intensity(path):
    """Sum of the weights / initial_ray_count (n/s for a McStas file); reads all particles."""
    import mcpl
    meta = read(path)
    if meta.initial_ray_count is None:
        raise ValueError(f'{path}: no initial_ray_count (not a McStas weight_mode=1 file, or an unfinished run)')
    mcplfile = mcpl.MCPLFile(str(path))
    total = sum(float(block.weight.sum()) for block in mcplfile.particle_blocks)
    return total / meta.initial_ray_count


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(prog='python3 -m LOKI.McplMetadata', description=__doc__.split('\n')[0])
    parser.add_argument('file')
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--get', nargs='+', metavar='NAME', help='print these values, one per line')
    output.add_argument('--json', action='store_true', help='print all metadata as JSON')
    output.add_argument('--intensity', action='store_true', help='also read the particles: sum of weights / ncount')
    args = parser.parse_args(argv)
    meta = read(args.file)
    if args.get:
        try:
            print('\n'.join(str(meta.get(name)) for name in args.get))
        except KeyError as error:
            parser.exit(1, f'ERROR: {error.args[0]}\n')
        return
    if args.json:
        print(json.dumps(dataclasses.asdict(meta), indent=1))
    else:
        print(meta.summary())
    if args.intensity:
        print(f'  intensity: {intensity(args.file):.6g} n/s')


if __name__ == '__main__':
    sys.exit(main())
