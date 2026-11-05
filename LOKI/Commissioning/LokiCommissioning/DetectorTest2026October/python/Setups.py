"""The setups of the LoKI detector test of October 2026 (data/setups.toml) as parameters of sb_lokiambe_sim.

The parameters of a setup: those of [common], then absorber_strips and absorber_plates (the items that were on the
detector in that setup, joined with ';' in the order of the file), then the parameters of the setup itself (they
replace the common ones).
"""
import hashlib
import tomllib

from Core import FindData

CATALOGUE = 'setups.toml'
ABSORBER_PARAMETERS = ('absorber_strips', 'absorber_plates')
SETUP_INFO = ('runs', 'description')  # setup entries that are not simulation parameters


def catalogue_path(path=None):
    return path or FindData('DetectorTest2026October', CATALOGUE)


def load(path=None):
    """The catalogue (checked)."""
    with open(catalogue_path(path), 'rb') as f:
        catalogue = tomllib.load(f)
    setups = catalogue.get('setups', {})
    if not setups:
        raise ValueError(f'{CATALOGUE}: no setups')
    for parameter in ABSORBER_PARAMETERS:
        for item in catalogue.get(parameter, []):
            unknown = set(item.get('setups', [])) - set(setups)
            if unknown or 'value' not in item:
                raise ValueError(f'{CATALOGUE}: {parameter} item {item.get("name")!r}: no value or unknown setups {sorted(unknown)}')
    return catalogue


def checksum(path=None):
    """The first 12 hex digits of the SHA-256 of the catalogue file (recorded with the simulation)."""
    with open(catalogue_path(path), 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()[:12]


def names(catalogue=None):
    return list((catalogue or load())['setups'])


def _text(value):
    return value if isinstance(value, str) else repr(value)


def parameters(name, catalogue=None):
    """{parameter: value} (strings) of the setup, in the order of the catalogue."""
    catalogue = catalogue or load()
    if name not in catalogue['setups']:
        raise ValueError(f'unknown setup {name!r} (setups: {", ".join(catalogue["setups"])})')
    out = {key: _text(value) for key, value in catalogue.get('common', {}).items()}
    for parameter in ABSORBER_PARAMETERS:
        items = [item['value'] for item in catalogue.get(parameter, []) if name in item['setups']]
        if items:
            out[parameter] = ';'.join(items)
    for key, value in catalogue['setups'][name].items():
        if key not in SETUP_INFO:
            out.pop(key, None)  # a setup parameter replaces a common one (and takes the place after the common ones)
            out[key] = _text(value)
    return out


def expand(args, catalogue=None):
    """(setup name or None, args): setup=NAME in the command line arguments replaced by the parameters of the setup
    (key=value); a parameter given on the command line replaces the one of the setup."""
    setups = [a for a in args if a.startswith('setup=')]
    if not setups:
        return None, list(args)
    if len(setups) > 1:
        raise ValueError('setup= given more than once')
    name = setups[0][len('setup='):]
    given = {a.split('=', 1)[0] for a in args if '=' in a and not a.startswith('-')}
    expanded = [f'{key}={value}' for key, value in parameters(name, catalogue).items() if key not in given]
    return name, expanded + [a for a in args if a != setups[0]]
