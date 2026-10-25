"""Stages B and C (Geant4, dgcode env).

B: coverage of CAD solids by Geant4 volumes: each CAD point is located in the Geant4 geometry, also shifted along the
   plate normal by up to +-SHIFT mm (thin plates: a sub-mm offset in thickness is irrelevant for the absorption); covered if
   the deepest volume's name starts with one of PREFIXES.
C: random points inside every Geant4 volume whose name starts with PREFIXES (sampled in the bounding box of the volume and
   kept when the navigator puts them in that volume), written for stage D (FreeCAD).

usage: python pts_g4.py CADPOINTS.json OUT_coverage.json OUT_g4points.json [PREFIXES=BoronMask] [SHIFT=3] [geometry par=value ...]"""
import sys, json, random
import numpy as np
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output

cadfile, outcov, outg4 = sys.argv[1:4]
opts = dict(a.split('=', 1) for a in sys.argv[4:])
prefixes = tuple(opts.pop('PREFIXES', 'BoronMask').split(','))
shift = float(opts.pop('SHIFT', '3'))
npts = int(opts.pop('NG4', '400'))
pars = dict(bank_calibration='nominal-geant4-geometry'); pars.update(opts)
geo = GEO.create()
for k, v in pars.items():
    t = type(getattr(geo, k))
    setattr(geo, k, (v.lower() in ('1', 'true', 'yes')) if t is bool else t(v))
with captured_fd_output():
    TU.construct_world(geo)


def deepest(p):
    return TU.locate(*p)[0][0]


cov = []
for s in json.load(open(cadfile)):
    n = np.array(s['normal']); hits = 0; names = {}
    for p in s['points']:
        p = np.array(p)
        for d in sorted(np.arange(-shift, shift + 1e-9, 0.5), key=abs):
            nm = deepest(p + d * n)
            if nm[0].startswith(prefixes):
                hits += 1; names[nm[0]] = names.get(nm[0], 0) + 1
                break
    cov.append(dict(label=s['label'], idx=s['idx'], coverage=hits / max(len(s['points']), 1), g4=names))
json.dump(cov, open(outcov, 'w'))

# stage C: points inside the Geant4 volumes
random.seed(2)
tree = TU.world_tree_solids(-1)
g4pts = []
for depth, pv, copy, lv, t, R, stype, p, mat in tree:
    if not pv.startswith(prefixes):
        continue
    R = np.array(R); t = np.array(t)
    half = np.array(p) / (1 if stype == 'G4Box' else 2)      # G4Box: half lengths; others: full bbox extent
    pts = []; tries = 0
    while len(pts) < npts and tries < 100 * npts:
        tries += 1
        q = t + R @ (np.array([random.uniform(-1, 1) for _ in range(3)]) * half)
        h = TU.locate(*q)[0][0]
        if h[0] == pv and h[1] == copy:
            pts.append(q.tolist())
    g4pts.append(dict(name=f'{pv}:{copy}', mat=mat, points=pts))
json.dump(g4pts, open(outg4, 'w'))
print(len(cov), 'CAD solids checked;', len(g4pts), 'Geant4 volumes sampled')
