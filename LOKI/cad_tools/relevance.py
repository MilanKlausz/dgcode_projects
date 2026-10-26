"""Relevance of CAD solids for the detector: from random points inside each solid (pts_cad.py), the fraction of
directions (isotropic, Fibonacci sphere) along which a straight line reaches a detector tube (TubeWall) before an
absorber (B4C, Cd, Gd) - i.e. the chance that a neutron scattered there flies straight into a tube. The relevance of a
solid is its volume times that fraction (cm3): passive material hidden behind the B4C panels, masks or other banks
scores low. (Geometry only: the neutron flux at the solid is not included.)

usage: python relevance.py OUT.json PTS.json [PTS.json ...] [NDIR=400] [geometry par=value ...]"""
import sys, json
import numpy as np
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output

out = sys.argv[1]
files = [a for a in sys.argv[2:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[2:] if '=' in a)
ndir = int(opts.pop('NDIR', '400'))
pars = dict(bank_calibration='nominal-geant4-geometry', rear_detector_distance_m='5.0')   # the frame of the CAD alignment
pars.update(opts)
geo = GEO.create()
for k, v in pars.items():
    t = type(getattr(geo, k))
    setattr(geo, k, (v.lower() in ('1', 'true', 'yes')) if t is bool else t(v))
with captured_fd_output():
    TU.construct_world(geo)
# isotropic directions: Fibonacci sphere
i = np.arange(ndir) + 0.5
cz = 1 - 2 * i / ndir; phi = np.pi * (1 + 5 ** 0.5) * i; sz = np.sqrt(1 - cz ** 2)
dirs = np.stack([sz * np.cos(phi), sz * np.sin(phi), cz], 1).tolist()
res = []; stopmats = set()
for fn in files:
    for s in json.load(open(fn)):
        if not s['points']:
            continue
        counts, mats = TU.ray_visibility(s['points'], dirs, 8000.0, ['TubeWall'], ['B4C', 'Cd', 'Gd'], ['VacuumTank'])
        stopmats.update(mats)
        vis = float(np.mean(counts)) / ndir
        c = np.mean(s['points'], 0)
        res.append(dict(file=fn.split('/')[-1], label=s['label'], idx=s['idx'], mat=s['mat'], volume_cm3=s['volume'] / 1e3,
                        visibility=vis, relevance_cm3=s['volume'] / 1e3 * vis, centre=c.tolist()))
        print(f"{s['label']:18s} {s['mat']:6s} V {s['volume'] / 1e3:8.0f} cm3  visible tubes {vis:6.3f}  -> {s['volume'] / 1e3 * vis:7.1f} cm3"
              f"  at {np.round(c).astype(int)}", flush=True)
print('absorber materials met:', sorted(stopmats))
json.dump(dict(pars=pars, ndir=ndir, solids=res), open(out, 'w'), indent=0)
