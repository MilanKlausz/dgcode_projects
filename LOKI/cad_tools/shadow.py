"""Impact metric for (unmodelled) CAD solids: the fraction of a solid's points that lie in front of a detector tube as
seen from a source point (default the sample position, the origin): the ray source -> point, continued beyond the point,
enters a TubeWall within MAXDIST mm. For each solid also the tubes' banks.

usage: python shadow.py CADPOINTS.json OUT.json [SRC=x,y,z mm] [MAXDIST=800] [STEP=4] [geometry par=value ...]"""
import sys, json
import numpy as np
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output

cadfile, out = sys.argv[1:3]
opts = dict(a.split('=', 1) for a in sys.argv[3:])
src = np.array([float(x) for x in opts.pop('SRC', '0,0,0').split(',')])
maxd = float(opts.pop('MAXDIST', '800')); step = float(opts.pop('STEP', '4'))
pars = dict(bank_calibration='nominal-geant4-geometry'); pars.update(opts)
geo = GEO.create()
for k, v in pars.items():
    t = type(getattr(geo, k))
    setattr(geo, k, (v.lower() in ('1', 'true', 'yes')) if t is bool else t(v))
with captured_fd_output():
    TU.construct_world(geo)
TUBE = ('TubeWall', 'EmptyTube', 'StrawWall', 'Converter', 'CountingGas')
res = []
for s in json.load(open(cadfile)):
    nhit = 0; nfree = 0; banks = {}
    for p in s['points']:
        p = np.array(p); L = np.linalg.norm(p - src); u = (p - src) / L
        for d in np.arange(step, maxd, step):
            h = TU.locate(*(p + d * u))[0]
            if h[0][0].startswith(TUBE):
                b = [x[1] for x in h if x[0].startswith('Bank') or x[2] == 'Bank']
                nhit += 1; banks[str(b[0] if b else '?')] = banks.get(str(b[0] if b else '?'), 0) + 1
                # not shielded: no Geant4 boron mask on the way source -> point (stop 6 mm before the point)
                if not any(TU.locate(*(src + x * u))[0][0][0].startswith('BoronMask') for x in np.arange(0, L - 6, 1.0)):
                    nfree += 1
                break
    res.append(dict(label=s['label'], idx=s['idx'], mat=s.get('mat'), volume=s.get('volume'),
                    front_fraction=nhit / max(len(s['points']), 1), unshielded_fraction=nfree / max(len(s['points']), 1), banks=banks))
    print(f"{s['label']:16s} {s.get('mat',''):6s} V {s.get('volume', 0)/1e3:8.1f} cm3  in front of tubes: {nhit / max(len(s['points']), 1):5.2f}  (not behind a boron mask: {nfree / max(len(s['points']), 1):5.2f})  banks {banks}", flush=True)
json.dump(res, open(out, 'w'))
