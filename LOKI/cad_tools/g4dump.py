"""Dump the Geant4 LOKI geometry (G4GeoLoki GeoBCSBanks) with solids, transforms and materials to JSON.
usage: python g4dump.py OUT.json [par=value ...]   (default: nominal bank placement, all masks)"""
import sys, json
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output
out = sys.argv[1]
pars = dict(bank_calibration='nominal-geant4-geometry')
for a in sys.argv[2:]:
    k, v = a.split('=', 1)
    pars[k] = v
geo = GEO.create()
for k, v in pars.items():
    t = type(getattr(geo, k))
    setattr(geo, k, (v.lower() in ('1', 'true', 'yes')) if t is bool else t(v))
with captured_fd_output() as lines:
    TU.construct_world(geo)
bad = [l for l in lines if l.startswith('Checking overlaps') and not l.rstrip().endswith('OK!')]
rows = []
for depth, pv, copy, lv, t, R, stype, p, mat in TU.world_tree_solids(-1):
    if lv in ('TubeWall', 'EmptyTube', 'StrawWall', 'Converter', 'CountingGas'):
        continue
    rows.append(dict(depth=depth, pv=pv, copy=copy, lv=lv, t=t, R=R, solid=stype, pars=p, mat=mat))
json.dump(dict(pars=pars, overlaps=bad, volumes=rows), open(out, 'w'))
print(len(rows), 'volumes;', len(bad), 'overlap problems')
