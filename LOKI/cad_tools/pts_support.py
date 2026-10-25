"""Stage D (FreeCAD): support of Geant4 volumes by CAD solids: the fraction of each Geant4 volume's points within TOL mm
of a CAD solid of the material(s) MAT (TOL: thin-plate offsets are irrelevant).
env: TAG, STEP, ALIGN, MAT, G4PTS (stage C json), TOL (default 3), OUT"""
import os, sys, json
sys.path.insert(0, os.environ.get('CAD_TOOLS', '.'))   # the directory of these scripts
import numpy as np
import FreeCAD, Part
from FreeCAD import Vector, Matrix
import cadlib

tag = os.environ['TAG']; mats = os.environ.get('MAT', 'B4C').split(','); tol = float(os.environ.get('TOL', '3'))
al = np.load(os.environ['ALIGN']); R, t = al['R'], al['t']
M = Matrix(R[0, 0], R[0, 1], R[0, 2], t[0], R[1, 0], R[1, 1], R[1, 2], t[1], R[2, 0], R[2, 1], R[2, 2], t[2], 0, 0, 0, 1)
allshape = Part.read(os.environ['STEP'])
solids = []
REG = [float(x) for x in os.environ['REGION'].split(',')] if os.environ.get('REGION') else None   # CAD x0,x1,y0,y1
def inreg(r):
    return REG is None or (REG[0] < r['centre'][0] < REG[1] and REG[2] < r['centre'][1] < REG[3])
for r in cadlib.load(tag):
    if r['mat'] in mats and inreg(r):
        so = allshape.Solids[r['solid']].copy(); so.transformShape(M)
        bb = so.BoundBox; bb.enlarge(tol)
        solids.append((r['label'], so, bb))
res = []
for g in json.load(open(os.environ['G4PTS'])):
    ok = 0; who = {}
    for p in g['points']:
        v = Vector(*p)
        for lab, so, bb in solids:
            if bb.isInside(v) and (so.isInside(v, 1e-6, True) or so.distToShape(Part.Vertex(v))[0] <= tol):
                ok += 1; who[lab] = who.get(lab, 0) + 1
                break
    res.append(dict(name=g['name'], support=ok / max(len(g['points']), 1), cad=who, n=len(g['points'])))
    print(f"{g['name']:30s} support {ok / max(len(g['points']), 1):5.3f}  by {who}")
json.dump(res, open(os.environ['OUT'], 'w'))
