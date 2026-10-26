"""Stage A (FreeCAD): random points inside every CAD solid of a material, in the Geant4 world frame.

env: TAG (cadlib inventory tag), STEP, ALIGN (npz: world = R cad + t), MAT (default B4C), N (points per solid, default 400),
     OUT (json: [{label, idx, normal, points}]), REGION (CAD x0,x1,y0,y1: solids with their centre inside), MINVOL (cm3)"""
import os, sys, json, random
sys.path.insert(0, os.environ.get('CAD_TOOLS', '.'))   # the directory of these scripts
import numpy as np
import FreeCAD, Part
from FreeCAD import Vector, Matrix
import cadlib

tag = os.environ['TAG']; mats = os.environ.get('MAT', 'B4C').split(','); N = int(os.environ.get('N', '400'))
al = np.load(os.environ['ALIGN']); R, t = al['R'], al['t']
M = Matrix(R[0, 0], R[0, 1], R[0, 2], t[0], R[1, 0], R[1, 1], R[1, 2], t[1], R[2, 0], R[2, 1], R[2, 2], t[2], 0, 0, 0, 1)
REG = [float(x) for x in os.environ['REGION'].split(',')] if os.environ.get('REGION') else None   # CAD x0,x1,y0,y1
def inreg(r):
    return REG is None or (REG[0] < r['centre'][0] < REG[1] and REG[2] < r['centre'][1] < REG[3])
rows = [r for r in cadlib.load(tag) if r['mat'] in mats and inreg(r) and r['vol'] >= float(os.environ.get('MINVOL', '0')) * 1e3]
allshape = Part.read(os.environ['STEP'])
random.seed(1)
out = []
for r in rows:
    so = allshape.Solids[r['solid']].copy(); so.transformShape(M)
    bb = so.BoundBox
    pts = []
    tries = 0
    while len(pts) < N and tries < 200 * N:
        tries += 1
        p = Vector(random.uniform(bb.XMin, bb.XMax), random.uniform(bb.YMin, bb.YMax), random.uniform(bb.ZMin, bb.ZMax))
        if so.isInside(p, 1e-6, True):
            pts.append([p.x, p.y, p.z])
    thin = int(np.argmin(r['size']))
    normal = R @ np.array(r['axes'][thin])
    out.append(dict(label=r['label'], idx=r['solid'], mat=r['mat'], normal=normal.tolist(), thickness=float(r['size'][thin]),
                    volume=r['vol'], points=pts))
json.dump(out, open(os.environ['OUT'], 'w'))
print(len(out), 'solids,', sum(len(o['points']) for o in out), 'points')
