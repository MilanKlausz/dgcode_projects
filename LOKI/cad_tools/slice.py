"""FreeCAD: planar section of selected solids (Part.read order index) at a plane z = Z0 (CAD frame): polygons per solid.
env: STEP, IDX (comma list), Z0, OUT (json: {idx: [[x,y],...] per wire})"""
import os, json, Part, FreeCAD
from FreeCAD import Vector
sh = Part.read(os.environ['STEP']); z0 = float(os.environ['Z0'])
out = {}
for i in [int(x) for x in os.environ['IDX'].split(',')]:
    so = sh.Solids[i]
    wires = so.slice(Vector(0, 0, 1), z0)
    out[i] = [[[p.x, p.y] for p in w.discretize(Distance=0.5)] for w in wires]
json.dump(out, open(os.environ['OUT'], 'w')); print(len(out), 'sections')
