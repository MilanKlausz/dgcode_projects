"""FreeCAD: coarse tessellation of every solid of a STEP file (Part.read order = the inventory3 'solid' index).
env: STEP, OUT (npz: v<i> vertices, f<i> triangles), TOL (mm, default 3)"""
import os, numpy as np, Part
sh = Part.read(os.environ['STEP']); tol = float(os.environ.get('TOL', '3'))
d = {}
for i, so in enumerate(sh.Solids):
    try:
        v, f = so.tessellate(tol)
        d[f'v{i}'] = np.array([[p.x, p.y, p.z] for p in v], np.float32); d[f'f{i}'] = np.array(f, np.int32)
    except Exception:
        pass
np.savez_compressed(os.environ['OUT'], **d); print(len(d) // 2, 'solids tessellated')
