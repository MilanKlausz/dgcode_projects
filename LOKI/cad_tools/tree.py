"""Dump the assembly tree of a STEP file (FreeCAD): indentation = depth, each leaf with its global bounding box (mm),
volume (cm3), boxness, and face count. usage: env STEP=... MAXDEPTH=... freecadcmd -c "exec(open('tree.py').read())" """
import os
import FreeCAD, Import
STEP = os.environ['STEP']; MAXD = int(os.environ.get('MAXDEPTH', '9'))
doc = FreeCAD.newDocument('t'); Import.insert(STEP, doc.Name)
roots = [o for o in doc.Objects if not o.InList]
def show(o, d):
    if d > MAXD: return
    kids = getattr(o, 'Group', []) or []
    s = getattr(o, 'Shape', None)
    info = ''
    if s is not None and not s.isNull() and s.Solids:
        try:
            sh = s.copy(); sh.Placement = o.getGlobalPlacement() if not kids else sh.Placement
        except Exception:
            sh = s
        bb = sh.BoundBox; v = sum(x.Volume for x in sh.Solids)
        info = f'  bb x {bb.XMin:8.1f}..{bb.XMax:8.1f} y {bb.YMin:8.1f}..{bb.YMax:8.1f} z {bb.ZMin:8.1f}..{bb.ZMax:8.1f}  V {v/1e3:8.1f} cm3  box {v/max(bb.XLength*bb.YLength*bb.ZLength,1e-9):.2f} n{len(sh.Solids)} f{len(sh.Faces)}'
    print('  ' * d + o.Label + (f' [{len(kids)}]' if kids else '') + info)
    for k in kids:
        show(k, d + 1)
for r in roots:
    show(r, 0)
