"""FreeCAD: cross-sections of elongated, non-box solids (candidates for extrusions: beams, channels, tubes).

For every solid of the inventory (OBB_JSON) with obb_fill < MAXFILL, not a plate, and longest side >= 3 x the middle
side: sections perpendicular to its longest OBB axis at 25, 50 and 75% of the length. Output (JSON, per solid index):
long axis, in-plane axes, OBB centre and extents, and per section its area and polygons in the in-plane coordinates
(relative to the OBB centre). An extrusion has three (nearly) equal sections.
env: STEP, OBB_JSON, OUT, MAXFILL (0.8)"""
import os, json
import Part
from FreeCAD import Vector

rows = json.load(open(os.environ['OBB_JSON'])); maxfill = float(os.environ.get('MAXFILL', '0.8'))
sh = Part.read(os.environ['STEP'])
out = {}
for r in rows:
    size = [e[1] - e[0] for e in r['ext']]
    s = sorted(size)
    if r['obb_fill'] >= maxfill or s[0] < 0.25 * s[1] or s[2] < 3 * s[1]:
        continue
    k = size.index(s[2]); others = [i for i in range(3) if i != k]
    L = Vector(*r['axes'][k]); A = Vector(*r['axes'][others[0]]); B = Vector(*r['axes'][others[1]])
    com = Vector(*r['com'])
    centre = com + L * (0.5 * (r['ext'][k][0] + r['ext'][k][1])) + A * (0.5 * (r['ext'][others[0]][0] + r['ext'][others[0]][1])) \
        + B * (0.5 * (r['ext'][others[1]][0] + r['ext'][others[1]][1]))
    so = sh.Solids[r['solid']]
    secs = []
    for f in (0.25, 0.5, 0.75):
        p = com + L * (r['ext'][k][0] + f * (r['ext'][k][1] - r['ext'][k][0]))
        d = L.dot(p)
        try:
            wires = so.slice(L, d)
            areas = sorted((abs(Part.Face(w).Area) for w in wires if w.isClosed()), reverse=True)
            area = (areas[0] - sum(areas[1:])) if areas else 0.0          # outer contour minus holes (one level)
            polys = [[[(q - centre).dot(A), (q - centre).dot(B)] for q in w.discretize(Distance=1.0)] for w in wires]
        except Exception:
            area, polys = 0.0, []
        secs.append(dict(area=area, polys=polys))
    out[r['solid']] = dict(label=r['label'], long=[L.x, L.y, L.z], a=[A.x, A.y, A.z], b=[B.x, B.y, B.z],
                           centre=[centre.x, centre.y, centre.z], length=s[2], vol=r['vol'], sections=secs)
json.dump(out, open(os.environ['OUT'], 'w'))
print(len(out), 'elongated non-box solids sectioned')
