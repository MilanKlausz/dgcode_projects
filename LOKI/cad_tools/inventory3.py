"""STEP inventory with ALL instances and oriented bounding boxes (FreeCAD).

FreeCAD's (legacy) Import keeps the product names but drops repeated instances of a product; Part.read keeps every
instance (with its placement) but no names. So: every solid of Part.read, named by its geometric signature (volume and
principal moments of inertia, which are the same for every instance of a product) from the named solids of Import.

For every solid: label, volume (mm3), centre of mass (mm), principal axes, extents along them (= an oriented bounding
box, OBB), obb_fill = volume / OBB volume. Fasteners are skipped. Output: JSON list (same format as inventory2.py).
usage: env STEP=... OUT=... freecadcmd -c "exec(open('inventory3.py').read())" """
import os
import json
import FreeCAD
import Import
import Part

STEP, OUT = os.environ['STEP'], os.environ['OUT']
# colours per product (step_colours.py): ties between a leaf part and a merged repeated sub-assembly (FreeCAD's legacy
# import merges repeated sub-assemblies into one compound named after the assembly) go to the single-colour product
COLS = json.load(open(os.environ['COL'])) if os.environ.get('COL') else {}


def ncol(label):
    n = len(COLS.get(label, {}))
    return n if n > 0 else 99
SKIP = ('SCREW', 'WASHER', 'NUT', 'BOLT')


def signature(so):
    pp = so.PrincipalProperties
    m = sorted(pp['Moments'])
    return (so.Volume, m[0], m[1], m[2])


def same(a, b, rtol=2e-4):
    return all(abs(x - y) <= rtol * max(abs(x), abs(y), 1.0) for x, y in zip(a, b))


doc = FreeCAD.newDocument('named')
Import.insert(STEP, doc.Name)
named = []
for o in doc.Objects:
    if not hasattr(o, 'Shape') or o.TypeId == 'App::Part':
        continue
    s = o.Shape
    if s.isNull():
        continue
    for so in s.Solids:
        try:
            named.append((signature(so), o.Label, len(s.Solids)))   # (assembly features also hold their parts' solids)
        except Exception:
            pass
named.sort(key=lambda x: x[0][0])
vols = [n[0][0] for n in named]

import bisect
rows = []
allshape = Part.read(STEP)
for k, so in enumerate(allshape.Solids):
    try:
        sig = signature(so)
    except Exception:
        continue
    label, nbest = '?', (10 ** 9, 99)
    i = bisect.bisect_left(vols, sig[0] * (1 - 2e-4))
    while i < len(named) and vols[i] <= sig[0] * (1 + 2e-4) + 1e-9:
        score = (named[i][2], ncol(named[i][1]))
        if same(named[i][0], sig) and score < nbest:     # prefer the leaf part (a single solid, one colour)
            label, nbest = named[i][1], score
        i += 1
    if any(w in label for w in SKIP):
        continue
    pp = so.PrincipalProperties
    axes = [pp['FirstAxisOfInertia'], pp['SecondAxisOfInertia'], pp['ThirdAxisOfInertia']]
    com = so.CenterOfMass
    pts = [v.Point for v in so.Vertexes]
    try:
        pts += so.tessellate(2.0)[0]
    except Exception:
        pass
    ext = []
    for a in axes:
        proj = [(p - com).dot(a) for p in pts]
        ext.append((min(proj), max(proj)))
    vobb = 1.0
    for lo, hi in ext:
        vobb *= max(hi - lo, 1e-9)
    rows.append(dict(path=label, label=label, solid=k, vol=so.Volume, com=[com.x, com.y, com.z],
                     axes=[[a.x, a.y, a.z] for a in axes], ext=ext, obb_fill=so.Volume / vobb, nfaces=len(so.Faces)))
json.dump(rows, open(OUT, 'w'))
print('solids', len(rows), 'unnamed', sum(r['label'] == '?' for r in rows))
