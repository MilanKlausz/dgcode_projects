"""FreeCAD: split a solid made of straight beams (a welded frame, e.g. the carriage ring SI-7615-178_01) into beams, each
described as an extrusion (cross-section at mid-length), for gen_components.py ('extra_solids').

1. the long planar faces of the solid (length >= MINLEN, length / width >= 4) are the walls of the beams;
2. faces with parallel long directions whose axis lines are closer than MAXSEP and whose extents overlap are one beam;
3. each beam: the long direction d, the extent along d (union of its faces), the cross-section box (its faces projected
   on the plane perpendicular to d, plus MARGIN); the part of the solid inside that box is the beam; its section at
   mid-length gives the polygons (even-odd), its volume the check.
env: STEP, OBB_JSON, LABEL (prefix of the solid label), OUT (json)"""
import os, json, math
import Part
from FreeCAD import Vector

MINLEN, MAXSEP, MARGIN = 200.0, 160.0, 3.0
rows = json.load(open(os.environ['OBB_JSON'])); label = os.environ['LABEL']
sh = Part.read(os.environ['STEP'])
out = {}
for r in rows:
    if not r['label'].startswith(label):
        continue
    so = sh.Solids[r['solid']]
    walls = []
    for f in so.Faces:
        if f.Surface.TypeId != 'Part::GeomPlane' or f.Area < 1000:
            continue
        n = f.normalAt(0, 0)
        pts = [v.Point for v in f.Vertexes]
        if len(pts) < 3:
            continue
        # principal in-plane direction of the face from its vertices
        c = sum(pts, Vector()) * (1.0 / len(pts))
        best = None
        for p in pts:
            for q in pts:
                d = q - p
                if d.Length > 1e-6 and (best is None or d.Length > best.Length):
                    best = d
        d = best.normalize()
        d = d - n * d.dot(n); d.normalize()
        e = n.cross(d)
        s = [(p - c).dot(d) for p in pts]; w = [(p - c).dot(e) for p in pts]
        length, width = max(s) - min(s), max(w) - min(w)
        if length >= MINLEN and length >= 4 * max(width, 1e-6):
            walls.append(dict(d=d, c=c, lo=min(s), hi=max(s), pts=pts))
    # cluster
    beams = []
    for wl in walls:
        for b in beams:
            if abs(b['d'].dot(wl['d'])) < 0.999:
                continue
            dd = b['d']; delta = wl['c'] - b['c']; perp = delta - dd * delta.dot(dd)
            s0 = delta.dot(dd)
            if perp.Length < MAXSEP and (s0 + wl['lo'] * math.copysign(1, dd.dot(wl['d'])) <= b['hi'] + 50 or True):
                b['walls'].append(wl); break
        else:
            beams.append(dict(d=wl['d'], c=wl['c'], hi=wl['hi'], walls=[wl]))
    res = []
    for k, b in enumerate(beams):
        d = b['d']
        a = d.cross(Vector(0, 0, 1)) if abs(d.z) < 0.9 else d.cross(Vector(1, 0, 0)); a.normalize()
        bb = d.cross(a); bb.normalize()
        P = [p for wl in b['walls'] for p in wl['pts']]
        o = b['c']
        s = [(p - o).dot(d) for p in P]; u = [(p - o).dot(a) for p in P]; v = [(p - o).dot(bb) for p in P]
        lo = Vector(min(s), min(u) - MARGIN, min(v) - MARGIN); hi = Vector(max(s), max(u) + MARGIN, max(v) + MARGIN)
        L = hi.x - lo.x
        if L < MINLEN:
            continue
        box = Part.makeBox(L, hi.y - lo.y, hi.z - lo.z)
        # box frame: x along d, y along a, z along bb, origin at lo
        from FreeCAD import Matrix
        m = Matrix(d.x, a.x, bb.x, 0, d.y, a.y, bb.y, 0, d.z, a.z, bb.z, 0, 0, 0, 0, 1)
        box.transformShape(m)
        box.translate(o + d * lo.x + a * lo.y + bb * lo.z)
        piece = so.common(box)
        if piece.Volume < 1000:
            continue
        mid = o + d * (0.5 * (lo.x + hi.x))
        wires = piece.slice(d, d.dot(mid))
        areas = sorted((abs(Part.Face(w).Area) for w in wires if w.isClosed()), reverse=True)
        area = (areas[0] - sum(areas[1:])) if areas else 0.0
        centre = mid + a * (0.5 * (lo.y + hi.y)) + bb * (0.5 * (lo.z + hi.z))
        polys = [[[(q - centre).dot(a), (q - centre).dot(bb)] for q in w.discretize(Distance=1.0)] for w in wires]
        res.append(dict(label=f"{r['label']}/beam{k}", long=[d.x, d.y, d.z], a=[a.x, a.y, a.z], b=[bb.x, bb.y, bb.z],
                        centre=[centre.x, centre.y, centre.z], length=L, vol=piece.Volume,
                        sections=[dict(area=area, polys=polys)] * 3))
    out[r['solid']] = dict(label=r['label'], volume=so.Volume, beams=res)
    print(r['label'], 'volume', round(so.Volume / 1e3), 'cm3:', len(walls), 'wall faces,', len(res), 'beams, beam volume',
          round(sum(x['vol'] for x in res) / 1e3), 'cm3; section x length / beam volume:',
          [round(x['sections'][0]['area'] * x['length'] / max(x['vol'], 1), 2) for x in res])
json.dump(out, open(os.environ['OUT'], 'w'))
