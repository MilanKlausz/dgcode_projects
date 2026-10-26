"""Generate a CAD pack-component set: boxes in the frame of the Geant4 detector pack (G4Para EmptyPackBox; x = depth,
y = across the tubes, z = along the tubes), from the mid-length section of CAD pack parts (SI-7615-131, slice.py).

The CAD pack SI-7615-133 is mapped to the Geant4 pack frame by G4 (x, y, z) = (-y_CAD - 18.84, x_CAD, z_CAD + 658)
(from the tube centres; the B4C and Al panels of the Geant4 pack coincide with the CAD parts). The section of each part
is rasterised and covered by disjoint rectangles that stay inside the pack's G4Para and clear of the tubes and of the
B4C/Al panels of the Geant4 pack. The boxes extend along the tubes as 'length' boxes: z in [-(L/2 + a), L/2 + a] with
L the straw length of the bank (the CAD parts are L + 2a long for the 1200 mm straws of the CAD pack).

usage: python gen_pack.py CONFIG.json   (config: name, output, description, parts: [{section, label, material, a}])"""
import sys
import json
import numpy as np
from matplotlib.path import Path

M = np.array([[0, -1], [1, 0.]]); T2 = np.array([-18.84, 0])
DX, DY, ALPHA = 0.799 * 0.5 * 265.7, 0.5 * 55.2, np.radians(13.45)     # the G4Para of the pack
TUBE_R = 12.7
TUBES = [(-0.5 * 265.7 + 36.5 + (k % 4) * 27.0 + (28.4 * np.sin(np.radians(13.45)) if k < 4 else 0),
          0.5 * 28.4 * np.cos(np.radians(13.45)) * (1 if k < 4 else -1)) for k in range(8)]
PANELS = [(7.65, 12.65, -24.0, 27.2), (-3.15, 7.65, -24.0, -19.0), (-3.15, -0.15, -27.2, -24.0),   # B4C (x0, x1, y0, y1)
          (12.65, 15.85, -27.2, 27.2), (-0.15, 12.65, -27.2, -24.0)]                                # Al
CELL = 0.25
GAP = 0.05
MAXBOX = 10


def free_cells(X, Y, beyond_tubes=False):
    """cell centres that are inside the Para and clear of tubes and panels (with GAP + half a cell); beyond the tube ends
    (end slabs, z > L/2 + 12 mm) only the Para limits"""
    m = CELL / 2 + GAP
    ok = (np.abs(Y) <= DY - m) & (np.abs(X - Y * np.tan(ALPHA)) <= DX - m / np.cos(ALPHA))
    if beyond_tubes:
        return ok
    for (tx, ty) in TUBES:
        ok &= np.hypot(X - tx, Y - ty) > TUBE_R + m * 1.5
    for (x0, x1, y0, y1) in PANELS:
        ok &= ~((X > x0 - m) & (X < x1 + m) & (Y > y0 - m) & (Y < y1 + m))
    return ok


def rectangles(mask):
    left = mask.copy(); out = []; nx, ny = mask.shape; want = left.sum()
    while left.sum() > 0.01 * want and len(out) < MAXBOX:
        best = (0, None); h = np.zeros(ny, int)
        for i in range(nx):
            h = np.where(left[i], h + 1, 0); st = []
            for j in range(ny + 1):
                hj = h[j] if j < ny else 0; s = j
                while st and st[-1][1] >= hj:
                    s, hs = st.pop()
                    if hs * (j - s) > best[0]:
                        best = (hs * (j - s), (i - hs + 1, i + 1, s, j))
                st.append((s, hj))
        if best[1] is None or best[0] < 0.005 * want:
            break
        i0, i1, j0, j1 = best[1]; left[i0:i1, j0:j1] = False; out.append(best[1])
    return out, 1 - left.sum() / max(want, 1)


def main(cfg):
    C = json.load(open(cfg))
    lines = []
    for part in C['parts']:
        for slab in part.get('slabs', []):              # an end part in slabs along the tubes: 'ends' boxes per slab
            S = json.load(open(slab['section']))
            for idx, wires in S.items():
                polys = [(M @ np.array(w).T).T + T2 for w in wires if len(w) > 2]
                if not polys:
                    continue
                P = np.vstack(polys); lo, hi = P.min(0), P.max(0)
                gx = np.arange(lo[0], hi[0] + CELL, CELL) + CELL / 2; gy = np.arange(lo[1], hi[1] + CELL, CELL) + CELL / 2
                X, Y = np.meshgrid(gx, gy, indexing='ij'); pts = np.stack([X.ravel(), Y.ravel()], 1)
                inside = np.zeros(len(pts), bool)
                for poly in polys:
                    inside ^= Path(poly).contains_points(pts)
                inside = inside.reshape(X.shape)
                mask = inside & free_cells(X, Y, beyond_tubes=True)
                rects, cov = rectangles(mask)
                area = inside.sum() * CELL ** 2; used = sum((i1 - i0) * (j1 - j0) for i0, i1, j0, j1 in rects) * CELL ** 2
                print(f"{part['label']} slab {slab['a']:.1f}-{slab['b']:.1f}: section {area:.0f} mm2, {len(rects)} boxes cover "
                      f"{used:.0f} mm2 ({used / max(area, 1):.1%}; in the pack volume: {mask.sum() / max(inside.sum(), 1):.1%})")
                for k, (i0, i1, j0, j1) in enumerate(rects):
                    x0, x1 = gx[i0] - CELL / 2, gx[i1 - 1] + CELL / 2; y0, y1 = gy[j0] - CELL / 2, gy[j1 - 1] + CELL / 2
                    lines.append(f"{part['label']}.z{slab['a']:.0f}.{k} {part['material']} {x0:.3f} {x1:.3f} {y0:.3f} {y1:.3f} ends "
                                 f"{slab['a']:.3f} {slab['b']:.3f}")
        if 'section' not in part:
            continue
        S = json.load(open(part['section']))
        for idx, wires in S.items():
            if part.get('index') is not None and int(idx) not in part['index']:
                continue
            polys = [(M @ np.array(w).T).T + T2 for w in wires]
            P = np.vstack(polys); lo, hi = P.min(0), P.max(0)
            gx = np.arange(lo[0], hi[0] + CELL, CELL) + CELL / 2; gy = np.arange(lo[1], hi[1] + CELL, CELL) + CELL / 2
            X, Y = np.meshgrid(gx, gy, indexing='ij'); pts = np.stack([X.ravel(), Y.ravel()], 1)
            inside = np.zeros(len(pts), bool)
            for poly in polys:                          # (even-odd: holes in a profile are separate wires)
                inside ^= Path(poly).contains_points(pts)
            inside = inside.reshape(X.shape)
            mask = inside & free_cells(X, Y)
            rects, cov = rectangles(mask)
            area = inside.sum() * CELL ** 2; used = sum((i1 - i0) * (j1 - j0) for i0, i1, j0, j1 in rects) * CELL ** 2
            print(f"{part['label']} (solid {idx}): section {area:.0f} mm2, {len(rects)} boxes cover {used:.0f} mm2 "
                  f"({used / max(area, 1):.1%}; {mask.sum() / max(inside.sum(), 1):.1%} of the section is free)")
            for k, (i0, i1, j0, j1) in enumerate(rects):
                x0, x1 = gx[i0] - CELL / 2, gx[i1 - 1] + CELL / 2; y0, y1 = gy[j0] - CELL / 2, gy[j1 - 1] + CELL / 2
                lines.append(f"{part['label']}.{k} {part['material']} {x0:.3f} {x1:.3f} {y0:.3f} {y1:.3f} length {part['a']:.3f} 0")
    # merge 'ends' boxes of consecutive slabs with the same rectangle (one box spanning the slabs)
    byrect = {}
    out_lines = []
    for ln in lines:
        f = ln.split()
        if f[6] != 'ends':
            out_lines.append(ln); continue
        key = (f[0].split('.z')[0], f[1], f[2], f[3], f[4], f[5])
        a, b = float(f[7]), float(f[8])
        if key in byrect and abs(byrect[key][1] - a) < 0.2:          # contiguous: extend
            byrect[key][1] = b
        else:
            if key in byrect:
                out_lines.append(('ENDS', key, byrect[key][0], byrect[key][1]))
            byrect[key] = [a, b]
    for key, (a, b) in byrect.items():
        out_lines.append(('ENDS', key, a, b))
    lines = []
    for k, ln in enumerate(out_lines):
        if isinstance(ln, tuple):
            _, key, a, b = ln
            lines.append(f"{key[0]}.e{k} {key[1]} {key[2]} {key[3]} {key[4]} {key[5]} ends {a:.3f} {b:.3f}")
        else:
            lines.append(ln)
    with open(C['output'], 'w') as f:
        f.write(f"# LOKI CAD pack component set '{C['name']}': {C.get('description', '')}\n")
        f.write("# generated by gen_pack.py; boxes in the frame of the detector pack (EmptyPackBox: x = depth, y = across the tubes)\n")
        f.write("# one box per line: name material x0 x1 y0 y1 [mm] zkind a b; zkind length: z in [-(L/2+a), L/2+a] (L: straw length),\n")
        f.write("#   ends: z in [L/2+a, L/2+b] and its mirror image\n")
        f.write('\n'.join(lines) + '\n')
    print(f"{C['name']}: {len(lines)} boxes written to {C['output']}")


if __name__ == '__main__':
    main(sys.argv[1])
