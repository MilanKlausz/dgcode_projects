"""Generate a CAD component set for the Geant4 LOKI geometry: simplified boxes in the frame of their bank.

For each selected CAD solid (inventory3 + cadlib material):
  1. its section: a plate (thin side < 0.25 x middle side) is rasterised in its own plane (the projection of the solid,
     RASTER mm cells, full thickness); a compact solid with OBB fill >= MINFILL is its full OBB section (thickness =
     its smallest side); other solids are skipped (reported);
  2. every cell (a small box: cell x cell x thickness, in the solid's frame, transformed to the Geant4 world with the
     alignment of its CAD file) is classified with a separating-axis test as
       - inside bank b: inside the box of bank b and clear of the bank's content (packs, masks, plates, and the boxes of
         this set already accepted in bank b),
       - outside: outside every bank box and clear of the world volumes (triangular masks, strips, tank, ... and the
         boxes of this set already accepted in the world),
       - conflict: anything else (the material is lost; reported);
  3. each class is covered by disjoint rectangles of cells (largest first, at most MAXBOX per class): the boxes.
     Boxes inside bank b are stored in the frame of bank b with mother "bank"; boxes outside are stored in the frame of
     the nearest bank with mother "world" (they follow the bank placement: calibration, rear detector distance).
Boxes never cross a bank box face or an existing volume, and the set does not overlap itself (accepted boxes become
obstacles). Every box is written SHRINK smaller per side.

The output data file (one box per line) is read by the geometry module (parameter cad_components).
usage: python gen_components.py CONFIG.json"""
import sys
import json
import collections
import numpy as np
from cadlib import load
from obb import Box, overlap_depth, separating_axis

RASTER = 2.0           # mm, the finest cell size (larger for large plates: at most MAXCELLS cells)
MAXCELLS = 120000
COVER = 0.99           # stop the rectangle cover of a class at this fraction of its cells
MAXBOX = 12            # at most this many boxes per class (bank / world) and solid
MINAREA = 0.003        # rectangles smaller than this fraction of the section are not used
MINFILL = 0.80
GAP = 0.05             # mm, clearance kept from other volumes
SHIFT_BANK = 2.0       # mm, a solid entering a bank box by at most this (through its thickness) is shifted out first
SHIFT_CONTENT = 3.5    # mm, the same for the content of a bank (masks: CAD vs Geant4 positions agree within about 3 mm)
CAL_MARGIN = 2.0       # mm, clearance of world boxes from the boxes of other banks (and of components attached to other
                       # banks) for calibrations not listed in the configuration ('calibrations': exact checks)
SHRINK = 0.01          # mm, every box is written this much smaller per side (touching boxes must not overlap after rounding)


def sat_depth_many(C, A, h, B):
    """separating-axis penetration depth of many boxes (centres C (n,3), common axes A (rows), half sizes h) with box B;
    < 0: separated by at least -depth"""
    axes = list(A) + list(B.A)
    for i in range(3):
        for j in range(3):
            c = np.cross(A[i], B.A[j])
            if np.linalg.norm(c) > 1e-9:
                axes.append(c / np.linalg.norm(c))
    depth = np.full(len(C), np.inf)
    D = B.c - C
    for L in axes:
        ra = np.sum(h * np.abs(A @ L)); rb = np.sum(B.h * np.abs(B.A @ L))
        depth = np.minimum(depth, ra + rb - np.abs(D @ L))
    return depth


def contained_many(C, A, h, B, margin):
    """True for the boxes (centres C, axes A, half sizes h) that are inside B (with a margin)"""
    off = (np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float) * h) @ A
    ok = np.ones(len(C), bool)
    for o in off:
        loc = (C + o - B.c) @ B.A.T
        ok &= np.all(np.abs(loc) <= B.h - margin, axis=1)
    return ok


def section_grid(r, tess, raster):
    """the solid's in-plane axes a, b and thin axis n, and its section raster"""
    order = np.argsort(r['size'])[::-1]
    A = r['axes']; a, b, n = A[order[0]], A[order[1]], A[order[2]]
    s = np.sort(r['size'])
    plate = s[0] < 0.25 * s[1]
    V, F = tess
    if plate and r['obb_fill'] < 0.97 and V is not None:
        loc = np.stack([(V - r['com']) @ a, (V - r['com']) @ b, (V - r['com']) @ n], 1)
        lo, hi = loc.min(0), loc.max(0)
    else:                                              # a box: its OBB
        e = np.array(r['ext'])[order]; lo, hi = e[:, 0], e[:, 1]; loc = None
    cell = max(raster, np.sqrt((hi[0] - lo[0]) * (hi[1] - lo[1]) / MAXCELLS))
    nx = max(1, int(np.ceil((hi[0] - lo[0]) / cell))); ny = max(1, int(np.ceil((hi[1] - lo[1]) / cell)))
    cx = (hi[0] - lo[0]) / nx; cy = (hi[1] - lo[1]) / ny
    gx = lo[0] + (np.arange(nx) + 0.5) * cx; gy = lo[1] + (np.arange(ny) + 0.5) * cy
    if loc is None:
        grid = np.ones((nx, ny), bool)
    else:
        grid = np.zeros((nx, ny), bool)
        for tri in loc[F][:, :, :2]:
            x0, y0 = tri.min(0); x1, y1 = tri.max(0)
            i0, i1 = np.searchsorted(gx, [x0, x1]); j0, j1 = np.searchsorted(gy, [y0, y1])
            if i1 <= i0 or j1 <= j0:
                continue
            X, Y = np.meshgrid(gx[i0:i1], gy[j0:j1], indexing='ij')
            (xa, ya), (xb, yb), (xc, yc) = tri
            d = (yb - yc) * (xa - xc) + (xc - xb) * (ya - yc)
            if abs(d) < 1e-12:
                continue
            l1 = ((yb - yc) * (X - xc) + (xc - xb) * (Y - yc)) / d
            l2 = ((yc - ya) * (X - xc) + (xa - xc) * (Y - yc)) / d
            grid[i0:i1, j0:j1] |= (l1 >= -1e-9) & (l2 >= -1e-9) & (1 - l1 - l2 >= -1e-9)
    return a, b, n, (lo, cx, cy), grid, hi[2] - lo[2], (lo[2] + hi[2]) / 2


def rectangles(mask, total):
    """disjoint maximal rectangles covering the True cells of mask (largest first)"""
    left = mask.copy(); out = []
    nx, ny = mask.shape
    want = left.sum()
    while left.sum() > (1 - COVER) * want and len(out) < MAXBOX:
        best = (0, None); h = np.zeros(ny, int)
        for i in range(nx):
            h = np.where(left[i], h + 1, 0)
            st = []
            for j in range(ny + 1):
                hj = h[j] if j < ny else 0
                s = j
                while st and st[-1][1] >= hj:
                    s, hs = st.pop()
                    if hs * (j - s) > best[0]:
                        best = (hs * (j - s), (i - hs + 1, i + 1, s, j))
                st.append((s, hj))
        if best[1] is None or best[0] < MINAREA * total:
            break
        i0, i1, j0, j1 = best[1]; left[i0:i1, j0:j1] = False; out.append(best[1])
    return out


def main(cfg):
    C = json.load(open(cfg))
    G = json.load(open(C.get('g4dump', 'g4_nominal.json')))['volumes']

    def g4box(v):
        R = np.array(v['R']); h = np.array(v['pars']) * (1 if v['solid'] == 'G4Box' else 0.5)
        return Box(v['t'], R.T, h, name=f"{v['pv']}:{v['copy']}")
    bankpv = {v['copy']: v for v in G if v['lv'] == 'Bank'}
    banks = {k: g4box(v) for k, v in bankpv.items()}
    content = {k: [] for k in banks}
    cur = None
    for v in G:
        if v['depth'] == 1:
            cur = v['copy'] if v['lv'] == 'Bank' else None
        elif v['depth'] == 2 and cur is not None:
            content[cur].append(g4box(v))
    world_other = [g4box(v) for v in G if v['depth'] == 1 and v['lv'] != 'Bank']
    # other bank placements (calibrations) in which a world box, moved with its own bank, must also be clear
    cals = []
    for fn in C.get('calibrations', []):
        Gc = json.load(open(fn))['volumes']
        bpv = {v['copy']: v for v in Gc if v['lv'] == 'Bank'}
        move = {}                                      # nominal world -> calibrated world, per bank
        for k in bankpv:
            Rn, tn = np.array(bankpv[k]['R']), np.array(bankpv[k]['t']); Rc, tc = np.array(bpv[k]['R']), np.array(bpv[k]['t'])
            Rm = Rc @ Rn.T; move[k] = (Rm, tc - Rm @ tn)
        cals.append(dict(name=fn, banks={k: g4box(v) for k, v in bpv.items()}, move=move,
                         world=[g4box(v) for v in Gc if v['depth'] == 1 and v['lv'] != 'Bank']))
    out = []; report = []; lost = collections.defaultdict(float); kept = collections.defaultdict(float)
    avoid = []                                         # boxes of other sets that may be used together with this one
    for fn in C.get('avoid', []):
        for line in open(fn):
            if line.startswith('#'):
                continue
            f = line.split(); bk = int(f[2]); v = np.array(f[4:], float)
            Rb = np.array(bankpv[bk]['R']); tb = np.array(bankpv[bk]['t'])
            loc = Box(v[:3], v[3:12].reshape(3, 3), v[12:15] + SHRINK, f[0], f[1])
            avoid.append(dict(world=loc.transformed(Rb, tb), bank=bk, mother=f[3], name=f[0]))
    perpart = collections.defaultdict(lambda: [0.0, 0.0])
    for sel in C['selections']:
        al = np.load(sel['align']); R, t = al['R'], al['t']
        tess = np.load(f"tess_{sel['tag']}.npz"); reg = sel.get('region')
        for r in load(sel['tag']):
            if r['mat'] not in sel['materials']:
                continue
            if sel.get('labels') and not any(r['label'].startswith(l) for l in sel['labels']):
                continue
            if sel.get('exclude') and any(r['label'].startswith(l) for l in sel['exclude']):
                continue
            if reg and not (reg[0] < r['centre'][0] < reg[1] and reg[2] < r['centre'][1] < reg[3]):
                continue
            if r['vol'] < sel.get('min_volume_cm3', 0) * 1e3:
                continue
            s = np.sort(r['size'])
            if not (s[0] < 0.25 * s[1] or r['obb_fill'] >= MINFILL):
                report.append(f"SKIPPED {r['label']} ({r['mat']}): not a plate and OBB fill {r['obb_fill']:.2f} < {MINFILL}")
                continue
            # small penetrations of the whole solid (along its thin axis): shift it out first
            shift = np.zeros(3)
            wsol = Box(r['centre'], r['axes'], r['size'] / 2).transformed(R, t)
            thin_w = R @ r['axes'][int(np.argmin(r['size']))]
            for O, lim in [(B, SHIFT_BANK) for B in banks.values()] + [(c, SHIFT_CONTENT) for k in banks for c in content[k]]:
                ws = Box(wsol.c + shift, wsol.A, wsol.h)
                d = overlap_depth(ws, O)
                if 0 < d + GAP <= lim + GAP:
                    L, pen = separating_axis(O, ws)
                    if abs(L @ thin_w) > 0.99:                # only through the thickness
                        shift += L * (pen + 2 * GAP)
            # nearly contained in a bank box (protruding by at most SHIFT_BANK along the thin axis): shift it in
            for B in banks.values():
                ws = Box(wsol.c + shift, wsol.A, wsol.h)
                if overlap_depth(ws, B) > wsol.h.min() and not B.contains(ws, GAP):
                    for step in np.arange(0.05, SHIFT_BANK + 1e-9, 0.05):
                        hit = [sg for sg in (1, -1) if B.contains(Box(ws.c + sg * step * thin_w, ws.A, ws.h), GAP)]
                        if hit:
                            shift += hit[0] * step * thin_w
                            break
            if np.linalg.norm(shift) > 0:
                report.append(f"shifted {r['label']}#{r['solid']} by {np.linalg.norm(shift):.2f} mm")
            key = f"v{r['solid']}"
            a, b, n, (lo, cx, cy), grid, thick, zmid = section_grid(
                r, (tess[key], tess[f"f{r['solid']}"]) if key in tess.files else (None, None), sel.get('raster', RASTER))
            Aw = (R @ np.array([a, b, n]).T).T              # the cells' axes in the world
            ii, jj = np.nonzero(grid)
            local = r['com'] + np.outer(lo[0] + (ii + 0.5) * cx, a) + np.outer(lo[1] + (jj + 0.5) * cy, b) + zmid * n
            Cw = (R @ local.T).T + t + shift
            hc = np.array([cx / 2, cy / 2, thick / 2])
            cls = np.full(len(Cw), -2)                     # -2 conflict, -1 world, k bank k
            wsh = Box(wsol.c + shift, wsol.A, wsol.h)
            ovb = {k: overlap_depth(wsh, B) for k, B in banks.items()}
            own = max(ovb, key=ovb.get) if max(ovb.values()) > -GAP else min(banks, key=lambda q: np.linalg.norm(banks[q].c - wsh.c))
            out_all = np.ones(len(Cw), bool)
            for k, B in banks.items():
                out_all &= sat_depth_many(Cw, Aw, hc, B) < -(GAP if k == own else CAL_MARGIN)
                free = contained_many(Cw, Aw, hc, B, GAP)
                if free.any():
                    for O in content[k] + [o['world'] for o in out + avoid if o['mother'] == 'bank' and o['bank'] == k]:
                        idx = np.nonzero(free)[0]
                        if not len(idx):
                            break
                        free[idx] &= sat_depth_many(Cw[idx], Aw, hc, O) < -GAP
                    cls[free] = k
            idx = np.nonzero(out_all)[0]
            for O, m in [(o, GAP) for o in world_other] + [(o['world'], GAP if o['bank'] == own else CAL_MARGIN)
                                                          for o in out + avoid if o['mother'] == 'world']:
                if len(idx):
                    idx = idx[sat_depth_many(Cw[idx], Aw, hc, O) < -m]
            for cal in cals:                               # the world cells moved with their own bank
                if not len(idx):
                    break
                Rm, tm = cal['move'][own]
                Cc = (Rm @ Cw[idx].T).T + tm; Ac = (Rm @ Aw.T).T
                keep = np.ones(len(idx), bool)
                for k, B in cal['banks'].items():
                    keep &= sat_depth_many(Cc, Ac, hc, B) < -GAP
                for O in cal['world'] + [o['world'].transformed(*cal['move'][o['bank']]) for o in out + avoid if o['mother'] == 'world']:
                    keep &= sat_depth_many(Cc, Ac, hc, O) < -GAP
                idx = idx[keep]
            cls[idx] = -1
            vcell = r['vol'] / max(grid.sum(), 1)            # the solid's volume per section cell
            lost[r['mat']] += (cls == -2).sum() * vcell; kept[r['mat']] += (cls != -2).sum() * vcell
            perpart[r['label']][0] += r['vol']; perpart[r['label']][1] += (cls != -2).sum() * vcell
            if (cls == -2).mean() > 0.05:
                report.append(f"{r['label']}#{r['solid']}: {(cls == -2).mean():.0%} of the section in conflict (not used)")
            nbox = 0
            for c in sorted(set(cls.tolist()) - {-2}):
                mask = np.zeros_like(grid); mask[ii[cls == c], jj[cls == c]] = True
                bk = c if c >= 0 else own
                Rb = np.array(bankpv[bk]['R']); tb = np.array(bankpv[bk]['t'])
                for (i0, i1, j0, j1) in rectangles(mask, grid.sum()):
                    ca = lo[0] + (i0 + i1) / 2 * cx; cb = lo[1] + (j0 + j1) / 2 * cy
                    hloc = np.array([(i1 - i0) * cx / 2, (j1 - j0) * cy / 2, thick / 2])
                    cw = R @ (r['com'] + ca * a + cb * b + zmid * n) + t + shift
                    wbox = Box(cw, Aw, hloc, f"{r['label'].replace(' ', '_')}#{r['solid']}.{nbox}", r['mat'])
                    nbox += 1
                    loc = wbox.transformed(Rb.T, -Rb.T @ tb)
                    out.append(dict(name=wbox.name, label=r['label'], mat=r['mat'], bank=bk, mother='bank' if c >= 0 else 'world',
                                    centre=loc.c.tolist(), axes=loc.A.tolist(), half=loc.h.tolist(), world=wbox))
    with open(C['output'], 'w') as f:
        f.write(f"# LOKI CAD component set '{C['name']}': {C.get('description', '')}\n")
        f.write(f"# generated by gen_components.py from {', '.join(sorted({s['step'].split('/')[-1] for s in C['selections']}))}\n")
        f.write("# one box per line: name material bankId mother(bank|world) cx cy cz [mm, bank frame] "
                "ax ay az bx by bz nx ny nz [the box axes in the bank frame] hx hy hz [half lengths along them, mm]\n")
        for o in out:
            A = np.array(o['axes']); h = np.array(o['half']) - SHRINK
            f.write(f"{o['name']} {o['mat']} {o['bank']} {o['mother']} " + ' '.join(f'{x:.4f}' for x in o['centre']) + ' '
                    + ' '.join(f'{x:.8f}' for x in A.ravel()) + ' ' + ' '.join(f'{x:.4f}' for x in h) + '\n')
    print(f"{C['name']}: {len(out)} boxes written to {C['output']}")
    for line in report:
        print('  ' + line)
    worst = sorted(perpart.items(), key=lambda kv: -(kv[1][0] - kv[1][1]))[:8]
    print('  largest losses (part: CAD cm3 -> kept cm3):', ', '.join(f'{k}: {v[0]/1e3:.0f} -> {v[1]/1e3:.0f}' for k, v in worst))
    print('  by (material, mother, bank):', dict(collections.Counter((o['mat'], o['mother'], o['bank']) for o in out)))
    vol = collections.defaultdict(float)
    for o in out:
        vol[o['mat']] += 8 * np.prod(o['half']) / 1e3
    print('  box volume [cm3]:', {k: round(v) for k, v in vol.items()},
          '; CAD volume kept / in conflict [cm3]:', {k: (round(kept[k] / 1e3), round(lost[k] / 1e3)) for k in kept})


if __name__ == '__main__':
    main(sys.argv[1])
