"""Generate SIMPLIFIED CAD component sets: few homogenised boxes per CAD solid, optionally grouped in envelopes.

Compared with gen_components.py (which covers the exact section of each solid with up to 12 boxes per region) every
selected solid becomes at most MAXBOX boxes per region (inside bank b / in the world), and a box may include the
small enclosed holes of the section (bolt holes, the windows of a truss, the inside of a hollow beam: holes up to
hole_max_cm2 of the selection, default HOLE_MAX_CM2 for aluminium, HOLE_MAX_OTHER_CM2 for other materials): its
material is the solid's material with the density scaled by the filled fraction ("fill"), so that the mass of the
solid is kept. Large openings (the windows of U- and Pi-shaped frame plates in front of the tubes) stay open. Steps:
  1. its section (gen_components.section_grid / extrusion_grid) on a coarse raster (cells of at least CELL mm, at most
     MAXCELLS cells), with the material fraction and the small-hole fraction of every cell; a solid that is neither a
     plate nor an extrusion is treated as a plate along its thinnest axis (its projection, the volume kept by the fill);
  2. every cell (material or not) is classified as inside bank b and clear of its content, outside all banks and clear
     of the world volumes, or conflict (as in gen_components.py, incl. the other calibrations and the tank radius);
  3. per class, rectangles of free cells are chosen greedily: the one with the largest sum of (material fraction -
     HOLE_COST x small-hole fraction - OPEN_COST x open fraction) first, until COVER of the class material is covered
     or MAXBOX rectangles;
  4. each rectangle is a box with fill = (the solid's volume in it) / (box volume).

The configuration defines several sets ("sets": [{name, output, description, selections}, ...]) generated together:
solids are taken set by set, in the order of the selections and, within one, by decreasing relevance (relevance.py:
volume x the fraction of directions in which a tube is visible) or volume; accepted boxes (of all sets) are obstacles
for later ones, so that any combination of the sets is free of overlaps. Selection options: tag, step, align, materials,
labels, exclude, region, beams, min_volume_cm3, min_relevance_cm3 (needs the 'relevance' files of the configuration),
hole_max_cm2.

Envelopes ("envelopes": true): the world boxes of each set and bank are grouped into boxes of the world material
(aligned with the bank frame) that are clear of all bank boxes, world volumes and the boxes of all sets (also in the
other calibrations), of the tank radius and of each other: greedy agglomeration (the pair whose common envelope adds the
least volume first).

Output format (G4GeoLoki/CadComponents.hh):
  box <name> <material> <fill> <bankId> <mother> cx cy cz ax ay az bx by bz nx ny nz hx hy hz
  envelope <name> <bankId> cx cy cz ax ay az bx by bz nx ny nz hx hy hz
(bank frame of bankId; mother: bank, world or the name of an envelope)
usage: python gen_simplified.py CONFIG.json"""
import sys
import json
import collections
import numpy as np
import gen_components as GC
from cadlib import load
from obb import Box, overlap_depth, separating_axis

CELL = 8.0            # mm, the smallest cell of the coarse raster
MAXCELLS = 6000       # cells per section at most
FINE = 40000          # cells of the fine raster (the material fraction of the coarse cells)
MAXBOX = 4            # boxes per solid and class at most
COVER = 0.97          # stop when this fraction of the material of a class is covered
HOLE_MAX_CM2 = 400.0  # aluminium: enclosed holes up to this area are homogenised (filled with the solid's material)
HOLE_MAX_OTHER_CM2 = 10.0   # other materials (hydrogenous Tufnol, B4C): only small holes
HOLE_COST = 0.15      # cost of a (small) hole cell in a box, relative to the gain of a full cell
OPEN_COST = 1.0       # cost of an open cell (outside the section or in a large opening)
MINPART = 0.02        # rectangles with less than this fraction of the solid's material are not used
GAP, SHRINK = GC.GAP, GC.SHRINK
CORNERS = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)], float)


def coarse_grid(fine, nx, ny):
    """material fraction of nx x ny coarse cells from a fine boolean grid"""
    fx, fy = fine.shape
    ix = (np.arange(fx) * nx) // fx; iy = (np.arange(fy) * ny) // fy
    I = np.repeat(ix[:, None], fy, 1); J = np.repeat(iy[None, :], fx, 0)
    m = np.zeros((nx, ny)); c = np.zeros((nx, ny))
    np.add.at(m, (I, J), fine); np.add.at(c, (I, J), 1)
    return m / np.maximum(c, 1)


def small_holes(fine, cell_area_mm2, hole_max_cm2):
    """the enclosed holes of a section (boolean grid) with an area up to hole_max_cm2"""
    from scipy import ndimage
    holes = ndimage.binary_fill_holes(fine) & ~fine
    lab, n = ndimage.label(holes)
    if not n:
        return holes
    area = ndimage.sum(holes, lab, range(1, n + 1)) * cell_area_mm2 / 100
    return np.isin(lab, 1 + np.nonzero(area <= hole_max_cm2)[0])


def best_rectangle(free, w):
    """the rectangle of free cells with the largest sum of the cell weights w: (score, (i0, i1, j0, j1))"""
    nx, ny = free.shape
    w = np.where(free, w, -1e9)
    best = (0.0, None)
    for i0 in range(nx):
        colfree = np.ones(ny, bool); colw = np.zeros(ny)
        for i1 in range(i0 + 1, nx + 1):
            colfree &= free[i1 - 1]; colw += w[i1 - 1]
            if not colfree.any():
                break
            run = 0.0; start = 0                 # Kadane over each run of free columns
            for j in range(ny):
                if not colfree[j]:
                    run = 0.0; start = j + 1
                    continue
                if run <= 0:
                    run = 0.0; start = j
                run += colw[j]
                if run > best[0]:
                    best = (run, (i0, i1, start, j + 1))
    return best


class Scene:
    """the Geant4 volumes (nominal placement and other calibrations) and the accepted boxes of all sets"""

    def __init__(self, C):
        G = json.load(open(C.get('g4dump', 'g4_nominal.json')))['volumes']
        self.bankpv = {v['copy']: v for v in G if v['lv'] == 'Bank'}
        self.banks = {k: self.g4box(v) for k, v in self.bankpv.items()}
        self.content = {k: [] for k in self.banks}
        cur = None
        for v in G:
            if v['depth'] == 1:
                cur = v['copy'] if v['lv'] == 'Bank' else None
            elif v['depth'] == 2 and cur is not None:
                self.content[cur].append(self.g4box(v))
        self.world_other = [self.g4box(v) for v in G if v['depth'] == 1 and v['lv'] != 'Bank']
        self.cals = []
        for fn in C.get('calibrations', []):
            Gc = json.load(open(fn))['volumes']
            bpv = {v['copy']: v for v in Gc if v['lv'] == 'Bank'}
            move = {}                              # nominal world -> calibrated world, per bank
            for k, v in self.bankpv.items():
                Rn, tn = np.array(v['R']), np.array(v['t']); Rc, tc = np.array(bpv[k]['R']), np.array(bpv[k]['t'])
                Rm = Rc @ Rn.T; move[k] = (Rm, tc - Rm @ tn)
            self.cals.append(dict(name=fn, banks={k: self.g4box(v) for k, v in bpv.items()}, move=move,
                                  world=[self.g4box(v) for v in Gc if v['depth'] == 1 and v['lv'] != 'Bank']))
        self.rmax = C.get('max_radius_mm')
        self.out = []                              # the accepted boxes of all sets

    @staticmethod
    def g4box(v):
        R = np.array(v['R']); h = np.array(v['pars']) * (1 if v['solid'] == 'G4Box' else 0.5)
        return Box(v['t'], R.T, h, name=f"{v['pv']}:{v['copy']}")

    def bank_frame(self, bk):
        return np.array(self.bankpv[bk]['R']), np.array(self.bankpv[bk]['t'])

    def bank_free(self, Cw, Aw, hc, k):
        """cells (centres Cw, common axes Aw, half sizes hc) inside bank k and clear of its content and boxes"""
        free = GC.contained_many(Cw, Aw, hc, self.banks[k], GAP)
        for O in self.content[k] + [o['world'] for o in self.out if o['mother'] == 'bank' and o['bank'] == k]:
            idx = np.nonzero(free)[0]
            if not len(idx):
                break
            free[idx] &= GC.sat_depth_many(Cw[idx], Aw, hc, O) < -GAP
        return free

    def world_free(self, Cw, Aw, hc, own):
        """indices of the cells clear of all banks, world volumes and world boxes (also in the other calibrations, where
        the cells move with bank `own`), inside the radius limit"""
        idx = np.arange(len(Cw))
        worldboxes = [o for o in self.out if o['mother'] != 'bank']
        for k, B in self.banks.items():
            idx = idx[GC.sat_depth_many(Cw[idx], Aw, hc, B) < -(GAP if k == own else GC.CAL_MARGIN)]
        for O, m in [(o, GAP) for o in self.world_other] + [(o['world'], GAP if o['bank'] == own else GC.CAL_MARGIN) for o in worldboxes]:
            if len(idx):
                idx = idx[GC.sat_depth_many(Cw[idx], Aw, hc, O) < -m]
        for cal in self.cals:
            if not len(idx):
                break
            Rm, tm = cal['move'][own]
            Cc = (Rm @ Cw[idx].T).T + tm; Ac = (Rm @ Aw.T).T
            keep = np.ones(len(idx), bool)
            for B in cal['banks'].values():
                keep &= GC.sat_depth_many(Cc, Ac, hc, B) < -GAP
            for O in cal['world'] + [o['world'].transformed(*cal['move'][o['bank']]) for o in worldboxes]:
                keep &= GC.sat_depth_many(Cc, Ac, hc, O) < -GAP
            idx = idx[keep]
        if self.rmax and len(idx):
            off = (CORNERS * hc) @ Aw
            r = np.max([np.hypot(*(Cw[idx] + o)[:, :2].T) for o in off], axis=0)
            idx = idx[r < self.rmax - GAP]
        return idx


def select_solids(sel, relevance, report):
    """the CAD solids of a selection (incl. the beams of welded frames), most relevant first"""
    rows = load(sel['tag'])
    try:
        xsec = json.load(open(f"xsec_{sel['tag']}.json"))
    except FileNotFoundError:
        xsec = {}
    if sel.get('beams'):                           # welded frames (beam_split.py): only their straight beams
        mats = {str(r['solid']): r['mat'] for r in rows}
        for sid, ent in json.load(open(sel['beams'])).items():
            parent = next(r for r in rows if str(r['solid']) == sid)
            for k, bm in enumerate(ent['beams']):
                pid = f'{sid}b{k}'; xsec[pid] = bm; L = bm['length']
                rows.append(dict(label=bm['label'], solid=pid, mat=mats[sid], vol=bm['vol'], obb_fill=0.3,
                                 size=np.array([L, 0.3 * L, 0.3 * L]), axes=np.array([bm['long'], bm['a'], bm['b']]),
                                 centre=np.array(bm['centre']), com=np.array(bm['centre']),
                                 ext=np.array([[-L / 2, L / 2]] * 3), colours=[], base=bm['label'],
                                 rel_key=(parent['label'], parent['solid']), rel_share=bm['vol'] / parent['vol']))
            rows = [r for r in rows if str(r['solid']) != sid]
    reg = sel.get('region')
    chosen = []
    for r in rows:
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
        rel = relevance.get(r.get('rel_key', (r['label'], r['solid'])))
        r['relevance'] = rel * r.get('rel_share', 1.0) if rel is not None else None
        if sel.get('min_relevance_cm3') is not None:
            if r['relevance'] is None:
                report.append(f"no relevance for {r['label']}#{r['solid']}: not used")
                continue
            if r['relevance'] < sel['min_relevance_cm3']:
                continue
        chosen.append(r)
    chosen.sort(key=lambda r: -(r['relevance'] if r['relevance'] is not None else r['vol'] / 1e3))
    return chosen, xsec


def add_solid(scene, setname, r, R, t, tess, xsec, perpart, report, hole_max_cm2=None):
    """the boxes of one CAD solid (CAD frame -> Geant4 world: R, t)"""
    banks, content, out = scene.banks, scene.content, scene.out
    s = np.sort(r['size'])
    xs = xsec.get(str(r['solid']))
    # the beams of a welded frame always by their (middle) section; other solids if they are extrusions
    extr = xs is not None and ('rel_key' in r or GC.is_extrusion(xs) and not (s[0] < 0.25 * s[1] or r['obb_fill'] >= GC.MINFILL))
    if not extr and f"v{r['solid']}" not in tess.files and r['obb_fill'] < 0.97:
        report.append(f"{r['label']}#{r['solid']}: no section (no tessellation): not used")
        return
    # small penetrations of the whole solid along its thin axis: shift it out first (as gen_components.py)
    shift = np.zeros(3)
    wsol = Box(r['centre'], r['axes'], r['size'] / 2).transformed(R, t)
    thin_w = R @ r['axes'][int(np.argmin(r['size']))]
    for O, lim in ([(B, GC.SHIFT_BANK) for B in banks.values()] + [(c, GC.SHIFT_CONTENT) for k in banks for c in content[k]]
                   + [(o['world'], GC.SHIFT_CONTENT) for o in out] + [(o, GC.SHIFT_CONTENT) for o in scene.world_other]):
        ws = Box(wsol.c + shift, wsol.A, wsol.h)
        d = overlap_depth(ws, O)
        if 0 < d + GAP <= lim + GAP:
            L, pen = separating_axis(O, ws)
            if abs(L @ thin_w) > 0.99:
                shift += L * (pen + 2 * GAP)
    for B in banks.values():                      # nearly inside a bank box: shift it in
        ws = Box(wsol.c + shift, wsol.A, wsol.h)
        if overlap_depth(ws, B) > wsol.h.min() and not B.contains(ws, GAP):
            for step in np.arange(0.05, GC.SHIFT_BANK + 1e-9, 0.05):
                hit = [sg for sg in (1, -1) if B.contains(Box(ws.c + sg * step * thin_w, ws.A, ws.h), GAP)]
                if hit:
                    shift += hit[0] * step * thin_w
                    break
    # the section: fine raster (material), coarse raster (cells)
    GC.MAXCELLS = FINE
    key = f"v{r['solid']}"
    if extr:
        a, b, n, base_pt, (lo, cx, cy), fine, thick = GC.extrusion_grid(xs, 1.0)
        zmid = 0.0
    else:
        a, b, n, (lo, cx, cy), fine, thick, zmid = GC.section_grid(
            r, (tess[key], tess[f"f{r['solid']}"]) if key in tess.files else (None, None), 1.0)
        base_pt = r['com']
    Lx, Ly = fine.shape[0] * cx, fine.shape[1] * cy
    cell = max(CELL, np.sqrt(Lx * Ly / MAXCELLS))
    nx = max(1, int(round(Lx / cell))); ny = max(1, int(round(Ly / cell)))
    hole_max = hole_max_cm2 if hole_max_cm2 is not None else (HOLE_MAX_CM2 if r['mat'] == 'Al' else HOLE_MAX_OTHER_CM2)
    mat = coarse_grid(fine, nx, ny); ccx, ccy = Lx / nx, Ly / ny
    hole = coarse_grid(small_holes(fine, cx * cy, hole_max), nx, ny)
    weight = mat - HOLE_COST * hole - OPEN_COST * (1 - mat - hole)
    Aw = (R @ np.array([a, b, n]).T).T
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij'); ii = ii.ravel(); jj = jj.ravel()
    local = base_pt + np.outer(lo[0] + (ii + 0.5) * ccx, a) + np.outer(lo[1] + (jj + 0.5) * ccy, b) + zmid * n
    Cw = (R @ local.T).T + t + shift
    hc = np.array([ccx / 2, ccy / 2, thick / 2])
    # classes: k inside bank k, -1 world, -2 conflict
    cls = np.full(len(Cw), -2)
    wsh = Box(wsol.c + shift, wsol.A, wsol.h)
    ovb = {k: overlap_depth(wsh, B) for k, B in banks.items()}
    own = max(ovb, key=ovb.get) if max(ovb.values()) > -GAP else min(banks, key=lambda q: np.linalg.norm(banks[q].c - wsh.c))
    for k in banks:
        cls[scene.bank_free(Cw, Aw, hc, k)] = k
    cls[scene.world_free(Cw, Aw, hc, own)] = -1
    mflat = mat.ravel(); mtot = mflat.sum()
    perpart[r['label']][0] += r['vol']
    nbox = 0
    for c in sorted(set(cls.tolist()) - {-2}):
        freeg = np.zeros((nx, ny), bool); freeg[ii[cls == c], jj[cls == c]] = True
        matc = np.where(freeg, mat, 0.0); want = matc.sum(); got = 0.0
        for _ in range(MAXBOX):
            if got >= COVER * want:
                break
            _, rect = best_rectangle(freeg, weight)
            if rect is None:
                break
            i0, i1, j0, j1 = rect
            m = matc[i0:i1, j0:j1].sum()
            if m < MINPART * mtot:
                break
            freeg[i0:i1, j0:j1] = False; matc[i0:i1, j0:j1] = 0; got += m
            ca = lo[0] + (i0 + i1) / 2 * ccx; cb = lo[1] + (j0 + j1) / 2 * ccy
            hloc = np.array([(i1 - i0) * ccx / 2, (j1 - j0) * ccy / 2, thick / 2])
            vol = r['vol'] * m / mtot
            cw = R @ (base_pt + ca * a + cb * b + zmid * n) + t + shift
            bk = c if c >= 0 else own
            wbox = Box(cw, Aw, hloc, f"{r['label'].replace(' ', '_')}#{r['solid']}.{nbox}", r['mat'])
            nbox += 1
            Rb, tb = scene.bank_frame(bk)
            out.append(dict(set=setname, name=wbox.name, label=r['label'], mat=r['mat'], fill=min(1.0, vol / (8 * np.prod(hloc))),
                            bank=bk, mother='bank' if c >= 0 else 'world', world=wbox, loc=wbox.transformed(Rb.T, -Rb.T @ tb), vol=vol))
            perpart[r['label']][1] += vol; perpart[r['label']][2] += 1
    lost = 1 - perpart_kept(r, out, setname) / r['vol']
    if lost > 0.05:
        report.append(f"{r['label']}#{r['solid']}: {lost:.0%} of the material not placed (conflicts)")


def perpart_kept(r, out, setname):
    prefix = f"{r['label'].replace(' ', '_')}#{r['solid']}."
    return sum(o['vol'] for o in out if o['set'] == setname and o['name'].startswith(prefix))


def make_envelopes(scene, report):
    """group the world boxes of each set and bank into envelope boxes clear of everything else. An envelope is oriented
    along the bank frame, the world frame or the axes of its largest member (the smallest clear one is taken)."""
    out = scene.out; envs = []
    for setname, bk in sorted({(o['set'], o['bank']) for o in out if o['mother'] == 'world'}):
        Rb, tb = scene.bank_frame(bk)
        members = [o for o in out if o['mother'] == 'world' and o['bank'] == bk and o['set'] == setname]
        blocked = collections.Counter()

        def envelope(group, F):
            """the box (world) around the members of the group with the axes F (rows, world frame), 0.5 mm margin"""
            Q = np.vstack([o['world'].corners() for o in group]) @ F.T
            lo, hi = Q.min(0) - 0.5, Q.max(0) + 0.5
            return Box(F.T @ ((lo + hi) / 2), F, (hi - lo) / 2)

        def clear(E, group):
            """'' if the envelope E of the group is clear of everything else, else the reason"""
            ingroup = {id(o) for o in group}
            others = [(o['world'], o['name']) for o in out if o['mother'] != 'bank' and id(o) not in ingroup]
            # the envelopes-to-be of the other groups (of 2 or more boxes) and the envelopes already made
            others += [(g['env'], 'group') for g in groups if len(g['members']) > 1 and not any(id(o) in ingroup for o in g['members'])]
            others += [(e['world'], e['name']) for e in envs]
            for k, B in scene.banks.items():
                if overlap_depth(E, B) > -GAP:
                    return f'bank {k}'
            for O, name in [(O, O.name) for O in scene.world_other] + others:
                if overlap_depth(E, O) > -GAP:
                    return name
            for cal in scene.cals:
                Ec = E.transformed(*cal['move'][bk])
                for k, B in cal['banks'].items():
                    if overlap_depth(Ec, B) > -GAP:
                        return f'bank {k} ({cal["name"]})'
                if any(overlap_depth(Ec, O) > -GAP for O in cal['world']):
                    return f'world volume ({cal["name"]})'
                for o in out:
                    if o['mother'] != 'bank' and id(o) not in ingroup and overlap_depth(Ec, o['world'].transformed(*cal['move'][o['bank']])) > -GAP:
                        return f'{o["name"]} ({cal["name"]})'
                for e in envs:
                    if overlap_depth(Ec, e['world'].transformed(*cal['move'][e['bank']])) > -GAP:
                        return f'{e["name"]} ({cal["name"]})'
            if scene.rmax and np.max(np.hypot(*E.corners()[:, :2].T)) > scene.rmax - GAP:
                return 'radius'
            return ''

        def best_envelope(group):
            """the smallest clear envelope of the group over the candidate orientations: (box, reason if none)"""
            largest = max(group, key=lambda o: np.prod(o['world'].h))
            frames = [Rb.T, np.eye(3), largest['world'].A]
            cands = sorted((envelope(group, F) for F in frames), key=lambda E: np.prod(E.h))
            why = ''
            for E in cands:
                why = clear(E, group)
                if not why:
                    return E, ''
            return None, why

        groups = [dict(members=[o], env=envelope([o], o['world'].A)) for o in members]
        while True:                                # merge the pair adding the least volume, as long as it is clear
            pairs = []
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    g = groups[i]['members'] + groups[j]['members']
                    tight = min(np.prod(envelope(g, F).h) for F in (Rb.T, np.eye(3)))
                    pairs.append((8 * tight - 8 * np.prod(groups[i]['env'].h) - 8 * np.prod(groups[j]['env'].h), i, j))
            for _, i, j in sorted(pairs):
                g = groups[i]['members'] + groups[j]['members']
                E, why = best_envelope(g)
                if E is not None:
                    groups[i] = dict(members=g, env=E); del groups[j]
                    break
                blocked[why] += 1
            else:
                break
        for g in groups:
            if len(g['members']) < 2:
                continue                           # a single box needs no envelope
            E = g['env']
            name = f"env{bk}_{sum(1 for e in envs if e['bank'] == bk and e['set'] == setname)}"
            envs.append(dict(set=setname, name=name, bank=bk, world=E, loc=E.transformed(Rb.T, -Rb.T @ tb)))
            for o in g['members']:
                o['mother'] = name
        report.append(f"{setname}, bank {bk}: {len(members)} world boxes -> {sum(1 for e in envs if e['bank'] == bk and e['set'] == setname)} "
                      f"envelopes, {sum(1 for o in members if o['mother'] == 'world')} boxes left in the world; merges blocked by "
                      f"{dict(blocked.most_common(4))}")
    return envs


def write(S, out, envelopes):
    with open(S['output'], 'w') as f:
        f.write(f"# LOKI CAD component set '{S['name']}' (simplified, {len(out)} boxes, {len(envelopes)} envelopes): {S.get('description', '')}\n")
        f.write(f"# generated by gen_simplified.py from {', '.join(sorted({s['step'].split('/')[-1] for s in S['selections']}))}\n")
        f.write("# box <name> <material> <fill> <bankId> <mother> cx cy cz ax ay az bx by bz nx ny nz hx hy hz\n"
                "# envelope <name> <bankId> cx cy cz ax ay az bx by bz nx ny nz hx hy hz\n"
                "# (bank frame of bankId [mm]; fill: the density factor of the homogenised box; mother: bank, world or an envelope)\n")
        for e in envelopes:
            loc = e['loc']
            f.write(f"envelope {e['name']} {e['bank']} " + ' '.join(f'{x:.4f}' for x in loc.c) + ' '
                    + ' '.join(f'{x:.8f}' for x in loc.A.ravel()) + ' ' + ' '.join(f'{x:.4f}' for x in loc.h) + '\n')
        for o in out:
            loc = o['loc']; h = loc.h - SHRINK
            f.write(f"box {o['name']} {o['mat']} {o['fill']:.4f} {o['bank']} {o['mother']} " + ' '.join(f'{x:.4f}' for x in loc.c) + ' '
                    + ' '.join(f'{x:.8f}' for x in loc.A.ravel()) + ' ' + ' '.join(f'{x:.4f}' for x in h) + '\n')


def main(cfg):
    C = json.load(open(cfg))
    scene = Scene(C)
    relevance = {}
    for fn in C.get('relevance', []):
        for s in json.load(open(fn))['solids']:
            relevance[(s['label'], s['idx'])] = s['relevance_cm3']
    report = []; perpart = collections.defaultdict(lambda: [0.0, 0.0, 0])
    for S in C['sets']:
        for sel in S['selections']:
            al = np.load(sel['align'])
            tess = np.load(f"tess_{sel['tag']}.npz")
            chosen, xsec = select_solids(sel, relevance, report)
            for r in chosen:
                add_solid(scene, S['name'], r, al['R'], al['t'], tess, xsec, perpart, report, sel.get('hole_max_cm2'))
    envelopes = make_envelopes(scene, report) if C.get('envelopes') else []
    for S in C['sets']:
        so = [o for o in scene.out if o['set'] == S['name']]; se = [e for e in envelopes if e['set'] == S['name']]
        write(S, so, se)
        print(f"{S['name']}: {len(so)} boxes, {len(se)} envelopes written to {S['output']}")
        print('  by (material, mother, bank):', dict(collections.Counter(
            (o['mat'], o['mother'] if o['mother'] in ('bank', 'world') else 'envelope', o['bank']) for o in so)))
    for line in report:
        print('  ' + line)
    print('  per part (CAD cm3 -> kept cm3, boxes):')
    for k, v in sorted(perpart.items(), key=lambda kv: -kv[1][0]):
        print(f'    {k:20s} {v[0] / 1e3:8.0f} -> {v[1] / 1e3:8.0f}  {v[2]}')


if __name__ == '__main__':
    main(sys.argv[1])
