"""Rigid transform CAD (one STEP file) -> Geant4 world, from the B4C plates that match Geant4 boron masks.

Candidate rotations: the 24 proper signed permutation matrices (CAD axes along world axes); the translation from each
CAD plate / G4 mask pair of equal dimensions; scored by the number of CAD plates whose centre lands within TOL of a
G4 mask of the same dimensions. The best is refined with a Kabsch fit over the matched centres.
Returns R, t (world = R @ cad + t), and the matches."""
import itertools
import json
import numpy as np
from cadlib import load

TOL = 15.0


def g4_boxes(fn, prefix=('BoronMask',)):
    G = json.load(open(fn))['volumes']
    out = []
    for v in G:
        if v['solid'] == 'G4Box' and v['pv'].startswith(prefix):
            R = np.array(v['R']); half = np.array(v['pars'])
            out.append(dict(name=v['pv'], centre=np.array(v['t']), axes=R.T, size=2 * half, dims=np.sort(2 * half)))
    return out


def rotations():
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, (p, s) in enumerate(zip(perm, signs)):
                R[i, p] = s
            if np.linalg.det(R) > 0:
                yield R


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def align(tag, g4file='g4_nominal.json', dimtol=3.0):
    cad = [r for r in load(tag) if r['mat'] == 'B4C' and np.sort(r['size'])[0] < 7]
    g4 = g4_boxes(g4file)
    pairs = [(c, g) for c in cad for g in g4 if np.all(np.abs(np.sort(c['size']) - g['dims']) < dimtol)]
    best = (0, None)
    for R in rotations():
        for c, g in pairs:
            t = g['centre'] - R @ c['centre']
            n = 0
            for c2 in cad:
                p = R @ c2['centre'] + t
                if any(np.all(np.abs(np.sort(c2['size']) - g2['dims']) < dimtol) and np.linalg.norm(p - g2['centre']) < TOL for g2 in g4):
                    n += 1
            if n > best[0]:
                best = (n, (R, t))
    R, t = best[1]
    for it in range(3):
        m = []
        for c in cad:
            p = R @ c['centre'] + t
            cands = [g for g in g4 if np.all(np.abs(np.sort(c['size']) - g['dims']) < dimtol)]
            if cands:
                g = min(cands, key=lambda g: np.linalg.norm(p - g['centre']))
                if np.linalg.norm(p - g['centre']) < 3 * TOL:
                    m.append((c, g))
        R, t = kabsch(np.array([c['centre'] for c, g in m]), np.array([g['centre'] for c, g in m]))
    return R, t, m, cad, g4


if __name__ == '__main__':
    import sys
    tag = sys.argv[1]
    R, t, m, cad, g4 = align(tag)
    print(f'{tag}: R =\n{np.round(R, 4)}\nt = {np.round(t, 1)}; {len(m)} matched of {len(cad)} CAD plates')
    for c, g in m:
        p = R @ c['centre'] + t
        print(f"  {c['label']:16s} -> {g['name']:16s} centre diff {np.round(p - g['centre'], 1)}  |d| {np.linalg.norm(p - g['centre']):5.1f} mm")
    np.savez(f'align_{tag}.npz', R=R, t=t)


def align_subset(tag, xr, yr, banks, g4file='g4_nominal.json', dimtol=3.0, out=None):
    """as align(), for the CAD plates with centre x in xr and y in yr (one frame assembly), to the masks of the given banks"""
    global TOL
    cad = [r for r in load(tag) if r['mat'] == 'B4C' and np.sort(r['size'])[0] < 7
           and xr[0] < r['centre'][0] < xr[1] and yr[0] < r['centre'][1] < yr[1]]
    g4 = [g for g in g4_boxes(g4file) if int(g['name'].split('-')[1]) in banks]
    pairs = [(c, g) for c in cad for g in g4 if np.all(np.abs(np.sort(c['size']) - g['dims']) < dimtol)]
    # rotation candidates: from three non-collinear matched pairs of plates (centres), all combinations
    best = (0, None, 1e9)
    import itertools as it
    for (c1, g1), (c2, g2), (c3, g3) in it.combinations(pairs, 3):
        if len({id(c1), id(c2), id(c3)}) < 3 or len({g1['name'], g2['name'], g3['name']}) < 3:
            continue
        P = np.array([c1['centre'], c2['centre'], c3['centre']]); Q = np.array([g1['centre'], g2['centre'], g3['centre']])
        if np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])) < 1e4:
            continue
        R, t = kabsch(P, Q)
        if np.abs(np.linalg.norm(P[1] - P[0]) - np.linalg.norm(Q[1] - Q[0])) > 5:
            continue
        n = 0; err = 0
        for c in cad:
            p = R @ c['centre'] + t
            d = [np.linalg.norm(p - g['centre']) for g in g4 if np.all(np.abs(np.sort(c['size']) - g['dims']) < dimtol)]
            if d and min(d) < TOL:
                n += 1; err += min(d)
        if n > best[0] or (n == best[0] and err < best[2]):
            best = (n, (R, t), err)
    R, t = best[1]
    m = []
    for c in cad:
        p = R @ c['centre'] + t
        cands = [g for g in g4 if np.all(np.abs(np.sort(c['size']) - g['dims']) < dimtol)]
        if cands:
            g = min(cands, key=lambda g: np.linalg.norm(p - g['centre']))
            if np.linalg.norm(p - g['centre']) < TOL:
                m.append((c, g))
    if len(m) >= 3:
        R, t = kabsch(np.array([c['centre'] for c, g in m]), np.array([g['centre'] for c, g in m]))
    print(f'{tag} subset x{xr} y{yr} -> banks {banks}: {len(m)} matched of {len(cad)} plates')
    for c, g in m:
        p = R @ c['centre'] + t
        print(f"  {c['label']:16s} -> {g['name']:16s} |d| {np.linalg.norm(p - g['centre']):5.1f} mm")
    if out:
        np.savez(out, R=R, t=t)
    return R, t, m
