"""Refined alignment of a separately placed frame (SI-7615-451) to the Geant4 masks, with every B4C plate:
plates matched to Geant4 masks of equal size, and plates represented in Geant4 by several boxes or by part of a box
(SI-7615-221: one half of a 1740 mm box; SI-7615-490 / -103: three boxes) matched to the bounding box of the union / to
the half box, in the plate frame. Fit: Kabsch on the 8 OBB corners of every pair (corners matched by position after
the initial alignment), which constrains rotations much better than plate centres alone."""
import sys, json, numpy as np
from cadlib import load
from obb import Box
from align import kabsch, g4_boxes


def corners_sorted(B, ref):
    """the 8 corners of box B ordered to match the 8 corners of ref (nearest neighbour after alignment)"""
    cb = B.corners(); cr = ref.corners()
    order = [int(np.argmin(np.linalg.norm(cb - c, axis=1))) for c in cr]
    return cb[order]


def union_box(boxes, axes):
    """axis-aligned bounding box (in the frame 'axes') of a set of boxes"""
    P = np.vstack([b.corners() for b in boxes]) @ axes.T
    lo, hi = P.min(0), P.max(0)
    return Box(((lo + hi) / 2) @ axes, axes, (hi - lo) / 2)


def main(tag, region, banks, init, out, special):
    a = np.load(init); R, t = a['R'], a['t']
    cad = [r for r in load(tag) if r['mat'] == 'B4C' and np.sort(r['size'])[0] < 7
           and region[0] < r['centre'][0] < region[1] and region[2] < r['centre'][1] < region[3]]
    g4 = [g for g in g4_boxes('g4_nominal.json') if int(g['name'].split('-')[1]) in banks]
    G = {g['name']: Box(g['centre'], g['axes'], g['size'] / 2, g['name']) for g in g4}
    for it in range(4):
        P, Q, used = [], [], []
        for r in cad:
            w = Box(r['centre'], r['axes'], r['size'] / 2).transformed(R, t)
            target = None
            sp = special.get(r['base'])
            if sp == 'half':                               # half of a Geant4 box twice as long
                cands = [b for b in G.values() if np.allclose(np.sort(2 * b.h), np.sort(r['size']) * [1, 1, 2], atol=12)
                         or np.allclose(np.sort(2 * b.h), np.sort([r['size'][0], r['size'][1], 2 * r['size'][2]]), atol=12)]
                best = None
                for b in G.values():
                    k = int(np.argmax(b.h))
                    for s in (-1, 1):
                        hb = b.h.copy(); hb[k] /= 2
                        half = Box(b.c + s * b.A[k] * hb[k], b.A, hb, b.name)
                        if np.allclose(np.sort(2 * half.h), np.sort(r['size']), atol=8):
                            d = np.linalg.norm(half.c - w.c)
                            if best is None or d < best[0]:
                                best = (d, half)
                target = best[1] if best and best[0] < 40 else None
            elif sp:                                       # a union of Geant4 boxes
                parts = [G[n] for n in G if n in sp]
                cand = []
                for grp in sp:                             # sp: list of name groups; choose the group nearest
                    bs = [G[n] for n in grp if n in G]
                    if bs:
                        u = union_box(bs, bs[0].A); cand.append((np.linalg.norm(u.c - w.c), u))
                if cand:
                    d, u = min(cand, key=lambda x: x[0])
                    target = u if d < 40 else None
            else:
                cands = [b for b in G.values() if np.allclose(np.sort(2 * b.h), np.sort(r['size']), atol=3)]
                if cands:
                    b = min(cands, key=lambda b: np.linalg.norm(b.c - w.c))
                    target = b if np.linalg.norm(b.c - w.c) < 40 else None
            if target is None:
                continue
            cw = corners_sorted(w, target)
            # back to the CAD frame for the fit
            P.append((cw - t) @ R); Q.append(target.corners()); used.append((r['label'], target.name))
        P = np.vstack(P); Q = np.vstack(Q)
        R, t = kabsch(P, Q)
    res = np.linalg.norm((P @ R.T + t) - Q, axis=1).reshape(-1, 8)
    print(f'{out}: {len(used)} plates, corner residual rms {np.sqrt((res ** 2).mean()):.2f} mm, max {res.max():.2f} mm')
    for (lab, name), rr in zip(used, res):
        print(f'   {lab:16s} -> {name:16s} corners {rr.mean():5.2f} mm (max {rr.max():.2f})')
    d = np.degrees(np.arccos(np.clip((np.trace(R.T @ a['R']) - 1) / 2, -1, 1)))
    print(f'   rotation change w.r.t. the first alignment: {d:.3f} deg; translation change {np.round(t - a["t"], 2)} mm')
    np.savez(out, R=R, t=t)


if __name__ == '__main__':
    main('frames', (12000, 15000, 30000, 32000), [0], 'align_frames_rear.npz', 'align_frames_rear_v2.npz',
         {'SI-7615-221': 'half'})
    groups13 = [['BoronMask-1-5', 'BoronMask-1-6', 'BoronMask-1-7'], ['BoronMask-3-5', 'BoronMask-3-6', 'BoronMask-3-7']]
    main('frames', (9500, 11500, 24000, 25300), [1, 2, 3, 4], 'align_frames_mid13.npz', 'align_frames_mid_v2.npz',
         {'SI-7615-490': groups13})
