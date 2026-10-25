"""Projections of CAD solids (OBB outlines, coloured by material) in the Geant4 world frame, with the Geant4 bank boxes.
usage: python plot_cad.py TAG ALIGN.npz OUT.png [x0,x1,y0,y1 (CAD region)]"""
import sys, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from cadlib import load, obb_corners
from matplotlib.collections import PolyCollection
tag, aln, out = sys.argv[1:4]
reg = [float(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
a = np.load(aln); R, t = a['R'], a['t']
COL = {'B4C': ('k', 0.75), 'Tufnol': ('C1', 0.8), 'Al': ('C0', 0.12), 'dark': ('C3', 0.6), 'tube': ('C2', 0.3)}
T = np.load(f'tess_{tag}.npz')
rows = [r for r in load(tag) if reg is None or (reg[0] < r['centre'][0] < reg[1] and reg[2] < r['centre'][1] < reg[3])]
G = json.load(open('g4_nominal.json'))['volumes']
banks = [v for v in G if v['lv'] == 'Bank']
views = (('front view (along the beam)', 0, 1), ('side view', 2, 1), ('top view', 2, 0))
fig, axs = plt.subplots(1, 3, figsize=(21, 7.5))
for ax, (title, i, j) in zip(axs, views):
    for mat in ('Al', 'tube', 'dark', 'Tufnol', 'B4C'):
        c, al = COL[mat]
        polys = []
        for r in rows:
            if r['mat'] != mat or f"v{r['solid']}" not in T: continue
            V = (R @ T[f"v{r['solid']}"].T).T + t; F = T[f"f{r['solid']}"]
            polys.extend(V[F][:, :, [i, j]])
        if polys:
            ax.add_collection(PolyCollection(polys, fc=c, ec='none', alpha=min(1, al * 0.25) if mat == 'Al' else al, lw=0))
    for b in banks:
        Rb = np.array(b['R']); h = np.array(b['pars']) * (1 if b['solid'] == 'G4Box' else 0.5)
        C = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]) * h
        P = (Rb @ C.T).T + np.array(b['t']); q = P[:, [i, j]]
        hh = ConvexHull(q); ax.add_patch(Polygon(q[hh.vertices], closed=True, fc='none', ec='m', lw=1.0, ls='--'))
        ax.text(*np.array(b['t'])[[i, j]], str(b['copy']), color='m', fontsize=12, ha='center')
    ax.set_title(title); ax.set_xlabel('xyz'[i] + ' [mm]'); ax.set_ylabel('xyz'[j] + ' [mm]'); ax.set_aspect('equal'); ax.autoscale_view()
from matplotlib.lines import Line2D
axs[0].legend([Line2D([], [], color=COL[m][0], lw=6, alpha=max(COL[m][1], 0.4)) for m in COL] + [Line2D([], [], color='m', ls='--')],
              list(COL) + ['Geant4 bank boxes'], fontsize=9, loc='upper right')
fig.suptitle(f'CAD {tag} in the Geant4 frame (oriented bounding boxes), coloured by material', fontsize=12)
fig.tight_layout(); fig.savefig(out, dpi=90)
