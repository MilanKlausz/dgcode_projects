"""CAD silhouettes (by material, light) with the Geant4 CAD-component boxes (outlines) and masks, three projections.
usage: python plot_components.py TAG ALIGN.npz G4DUMP.json OUT.png [x0,x1,y0,y1 (CAD region)] [zmin,zmax (world z window)]"""
import sys, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from cadlib import load
tag, aln, g4f, out = sys.argv[1:5]
reg = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 and sys.argv[5] != '-' else None
zw = [float(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else None
a = np.load(aln); R, t = a['R'], a['t']; T = np.load(f'tess_{tag}.npz')
COL = {'B4C': 'k', 'Tufnol': 'C1', 'Al': 'C0'}
rows = [r for r in load(tag) if r['mat'] in COL and (reg is None or (reg[0] < r['centre'][0] < reg[1] and reg[2] < r['centre'][1] < reg[3]))]
G = json.load(open(g4f))['volumes']
def corners(v):
    Rb = np.array(v['R']); h = np.array(v['pars']) * (1 if v['solid'] == 'G4Box' else 0.5)
    C = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)]) * h
    return (Rb @ C.T).T + np.array(v['t'])
comp = [v for v in G if v['pv'].startswith('CadComponent')]
masks = [v for v in G if v['pv'].startswith('BoronMask')]
banks = [v for v in G if v['lv'] == 'Bank']
inwin = lambda P: zw is None or (zw[0] < P[:, 2].mean() < zw[1])
views = (('front view (along the beam)', 0, 1), ('side view', 2, 1), ('top view', 2, 0))
fig, axs = plt.subplots(1, 3, figsize=(24, 8))
for ax, (title, i, j) in zip(axs, views):
    for mat, c in COL.items():
        polys = []
        for r in rows:
            if r['mat'] != mat or f"v{r['solid']}" not in T: continue
            V = (R @ T[f"v{r['solid']}"].T).T + t
            if not inwin(V): continue
            polys.extend(V[T[f"f{r['solid']}"]][:, :, [i, j]])
        if polys:
            ax.add_collection(PolyCollection(polys, fc=c, ec='none', alpha=0.08 if mat == 'Al' else 0.25))
    for vs, ec, lw in ((banks, 'm', 0.8), (masks, '0.3', 0.5), (comp, None, 0.6)):
        for v in vs:
            P = corners(v)
            if not inwin(P): continue
            q = P[:, [i, j]]
            try:
                h = ConvexHull(q)
            except Exception:
                continue
            col = ec or COL[v['pv'].split('-')[1]]
            ax.add_patch(Polygon(q[h.vertices], closed=True, fill=False, ec=col, lw=lw, ls='--' if vs is banks else '-'))
    ax.set_title(title); ax.set_xlabel('xyz'[i] + ' [mm]'); ax.set_ylabel('xyz'[j] + ' [mm]'); ax.set_aspect('equal'); ax.autoscale_view()
from matplotlib.lines import Line2D
axs[0].legend([Line2D([], [], color=c, lw=6, alpha=0.4) for c in COL.values()] + [Line2D([], [], color=c, lw=1) for c in COL.values()]
              + [Line2D([], [], color='0.3', lw=1), Line2D([], [], color='m', ls='--')],
              [f'CAD {m}' for m in COL] + [f'Geant4 {m} component boxes' for m in COL] + ['Geant4 boron masks', 'Geant4 bank boxes'], fontsize=8, loc='upper right')
fig.suptitle(f'CAD {tag} (filled) vs the Geant4 CAD-component boxes (outlines), Geant4 frame', fontsize=13)
fig.tight_layout(); fig.savefig(out, dpi=80)
