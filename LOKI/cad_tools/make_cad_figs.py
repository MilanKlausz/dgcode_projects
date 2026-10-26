"""Figures of the CAD-components report: the detector pack section (CAD vs Geant4, with the rail boxes)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Circle
from cadlib import load


def pack_section(out='fig_pack_section.png'):
    S = json.load(open('sec_131.json')); lab = {str(r['solid']): r['base'] for r in load('131')}
    M = np.array([[0, -1], [1, 0.]]); t = np.array([-18.84, 0])
    fig, ax = plt.subplots(figsize=(12, 4.6))
    for i, ws in S.items():
        col = {'SI-7615-141': 'C0', 'SI-7615-143': 'k', 'SI-7615-148': 'C2', 'SI-7615-061': 'C1'}[lab[i]]
        for w in ws:
            ax.add_patch(Polygon((M @ np.array(w).T).T + t, fc=col, ec='none', alpha=0.35))
    dx, dy, al = 106.15, 27.6, np.radians(13.45)
    ax.add_patch(Polygon(np.array([[-dx - dy * np.tan(al), -dy], [dx - dy * np.tan(al), -dy], [dx + dy * np.tan(al), dy],
                                   [-dx + dy * np.tan(al), dy]]), fill=False, ec='m', lw=1.5, label='Geant4 pack volume (G4Para)'))
    for k, (x0, x1, y0, y1, c) in enumerate(((7.65, 12.65, -24, 27.2, 'k'), (-3.15, 7.65, -24, -19, 'k'), (-3.15, -0.15, -27.2, -24, 'k'),
                                              (12.65, 15.85, -27.2, 27.2, 'g'), (-0.15, 12.65, -27.2, -24, 'g'))):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=c, lw=1.2, ls='--',
                               label=('Geant4 B4C panel' if k == 0 else 'Geant4 Al panel' if k == 3 else None)))
    for k in range(8):
        ax.add_patch(Circle((-96.35 + (k % 4) * 27 + (6.61 if k < 4 else 0), 13.81 if k < 4 else -13.81), 12.7, fill=False, ec='C1',
                            ls='--', label='Geant4 tubes' if k == 0 else None))
    for j, line in enumerate(l for l in open('cad_pack_rail.txt') if not l.startswith('#')):
        x0, x1, y0, y1 = map(float, line.split()[2:6])
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec='C0', lw=1.4,
                               label='new: rail boxes (cad_pack_components=rail)' if j == 0 else None))
    ax.set_aspect('equal'); ax.set_xlim(-120, 125); ax.set_ylim(-34, 34)
    ax.set_xlabel('pack x (depth) [mm]'); ax.set_ylabel('pack y [mm]')
    ax.legend(fontsize=7.5, loc='center right', bbox_to_anchor=(1.0, 0.45))
    ax.set_title('Detector pack at mid-length: CAD SI-7615-131 (filled: tubes orange, B4C black, Al panel green, C-channel '
                 'SI-7615-141 blue) vs Geant4 (outlines)', fontsize=10)
    fig.tight_layout(); fig.savefig(out, dpi=100)


if __name__ == '__main__':
    pack_section()
