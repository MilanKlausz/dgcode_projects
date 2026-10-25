"""Shared helpers for the CAD analysis: load the OBB inventory of a STEP file and classify each solid by material.

Material from the STEP colours (step_colours.py) and the B4C part list of 'Loki detector frame B4C component matrix.pdf':
  B4C     dark grey (0.21, 0.2, 0.189) or a part number of the B4C matrix (S-DOUGH: 82% B4C, 18% epoxy)
  Tufnol  brown (0.61, 0.209, 0.121) (e.g. SI-7615-180 'detector frame wedge plate', TUFNOL CARP)
  dark    (0.2, 0.2, 0.2) / (0.3, 0.226, 0.121): unknown (electronics, PCB, rubber?)
  tube    the detector tubes (SI-7615-061*)
  Al      everything else (Al 6082 frame plates, pack parts)"""
import os
import re
import json
import numpy as np

CAD = os.environ.get('CAD_WORK', '.') + '/'   # the directory of the inventories (obb_*.json, col_*.json, tess_*.npz)
B4C_PARTS = {f'SI-7615-{n}' for n in (476, 477, 104, 103, 519, 479, 510, 511, 512, 513, 485, 486, 117, 488, 490, 492, 493,
                                      559, 555, 556, 221, 222, 223, 192, 143, 514, 518, 480, 489, 554,
                                      691, 692, 693, 694, 695, 696, 697, 702, 703, 704)}
COL = {(0.21, 0.2, 0.189): 'B4C', (0.61, 0.209, 0.121): 'Tufnol', (0.2, 0.2, 0.2): 'dark', (0.3, 0.226, 0.121): 'dark',
       (0.898, 0.0, 0.0): 'B4C', (1.0, 0.498, 0.0): 'B4C'}


def base(label):
    return re.sub(r'(_\d+)+$', '', label)


def load(tag):
    """tag: 131, 243, frames, 097, 678"""
    rows = json.load(open(f'{CAD}obb_{tag}.json'))
    cols = json.load(open(f'{CAD}col_{tag}.json'))
    for r in rows:
        b = base(r['label'])
        cs = [tuple(float(x) for x in c.split(',')) for c in cols.get(r['label'], cols.get(b, {}))]
        mat = 'Al'
        if b.startswith('SI-7615-061'):
            mat = 'tube'
        elif any(COL.get(c) == 'B4C' for c in cs) or (not cs and b in B4C_PARTS):
            mat = 'B4C'          # the part list only decides for parts without colour (colour copies are lost in some files)
        elif any(COL.get(c) == 'Tufnol' for c in cs):
            mat = 'Tufnol'
        elif cs and all(COL.get(c) == 'dark' for c in cs):
            mat = 'dark'
        r['base'] = b; r['mat'] = mat; r['colours'] = cs
        r['com'] = np.array(r['com']); r['axes'] = np.array(r['axes']); r['ext'] = np.array(r['ext'])
        r['size'] = r['ext'][:, 1] - r['ext'][:, 0]
        # OBB centre (the COM is not the box centre for asymmetric solids)
        r['centre'] = r['com'] + r['axes'].T @ r['ext'].mean(1)
    return rows


def obb_corners(r):
    c = []
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                c.append(r['com'] + r['axes'][0] * r['ext'][0][i] + r['axes'][1] * r['ext'][1][j] + r['axes'][2] * r['ext'][2][k])
    return np.array(c)
