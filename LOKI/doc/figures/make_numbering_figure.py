#!/usr/bin/env python3
"""Writes numbering.svg: the tube and straw numbering of the rear bank (bank 0), side view, from the Geant4
geometry (G4GeoLoki AimHelper, nominal placement). Run in the dgcode environment after building (sb):

    python3 LOKI/doc/figures/make_numbering_figure.py
"""

import pathlib
import G4GeoLoki.LokiAimHelper as LokiAim
import G4GeoLokiTests.TestUtils as TU

NSTRAWS = 7
BANK = 0
LAYERS, ROWS = 4, 8  # the top 8 rows of the 4 layers
TUBE_RADIUS, STRAW_RADIUS = 12.7, 3.75  # [mm]
SCALE = 3.2  # px per mm (left part)
ZOOM = 14.0  # px per mm (straw inset)

aim = LokiAim.AimHelper(5000.0, 1, 9, LokiAim.NOMINAL_BANK_CALIBRATION)  # 1 pixel per straw: pixel centre = straw centre
npacks = TU.getNumberOfPacksByBankId(BANK)


def straw_centres(tube):
    """(z, y) [mm] of the 7 straw centres (side view: z along the beam, y up)."""
    return [aim.getPixelCentreCoordinates(aim.getBankPixelOffset(BANK) + tube * NSTRAWS + s)[1:][::-1]
            for s in range(NSTRAWS)]


def tube_centre(tube):
    zs, ys = zip(*straw_centres(tube))
    return sum(zs) / NSTRAWS, sum(ys) / NSTRAWS


tubes = {(layer, row): layer * 2 * npacks + row for layer in range(LAYERS) for row in range(ROWS)}
centres = {key: tube_centre(t) for key, t in tubes.items()}
z0 = min(z for z, _ in centres.values()) - TUBE_RADIUS
y1 = max(y for _, y in centres.values()) + TUBE_RADIUS
X0, Y0 = 190.0, 60.0  # left-top margin [px]


def to_px(z, y, scale=SCALE, x0=X0, y0=Y0, zref=z0, yref=y1):
    return x0 + (z - zref) * scale, y0 + (yref - y) * scale


parts = []
for (layer, row), tube in tubes.items():
    x, y = to_px(*centres[(layer, row)])
    parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{TUBE_RADIUS * SCALE:.1f}" class="tube"/>')
    parts.append(f'<text x="{x:.1f}" y="{y + 5:.1f}" class="id">{tube}</text>')
for layer in range(LAYERS):
    x, _ = to_px(*centres[(layer, 0)])
    parts.append(f'<text x="{x:.1f}" y="{Y0 - 14:.1f}" class="label">layer {layer}</text>')
_, ytop = to_px(*centres[(0, 0)])
_, ybottom = to_px(*centres[(0, ROWS - 1)])
xbeam = X0 - 25
parts.append(f'<line x1="{X0 - 135}" y1="{(ytop + ybottom) / 2:.1f}" x2="{xbeam:.1f}" y2="{(ytop + ybottom) / 2:.1f}" '
             'class="beam" marker-end="url(#arrow)"/>')
parts.append(f'<text x="{X0 - 135}" y="{(ytop + ybottom) / 2 - 10:.1f}" class="label" text-anchor="start">beam</text>')
parts.append(f'<text x="{X0 - 135}" y="{ybottom + TUBE_RADIUS * SCALE + 30:.1f}" class="note" text-anchor="start">rows {ROWS}-{2 * npacks - 1} '
             f'continue downwards ({2 * npacks} rows, {npacks} packs)</text>')

# inset: the straws of tube 0
zc, yc = centres[(0, 0)]
ix0, iy0 = to_px(*centres[(LAYERS - 1, 0)])[0] + 120, Y0 + 20
cx, cy = ix0 + TUBE_RADIUS * ZOOM, iy0 + TUBE_RADIUS * ZOOM
parts.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{TUBE_RADIUS * ZOOM:.1f}" class="tube"/>')
for s, (z, y) in enumerate(straw_centres(tubes[(0, 0)])):
    x, yy = cx + (z - zc) * ZOOM, cy - (y - yc) * ZOOM
    parts.append(f'<circle cx="{x:.1f}" cy="{yy:.1f}" r="{STRAW_RADIUS * ZOOM:.1f}" class="straw"/>')
    parts.append(f'<text x="{x:.1f}" y="{yy + 6:.1f}" class="id big">{s}</text>')
parts.append(f'<text x="{cx:.1f}" y="{iy0 - 8:.1f}" class="label">straws of tube 0</text>')

width = cx + TUBE_RADIUS * ZOOM + 30
height = max(ybottom + TUBE_RADIUS * SCALE + 50, cy + TUBE_RADIUS * ZOOM + 30)
svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" font-family="sans-serif">
<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
<path d="M0,0 L10,5 L0,10 z" fill="#c0392b"/></marker></defs>
<style>
.tube {{ fill: #dbe6f3; stroke: #3d6a99; stroke-width: 1.5; }}
.straw {{ fill: #f6e3c8; stroke: #b07a2a; stroke-width: 1.5; }}
.id {{ font-size: 13px; text-anchor: middle; fill: #1d2b3a; }}
.big {{ font-size: 17px; }}
.label {{ font-size: 14px; text-anchor: middle; fill: #1d2b3a; }}
.note {{ font-size: 12px; fill: #555; }}
.beam {{ stroke: #c0392b; stroke-width: 2.5; }}
</style>
<rect width="100%" height="100%" fill="white"/>
{chr(10).join(parts)}
</svg>
'''
out = pathlib.Path(__file__).parent.parent / 'images' / 'numbering.svg'
out.parent.mkdir(exist_ok=True)
out.write_text(svg)
print(f'written {out}')
