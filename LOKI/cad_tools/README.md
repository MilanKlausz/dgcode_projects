# CAD-based optional components for the LOKI Geant4 geometry

Tools that turn the STEP models of the LOKI detector system (`cad_models/`) into simple, overlap-free boxes for the
Geant4 geometry `G4GeoLoki/GeoBCSBanks`, and that check the existing geometry against the CAD models. Nothing here is
built by simplebuild (no `pkg.info`); the scripts need FreeCAD 0.21 (`freecadcmd`) for the STEP files and the dgcode
environment (G4GeoLoki, G4GeoLokiTests) for the Geant4 side.

## Using the components in a simulation

All components are optional and off by default:

| parameter | meaning |
|---|---|
| `cad_components` | comma separated component sets (data files `G4GeoLoki/data/cad_components_<set>.txt`, or file paths) |
| `cad_pack_components` | sets of boxes added to every detector pack (`G4GeoLoki/data/cad_pack_<set>.txt`) |
| `cad_B4C_material` | material of the B4C parts (default `MAT_B4C:b10_enrichment=0.95`, as the boron masks) |
| `cad_Tufnol_material` | Tufnol (phenolic laminate) parts (default `G4_BAKELITE`) |
| `cad_Al_material` | aluminium alloy parts (default `NCrystal:cfg=Al_sg225.ncmat`, as the pack Al panel) |

e.g. `sb_loki_sim cad_components=tufnol,frames_al cad_pack_components=rail ...`

Volumes are named `CadComponent-<material>-<set>:<CAD part>#<solid>.<box>` and `CadPack-<material>-<set>:<part>.<box>`
(copy numbers -11 and -12: passive). Pack boxes are siblings of the tubes, so the copy-number depths used by the
analyses are unchanged.

### The sets (October 2026)

| set | content | boxes | volume (CAD -> boxes) |
|---|---|---|---|
| `front_b4c_2021` | the B4C plates added to the front carriage with the updated boron masks (SI-7615-678): SI-7615-695/-696/-697/-703 (side plates, 300 x 600 x 5 mm, outside banks 6/8) and -702/-704 (bottom plates, 250 x 600 x 5 mm, in front of bank 7) | 8 | 5.1 -> 5.0 l |
| `tufnol` | the Tufnol parts of the detector frames: front carriage (wedge pads -229, plates -217), the Pi-shaped 16 mm frames -558 behind the masks of banks 1/3, the 4 mm side plates of the rear bank (-195, -219, -197, -208); not the small blocks -196 and the wedges -180 (not plate-like) | 67 | 13.2 -> 12.2 l |
| `frames_al` | aluminium plates (and compact parts) >= 200 cm3 of the front carriage (SI-7615-678) and of the middle and rear frames (SI-7615-451); hollow / complex parts (beams, extrusions) are not included | 422 | 258 -> 178 l |
| `rail` (pack) | the aluminium C-channel behind the B4C/Al panels of every pack (SI-7615-141/-162), 6 boxes per pack, L + 29 mm long | 6 per pack | 97% of the section |

`tufnol` is generated to be compatible with `front_b4c_2021`, `frames_al` with both (a later set yields to the earlier
ones where the CAD parts touch). All sets are free of overlaps with the nominal and with the SAM-606 bank placement
(construction-time check and the dense checks of `dense_check.py`: rear bank at 5 m and 10 m, calibration slits, vacuum
tank).

What is lost (`frames_al`): material where a CAD plate collides with an existing Geant4 volume (a mask passing through
the footprint of a frame plate; the CAD and Geant4 masks agree within about 3 mm), where it would cross a bank-box face
in a way a box cannot follow, or where it would collide with a neighbouring bank in a calibrated placement. The largest
losses are the end walls of the middle frame (SI-7615-112, -401), the front frame plate -108, the perforated plates -246
and the rear frame plates -119/-121/-123.

## Validation of the existing geometry against the CAD models

* Pack (SI-7615-131): tube positions, the L-shaped B4C panel (cross-section 320 mm2 in Geant4, 314 mm2 in the CAD) and
  the Al panel (215 vs 214 mm2) agree; see the mid-length section `fig_pack_section.png`.
* Boron masks: every rectangular mask of the front banks (SI-7615-678), of the middle banks and of the rear bank
  (SI-7615-451) is within about 1 mm of its CAD plate; the Geant4 boxes cover 98-100% of the CAD plates and every
  Geant4 mask lies on CAD B4C (point tests with a +-3 mm tolerance across the 5 mm plates). Non-rectangular plates are
  decomposed into boxes of the right area (SI-7615-103: 250 x 1119 + 2 x 173 x 305; -490 likewise). Simplifications
  found: the rear masks SI-7615-221 are two 870 mm halves staggered by about +-4 mm in depth (one 1740 mm box in
  Geant4); the boxes 1-5/1-6/3-5/3-6 (the decomposition of -490) are about 3 mm deeper than the CAD plate.
* The triangular masks (SI-7615-691 to -694) agree with the CAD parts.
* The front masks follow the updated design of SI-7615-678 (2021: -476/-511 210 mm wide, -477/-510 1045 mm long).

## Workflow (to regenerate or add a set)

Environment: `CAD_MODELS` (the `cad_models` directory), `CAD_WORK` (a work directory), `CAD_TOOLS` (this directory).

1. Inventories of the STEP files (FreeCAD; every solid with all instances, oriented bounding box, name):
   `python step_colours.py FILE.stp col_<tag>.json` (material from the colours: B4C dark grey, Tufnol brown),
   `STEP=FILE.stp OUT=obb_<tag>.json COL=col_<tag>.json freecadcmd -c "exec(open('inventory3.py').read())"`,
   `STEP=FILE.stp OUT=tess_<tag>.npz freecadcmd -c "exec(open('tess.py').read())"`.
   Materials (`cadlib.py`): colours plus the B4C part list of `Loki detector frame B4C component matrix.pdf`.
2. The Geant4 geometry: `python g4dump.py g4_nominal.json` (and `bank_calibration=2026-September-SAM-606` for
   `g4_sam606.json`).
3. Alignment CAD -> Geant4 world (`align.py`): from the B4C plates that match Geant4 masks. SI-7615-678 is already in the
   Geant4 frame (axes permuted, offset 0); the frames of SI-7615-451 are placed one by one (`align_subset`).
4. Checks: `pts_cad.py` (points in CAD solids) + `pts_g4.py` (located in Geant4: coverage; points in Geant4 volumes) +
   `pts_support.py` (support by CAD solids); `shadow.py` (fraction of a solid in front of the tubes as seen from a
   source, and whether a boron mask shields it); `plot_cad.py` (projections with the bank boxes).
5. A set: `python gen_components.py configs/cfg_<set>.json` (plates are rasterised in their plane, every cell is
   classified as inside bank b / outside the banks / conflict with a separating-axis test against the bank boxes, their
   content, the world volumes, the other sets (`avoid`) and the calibrated placements (`calibrations`); each class is
   covered by disjoint rectangles). Pack sets: `slice.py` + `gen_pack.py`.
6. Copy the data file to `G4GeoLoki/data/`, rebuild, run `dense_check.py` (all sets together).
