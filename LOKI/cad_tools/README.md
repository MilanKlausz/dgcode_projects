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

e.g. `sb_loki_sim cad_components=simple_tufnol,simple_al cad_pack_components=simple_ends ...`

Volumes are named `CadComponent-<material>-<set>:<CAD part>#<solid>.<box>` and `CadPack-<material>-<set>:<part>.<box>`
(copy numbers -11 and -12: passive), envelopes `CadEnvelope-<set>:<name>` (copy number -13). Pack boxes are siblings of
the tubes and envelopes are placed in the world, so the copy-number depths used by the analyses are unchanged.

### The simplified sets (branch cad-simplified, recommended)

Few volumes per CAD part: each selected solid is at most a few boxes per region (inside a bank volume / in the world),
**homogenised**: a box may contain the small enclosed holes of the part (bolt holes, the windows of a truss, the inside of
a hollow beam: aluminium holes up to 400 cm2, Tufnol and B4C up to 10 cm2) and has the part's material with the density
scaled by its filled fraction, so that the mass is kept. Large openings (the windows of the U- and Pi-shaped frame plates
in front of the tubes) stay open. The density factors are rounded to steps of 5% (relative), i.e. a few materials per set.
Aluminium parts are selected by their **relevance** (`relevance.py`): the volume of a part times the fraction of directions
in which a straight line from a point of the part reaches a detector tube before any B4C/Cd/Gd (the chance that a neutron
scattered there flies straight into a tube). Parts hidden behind the B4C panels, masks or other banks, and the small
brackets, score low and are left out. World boxes of a set are grouped into **envelopes** (boxes of the world material)
where that is possible without touching a bank box or another volume (most frame parts wrap around a bank box, so many
stay directly in the world).

| set | content | volumes |
|---|---|---|
| `simple_b4c2021` | the B4C plates added with the updated boron masks of the front carriage (SI-7615-695/-696/-697/-703, -702/-704; as `front_b4c_2021`) | 8 boxes, 3 envelopes |
| `simple_tufnol` | the Tufnol parts with a relevance >= 1.5 cm3 (Tufnol scatters ~30x more than aluminium per cm): the Pi-shaped frames of banks 1/3 (3 boxes each), the side plates of the rear bank, the pads and plates of the front carriage, the wedges -180; 12.8 l | 18 boxes |
| `simple_al` | the aluminium of the front carriage, middle and rear frames with a relevance >= 40 cm3 (97% of the total relevance): plates, hollow beams (one homogenised box each) and the beams of the carriage ring, without the thin carriage-ring bar SI-7615-178_01/beam18 (8 x 8 mm, 40 cm3 of Al); 213 l (the detailed `frames_al`: 228 l) | 136 boxes, 8 envelopes |
| `simple_ends` (pack) | the end blocks and end plates beyond both straw ends of every pack (548 cm3 per end, the whole CAD volume) as two homogenised boxes per end (upper / lower half of the pack section, fill 0.50) | 4 per pack |

The C-channel of the packs (`rail`) has no simplified version: it had no significant effect on any setup. All four sets
together, and any subset, are free of overlaps (construction-time check; `dense_check.py simple_b4c2021,simple_tufnol,simple_al
simple_ends`: nominal and SAM-606 placement, rear bank at 5 m and 10 m, calibration slits, vacuum tank; the setup-B
strips and plates). Generated with `gen_simplified.py configs/cfg_simplified.json` and `gen_pack_simplified.py
configs/cfg_pack_ends.json cad_pack_simple_ends.txt simple_ends` (see the workflow below).

**Simulation time** (CPU time per event of `sb_loki_sim` with GRIFF output, AmBe source with the reference setup: air,
tank, agreed masks; 200k events x 6 runs per variant, each run pinned to its own core, the runs spread over both sockets
of one DMSC quark node per source position; relative to the geometry without CAD components, uncertainty about +-2%):

| variant | nominal source position (0, 0, 0) | setup B (source among the front banks) | setup E (in front of the rear bank) |
|---|---|---|---|
| detailed sets `front_b4c_2021,tufnol,frames_al` + pack `rail,ends` | +13% | +17% | +11% |
| simplified `simple_b4c2021,simple_tufnol,simple_al` + pack `simple_ends` | +8% | +10% | +4% |
| the same, all component materials replaced by air (geometry only) | +3% | +6% | +3% |
| the same without envelopes | +7% | +12% | +6% |
| `simple_tufnol` / `simple_al` / pack `simple_ends` alone | 0% / +4% / +2% | 0% / +5% / +3% | 0% / +1% / +1% |

The geometry itself (navigation through 163 boxes, 12 envelopes and 4 boxes per pack; measured before beam18 was
excluded, now 162 boxes and 11 envelopes) costs 3-6%; the rest is the
physics of the added material (more neutron scattering, longer histories: the GRIFF output grows like the CPU time). The
envelopes change the time by about 1-2% (within the uncertainty): the world has 91 daughters with the simplified sets
(17 without CAD components). Memory: +20 MB per process (detailed sets: +55 MB). Initialisation: unchanged (~11 s).
The geometry-only cost of the detailed sets is +5-9%. `cad_Al_material=G4_Al` (instead of NCrystal) is not faster.

**Physics: the simplified sets give the same results as the detailed ones** (AmBe reference setups, the same seeds,
rates relative to the geometry without CAD components; statistical errors 0.05-0.12 for the back-lit banks of setup E,
<= 0.01 for setup B):

| | bank 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| E (100M): detailed `tufnol,frames_al` + `rail,ends` | 1.025 | 1.48 | 1.39 | 1.56 | 1.30 | 0.93 | 1.19 | 0.92 | 1.14 |
| E (100M): simplified `simple_tufnol,simple_al` + `simple_ends` | 1.024 | 1.43 | 1.49 | 1.55 | 1.41 | 1.06 | 1.20 | 0.98 | 1.11 |
| B (50M): detailed `tufnol,frames_al` | 0.94(4) | 1.021 | 1.008 | 1.018 | 1.001 | 1.073 | 0.967 | 0.986 | 0.953 |
| B (50M): simplified `simple_tufnol,simple_al` | 1.01(4) | 1.012 | 1.019 | 1.029 | 1.014 | 1.069 | 0.970 | 0.984 | 0.956 |

The excess of the data at the tube ends of the rear bank (setup E, data/sim at the ends relative to the centre, layers
0-3): 1.56 / 1.48 / 1.44 / 1.46 without components, 1.46 / 1.33 / 1.30 / 1.32 with the detailed sets, 1.42 / 1.33 / 1.30 /
1.29 with the simplified sets.

### The detailed sets (October 2026)

| set | content | boxes | volume (CAD -> boxes) |
|---|---|---|---|
| `front_b4c_2021` | the B4C plates added to the front carriage with the updated boron masks (SI-7615-678): SI-7615-695/-696/-697/-703 (side plates, 300 x 600 x 5 mm, outside banks 6/8) and -702/-704 (bottom plates, 250 x 600 x 5 mm, in front of bank 7) | 8 | 5.1 -> 5.0 l |
| `tufnol` | the Tufnol parts of the detector frames: front carriage (wedge pads -229, plates -217), the Pi-shaped 16 mm frames -558 behind the masks of banks 1/3, the 4 mm side plates of the rear bank (-195, -219, -197, -208); not the small blocks -196 and the wedges -180 (not plate-like) | 67 | 13.2 -> 12.2 l |
| `frames_al` | aluminium plates (and compact parts) >= 200 cm3 of the front carriage (SI-7615-678) and of the middle and rear frames (SI-7615-451), and the hollow beams as extrusions (SI-7615-223_01, -209_01 and the beams of the welded carriage ring SI-7615-178_01, split by `beam_split.py`); other complex parts (brackets) are not included | 613 | 330 -> 225 l |
| `rail` (pack) | the aluminium C-channel behind the B4C/Al panels of every pack (SI-7615-141/-162), 6 boxes per pack, L + 29 mm long | 6 per pack | 97% of the section |
| `ends` (pack) | the aluminium end blocks (SI-7615-134/-135, in slabs) and end plates (SI-7615-073) beyond both straw ends of every pack, clipped to the pack volume (the real pack reaches 24 mm further forward and beyond the back) | 2 x 34 per pack | 586 -> 325 cm3 per end |

`tufnol` is generated to be compatible with `front_b4c_2021`, `frames_al` with both (a later set yields to the earlier
ones where the CAD parts touch). All sets are free of overlaps with the nominal and with the SAM-606 bank placement
(construction-time check and the dense checks of `dense_check.py`: rear bank at 5 m and 10 m, calibration slits, vacuum
tank; the boxes stay within 1620 mm of the beam axis, inside the default tank radius of 1625 mm).

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
  decomposed into boxes of the right area (SI-7615-103: 250 x 1119 + 2 x 173 x 305; -490 likewise). The rear masks
  0-0/0-2 (one 1740 mm box each) are the two 870 mm plates SI-7615-221 of each side, at the same depth (59.00 mm in front of
  -222/-223 in the CAD, 59 mm in the Geant4 table); the boxes 1-5/1-6/1-7 coincide with SI-7615-490 (plane offsets to the
  neighbouring masks equal within 0.05 mm). Note: an alignment fitted only to plates near the centre line of a frame can be
  rotated by a fraction of a degree (the first rear-frame fit was off by 0.44 deg, which made the -221 plates look staggered
  by +-4 mm); `align_full.py` fits all plates by their corners (residuals <= 1.4 mm).
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
3. Alignment CAD -> Geant4 world (`align.py`, refined with `align_full.py`): from the B4C plates that match Geant4 masks. SI-7615-678 is already in the
   Geant4 frame (axes permuted, offset 0); the frames of SI-7615-451 are placed one by one (`align_subset`).
4. Checks: `pts_cad.py` (points in CAD solids) + `pts_g4.py` (located in Geant4: coverage; points in Geant4 volumes) +
   `pts_support.py` (support by CAD solids); `shadow.py` (fraction of a solid in front of the tubes as seen from a
   source, and whether a boron mask shields it); `plot_cad.py` (projections with the bank boxes).
5. Hollow beams: `extrude_sections.py` (sections of elongated non-box solids at 25/50/75% of their length: an extrusion
   has equal sections) and `beam_split.py` (a welded frame of straight beams split into beams, each an extrusion).
6. A set: `python gen_components.py configs/cfg_<set>.json` (plates are rasterised in their plane, every cell is
   classified as inside bank b / outside the banks / conflict with a separating-axis test against the bank boxes, their
   content, the world volumes, the other sets (`avoid`) and the calibrated placements (`calibrations`); each class is
   covered by disjoint rectangles). Pack sets: `slice.py` + `gen_pack.py` (parts along the tubes from one section; end parts
   in slabs along the tubes, one section per slab, merged where consecutive slabs are equal).
7. Copy the data file to `G4GeoLoki/data/`, rebuild, run `dense_check.py` (all sets together); `plot_components.py`
   draws the Geant4 boxes over the CAD silhouettes.
8. Simplified sets: points in all solids (`pts_cad.py` with `MAT=Al,Tufnol,B4C,dark MINVOL=50 N=60`), their relevance
   (`python relevance.py relevance_nominal.json pts_*.json`: rays through the Geant4 geometry with
   `G4GeoLokiTests.TestUtils.ray_visibility`), then `python gen_simplified.py configs/cfg_simplified.json` (all simplified
   sets in one run, so that they are compatible; `"envelopes": false` for a version without envelopes) and
   `python gen_pack_simplified.py configs/cfg_pack_ends.json cad_pack_simple_ends.txt simple_ends`.
