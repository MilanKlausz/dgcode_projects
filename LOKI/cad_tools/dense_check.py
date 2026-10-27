"""Dense overlap checks (as G4GeoLokiTests/test_geo_dense_overlaps) of the geometry with CAD component sets: world
daughters (incl. CAD world boxes and envelopes) against each other, the daughters of the banks and of the envelopes
inside their mother. Construction-time checks (Geant4) for everything else (pack boxes, siblings in a mother). With
bank mothers (with_bank_mothers=yes), also the daughters of each bank mother against each other and inside it.

usage: python dense_check.py [cad_components] [cad_pack_components] [par=value ...]
(par=value: geometry parameters for all configurations, e.g. with_bank_mothers=yes)"""
import sys
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output
sets = sys.argv[1] if len(sys.argv) > 1 else 'front_b4c_2021,tufnol'
packsets = sys.argv[2] if len(sys.argv) > 2 else 'rail,ends'
value = lambda v: {'yes': True, 'no': False}.get(v, v)
common = {k: value(v) for k, v in (a.split('=', 1) for a in sys.argv[3:])}
B_STRIPS = ('1,-124.8,695.6,2846.6,-593.3,262.6,2888.4,61.2,1;2,528.6,358.2,3306.7,199.2,-225.2,3331.0,38.2,1;'
            '3,-566.8,-258.6,2888.9,-2.3,-691.2,2846.4,38,1;4,-208.3,-362.0,3328.9,-537.7,376.3,3302.7,55,1')  # setup B
for cal in ('nominal-geant4-geometry', '2026-September-SAM-606'):
    for extra in ({}, {'rear_detector_distance_m': 10.0, 'with_rear_calibration_mask': True},
                  {'vacuum_tank_wall_mm': 10.0, 'rear_detector_distance_m': 9.985615, 'world_material': 'G4_AIR', 'pack_fill_material': 'G4_AIR'},
                  {'vacuum_tank_wall_mm': 10.0, 'vacuum_tank_rear_m': 5.6, 'with_rear_calibration_mask': True, 'absorber_strips': B_STRIPS}):
        geo = GEO.create(); geo.bank_calibration = cal; geo.cad_components = sets; geo.cad_pack_components = packsets
        for k, v in {**common, **extra}.items():
            setattr(geo, k, v)
        with captured_fd_output() as lines:
            TU.construct_world(geo)
        bad = [l for l in lines if l.startswith('Checking overlaps') and not l.rstrip().endswith('OK!')]
        ov = TU.dense_world_overlaps(200000, 1e-6)
        out = TU.dense_bank_containment(50000, 1e-6, ['Bank', 'CadEnvelope-', 'BankMother-'])
        mothers = [l.strip() for l in lines if l.startswith('GeoBCSBanks: bank mother')]
        print(f'{cal:26s} {str(extra)[:95]:95s}: construction problems {len(bad)}, dense world overlaps {len(ov)}, outside bank {len(out)}', flush=True)
        for o in list(ov)[:5] + list(out)[:5] + [l for l in lines if l.startswith('ERROR') or 'xception' in l][:5]:
            print('    ', o)
        for m in mothers:
            print('    ', m)
