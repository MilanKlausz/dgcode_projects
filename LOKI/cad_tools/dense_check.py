"""Dense overlap checks (as G4GeoLokiTests/test_geo_dense_overlaps) of the geometry with CAD component sets."""
import sys
import G4GeoLoki.GeoBCSBanks as GEO
import G4GeoLokiTests.TestUtils as TU
from G4GeoLokiTests.testhelpers import captured_fd_output
sets = sys.argv[1] if len(sys.argv) > 1 else 'front_b4c_2021,tufnol'
for cal in ('nominal-geant4-geometry', '2026-September-SAM-606'):
    for extra in ({}, {'rear_detector_distance_m': 10.0, 'with_calibration_slits': True},
                  {'vacuum_tank_wall_mm': 10.0, 'rear_detector_distance_m': 9.985615, 'world_material': 'G4_AIR', 'pack_fill_material': 'G4_AIR'}):
        geo = GEO.create(); geo.bank_calibration = cal; geo.cad_components = sets; geo.cad_pack_components = 'rail,ends'
        for k, v in extra.items():
            setattr(geo, k, v)
        with captured_fd_output() as lines:
            TU.construct_world(geo)
        bad = [l for l in lines if l.startswith('Checking overlaps') and not l.rstrip().endswith('OK!')]
        ov = TU.dense_world_overlaps(200000, 1e-6)
        out = TU.dense_bank_containment(50000, 1e-6)
        print(f'{cal:26s} {str(extra):95s}: construction problems {len(bad)}, dense world overlaps {len(ov)}, outside bank {len(out)}', flush=True)
        for o in list(ov)[:5] + list(out)[:5]:
            print('    ', o)
