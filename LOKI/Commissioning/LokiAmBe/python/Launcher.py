"""The LOKI simulation (LOKI/python/Launcher.py) with the AmBe source in addition (event_gen=ambe)."""


def launch(user_data=None):
    import G4GeoLoki.GeoBCSBanks as GEO
    import LOKI.Launcher
    from LokiAmBe.AmBeSourceGen import AmBeSourceGen
    LOKI.Launcher.launch(GEO.create(), extra_generators={'ambe': lambda launcher: AmBeSourceGen()}, user_data=user_data)
