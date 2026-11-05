import G4GeoLoki.LokiAimHelper as LokiAim

def launch(geo):
    import G4Launcher
    launcher = G4Launcher()

    launcher.addParameterInt("analysis_straw_pixel_number", 256)
    launcher.addParameterDouble("gen_x_offset_meters", 0.0)
    launcher.addParameterBoolean('masking_only',False)
    launcher.addParameterInt('aiming_bank_id',-1) #aim only at a certain bank

    #geometry:
    launcher.setGeo(geo)

    #generator:
    from LokiMasking.MaskingSourceGen import MaskingSourceGen as Gen
    gen = Gen()
    gen.exposeParameter("rear_detector_distance_m",geo,"geo_rear_detector_distance_m")
    gen.gen_x_offset_meters = launcher.getParameterDouble('gen_x_offset_meters')
    gen.aiming_bank_id = launcher.getParameterInt('aiming_bank_id')
    if geo.hasParameterString("bank_calibration"):
        gen.exposeParameter("bank_calibration",geo,"geo_bank_calibration") #aim at the pixels of the calibrated banks
    launcher.setGen(gen)

    def addUserData():
      launcher.setUserData("analysis_straw_pixel_number", str(launcher.getParameterInt('analysis_straw_pixel_number')))
      launcher.setUserData("rear_detector_distance_m", str(launcher.getGeo().getParameterDouble("rear_detector_distance_m")))
      launcher.setUserData("aiming_bank_id", str(launcher.getParameterInt('aiming_bank_id')))
      for key, value in LokiAim.bankPlacementUserData(launcher.getGeo()).items():
        launcher.setUserData(key, value)

    launcher.addPrePreInitHook(addUserData) #add userdata when all parameters are available

    #filter:
    if not launcher.getParameterBoolean('masking_only'):#TODO reduced material list should be the default
        launcher.setOutput('loki_masking','REDUCED')
    else:
        # the converters and every absorbing volume (masks, beamstop, pack panels, rear calibration paint; the names
        # from G4GeoLoki: BcsBanks::getAbsorberVolumeNames)
        griff_output_volumes = ["Converter", *LokiAim.absorberVolumeNames()]
        import G4CollectFilters.StepFilterVolume
        f = G4CollectFilters.StepFilterVolume.create()
        f.volumeList = griff_output_volumes
        launcher.setFilter(f)
        launcher.setOutput('loki_masking','REDUCED')

    #launch:
    launcher.go()
