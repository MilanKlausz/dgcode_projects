"""Pixel positions and ids of the LoKI detector banks (see G4GeoLoki/AimHelper.hh)"""

from . import _LokiAimHelper as _impl

# the bank calibration names (see G4GeoLoki/BankCalibration.hh)
NOMINAL_BANK_CALIBRATION = _impl.nominal_bank_calibration()
DEFAULT_BANK_CALIBRATION = _impl.default_bank_calibration()
DEFAULT_NUMBER_OF_PIXELS_IN_STRAW = _impl.default_number_of_pixels_in_straw()

# the text of a bank calibration file ("" for the nominal geometry), with the banks shifted by bank_shifts_mm
# ("bank:dx,dy,dz;..." [mm], see BankCalibration::shifted); throws for an unknown name or a wrong shift
def bankCalibrationText( name, shifts = "" ):
    return _impl.bankCalibrationText( name, shifts )

# the user data of the bank placements of a GeoBCSBanks geometry: the analysis uses the recorded text, so that it has the
# placements of the simulation (with the bank shifts) even if the calibration file changes or is not there
def bankPlacementUserData( geo ):
    name, shifts = geo.getParameterString( "bank_calibration" ), geo.getParameterString( "bank_shifts_mm" )
    return { "bank_calibration": name, "bank_shifts_mm": shifts, "bank_calibration_text": bankCalibrationText( name, shifts ) }

# the volume names of the geometry that analyses look for (see BcsBanks::getAbsorberVolumeNames)
def absorberVolumeNames():
    return [ _impl.absorber_volume_name( i ) for i in range( _impl.number_of_absorber_volume_names() ) ]

TUBE_WALL_NAME = _impl.tube_wall_name()

def AimHelper( rearBankDistance, strawPixelNumber = DEFAULT_NUMBER_OF_PIXELS_IN_STRAW,
               numberOfBanks = 9, bankCalibration = NOMINAL_BANK_CALIBRATION ):
    """Create an AimHelper object (the defaults are those of the C++ class)"""
    return _impl.AimHelper( float( rearBankDistance ), int( strawPixelNumber ),
                            int( numberOfBanks ), bankCalibration )
