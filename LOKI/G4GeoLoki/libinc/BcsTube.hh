#ifndef G4GeoLoki_BcsTube_hh
#define G4GeoLoki_BcsTube_hh

#include "G4Materials/NamedMaterialProvider.hh"
#include <string>

class BcsTube {
public:
  /// The tube (CAD part SI-7615-061; variants _02: 1.2 m, _04: 0.5 m, the 1.0 m variant not in the CAD files at hand)
  static const std::string tubeWallName; // "TubeWall-SI-7615-061"

    /// materials ///
  static G4Material* tubeWallMaterial;
  static G4Material* tubeInnerGas;
  static G4Material* strawWallMaterial;
  static G4Material* converterMaterial;
  static G4Material* countingGas;

    /// straw ///
  static double getStrawOuterRadius();
  static double getStrawInnerRadius();
  static double getStrawWallThickness();

  static double getFrontTubeConverterThickness();
  static double getBackTubeConverterThickness();

    /// tube ///
  static double getTubeOuterRadius();
  static double getTubeInnerRadius();

  static double getStrawPositionX(const int strawId);
  static double getStrawPositionY(const int strawId);

private:
  const static double strawOuterRadius;
  const static double strawWallThickness;

  const static double frontTubeConverterThickness;
  const static double backTubeConverterThickness;

  const static double tubeOuterRadius;
  const static double tubeWallThickness;

  const static double strawPositionsInTube[7][2];
  const static double strawStrawDistance;
};

#endif
