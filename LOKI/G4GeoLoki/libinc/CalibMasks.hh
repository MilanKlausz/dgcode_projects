#ifndef G4GeoLoki_CalibMasks_hh
#define G4GeoLoki_CalibMasks_hh

#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include "G4Materials/NamedMaterialProvider.hh"

/// Calibration masks: CalibMasksBase, a slit mask of absorber strips (the Cd mask of the Larmor geometry), and
/// RearPlate, the rear calibration mask plate of LOKI.
class CalibMasks {
public:
  /// The REAR CALIBRATION MASK PLATE of LOKI (drawing SI-7615-993, ESS-5767239): a 2 mm aluminium sheet, its front face
  /// painted with Gd2O3 paint except ten 6 mm bare slits; the face 1742 mm tall and 996 mm wide, the paint between its side
  /// walls (992 mm); mounted with its lips on the front of the bank's side plates, the face 51 mm deeper. Lengths in mm.
  class RearPlate {
  public:
    static constexpr double thickness = 2.0;
    static constexpr double faceHeight = 1742.0;
    static constexpr double faceWidth = 996.0;
    static constexpr double paintedWidth = 992.0;
    /// across the tubes over the painted width: paint, slit, paint, ..., paint (43, then (6, 94) x 9, 6, 43: slit
    /// centres at +-50, ..., +-450 mm from the centre)
    static std::vector<double> getPaintPattern();
    static int getNumberOfPaintStrips();
    static std::string getPlateName();             // "CalibPlate-0-SI-7615-993"
    static std::string getPaintName(const int k);  // "CalibPaint-0-SI-7615-993.<k>"
  };

  class CalibMasksBase {
    public:
      CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, double elevationFromTubeFront, std::vector<double> pattern);

      std::string getName() const;
      double getThickness() const;
      double getHeight() const;
      double getLeftTubeEndDistance() const;
      double getElevationFromTubeFront() const;
      std::vector<double> getPattern() const;
      double getWidth() const;
    private:
      const std::string m_name;
      const double m_thickness;
      const double m_height;
      const double m_leftTubeEndDistance;
      const double m_elevationFromTubeFront;
      const std::vector<double> m_pattern;
  };
  static G4Material* maskMaterial;
  static G4Material* maskBoxMaterial;
};

#endif
