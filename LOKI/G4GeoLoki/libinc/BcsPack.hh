#ifndef G4GeoLoki_BcsPack_hh
#define G4GeoLoki_BcsPack_hh

#include "G4Materials/NamedMaterialProvider.hh"
#include <string>

class BcsPack {
public:
  static double getHorizontalTubeDistanceInPack();
  static double getVerticalTubeDistanceInPack();
  static double getTubeGridParallelogramAngle();

  static double getTopRowOffsetInPack();
  static double getHorizontalTubeCentreOffsetInPack();

  static double getTubeGridParallelogramSide();
  static double getTubeCentreDistanceFromPackFront();

  static G4Material* packBoxFillMaterial;
  static double getPackBoxWidth();
  static double getPackBoxHeight();
  static double getPackBoxIdleLengthOnOneEnd();

  static double getTubeRotation();
  static double getHorizontalTubeOffset(const int inPackTubeId);
  static double getVerticalTubeOffset(const int inPackTubeId);

  /// B4C panel parts ///
  static G4Material* B4CPanelMaterial;
  static double getB4CLengthOverStrawOnOneEnd();

  /// The B4C panel (CAD part SI-7615-143_01, S-DOUGH) as boxes: "B4CPanel-SI-7615-143_01.<partId>"
  static int getNumberOfB4CParts();
  static std::string getB4CPartName(const int partId);
  static double getB4CPartThickness(const int partId);
  static double getB4CPartHeight(const int partId);
  static double getB4CPartHorizontalOffset(const int partId);
  static double getB4CPartVerticalOffset(const int partId);

  /// Al panel parts ///
  static G4Material* AlPanelMaterial;
  /// The Al panel (CAD part SI-7615-142) as boxes: "AlPanel-SI-7615-142.<partId>"
  static int getNumberOfAlParts();
  static std::string getAlPartName(const int partId);
  static double getAlPartThickness(const int partId);
  static double getAlPartHeight(const int partId);
  static double getAlPartHorizontalOffset(const int partId);
  static double getAlPartVerticalOffset(const int partId);

private:
  const static double tubeGridParallelogramBase;
  const static double tubeGridParallelogramSide;
  const static double tubeGridParallelogramAngle;

  const static double tubeRotationAngle;

  const static double packBoxWidth;
  const static double packBoxHeight;
  const static double packBoxIdleLengthOnOneEnd;
  const static double tubeCentreDistanceFromPackFront; //1st lower (closest) tube

  const static double packHolderDistanceFromPackTop;
  const static double packHolderDistanceFromPackFront;

  /// B4C panel parts ///
  const static double B4CLengthOverStrawOnOneEnd;
  const static double B4CDistanceFromLastTubeCentre;
  const static double B4CPanelPartThickness[3];
  const static double B4CPanelPartHeight[3];
};

#endif
