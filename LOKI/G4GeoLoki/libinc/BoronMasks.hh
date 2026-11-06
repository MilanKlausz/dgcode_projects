#ifndef G4GeoLoki_BoronMasks_hh
#define G4GeoLoki_BoronMasks_hh

#include "G4Materials/NamedMaterialProvider.hh"
#include <array>
#include <string>

class BoronMasks {
public:
  static G4Material* maskMaterial;
  /// regular rectangular masks ///
  static int getNumberOfBoronMasks(const int bankId);
  static double getSize(const int bankId, const int maskId, const int axisIndex); // 0: width, 1: height, 2: thickness
  /// The size in the bank frame (x = depth: the thickness, y = across the tubes, z = along the tubes).
  static std::array<double,3> getSizeInBankFrame(const int bankId, const int maskId);
  static double getPosition(const int bankId, const int maskId, const int axisIndex);
  static double getRotation(const int bankId, const int maskId);
  /// The CAD part of a mask: the ESS drawing number of the B4C part ("Loki detector frame B4C component matrix";
  /// several masks of a bank can be boxes of the same part, e.g. the three top parts of a front bank: SI-7615-103).
  static std::string getDrawing(const int bankId, const int maskId);
  /// The volume name: "BoronMask-<bank>-<drawing>", with ".<k>" if the drawing occurs more than once in the bank (the
  /// k-th such mask in the order of the mask ids, from 0).
  static std::string getName(const int bankId, const int maskId);
  /// The beamstop (bank 0, beamstop_id; no CAD part, sizes from ESS-1178830).
  static const std::string beamstopName;

  /// triangular masks ///
  static double getHalfSizeOfTriangularMask(const int maskId, const int axisIndex);
  static double getCutPointOfTriangularMask(const int maskId, const int axisIndex);
  static double getCutDirOfTriangularMask(const int maskId, const int axisIndex);
  static double getBankIdOfTriangularMask(const int maskId);
  static double getPosInBankOfTriangularMask(const int maskId, const int axisIndex);
  static std::string getDrawingOfTriangularMask(const int maskId);
  /// "BoronMask-triangular-<bank>-<drawing>"
  static std::string getNameOfTriangularMask(const int maskId);

private:
  /// regular rectangular masks ///
  static double getBoronMaskParameter(const int bankId, const int maskId, const int parameterIndex);
  const static std::array<std::array<double, 7>, 6> rearBoronMasks;
  const static std::array<std::array<double, 7>, 8> midVerticalBoronMasks;
  const static std::array<std::array<double, 7>, 8> midHorizontalBoronMasks;
  const static std::array<std::array<double, 7>, 8> frontTopBoronMasks;
  const static std::array<std::array<double, 7>, 8> frontBottomBoronMasks;
  const static std::array<std::array<double, 7>, 8> frontVerticalBoronMasks;
  const static std::array<const char*, 6> rearDrawings;
  const static std::array<const char*, 8> midVerticalDrawings;
  const static std::array<const char*, 8> midHorizontalDrawings;
  const static std::array<const char*, 8> frontTopDrawings;
  const static std::array<const char*, 8> frontBottomDrawings;
  const static std::array<const char*, 8> frontVerticalDrawings;

  /// triangular masks ///
  static double getTriangularBoronMaskParameter(const int maskId, const int parameterIndex);
  const static std::array<std::array<double, 10>, 4> triangularBoronMasks;
  const static std::array<const char*, 4> triangularDrawings;
};

#endif
