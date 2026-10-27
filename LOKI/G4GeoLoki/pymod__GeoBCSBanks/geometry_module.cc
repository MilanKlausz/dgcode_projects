/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructToPy.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"
#include "G4Vector3D.hh"
#include "G4SubtractionSolid.hh"
#include "G4AffineTransform.hh"
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <string>
#include <cassert>
#include <set>
#include <vector>
#include <sstream>

#include "G4GeoLoki/BcsBanks.hh"

class GeoBCS : public G4Interfaces::GeoConstructBase
{
public:
  GeoBCS();
  virtual ~GeoBCS(){}
  virtual G4VPhysicalVolume* Construct();

protected:
  virtual bool validateParameters();
private:
  BcsBanks* banks;
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(int bankId, int packNumber);
  G4LogicalVolume * createBankLV(int bankId);
  G4LogicalVolume * createTriangularMaskLV(int maskId);
  bool hasBoronMasks(int bankId);
  bool isExcludedBoronMask(int bankId, int maskId);
  bool isListedBoronMask(const std::string& parameter, int bankId, int maskId);
};

// this line is necessary to be able to declare the geometry in the python simulation script
SBLDPY_MODULE( _GeoBCSBanks ) { G4Interfaces::geoToPy<GeoBCS>( m ); }

namespace {
  // the banks of a bank list parameter: "all", "none" (or empty), or comma separated ids and ranges, e.g. "0,5-8"
  std::set<int> parseBankList(const std::string& text, int numberOfBanks = 9) {
    std::set<int> banks;
    if (text == "all") {
      for (int b = 0; b < numberOfBanks; b++)
        banks.insert(b);
      return banks;
    }
    if (text.empty() || text == "none")
      return banks;
    std::stringstream items(text);
    std::string item;
    while (std::getline(items, item, ',')) {
      const auto dash = item.find('-');
      try {
        size_t used = 0;
        const int first = std::stoi(item.substr(0, dash), &used);
        if (used != item.substr(0, dash).size())
          throw std::invalid_argument(item);
        int last = first;
        if (dash != std::string::npos) {
          last = std::stoi(item.substr(dash + 1), &used);
          if (used != item.substr(dash + 1).size())
            throw std::invalid_argument(item);
        }
        if (first < 0 || last >= numberOfBanks || first > last)
          throw std::out_of_range(item);
        for (int b = first; b <= last; b++)
          banks.insert(b);
      }
      catch (const std::logic_error&) {
        throw std::runtime_error("invalid bank list \"" + text + "\" (expected all, none, or ids and ranges such as 0,5-8)");
      }
    }
    return banks;
  }

  // a list of absorber items: ";"-separated items of n comma-separated numbers
  std::vector<std::vector<double>> parseNumberList(const std::string& text, size_t n, const std::string& what, size_t extra = 0) {
    std::vector<std::vector<double>> items;
    std::stringstream parts(text);
    std::string item;
    while (std::getline(parts, item, ';')) {
      if (item.find_first_not_of(" ") == std::string::npos)
        continue;
      std::vector<double> values;
      std::stringstream numbers(item);
      std::string number;
      while (std::getline(numbers, number, ','))
        values.push_back(std::stod(number));
      if ((values.size() != n && values.size() != n + extra) || values[0] < 0 || values[0] > 8)
        throw std::runtime_error("invalid " + what + ": " + item);
      items.push_back(values);
    }
    return items;
  }

  // absorber strips: "bank,x1,y1,z1,x2,y2,z2,width,thickness;..." -> 9 numbers per strip
  std::vector<std::array<double, 9>> parseAbsorberStrips(const std::string& text) {
    std::vector<std::array<double, 9>> strips;
    std::stringstream items(text);
    std::string item;
    while (std::getline(items, item, ';')) {
      if (item.find_first_not_of(" ") == std::string::npos)
        continue;
      std::array<double, 9> values;
      std::stringstream numbers(item);
      std::string number;
      size_t n = 0;
      while (std::getline(numbers, number, ',')) {
        if (n >= 9)
          throw std::runtime_error("absorber strip with more than 9 values: " + item);
        values[n++] = std::stod(number);
      }
      if (n != 9 || values[0] < 0 || values[0] > 8 || values[7] <= 0 || values[8] <= 0)
        throw std::runtime_error("invalid absorber strip (bank,x1,y1,z1,x2,y2,z2,width,thickness): " + item);
      strips.push_back(values);
    }
    return strips;
  }

  // object rotation of a bank (the columns are the bank-local axes in the world)
  G4RotationMatrix toG4Rotation(const BankTransform& transform) {
    const auto& R = transform.rotation;
    return G4RotationMatrix(G4ThreeVector(R[0][0], R[1][0], R[2][0]),
                            G4ThreeVector(R[0][1], R[1][1], R[2][1]),
                            G4ThreeVector(R[0][2], R[1][2], R[2][2]));
  }
}

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

GeoBCS::GeoBCS()
  : GeoConstructBase("G4GeoLoki/GeoBCSBanks"){
  // declare all parameters that can be used from the command line,
  addParameterDouble("rear_detector_distance_m", 5.0, 4.0, 10.0); // default, min, max
  addParameterInt("beamstop_id", 0, 0, 5); // id [1-5] from ESS-1178830 Table 4.6 'Selected beamstop sizes' (0 = no beamstop)
  // the calibration mask of the rear bank (the manufactured REAR CALIBRATION MASK PLATE, drawing ESS-5767239 /
  // SI-7615-993: a 2 mm aluminium sheet, its front face painted with Gd2O3 paint except ten 6 mm bare slits, the face
  // 1742 mm tall and 992 mm wide between its side walls, mounted with its lips on the front of the bank's side plates,
  // the face 51 mm deeper; see CalibMasks::RearPlate)
  addParameterBoolean("with_rear_calibration_mask", false);
  addParameterDouble("calibration_plate_gap_mm", 8.0, 0.0, 50.0);         // back of the 2 mm sheet - front of the tubes
  addParameterDouble("calibration_paint_thickness_mm", 0.15, 0.0, 5.0);   // the Gd2O3 paint layer (not specified; 0.15 mm fits run 360)
  // 50% polyurethane varnish (approximated as C25H42N2O6, 1.2 g/cm3) + 50% Gd2O3 (7.41 g/cm3) by weight
  addParameterString("calibration_paint_material", "solid::Gd10C100H167N8O39/2.07gcm3");
  addParameterString("calibration_plate_material", "G4_Al");
  // the B4C masks of the banks (the rectangular masks of each bank and the triangular masks of banks 5 and 7; not the
  // beamstop): "all" (the design), "none", or a bank list as above
  addParameterString("boron_mask_banks", "all");
  // the triangular B4C masks of banks 5 and 7 (part of their design masks; switch them off separately)
  addParameterBoolean("with_triangular_masks", true);
  // individual design masks not present: "bank:mask,mask;bank:..." (the mask ids of the volume names BoronMask-<bank>-<mask>)
  addParameterString("boron_mask_exclude", "");
  // design masks drawn in colour in the viewers (display only, same format as boron_mask_exclude): red, and orange
  addParameterString("highlight_boron_masks", "");
  addParameterString("highlight_boron_masks_2", "");
  // a cylinder (vertical axis) of the world material marking the source moderator, for the visualisation only (no effect
  // on the physics, but it must not overlap other volumes in a simulation; radius 0: none)
  addParameterDouble("show_source_cylinder_radius_cm", 0.0, 0.0, 100.0);
  addParameterDouble("show_source_cylinder_height_cm", 27.0, 0.1, 200.0);
  addParameterDouble("show_source_cylinder_x_meters", 0.0, -5.0, 5.0);
  addParameterDouble("show_source_cylinder_y_meters", 0.0, -5.0, 5.0);
  addParameterDouble("show_source_cylinder_z_meters", 0.0, -5.0, 15.0);
  // absorber strips on the front faces of banks (e.g. the B4C sheets taped on the banks in source tests):
  // "bank,x1,y1,z1,x2,y2,z2,width,thickness;..." [mm], the strip from the point 1 to the point 2 (world coordinates,
  // projected on the plane just in front of the tubes of the bank), width and thickness of the strip
  addParameterString("absorber_strips", "");
  addParameterString("absorber_strip_material", "MAT_B4C");
  // absorber plates on a side face of a bank (the face across the tubes nearest to the given points, e.g. the
  // inner face of a front bank towards a source between the banks): "bank,x1,y1,z1,x2,y2,z2,width,thickness;..."
  // [mm]: the centre line of the plate from point 1 to point 2 (world coordinates, projected on the plane just outside
  // the face: along the tubes and in depth), its width (in the plane) and thickness (strip material). An optional 10th
  // value (> 0) places the plate inside the bank volume instead, its centre at this distance from the points across the
  // tubes (e.g. a sheet lying directly on the tubes of the outer row)
  addParameterString("absorber_plates", "");

  // Optional vacuum tank around the detector system (background studies): a cylindrical shell around the beam axis
  // (z) with a closed rear end and a front end plate with the beam window, lined with B4C on the inside (0: no tank)
  addParameterDouble("vacuum_tank_wall_mm", 0.0, 0.0, 200.0);
  addParameterDouble("vacuum_tank_inner_radius_m", 1.625, 0.9, 5.0);
  addParameterDouble("vacuum_tank_front_m", 0.1, 0.0, 1.0); // the inner face of the front end plate (z)
  addParameterDouble("vacuum_tank_rear_m", 10.6, 5.0, 11.5); // the inner face of the rear end plate (z)
  addParameterDouble("vacuum_tank_window_radius_m", 0.15, 0.0, 1.0); // the opening of the front end plate
  addParameterDouble("vacuum_tank_front_B4C_mm", 5.0, 0.0, 50.0); // the B4C lining of the front end plate
  addParameterString("vacuum_tank_material", "G4_Al");
  addParameterString("vacuum_tank_lining_material", "MAT_B4C");

  // bank placements: "nominal-geant4-geometry" or the name of a calibration in G4GeoLoki/data/bank_calibration_<name>.txt
  // (or the path of such a file), see G4GeoLoki/BankCalibration.hh
  addParameterString("bank_calibration", BankCalibration::defaultName);

  addParameterString("world_material","G4_Vacuum");
  // the gas around the tubes inside the bank and pack volumes (the measurements at atmospheric pressure: G4_AIR)
  addParameterString("pack_fill_material","G4_Vacuum");
  addParameterString("B4C_panel_material","MAT_B4C:b10_enrichment=0.95"); // the B4C panels behind the tubes
}

G4LogicalVolume * GeoBCS::createTubeLV(double converterThickness, double strawLength){
  const double effectiveStrawLength = strawLength - BcsTube::getStrawWallThickness(); //This is only epsilon difference...

  auto lv_tube = new G4LogicalVolume(new G4Tubs("TubeWall",0,  BcsTube::getTubeOuterRadius(), 0.5*strawLength, 0., 2*M_PI),
                                     BcsTube::tubeWallMaterial, "TubeWall");

  auto lv_empty_tube = place(new G4Tubs("EmptyTube", 0., BcsTube::getTubeInnerRadius(), 0.5*strawLength, 0., 2*M_PI),
                             BcsTube::tubeInnerGas, 0,0,0, lv_tube, G4Colour(0,1,1),-2,0,0).logvol;

  for (int cpNo = 0; cpNo <= 6; cpNo++){
    auto lv_straw_wall = place(new G4Tubs("StrawWall", 0, BcsTube::getStrawOuterRadius(), 0.5*strawLength, 0., 2 * M_PI),
                               BcsTube::strawWallMaterial, BcsTube::getStrawPositionX(cpNo), BcsTube::getStrawPositionY(cpNo), 0, lv_empty_tube, ORANGE, cpNo, 0, 0).logvol;

    auto lv_converter = place(new G4Tubs("Converter", 0., BcsTube::getStrawInnerRadius(), 0.5*effectiveStrawLength, 0., 2 * M_PI),
                              BcsTube::converterMaterial, 0, 0, 0, lv_straw_wall, G4Colour(0, 1, 1), cpNo + 100, 0, 0).logvol;

    place(new G4Tubs("CountingGas", 0., BcsTube::getStrawInnerRadius() - converterThickness, 0.5*effectiveStrawLength, 0., 2 * M_PI),
          BcsTube::countingGas, 0, 0, 0, lv_converter, G4Colour(0, 0, 1), 0, 0, 0);
  }
  return lv_tube;
}

///////////  CREATE PACK BOX LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createPackBoxLV(int bankId, int packNumber){
  const double strawLength = banks->getStrawLengthByBankId(bankId);
  const double packRotation = banks->getPackRotation();
  // Instead of a rectangular box, a detector pack is encapsulated in parallelepiped, to avoid collision of the corners with the calibraion slits after applying the pack rotation.
  // The PackBoxWidth corresponds to the size of the volume encapsulating the electronics on the sides as well, not just the detectors, but that would cause collision with the calibration slits, so a multiplication factor of 0.799 is applied, to get a volume just large enough to fit in the detectors in the front.
  auto lv_pack_box = new G4LogicalVolume(
    new G4Para("EmptyPackBox", 0.799*0.5*BcsPack::getPackBoxWidth(), 0.5*BcsPack::getPackBoxHeight(), 0.5 * strawLength + BcsPack::getPackBoxIdleLengthOnOneEnd(), packRotation, 0, 0),
    getParameterMaterial("pack_fill_material"), "EmptyPackBox");

  /// Add 8 BCS detector tubes ///
  auto lv_front_tube = createTubeLV(BcsTube::getFrontTubeConverterThickness(), strawLength);
  auto lv_back_tube = createTubeLV(BcsTube::getBackTubeConverterThickness(), strawLength);
  G4RotationMatrix* tubeRotationMatrix = new G4RotationMatrix(0, 0, BcsPack::getTubeRotation());

  for (int inPackTubeId = 0; inPackTubeId < 8; inPackTubeId++) {
    place((inPackTubeId % 4 < 2) ? lv_front_tube : lv_back_tube,
          BcsPack::getHorizontalTubeOffset(inPackTubeId), BcsPack::getVerticalTubeOffset(inPackTubeId), 0,
          lv_pack_box, SILVER, banks->getTubeIdInBank(bankId, packNumber, inPackTubeId), 0, tubeRotationMatrix);
  }
  /// Add B4C panel behind detectors in 3 parts ///
  const double B4CLengthHalf = 0.5*strawLength + BcsPack::getB4CLengthOverStrawOnOneEnd();

  auto B4CPanelMaterial = getParameterMaterial("B4C_panel_material");
  for (int partId = 0; partId < 3; partId++){
    place(new G4Box("B4CPanel", 0.5*BcsPack::getB4CPartThickness(partId), 0.5*BcsPack::getB4CPartHeight(partId), B4CLengthHalf),
          B4CPanelMaterial,
          BcsPack::getB4CPartHorizontalOffset(partId), BcsPack::getB4CPartVerticalOffset(partId), 0,
          lv_pack_box, G4Colour(0, 1, 0), -2, 0, new G4RotationMatrix());
    }
  /// Add Al behing the B4C panel in 2 parts ///
  for (int partId = 0; partId < 2; partId++){
    place(new G4Box("AlPanel", 0.5*BcsPack::getAlPartThickness(partId), 0.5*BcsPack::getAlPartHeight(partId), B4CLengthHalf),
          BcsPack::AlPanelMaterial,
          BcsPack::getAlPartHorizontalOffset(partId),  BcsPack::getAlPartVerticalOffset(partId), 0,
          lv_pack_box, SILVER, -2, 0, new G4RotationMatrix());
    }
  return lv_pack_box;
}

///////////  CREATE DETECTOR BANK LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createBankLV(int bankId){
  const int numberOfPacks = banks->getNumberOfPacksByBankId(bankId);

  const double packRotation = banks->getPackRotation();

  // all positions and sizes in the bank volume are in the bank frame (x = depth, y = across the tubes, z = along
  // the tubes), see BcsBanks.hh
  const auto bankHalfSize = banks->getBankHalfSizeInBankFrame(bankId);

  auto lv_bank = new G4LogicalVolume(new G4Box("EmptyPanelBox", bankHalfSize[0], bankHalfSize[1], bankHalfSize[2]),
                                     getParameterMaterial("pack_fill_material"), "Bank");

  // override lv_bank to subtract some empty part of front left and right banks where front top and bottom banks would overlap with them
  // (the notch is computed from the bank placements in use, see BcsBanks::getBankNotch)
  if (banks->hasBankNotch(bankId)) {
    const auto notch = banks->getBankNotch(bankId);
    if (!notch.isNominal)
      printf("GeoBCSBanks: bank %d: notch enlarged for bank_calibration=%s to the half sizes (%g, %g, %g) mm at (%g, %g, %g) mm (nominal: (15, 70, 60) mm at (-145, -42, -425) mm)\n",
             bankId, banks->getBankCalibration().name().c_str(), notch.halfSize[0], notch.halfSize[1], notch.halfSize[2], notch.centre[0], notch.centre[1], notch.centre[2]);
    auto fullBankBox = new G4Box("EmptyPanelBox", bankHalfSize[0], bankHalfSize[1], bankHalfSize[2]);
    auto bankBoxCut = new G4Box("EmptyPanelBox", notch.halfSize[0], notch.halfSize[1], notch.halfSize[2]);
    auto bankBox = new G4SubtractionSolid("EmptyPanelBox", fullBankBox, bankBoxCut, 0, G4ThreeVector(notch.centre[0], notch.centre[1], notch.centre[2]));

    lv_bank = new G4LogicalVolume(bankBox, getParameterMaterial("pack_fill_material"), "Bank");
  }

  for (int packNumber = 0; packNumber < numberOfPacks; ++packNumber){
    auto lv_pack_box = createPackBoxLV(bankId, packNumber);
    const auto packPosition = banks->getPackPositionInBankFrame(bankId, packNumber);
    place(lv_pack_box, packPosition[0], packPosition[1], packPosition[2],
          lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  const int numberOfBoronMasks = hasBoronMasks(bankId) ? BoronMasks::getNumberOfBoronMasks(bankId) : 0;
  for (int maskId = 0; maskId < numberOfBoronMasks; ++maskId){
    if (isExcludedBoronMask(bankId, maskId))
      continue;
    const std::string maskName = "BoronMask-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    const auto maskSize = BoronMasks::getSizeInBankFrame(bankId, maskId);
    const auto maskPosition = banks->getBoronMaskPositionInBankFrame(bankId, maskId);
    const G4Colour maskColour = isListedBoronMask("highlight_boron_masks", bankId, maskId) ? G4Colour(1, 0, 0)
                              : isListedBoronMask("highlight_boron_masks_2", bankId, maskId) ? G4Colour(1, 0.5, 0) : BLACK;
    place(new G4Box(maskName, 0.5*maskSize[0], 0.5*maskSize[1], 0.5*maskSize[2]),
            BoronMasks::maskMaterial,
            maskPosition[0], maskPosition[1], maskPosition[2],
            lv_bank, maskColour, -2, 0, new G4RotationMatrix(0, 0, BoronMasks::getRotation(bankId, maskId)));
  }

  // Add Beamstop to the Rear Bank volume
  const int beamstopId = getParameterInt("beamstop_id");
  if (bankId == 0 && beamstopId) { // beamstopId==0 means no beamstop
    const std::string maskName = "BoronMask-Beamstop";
    const double detBankFrontDistance = banks->detectorSystemFrontDistanceFromBankFront(bankId);
    const double distanceFromDetectorFront = 5*Units::cm;

    const double width = banks->getBeamstopSize(beamstopId, 0);
    const double height = banks->getBeamstopSize(beamstopId, 1);
    const double thickness = banks->getBeamstopSize(beamstopId, 2);

    // on the beam axis (world x = y = 0), 5 cm in front of the detector front: the point of the beam axis at that
    // depth, in the bank frame (this compensates the bank elevation, also of a calibrated bank)
    const double depth = -bankHalfSize[0] + detBankFrontDistance - distanceFromDetectorFront;
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto& R = transform.rotation;
    const auto& t = transform.translation;
    const double z = (depth + R[0][0]*t[0] + R[1][0]*t[1] + R[2][0]*t[2]) / R[2][0]; // n.(P - t) = depth
    const auto local = transform.toLocal({0.0, 0.0, z});
    const G4ThreeVector position(local[0], local[1], local[2]);
    place(new G4Box(maskName, 0.5* thickness, 0.5* height, 0.5* width),
          BoronMasks::maskMaterial,
          position.x(), position.y(), position.z(),
          lv_bank, BLACK, -5, 0, new G4RotationMatrix());
  }

  if (bankId == 0 && getParameterBoolean("with_rear_calibration_mask")) {
    // the REAR CALIBRATION MASK PLATE (drawing SI-7615-993), bank frame: x = depth (front at -0.5*depth), y = across the
    // tubes, z = along the tubes; centred on the bank (as the design mask)
    const double tubeFront = BcsBanks::getTubeFrontPositionInBankFrame(bankId);
    const double gap = getParameterDouble("calibration_plate_gap_mm") * Units::mm;
    const double plateThickness = 2.0 * Units::mm, faceHeight = 1742.0 * Units::mm, faceWidth = 996.0 * Units::mm;
    const double paintThickness = getParameterDouble("calibration_paint_thickness_mm") * Units::mm;
    const double plateX = tubeFront - gap - 0.5 * plateThickness;
    place(new G4Box("CalibPlate-0", 0.5 * plateThickness, 0.5 * faceHeight, 0.5 * faceWidth),
          getParameterMaterial("calibration_plate_material"), plateX, 0., 0., lv_bank, G4Colour(0.75, 0.75, 0.75), -5, 0,
          new G4RotationMatrix());
    if (paintThickness > 0) {
      // the paint between the side walls (992 mm): 43 mm paint, then (6 mm slit, 94 mm paint) x 9, 6 mm slit, 43 mm paint
      // (slit centres at +-50, ..., +-450 mm from the centre)
      std::vector<double> pattern = {43.0};
      for (int k = 0; k < 9; ++k) { pattern.push_back(6.0); pattern.push_back(94.0); }
      pattern.push_back(6.0); pattern.push_back(43.0);
      const double paintX = plateX - 0.5 * plateThickness - 0.5 * paintThickness;
      auto paintMaterial = getParameterMaterial("calibration_paint_material");
      double offset = -0.5 * 992.0 * Units::mm;
      for (size_t i = 0; i < pattern.size(); ++i) {
        const double w = pattern[i] * Units::mm;
        if (i % 2 == 0)
          place(new G4Box("CalibPaint-0-" + std::to_string(i / 2), 0.5 * paintThickness, 0.5 * faceHeight, 0.5 * w),
                paintMaterial, paintX, 0., offset + 0.5 * w, lv_bank, DARKPURPLE, -5, 0, new G4RotationMatrix());
        offset += w;
      }
    }
    printf("GeoBCSBanks: rear calibration mask SI-7615-993: 2 mm %s sheet %g mm in front of the tubes, %g mm %s paint, "
           "slits 6 mm at +-50..+-450 mm\n", getParameterString("calibration_plate_material").c_str(), gap / Units::mm,
           paintThickness / Units::mm, getParameterString("calibration_paint_material").c_str());
  }

  // Absorber plates lying on the outer tube row (inside the bank volume; bank frame: x = depth, y = across the tubes)
  for (const auto& plate : parseNumberList(getParameterString("absorber_plates"), 9, "absorber plate", 1)) {
    if (static_cast<int>(plate[0]) != bankId || plate.size() != 10 || plate[9] <= 0)
      continue;
    const double width = plate[7], thickness = plate[8], distance = plate[9];
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto p1 = transform.toLocal({plate[1], plate[2], plate[3]});
    const auto p2 = transform.toLocal({plate[4], plate[5], plate[6]});
    const double side = (p1[1] + p2[1] > 0) ? 1.0 : -1.0;
    const double outer = side > 0 ? std::max(p1[1], p2[1]) : std::min(p1[1], p2[1]);  // the outermost of the points
    const double dx = p2[0] - p1[0], dz = p2[2] - p1[2];
    const double length = std::sqrt(dx * dx + dz * dz);
    auto rotation = new G4RotationMatrix();
    rotation->rotateY(-std::atan2(dx, dz));  // frame rotation: the inverse of the turn in the plane
    place(new G4Box("AbsorberPlate", 0.5 * width, 0.5 * thickness, 0.5 * length), getParameterMaterial("absorber_strip_material"),
          0.5 * (p1[0] + p2[0]), outer + side * distance, 0.5 * (p1[2] + p2[2]), lv_bank, BLACK, -7, 0, rotation);
    printf("GeoBCSBanks: absorber plate in bank %d, %g mm from the outer tube row: %g x %g x %g mm, turned by %.1f deg from the tube direction\n",
           bankId, distance, length, width, thickness, std::atan2(dx, dz) * 180 / M_PI);
  }

  // An enlarged notch must not cut into the content of the bank (packs, masks): it is only meant to remove empty
  // space of the bank volume. (Conservative check with the bounding boxes of the daughters in the bank frame.)
  if (banks->hasBankNotch(bankId) && !banks->getBankNotch(bankId).isNominal) {
    const auto notch = banks->getBankNotch(bankId);
    for (size_t i = 0; i < lv_bank->GetNoDaughters(); i++) {
      const auto daughter = lv_bank->GetDaughter(i);
      G4ThreeVector localMin, localMax;
      daughter->GetLogicalVolume()->GetSolid()->BoundingLimits(localMin, localMax);
      const G4RotationMatrix* frameRotation = daughter->GetRotation();
      const G4AffineTransform toBank = frameRotation ? G4AffineTransform(frameRotation, daughter->GetTranslation()) : G4AffineTransform(daughter->GetTranslation());
      G4ThreeVector low(1e99, 1e99, 1e99), high(-1e99, -1e99, -1e99);
      for (int corner = 0; corner < 8; corner++) {
        const G4ThreeVector p = toBank.TransformPoint(G4ThreeVector(corner & 1 ? localMax.x() : localMin.x(),
                                                                    corner & 2 ? localMax.y() : localMin.y(),
                                                                    corner & 4 ? localMax.z() : localMin.z()));
        for (int k = 0; k < 3; k++) {
          low[k] = std::min(low[k], p[k]);
          high[k] = std::max(high[k], p[k]);
        }
      }
      bool apart = false;
      for (int k = 0; k < 3; k++)
        if (high[k] <= notch.centre[k] - notch.halfSize[k] || low[k] >= notch.centre[k] + notch.halfSize[k])
          apart = true;
      if (!apart)
        throw std::runtime_error("GeoBCSBanks: bank_calibration=" + banks->getBankCalibration().name() + ": the notch of bank "
                                 + std::to_string(bankId) + ", enlarged because bank " + std::to_string(banks->getBankNotchNeighbour(bankId))
                                 + " overlaps it, would cut into " + std::string(daughter->GetName())
                                 + ": the two banks collide (check the calibration)");
    }
  }

  return lv_bank;
 }

G4LogicalVolume *GeoBCS::createTriangularMaskLV(int maskId){
    //creating special triangular volume by subtraction of 2 boxes
    const double xHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 0);
    const double yHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 1);
    const double zHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 2);
    const double xSideToCut = 2 * xHalf - BoronMasks::getCutPointOfTriangularMask(maskId, 0);
    const double ySideToCut = 2 * yHalf - BoronMasks::getCutPointOfTriangularMask(maskId, 1);
    const double xCutDir = BoronMasks::getCutDirOfTriangularMask(maskId, 0);
    const double yCutDir = BoronMasks::getCutDirOfTriangularMask(maskId, 1);

    auto triangularMaskBox = new G4Box("TriangularMaskBox", xHalf, yHalf, zHalf);
    const double xCutHalf = sqrt(pow(ySideToCut, 2) + pow(xSideToCut, 2)) / 2;
    const double yCutHalf = ySideToCut * xSideToCut / (2 * xCutHalf) / 2;
    const double alpha = atan(ySideToCut / xSideToCut); //rad

    auto triangularMaskCut = new G4Box("TriangularMaskSustractBox", xCutHalf, yCutHalf, zHalf * 1.1);
    G4RotationMatrix *cutRotationMatrix = new G4RotationMatrix(0, 0, -alpha * xCutDir * yCutDir);
    auto cutCentre = G4ThreeVector(xCutDir * (xHalf - xCutHalf * cos(alpha) + yCutHalf * sin(alpha)), yCutDir * (-(ySideToCut - yHalf) + xCutHalf * sin(alpha) + yCutHalf * cos(alpha)), 0);
    auto triangularMask = new G4SubtractionSolid("TriangularMask", triangularMaskBox, triangularMaskCut, cutRotationMatrix, cutCentre);

    const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId);
    const std::string maskName = "BoronMask-triangular-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    auto lv_triangularMask = new G4LogicalVolume(triangularMask, BoronMasks::maskMaterial, maskName);

    return lv_triangularMask;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* GeoBCS::Construct(){
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below
  const double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  banks = new BcsBanks(rear_detector_distance, 9, getParameterString("bank_calibration"));

  // calculate a value that is big enough to fit your world volume, the "super mother"
  double big_dimension = 1.1*( 1 *Units::m + rear_detector_distance);

  //World volume:
  auto world_material = getParameterMaterial("world_material");
  auto worldvols = place(new G4Box("World", big_dimension, big_dimension, big_dimension), world_material, 0, 0, 0, 0, INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  // Create and place detector banks
  for (int bankId = 0; bankId < banks->getNumberOfBanks(); bankId++){
    auto lv_bank = createBankLV(bankId);

    // bank placement (the same transform as used by AimHelper and PixelatedBanks)
    const BankTransform transform = banks->getBankTransform(bankId);
    // G4PVPlacement takes the frame rotation, i.e. the inverse of the bank rotation
    auto rotation = new G4RotationMatrix(toG4Rotation(transform).inverse());

    place(lv_bank, transform.translation[0], transform.translation[1], transform.translation[2], lvWorld, ORANGE, bankId, 0, rotation);
  }

  // Add 4 triangular boron masks (added to the World instead of the banks)
  for (int maskId = 0; maskId <= 3; maskId++) {
    if (!getParameterBoolean("with_triangular_masks") || !hasBoronMasks(BoronMasks::getBankIdOfTriangularMask(maskId)))
      continue;
    auto lv_triangularMask = createTriangularMaskLV(maskId);

    // mounted on the front of bank 5 or 7: its placement in the bank frame, moved with the bank transform
    const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId);
    const auto inBank = banks->getTriangularBoronMaskPlacementInBankFrame(maskId);
    const BankTransform bank = banks->getBankTransform(bankId);
    const auto position = bank.toGlobal(inBank.position);
    // G4PVPlacement takes the frame rotation: the inverse of the object rotation (bank rotation * rotation in bank)
    const G4RotationMatrix rotationInBank(G4ThreeVector(inBank.rotation[0][0], inBank.rotation[1][0], inBank.rotation[2][0]),
                                          G4ThreeVector(inBank.rotation[0][1], inBank.rotation[1][1], inBank.rotation[2][1]),
                                          G4ThreeVector(inBank.rotation[0][2], inBank.rotation[1][2], inBank.rotation[2][2]));
    auto rotation = new G4RotationMatrix((toG4Rotation(bank) * rotationInBank).inverse());

    place(lv_triangularMask, position[0], position[1], position[2], lvWorld, BLACK, -5, 0, rotation);
  }

  // Absorber strips in front of the front face of a bank volume (in the world: taped over the structural masks and
  // everything else inside the bank volume); bank frame: x = depth, y = across the tubes, z = along the tubes
  for (const auto& strip : parseAbsorberStrips(getParameterString("absorber_strips"))) {
    const int bankId = static_cast<int>(strip[0]);
    const double width = strip[7], thickness = strip[8];
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto half = banks->getBankHalfSizeInBankFrame(bankId);
    const auto p1 = transform.toLocal({strip[1], strip[2], strip[3]});
    const auto p2 = transform.toLocal({strip[4], strip[5], strip[6]});
    const double dy = p2[1] - p1[1], dz = p2[2] - p1[2];
    const double length = std::sqrt(dy * dy + dz * dz);
    const double depth = -half[0] - 0.5 * Units::mm - 0.5 * thickness;  // just in front of the bank volume
    const auto position = transform.toGlobal({depth, 0.5 * (p1[1] + p2[1]), 0.5 * (p1[2] + p2[2])});
    // object rotation: the bank rotation, then the box z axis turned in the y-z plane to the strip direction
    G4RotationMatrix inPlane;
    inPlane.rotateX(std::atan2(-dy, dz));
    auto rotation = new G4RotationMatrix((toG4Rotation(transform) * inPlane).inverse());
    place(new G4Box("AbsorberStrip", 0.5 * thickness, 0.5 * width, 0.5 * length), getParameterMaterial("absorber_strip_material"),
          position[0], position[1], position[2], lvWorld, BLACK, -6, 0, rotation);
    printf("GeoBCSBanks: absorber strip on bank %d: %g x %g x %g mm, %.1f mm in front of the tubes (in front of the bank volume)\n",
           bankId, length, width, thickness, -half[0] + banks->detectorSystemFrontDistanceFromBankFront(bankId) - depth);
  }

  // Absorber plates on a side face of a bank, outside the bank volume (in the world): a rectangle in the plane of the
  // face (bank frame: x = depth, y = across the tubes, z = along the tubes), its centre line from point 1 to point 2
  // (projected on the plane: along the tubes and in depth), turned in the plane accordingly
  for (const auto& plate : parseNumberList(getParameterString("absorber_plates"), 9, "absorber plate (bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance])", 1)) {
    if (plate.size() == 10 && plate[9] > 0)
      continue;  // inside the bank volume (createBankLV)
    const int bankId = static_cast<int>(plate[0]);
    const double width = plate[7], thickness = plate[8];
    if (width <= 0 || thickness <= 0)
      throw std::runtime_error("absorber plate with a width or thickness <= 0");
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto half = banks->getBankHalfSizeInBankFrame(bankId);
    const auto p1 = transform.toLocal({plate[1], plate[2], plate[3]});
    const auto p2 = transform.toLocal({plate[4], plate[5], plate[6]});
    const double side = (p1[1] + p2[1] > 0) ? 1.0 : -1.0;  // the face across the tubes nearest to the points
    const double dx = p2[0] - p1[0], dz = p2[2] - p1[2];
    const double length = std::sqrt(dx * dx + dz * dz);
    const std::array<double, 3> centre = {0.5 * (p1[0] + p2[0]), side * (half[1] + 0.5 * Units::mm + 0.5 * thickness), 0.5 * (p1[2] + p2[2])};
    const auto position = transform.toGlobal(centre);
    // object rotation: the bank rotation, then the box z axis turned in the face plane to the centre line
    G4RotationMatrix inPlane;
    inPlane.rotateY(std::atan2(dx, dz));
    auto rotation = new G4RotationMatrix((toG4Rotation(transform) * inPlane).inverse());
    place(new G4Box("AbsorberPlate", 0.5 * width, 0.5 * thickness, 0.5 * length),
          getParameterMaterial("absorber_strip_material"), position[0], position[1], position[2], lvWorld, BLACK, -7, 0, rotation);
    printf("GeoBCSBanks: absorber plate on bank %d (face %s of the bank volume, %.1f mm from the given points): %g x %g x %g mm, "
           "turned by %.1f deg from the tube direction, centre (%.1f, %.1f, %.1f) mm\n",
           bankId, side > 0 ? "+y" : "-y", side * half[1] - 0.5 * (p1[1] + p2[1]), length, width, thickness,
           std::atan2(dx, dz) * 180 / M_PI, position[0], position[1], position[2]);
  }

  // The source cylinder for the visualisation (world material)
  const double sourceRadius = getParameterDouble("show_source_cylinder_radius_cm") * Units::cm;
  if (sourceRadius > 0) {
    auto rotation = new G4RotationMatrix();
    rotation->rotateX(90 * Units::deg);  // the G4Tubs axis (z) -> vertical (y)
    place(new G4Tubs("SourceCylinder", 0, sourceRadius, 0.5 * getParameterDouble("show_source_cylinder_height_cm") * Units::cm, 0., 2 * M_PI),
          world_material, getParameterDouble("show_source_cylinder_x_meters") * Units::m, getParameterDouble("show_source_cylinder_y_meters") * Units::m,
          getParameterDouble("show_source_cylinder_z_meters") * Units::m, lvWorld, G4Colour(1, 1, 1), -9, 0, rotation);
  }

  // Optional vacuum tank (outside all banks: checked in validateParameters)
  const double tankWall = getParameterDouble("vacuum_tank_wall_mm")*Units::mm;
  if (tankWall > 0) {
    const double rIn = getParameterDouble("vacuum_tank_inner_radius_m")*Units::m;
    const double zFront = getParameterDouble("vacuum_tank_front_m")*Units::m;
    const double zRear = getParameterDouble("vacuum_tank_rear_m")*Units::m;
    const double rWindow = getParameterDouble("vacuum_tank_window_radius_m")*Units::m;
    const double lining = getParameterDouble("vacuum_tank_front_B4C_mm")*Units::mm;
    auto tankMaterial = getParameterMaterial("vacuum_tank_material");
    // the shell (from the outer face of the front plate to the outer face of the rear plate)
    const double zShellMin = zFront - lining - tankWall;
    const double zShellMax = zRear + tankWall;
    place(new G4Tubs("VacuumTankShell", rIn, rIn + tankWall, 0.5*(zShellMax - zShellMin), 0., 2*M_PI),
          tankMaterial, 0, 0, 0.5*(zShellMin + zShellMax), lvWorld, SILVER, -10, 0, 0);
    place(new G4Tubs("VacuumTankRear", 0, rIn, 0.5*tankWall, 0., 2*M_PI),
          tankMaterial, 0, 0, zRear + 0.5*tankWall, lvWorld, SILVER, -10, 0, 0);
    place(new G4Tubs("VacuumTankFront", rWindow, rIn, 0.5*tankWall, 0., 2*M_PI),
          tankMaterial, 0, 0, zFront - lining - 0.5*tankWall, lvWorld, SILVER, -10, 0, 0);
    if (lining > 0)
      place(new G4Tubs("VacuumTankFrontLining", rWindow, rIn, 0.5*lining, 0., 2*M_PI),
            getParameterMaterial("vacuum_tank_lining_material"), 0, 0, zFront - 0.5*lining, lvWorld, BLACK, -10, 0, 0);
    printf("GeoBCSBanks: vacuum tank: %s, wall %g mm, inner radius %g m, z %g - %g m (inner faces), window radius %g m, "
           "front lining %g mm\n", getParameterString("vacuum_tank_material").c_str(), tankWall/Units::mm, rIn/Units::m,
           zFront/Units::m, zRear/Units::m, rWindow/Units::m, lining/Units::mm);
  }

  delete banks;
  return pvWorld;
}


bool GeoBCS::isExcludedBoronMask(int bankId, int maskId) {
  return isListedBoronMask("boron_mask_exclude", bankId, maskId);
}

bool GeoBCS::isListedBoronMask(const std::string& parameter, int bankId, int maskId) {
  std::stringstream items(getParameterString(parameter));
  std::string item;
  while (std::getline(items, item, ';')) {
    const auto colon = item.find(':');
    if (colon == std::string::npos)
      throw std::runtime_error("invalid " + parameter + " item (bank:mask,mask): " + item);
    if (std::stoi(item.substr(0, colon)) != bankId)
      continue;
    std::stringstream masks(item.substr(colon + 1));
    std::string m;
    while (std::getline(masks, m, ','))
      if (std::stoi(m) == maskId)
        return true;
  }
  return false;
}

bool GeoBCS::hasBoronMasks(int bankId) {
  return parseBankList(getParameterString("boron_mask_banks")).count(bankId) > 0;
}

bool GeoBCS::validateParameters() {
  // you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc
  double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  const std::string bankCalibration = getParameterString("bank_calibration");
  try {
    BankCalibration::load(bankCalibration);
  }
  catch (const std::exception& error) {
    printf("ERROR: Wrong bank_calibration value: %s\n", error.what());
    return false;
  }
  try {
    parseBankList(getParameterString("boron_mask_banks"));
    parseAbsorberStrips(getParameterString("absorber_strips"));
    parseNumberList(getParameterString("absorber_plates"), 9, "absorber plate (bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance])", 1);
  }
  catch (const std::exception& error) {
    printf("ERROR: %s\n", error.what());
    return false;
  }
  if (getParameterDouble("vacuum_tank_wall_mm") > 0) {
    // the rear end plate behind the rear bank (its back is less than 0.3 m behind the rear detector distance)
    if (getParameterDouble("vacuum_tank_rear_m")*Units::m < rear_detector_distance + 0.3*Units::m) {
      printf("ERROR: vacuum_tank_rear_m must be at least 0.3 m behind rear_detector_distance_m\n");
      return false;
    }
    if (getParameterDouble("vacuum_tank_rear_m") + getParameterDouble("vacuum_tank_wall_mm")/1000.0 > 1.1*(1.0 + getParameterDouble("rear_detector_distance_m"))) {
      printf("ERROR: the vacuum tank does not fit in the world volume\n");
      return false;
    }
  }
  if(rear_detector_distance < 5.0 *Units::m) {
    printf("ERROR: Wrong rear_detector_distance_m value for LOKI! (It should be >=5.0 m)\n");
      return false;
  }
  return true;
}
