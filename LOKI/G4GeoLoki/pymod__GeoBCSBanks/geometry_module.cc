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

namespace {
  // an item of absorber_strips or absorber_plates: bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance[,on_pack]]
  struct AbsorberItem {
    int bankId;
    std::array<double,3> p1, p2;  // [mm] world frame, for the calibrated position of the bank
    double width, thickness;
    double distance;              // the 10th value (0 if not given); > 0: inside the bank volume
    bool onPack;                  // the 11th value (plates)
    bool insideBank() const { return distance > 0; }
  };
}

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
  std::map<int, std::array<double,3>> bankShifts;  // bank_shifts_mm
  // the end points of an absorber in the frame of its bank (a shifted bank carries its absorbers along)
  std::array<std::array<double,3>,2> pointsInBankFrame(const AbsorberItem& item) const;
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(int bankId, int packNumber);
  G4LogicalVolume * createBankLV(int bankId);
  G4LogicalVolume * createTriangularMaskLV(int maskId);
  bool hasBoronMasks(int bankId);
  // an id listed for the bank in a parameter of the format "bank:id,id;bank:..."
  bool isListed(const std::string& parameter, int bankId, int id);
  void placeRearCalibrationMask(G4LogicalVolume* lvBank);
  void placeAbsorberPlatesInBank(G4LogicalVolume* lvBank, int bankId);
  void placeAbsorberPlateOnPack(G4LogicalVolume* lvBank, const AbsorberItem& plate, const std::array<double,3>& p1,
                                const std::array<double,3>& p2, double side);
  void placeAbsorberStripsInBank(G4LogicalVolume* lvBank, int bankId);
  void placeAbsorberStripsOnBanks(G4LogicalVolume* lvWorld);
  void placeAbsorberPlatesOnBanks(G4LogicalVolume* lvWorld);
  void placeSourceCylinder(G4LogicalVolume* lvWorld, G4Material* material);
  void placeVacuumTank(G4LogicalVolume* lvWorld);
};

// this line is necessary to be able to declare the geometry in the python simulation script
SBLDPY_MODULE( _GeoBCSBanks ) { G4Interfaces::geoToPy<GeoBCS>( m ); }

namespace {
  const double absorberClearance = 0.5 * Units::mm;  // between an absorber in the world and the bank volume

  // the banks of a bank list parameter: "all", "none" (or empty), or comma separated ids and ranges, e.g. "0,5-8"
  std::set<int> parseBankList(const std::string& text) {
    const int numberOfBanks = 9;
    std::set<int> banks;
    if (text == "all") {
      for (int b = 0; b < numberOfBanks; b++)
        banks.insert(b);
      return banks;
    }
    if (text.empty() || text == "none")
      return banks;
    const auto toInt = [](const std::string& number) {
      size_t used = 0;
      const int value = std::stoi(number, &used);
      if (used != number.size())
        throw std::invalid_argument(number);
      return value;
    };
    std::stringstream items(text);
    std::string item;
    while (std::getline(items, item, ',')) {
      const auto dash = item.find('-');
      try {
        const int first = toInt(item.substr(0, dash));
        const int last = dash == std::string::npos ? first : toInt(item.substr(dash + 1));
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

  // the bank face across the tubes nearest to the two points of an absorber (bank frame): +1 or -1
  double sideOfBank(const std::array<double,3>& p1, const std::array<double,3>& p2) {
    return (p1[1] + p2[1] > 0) ? 1.0 : -1.0;
  }

  // the pack of a bank nearest to its face across the tubes on this side
  int outermostPack(int bankId, double side) {
    int outermost = 0;
    for (int packNumber = 1; packNumber < BcsBanks::getNumberOfPacksByBankId(bankId); ++packNumber)
      if (side * BcsBanks::getPackPositionInBankFrame(bankId, packNumber)[1] >
          side * BcsBanks::getPackPositionInBankFrame(bankId, outermost)[1])
        outermost = packNumber;
    return outermost;
  }

  // a list of absorber items: ";"-separated items of 9 to maxValues comma-separated numbers (see AbsorberItem)
  std::vector<AbsorberItem> parseAbsorbers(const std::string& text, size_t maxValues, const std::string& what) {
    std::vector<AbsorberItem> absorbers;
    std::stringstream items(text);
    std::string item;
    while (std::getline(items, item, ';')) {
      if (item.find_first_not_of(" ") == std::string::npos)
        continue;
      std::vector<double> values;
      std::stringstream numbers(item);
      std::string number;
      while (std::getline(numbers, number, ','))
        values.push_back(std::stod(number));
      if (values.size() < 9 || values.size() > maxValues || values[0] < 0 || values[0] > 8 || values[7] <= 0 || values[8] <= 0)
        throw std::runtime_error("invalid " + what + ": " + item);
      absorbers.push_back({static_cast<int>(values[0]), {values[1], values[2], values[3]}, {values[4], values[5], values[6]},
                           values[7], values[8], values.size() > 9 ? values[9] : 0.0, values.size() == 11 && values[10] > 0});
    }
    return absorbers;
  }

  std::vector<AbsorberItem> parseAbsorberStrips(const std::string& text) {
    return parseAbsorbers(text, 10, "absorber strip (bank,x1,y1,z1,x2,y2,z2,width,thickness[,inside])");
  }

  std::vector<AbsorberItem> parseAbsorberPlates(const std::string& text) {
    return parseAbsorbers(text, 11, "absorber plate (bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance[,on_pack]])");
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
  // the calibration mask of the rear bank (the manufactured plate, drawing SI-7615-993, see CalibMasks::RearPlate)
  addParameterBoolean("with_rear_calibration_mask", false);
  addParameterDouble("calibration_plate_gap_mm", 8.0, 0.0, 50.0);         // back of the 2 mm sheet - front of the tubes
  addParameterDouble("calibration_paint_thickness_mm", 0.15, 0.0, 5.0);   // the Gd2O3 paint layer (not specified in the drawing)
  // 50% polyurethane varnish (approximated as C25H42N2O6, 1.2 g/cm3) + 50% Gd2O3 (7.41 g/cm3) by weight
  addParameterString("calibration_paint_material", "solid::Gd10C100H167N8O39/2.07gcm3");
  addParameterString("calibration_plate_material", "G4_Al");
  // the B4C masks of the banks (the rectangular masks of each bank and the triangular masks of banks 5 and 7; not the
  // beamstop): "all" (the design), "none", or a bank list as above
  addParameterString("boron_mask_banks", "all");
  // the triangular B4C masks of banks 5 and 7 (part of their design masks; switch them off separately)
  addParameterBoolean("with_triangular_masks", true);
  // individual design masks not present: "bank:mask,mask;bank:..." (the mask ids: the indices of the mask tables of
  // BoronMasks.cc; the volume names are BoronMask-<bank>-<CAD drawing>[.<k>], see BoronMasks::getName)
  addParameterString("boron_mask_exclude", "");
  // detector packs not on the detector: "bank:pack,pack;bank:..." (pack numbers in the bank, 0 = the pack of tubes 0-1 of
  // each layer); their pixels get no hits
  addParameterString("absent_packs", "");
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
  // absorber strips on the front of a bank (e.g. the B4C sheets taped on the banks in source tests):
  // "bank,x1,y1,z1,x2,y2,z2,width,thickness[,inside];..." [mm], the strip from point 1 to point 2 (world coordinates):
  // - no inside value (or <= 0): just in front of the bank volume, the points projected on the plane of the strip
  // - inside > 0: inside the bank volume, from point 1 to point 2 as given (not projected: e.g. a strip running from the
  //   front of a mask down to a deeper mask), its width parallel to the bank front
  addParameterString("absorber_strips", "");
  addParameterString("absorber_strip_material", "MAT_B4C");
  // absorber plates on the side of a bank (the same kind of item as absorber_strips, in the strip material):
  // "bank,x1,y1,z1,x2,y2,z2,width,thickness[,distance[,on_pack]];..." [mm], the centre line from point 1 to point 2
  // (world coordinates), the width in the plane of the plate, on the side across the tubes nearest to the points (e.g.
  // the inner face of a front bank towards a source between the banks):
  // - no distance (or <= 0): just outside the bank volume, the points projected on that plane (on_pack is then ignored)
  // - distance > 0: inside the bank volume, the centre this far beyond the outermost point across the tubes (e.g. a
  //   sheet lying directly on the tubes of the outer row)
  // - distance > 0 and on_pack 1: on the outer face of the outermost pack on that side, parallel to the (tilted) packs,
  //   the points taken in the frame of that pack, the plate this far from the face
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
  // banks moved from their calibrated positions (e.g. on their rails for a measurement): "bank:dx,dy,dz;..." [mm, world
  // frame], e.g. "2:0,0,700;4:0,0,700"; the absorber strips and plates of a shifted bank (given for its calibrated
  // position) move with it; the analysis gets the shifted placements with the recorded calibration text. The z of the
  // rear bank (0) is rear_detector_distance_m, not shifted here (a dz of bank 0 only moves its absorbers)
  addParameterString("bank_shifts_mm", "");

  addParameterString("world_material","G4_Vacuum");
  // the gas around the tubes inside the bank and pack volumes (the measurements at atmospheric pressure: G4_AIR)
  addParameterString("pack_fill_material","G4_Vacuum");
  addParameterString("B4C_panel_material","MAT_B4C:b10_enrichment=0.95"); // the B4C panels behind the tubes
}

G4LogicalVolume * GeoBCS::createTubeLV(double converterThickness, double strawLength){
  const double effectiveStrawLength = strawLength - BcsTube::getStrawWallThickness(); //This is only epsilon difference...

  auto lv_tube = new G4LogicalVolume(new G4Tubs(BcsTube::tubeWallName, 0, BcsTube::getTubeOuterRadius(), 0.5*strawLength, 0., 2*M_PI),
                                     BcsTube::tubeWallMaterial, BcsTube::tubeWallName);

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
  for (int partId = 0; partId < BcsPack::getNumberOfB4CParts(); partId++){
    place(new G4Box(BcsPack::getB4CPartName(partId), 0.5*BcsPack::getB4CPartThickness(partId), 0.5*BcsPack::getB4CPartHeight(partId), B4CLengthHalf),
          B4CPanelMaterial,
          BcsPack::getB4CPartHorizontalOffset(partId), BcsPack::getB4CPartVerticalOffset(partId), 0,
          lv_pack_box, G4Colour(0, 1, 0), -2, 0, new G4RotationMatrix());
    }
  /// Add Al behing the B4C panel in 2 parts ///
  for (int partId = 0; partId < BcsPack::getNumberOfAlParts(); partId++){
    place(new G4Box(BcsPack::getAlPartName(partId), 0.5*BcsPack::getAlPartThickness(partId), 0.5*BcsPack::getAlPartHeight(partId), B4CLengthHalf),
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
    if (isListed("absent_packs", bankId, packNumber))
      continue;  // a pack not on the detector (e.g. removed for repair): no tubes, no material
    auto lv_pack_box = createPackBoxLV(bankId, packNumber);
    const auto packPosition = banks->getPackPositionInBankFrame(bankId, packNumber);
    place(lv_pack_box, packPosition[0], packPosition[1], packPosition[2],
          lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  const int numberOfBoronMasks = hasBoronMasks(bankId) ? BoronMasks::getNumberOfBoronMasks(bankId) : 0;
  for (int maskId = 0; maskId < numberOfBoronMasks; ++maskId){
    if (isListed("boron_mask_exclude", bankId, maskId))
      continue;
    const std::string maskName = BoronMasks::getName(bankId, maskId);
    const auto maskSize = BoronMasks::getSizeInBankFrame(bankId, maskId);
    const auto maskPosition = banks->getBoronMaskPositionInBankFrame(bankId, maskId);
    const G4Colour maskColour = isListed("highlight_boron_masks", bankId, maskId) ? G4Colour(1, 0, 0)
                              : isListed("highlight_boron_masks_2", bankId, maskId) ? G4Colour(1, 0.5, 0) : BLACK;
    place(new G4Box(maskName, 0.5*maskSize[0], 0.5*maskSize[1], 0.5*maskSize[2]),
            BoronMasks::maskMaterial,
            maskPosition[0], maskPosition[1], maskPosition[2],
            lv_bank, maskColour, -2, 0, new G4RotationMatrix(0, 0, BoronMasks::getRotation(bankId, maskId)));
  }

  // Add Beamstop to the Rear Bank volume
  const int beamstopId = getParameterInt("beamstop_id");
  if (bankId == 0 && beamstopId) { // beamstopId==0 means no beamstop
    const std::string maskName = BoronMasks::beamstopName;
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

  if (bankId == 0 && getParameterBoolean("with_rear_calibration_mask"))
    placeRearCalibrationMask(lv_bank);
  placeAbsorberPlatesInBank(lv_bank, bankId);
  placeAbsorberStripsInBank(lv_bank, bankId);

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

// The REAR CALIBRATION MASK PLATE (drawing SI-7615-993), bank frame: x = depth (front at -0.5*depth), y = across the
// tubes, z = along the tubes; centred on the bank
void GeoBCS::placeRearCalibrationMask(G4LogicalVolume* lvBank) {
  const double tubeFront = BcsBanks::getTubeFrontPositionInBankFrame(0);
  const double gap = getParameterDouble("calibration_plate_gap_mm") * Units::mm;
  using Plate = CalibMasks::RearPlate;
  const double plateThickness = Plate::thickness * Units::mm;
  const double faceHeight = Plate::faceHeight * Units::mm;
  const double faceWidth = Plate::faceWidth * Units::mm;
  const double paintThickness = getParameterDouble("calibration_paint_thickness_mm") * Units::mm;
  const double plateX = tubeFront - gap - 0.5 * plateThickness;
  place(new G4Box(Plate::getPlateName(), 0.5 * plateThickness, 0.5 * faceHeight, 0.5 * faceWidth),
        getParameterMaterial("calibration_plate_material"), plateX, 0., 0., lvBank, G4Colour(0.75, 0.75, 0.75), -5, 0,
        new G4RotationMatrix());
  if (paintThickness > 0) {
    // the paint between the side walls: paint strips and slits (Plate::getPaintPattern)
    const std::vector<double> pattern = Plate::getPaintPattern();
    const double paintX = plateX - 0.5 * plateThickness - 0.5 * paintThickness;
    auto paintMaterial = getParameterMaterial("calibration_paint_material");
    double offset = -0.5 * Plate::paintedWidth * Units::mm;
    for (size_t i = 0; i < pattern.size(); ++i) {
      const double w = pattern[i] * Units::mm;
      if (i % 2 == 0)
        place(new G4Box(Plate::getPaintName(static_cast<int>(i / 2)), 0.5 * paintThickness, 0.5 * faceHeight, 0.5 * w),
              paintMaterial, paintX, 0., offset + 0.5 * w, lvBank, DARKPURPLE, -5, 0, new G4RotationMatrix());
      offset += w;
    }
  }
  printf("GeoBCSBanks: rear calibration mask SI-7615-993: 2 mm %s sheet %g mm in front of the tubes, %g mm %s paint, "
         "slits 6 mm at +-50..+-450 mm\n", getParameterString("calibration_plate_material").c_str(), gap / Units::mm,
         paintThickness / Units::mm, getParameterString("calibration_paint_material").c_str());
}

// Absorber plates inside the bank volume, see absorber_plates (bank frame: x = depth, y = across the tubes)
void GeoBCS::placeAbsorberPlatesInBank(G4LogicalVolume* lvBank, int bankId) {
  for (const auto& plate : parseAbsorberPlates(getParameterString("absorber_plates"))) {
    if (plate.bankId != bankId || !plate.insideBank())
      continue;
    const auto [p1, p2] = pointsInBankFrame(plate);
    const double side = sideOfBank(p1, p2);
    if (plate.onPack) {
      placeAbsorberPlateOnPack(lvBank, plate, p1, p2, side);
      continue;
    }
    const double outer = side > 0 ? std::max(p1[1], p2[1]) : std::min(p1[1], p2[1]);  // the outermost of the points
    const double dx = p2[0] - p1[0], dz = p2[2] - p1[2];
    const double length = std::sqrt(dx * dx + dz * dz);
    auto rotation = new G4RotationMatrix();
    rotation->rotateY(-std::atan2(dx, dz));  // frame rotation: the inverse of the turn in the plane
    place(new G4Box("AbsorberPlate", 0.5 * plate.width, 0.5 * plate.thickness, 0.5 * length), getParameterMaterial("absorber_strip_material"),
          0.5 * (p1[0] + p2[0]), outer + side * plate.distance, 0.5 * (p1[2] + p2[2]), lvBank, BLACK, -7, 0, rotation);
    printf("GeoBCSBanks: absorber plate in bank %d, %g mm from the given points: %g x %g x %g mm, turned by %.1f deg from the tube direction\n",
           bankId, plate.distance, length, plate.width, plate.thickness, std::atan2(dx, dz) * 180 / M_PI);
  }
}

// A plate lying on the outer face of the outermost pack on this side, parallel to the (tilted) packs, this far from the
// face (the points are taken in the frame of that pack: x = depth, z = along the tubes)
void GeoBCS::placeAbsorberPlateOnPack(G4LogicalVolume* lvBank, const AbsorberItem& plate, const std::array<double,3>& p1,
                                      const std::array<double,3>& p2, double side) {
  const int bankId = plate.bankId;
  const int outermost = outermostPack(bankId, side);
  const auto pp = banks->getPackPositionInBankFrame(bankId, outermost);
  const G4ThreeVector packCentre(pp[0], pp[1], pp[2]);
  // packs are placed with this frame rotation (as in createBankLV): packFrame * (bank point - pack centre) = pack point
  const G4RotationMatrix packFrame(0, 0, banks->getPackRotation());
  const G4RotationMatrix packObject = packFrame.inverse();
  // the plate lies at a fixed height over the pack face: only x (depth) and z (along the tubes) of the points are used
  const G4ThreeVector q1 = packFrame * (G4ThreeVector(p1[0], p1[1], p1[2]) - packCentre);
  const G4ThreeVector q2 = packFrame * (G4ThreeVector(p2[0], p2[1], p2[2]) - packCentre);
  const double qx = q2.x() - q1.x(), qz = q2.z() - q1.z();
  const double length = std::sqrt(qx * qx + qz * qz);
  G4RotationMatrix inPlane;
  inPlane.rotateY(std::atan2(qx, qz));
  const G4ThreeVector centre = packCentre + packObject * G4ThreeVector(0.5 * (q1.x() + q2.x()),
      side * (0.5 * BcsPack::getPackBoxHeight() + plate.distance + 0.5 * plate.thickness), 0.5 * (q1.z() + q2.z()));
  place(new G4Box("AbsorberPlate", 0.5 * plate.width, 0.5 * plate.thickness, 0.5 * length), getParameterMaterial("absorber_strip_material"),
        centre.x(), centre.y(), centre.z(), lvBank, BLACK, -7, 0, new G4RotationMatrix((packObject * inPlane).inverse()));
  printf("GeoBCSBanks: absorber plate in bank %d on pack %d (parallel to the packs, %g mm from its face): %g x %g x %g mm, turned by %.1f deg from the tube direction\n",
         bankId, outermost, plate.distance, length, plate.width, plate.thickness, std::atan2(qx, qz) * 180 / M_PI);
}

// Absorber strips inside the bank volume (bank frame: x = depth): a box from the point 1 to the point 2 as given,
// its width across the centre line and parallel to the bank front, its thickness about along the depth
void GeoBCS::placeAbsorberStripsInBank(G4LogicalVolume* lvBank, int bankId) {
  for (const auto& strip : parseAbsorberStrips(getParameterString("absorber_strips"))) {
    if (strip.bankId != bankId || !strip.insideBank())
      continue;
    const auto [p1, p2] = pointsInBankFrame(strip);
    const G4ThreeVector a(p1[0], p1[1], p1[2]), b(p2[0], p2[1], p2[2]);
    const G4ThreeVector along = (b - a).unit();
    const G4ThreeVector across = G4ThreeVector(1, 0, 0).cross(along).unit();
    const G4ThreeVector normal = across.cross(along);
    const G4RotationMatrix objectRotation(normal, across, along);  // columns: the box axes in the bank frame (right-handed)
    const G4ThreeVector centre = 0.5 * (a + b);
    place(new G4Box("AbsorberStrip", 0.5 * strip.thickness, 0.5 * strip.width, 0.5 * (b - a).mag()), getParameterMaterial("absorber_strip_material"),
          centre.x(), centre.y(), centre.z(), lvBank, BLACK, -6, 0, new G4RotationMatrix(objectRotation.inverse()));
    printf("GeoBCSBanks: absorber strip in bank %d: %g x %g x %g mm, depth %.1f to %.1f mm in the bank frame (tilted by %.1f deg from the bank front)\n",
           bankId, (b - a).mag(), strip.width, strip.thickness, a.x(), b.x(), std::asin(std::fabs(along.x())) * 180 / M_PI);
  }
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

    const std::string maskName = BoronMasks::getNameOfTriangularMask(maskId);
    auto lv_triangularMask = new G4LogicalVolume(triangularMask, BoronMasks::maskMaterial, maskName);

    return lv_triangularMask;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* GeoBCS::Construct(){
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below
  const double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  bankShifts = BankCalibration::parseShifts(getParameterString("bank_shifts_mm"));
  banks = new BcsBanks(rear_detector_distance, 9,
                       BankCalibration::load(getParameterString("bank_calibration")).shifted(getParameterString("bank_shifts_mm")));
  for (const auto& [bankId, shift] : bankShifts)
    printf("GeoBCSBanks: bank %d shifted from its bank_calibration position by (%g, %g, %g) mm\n", bankId,
           shift[0] / Units::mm, shift[1] / Units::mm, shift[2] / Units::mm);

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

  placeAbsorberStripsOnBanks(lvWorld);
  placeAbsorberPlatesOnBanks(lvWorld);
  placeSourceCylinder(lvWorld, world_material);
  placeVacuumTank(lvWorld);

  delete banks;
  return pvWorld;
}

// Absorber strips in front of the front face of a bank volume (in the world: taped over the structural masks and
// everything else inside the bank volume); bank frame: x = depth, y = across the tubes, z = along the tubes
void GeoBCS::placeAbsorberStripsOnBanks(G4LogicalVolume* lvWorld) {
  for (const auto& strip : parseAbsorberStrips(getParameterString("absorber_strips"))) {
    if (strip.insideBank())
      continue;  // placed by placeAbsorberStripsInBank
    const int bankId = strip.bankId;
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto half = banks->getBankHalfSizeInBankFrame(bankId);
    const auto [p1, p2] = pointsInBankFrame(strip);
    const double dy = p2[1] - p1[1], dz = p2[2] - p1[2];
    const double length = std::sqrt(dy * dy + dz * dz);
    const double depth = -half[0] - absorberClearance - 0.5 * strip.thickness;  // just in front of the bank volume
    const auto position = transform.toGlobal({depth, 0.5 * (p1[1] + p2[1]), 0.5 * (p1[2] + p2[2])});
    // object rotation: the bank rotation, then the box z axis turned in the y-z plane to the strip direction
    G4RotationMatrix inPlane;
    inPlane.rotateX(std::atan2(-dy, dz));
    auto rotation = new G4RotationMatrix((toG4Rotation(transform) * inPlane).inverse());
    place(new G4Box("AbsorberStrip", 0.5 * strip.thickness, 0.5 * strip.width, 0.5 * length), getParameterMaterial("absorber_strip_material"),
          position[0], position[1], position[2], lvWorld, BLACK, -6, 0, rotation);
    printf("GeoBCSBanks: absorber strip on bank %d: %g x %g x %g mm, %.1f mm in front of the tubes (in front of the bank volume)\n",
           bankId, length, strip.width, strip.thickness, BcsBanks::getTubeFrontPositionInBankFrame(bankId) - depth);
  }
}

// Absorber plates on a side face of a bank, outside the bank volume (in the world): a rectangle in the plane of the
// face (bank frame: x = depth, y = across the tubes, z = along the tubes), its centre line from point 1 to point 2
// (projected on the plane: along the tubes and in depth), turned in the plane accordingly
void GeoBCS::placeAbsorberPlatesOnBanks(G4LogicalVolume* lvWorld) {
  for (const auto& plate : parseAbsorberPlates(getParameterString("absorber_plates"))) {
    if (plate.insideBank())
      continue;  // placed by placeAbsorberPlatesInBank
    const int bankId = plate.bankId;
    const BankTransform transform = banks->getBankTransform(bankId);
    const auto half = banks->getBankHalfSizeInBankFrame(bankId);
    const auto [p1, p2] = pointsInBankFrame(plate);
    const double side = sideOfBank(p1, p2);
    const double dx = p2[0] - p1[0], dz = p2[2] - p1[2];
    const double length = std::sqrt(dx * dx + dz * dz);
    const std::array<double, 3> centre = {0.5 * (p1[0] + p2[0]), side * (half[1] + absorberClearance + 0.5 * plate.thickness), 0.5 * (p1[2] + p2[2])};
    const auto position = transform.toGlobal(centre);
    // object rotation: the bank rotation, then the box z axis turned in the face plane to the centre line
    G4RotationMatrix inPlane;
    inPlane.rotateY(std::atan2(dx, dz));
    auto rotation = new G4RotationMatrix((toG4Rotation(transform) * inPlane).inverse());
    place(new G4Box("AbsorberPlate", 0.5 * plate.width, 0.5 * plate.thickness, 0.5 * length),
          getParameterMaterial("absorber_strip_material"), position[0], position[1], position[2], lvWorld, BLACK, -7, 0, rotation);
    printf("GeoBCSBanks: absorber plate on bank %d (face %s of the bank volume, %.1f mm from the given points): %g x %g x %g mm, "
           "turned by %.1f deg from the tube direction, centre (%.1f, %.1f, %.1f) mm\n",
           bankId, side > 0 ? "+y" : "-y", side * half[1] - 0.5 * (p1[1] + p2[1]), length, plate.width, plate.thickness,
           std::atan2(dx, dz) * 180 / M_PI, position[0], position[1], position[2]);
  }
}

// The source cylinder for the visualisation (world material)
void GeoBCS::placeSourceCylinder(G4LogicalVolume* lvWorld, G4Material* material) {
  const double sourceRadius = getParameterDouble("show_source_cylinder_radius_cm") * Units::cm;
  if (sourceRadius > 0) {
    auto rotation = new G4RotationMatrix();
    rotation->rotateX(90 * Units::deg);  // the G4Tubs axis (z) -> vertical (y)
    place(new G4Tubs("SourceCylinder", 0, sourceRadius, 0.5 * getParameterDouble("show_source_cylinder_height_cm") * Units::cm, 0., 2 * M_PI),
          material, getParameterDouble("show_source_cylinder_x_meters") * Units::m, getParameterDouble("show_source_cylinder_y_meters") * Units::m,
          getParameterDouble("show_source_cylinder_z_meters") * Units::m, lvWorld, G4Colour(1, 1, 1), -9, 0, rotation);
  }
}

// Optional vacuum tank (its rear end plate is checked in validateParameters, the radius only by the overlap check)
void GeoBCS::placeVacuumTank(G4LogicalVolume* lvWorld) {
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
}


bool GeoBCS::isListed(const std::string& parameter, int bankId, int id) {
  std::stringstream items(getParameterString(parameter));
  std::string item;
  while (std::getline(items, item, ';')) {
    const auto colon = item.find(':');
    if (colon == std::string::npos)
      throw std::runtime_error("invalid " + parameter + " item (bank:mask,mask): " + item);
    if (std::stoi(item.substr(0, colon)) != bankId)
      continue;
    std::stringstream ids(item.substr(colon + 1));
    std::string listed;
    while (std::getline(ids, listed, ','))
      if (std::stoi(listed) == id)
        return true;
  }
  return false;
}

std::array<std::array<double,3>,2> GeoBCS::pointsInBankFrame(const AbsorberItem& item) const {
  std::array<double,3> p1 = item.p1, p2 = item.p2;
  const auto shift = bankShifts.find(item.bankId);
  if (shift != bankShifts.end())
    for (int i = 0; i < 3; i++) {
      p1[i] += shift->second[i];
      p2[i] += shift->second[i];
    }
  const BankTransform transform = banks->getBankTransform(item.bankId);
  return {{transform.toLocal(p1), transform.toLocal(p2)}};
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
    BankCalibration::load(bankCalibration).shifted(getParameterString("bank_shifts_mm"));
  }
  catch (const std::exception& error) {
    printf("ERROR: Wrong bank_shifts_mm value: %s\n", error.what());
    return false;
  }
  try {
    parseBankList(getParameterString("boron_mask_banks"));
    parseAbsorberStrips(getParameterString("absorber_strips"));
    parseAbsorberPlates(getParameterString("absorber_plates"));
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
