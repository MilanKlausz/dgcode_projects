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
#include <array>
#include <stdexcept>
#include <cmath>
#include <string>
#include <cassert>
#include <set>
#include <vector>
#include <sstream>

#include "G4GeoLoki/BcsBanks.hh"
#include "G4GeoLoki/CadComponents.hh"
#include "G4Materials/NamedMaterialProvider.hh"
#include <map>

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
  std::vector<CadComponents::Box> cadBoxes; // the boxes of the CAD component sets (cad_components)
  std::vector<CadComponents::PackBox> cadPackBoxes; // the boxes of the CAD pack component sets (cad_pack_components)
  std::vector<CadComponents::Envelope> cadEnvelopes; // the envelopes of the CAD component sets (simplified sets)
  std::map<std::string, G4Material*> cadMaterials;  // the materials of the CAD components by (key, density factor)
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(int bankId, int packNumber);
  G4LogicalVolume * createBankLV(int bankId);
  G4LogicalVolume * createTriangularMaskLV(int maskId);
  bool hasBoronMasks(int bankId);
  bool isExcludedBoronMask(int bankId, int maskId);
  bool isListedBoronMask(const std::string& parameter, int bankId, int maskId);
  G4LogicalVolume * createCadBoxLV(const CadComponents::Box& box);
  G4Material * getCadMaterial(const std::string& key, double fill);
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

  // display colour of the CAD components by material: B4C black, Tufnol orange, aluminium silver
  G4Colour cadColour(const std::string& material) {
    if (material == "B4C")
      return G4Colour(0.1, 0.1, 0.1);
    if (material == "Tufnol")
      return G4Colour(1.0, 0.55, 0.0);
    return G4Colour(0.75, 0.75, 0.75);
  }

  // the distance of the origin from the convex hull of points in a plane (0: inside the hull)
  double distanceFromConvexHull(std::vector<std::array<double,2>> p) {
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    if (p.size() == 1)
      return std::hypot(p[0][0], p[0][1]);
    auto cross = [](const std::array<double,2>& o, const std::array<double,2>& a, const std::array<double,2>& b) {
      return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    std::vector<std::array<double,2>> hull(2 * p.size());   // monotone chain: counter-clockwise
    size_t k = 0;
    for (size_t i = 0; i < p.size(); i++) {
      while (k >= 2 && cross(hull[k - 2], hull[k - 1], p[i]) <= 0)
        k--;
      hull[k++] = p[i];
    }
    for (size_t i = p.size() - 1, lower = k + 1; i > 0; i--) {
      while (k >= lower && cross(hull[k - 2], hull[k - 1], p[i - 1]) <= 0)
        k--;
      hull[k++] = p[i - 1];
    }
    hull.resize(k - 1);
    bool inside = hull.size() >= 3;
    double distance = 1e99;
    const std::array<double,2> origin = {0.0, 0.0};
    for (size_t i = 0; i < hull.size(); i++) {
      const auto& a = hull[i];
      const auto& b = hull[(i + 1) % hull.size()];
      if (cross(a, b, origin) < 0)
        inside = false;
      const double dx = b[0] - a[0], dy = b[1] - a[1], length2 = dx * dx + dy * dy;
      const double t = length2 > 0 ? std::clamp(-(a[0] * dx + a[1] * dy) / length2, 0.0, 1.0) : 0.0;
      distance = std::min(distance, std::hypot(a[0] + t * dx, a[1] + t * dy));
    }
    return inside ? 0.0 : distance;
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
  // individual design masks not present: "bank:mask,mask;bank:..." (the mask ids: the indices of the mask tables of
  // BoronMasks.cc; the volume names are BoronMask-<bank>-<CAD drawing>[.<k>], see BoronMasks::getName)
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

  // Optional components from the CAD models, simplified to boxes (see G4GeoLoki/CadComponents.hh): a comma separated
  // list of component sets (data files G4GeoLoki/data/cad_components_<set>.txt, or file paths), e.g.
  // "front_b4c_2021,tufnol"; empty: none. The materials of the components (by their material key):
  addParameterString("cad_components", "");
  // the same for components of every detector pack (data files G4GeoLoki/data/cad_pack_<set>.txt), e.g. "rail"
  addParameterString("cad_pack_components", "");
  addParameterString("cad_B4C_material", "MAT_B4C:b10_enrichment=0.95"); // B4C S-DOUGH (82% B4C) parts, as the masks
  addParameterString("cad_Tufnol_material", "G4_BAKELITE");              // Tufnol (phenolic laminate)
  addParameterString("cad_Al_material", "NCrystal:cfg=Al_sg225.ncmat");  // aluminium alloy parts

  // Optional mother volumes of the three bank groups (navigation speed): front (banks 5-8), middle (banks 1-4) and rear
  // (bank 0). Each is a cylinder (a ring where the content leaves the beam axis free) of the world material around the
  // beam axis, sized to its content with a 1 mm gap: the banks and everything placed in the world for them (triangular
  // masks, CAD envelopes and world boxes, absorber strips and plates). The source cylinder goes into the mother that
  // contains its centre. The copy-number depths below the banks are unchanged.
  addParameterBoolean("with_bank_mothers", false);
  addParameterBoolean("show_bank_mothers", false);  // draw the bank mothers (transparent) in the viewers
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
  /// CAD pack components (cad_pack_components): boxes along the tubes or at both ends ///
  for (const auto& box : cadPackBoxes) {
    const std::string name = "CadPack-" + box.material + "-" + box.name;
    auto material = getCadMaterial(box.material, box.fill);
    const double hx = 0.5 * (box.x1 - box.x0), hy = 0.5 * (box.y1 - box.y0);
    const double cx = 0.5 * (box.x0 + box.x1), cy = 0.5 * (box.y0 + box.y1);
    if (box.alongLength) {
      place(new G4Box(name, hx, hy, 0.5 * strawLength + box.a), material, cx, cy, 0, lv_pack_box, cadColour(box.material), -12, 0, new G4RotationMatrix());
    } else {
      const double hz = 0.5 * (box.b - box.a), cz = 0.5 * strawLength + 0.5 * (box.a + box.b);
      for (int side : {1, -1})
        place(new G4Box(name, hx, hy, hz), material, cx, cy, side * cz, lv_pack_box, cadColour(box.material), -12, 0, new G4RotationMatrix());
    }
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
    const std::string maskName = BoronMasks::getName(bankId, maskId);
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

  if (bankId == 0 && getParameterBoolean("with_rear_calibration_mask")) {
    // the REAR CALIBRATION MASK PLATE (drawing SI-7615-993), bank frame: x = depth (front at -0.5*depth), y = across the
    // tubes, z = along the tubes; centred on the bank (as the design mask)
    const double tubeFront = BcsBanks::getTubeFrontPositionInBankFrame(bankId);
    const double gap = getParameterDouble("calibration_plate_gap_mm") * Units::mm;
    using Plate = CalibMasks::RearPlate;
    const double plateThickness = Plate::thickness * Units::mm, faceHeight = Plate::faceHeight * Units::mm;
    const double faceWidth = Plate::faceWidth * Units::mm;
    const double paintThickness = getParameterDouble("calibration_paint_thickness_mm") * Units::mm;
    const double plateX = tubeFront - gap - 0.5 * plateThickness;
    place(new G4Box(Plate::getPlateName(), 0.5 * plateThickness, 0.5 * faceHeight, 0.5 * faceWidth),
          getParameterMaterial("calibration_plate_material"), plateX, 0., 0., lv_bank, G4Colour(0.75, 0.75, 0.75), -5, 0,
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

  // CAD components inside this bank volume (bank frame)
  for (const auto& box : cadBoxes) {
    if (!box.inBank || box.bankId != bankId)
      continue;
    const auto& a = box.axes;
    const G4RotationMatrix inBank(G4ThreeVector(a[0][0], a[0][1], a[0][2]), G4ThreeVector(a[1][0], a[1][1], a[1][2]),
                                  G4ThreeVector(a[2][0], a[2][1], a[2][2]));   // columns: the box axes
    place(createCadBoxLV(box), box.centre[0], box.centre[1], box.centre[2], lv_bank,
          cadColour(box.material), -11, 0, new G4RotationMatrix(inBank.inverse()));
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

    const std::string maskName = BoronMasks::getNameOfTriangularMask(maskId);
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

  cadBoxes = CadComponents::load(getParameterString("cad_components"));
  cadEnvelopes = CadComponents::loadEnvelopes(getParameterString("cad_components"));
  cadPackBoxes = CadComponents::loadPack(getParameterString("cad_pack_components"));

  // The volumes of the world are collected first (with the bank they belong to, -1: none) and placed at the end: in the
  // world, or in the mother volume of their bank group (with_bank_mothers)
  struct WorldItem {
    G4LogicalVolume* lv;
    G4ThreeVector position;           // in the world
    G4RotationMatrix* rotation;       // the frame rotation (G4PVPlacement)
    G4Colour colour;
    int copyNumber;
    int bankId;
  };
  std::vector<WorldItem> worldItems;

  // Create and place detector banks
  for (int bankId = 0; bankId < banks->getNumberOfBanks(); bankId++){
    auto lv_bank = createBankLV(bankId);

    // bank placement (the same transform as used by AimHelper and PixelatedBanks)
    const BankTransform transform = banks->getBankTransform(bankId);
    // G4PVPlacement takes the frame rotation, i.e. the inverse of the bank rotation
    auto rotation = new G4RotationMatrix(toG4Rotation(transform).inverse());

    worldItems.push_back({lv_bank, G4ThreeVector(transform.translation[0], transform.translation[1], transform.translation[2]),
                          rotation, ORANGE, bankId, bankId});
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

    worldItems.push_back({lv_triangularMask, G4ThreeVector(position[0], position[1], position[2]), rotation, BLACK, -5, bankId});
  }

  // CAD components outside the bank volumes: given in the frame of their bank, placed with the bank transform; in the
  // world or in an envelope (a box of the world material, itself placed in the world with the bank transform)
  std::map<std::string, std::pair<G4LogicalVolume*, const CadComponents::Envelope*>> envelopeLVs;
  for (const auto& env : cadEnvelopes) {
    const BankTransform bank = banks->getBankTransform(env.bankId);
    const auto position = bank.toGlobal(env.centre);
    const auto& a = env.axes;
    const G4RotationMatrix inBank(G4ThreeVector(a[0][0], a[0][1], a[0][2]), G4ThreeVector(a[1][0], a[1][1], a[1][2]),
                                  G4ThreeVector(a[2][0], a[2][1], a[2][2]));   // columns: the envelope axes
    auto rotation = new G4RotationMatrix((toG4Rotation(bank) * inBank).inverse());
    const std::string name = "CadEnvelope-" + env.name;
    auto lv = new G4LogicalVolume(new G4Box(name, env.halfSize[0], env.halfSize[1], env.halfSize[2]), world_material, name);
    worldItems.push_back({lv, G4ThreeVector(position[0], position[1], position[2]), rotation, INVISIBLE, -13, env.bankId});
    envelopeLVs[env.name] = {lv, &env};
  }
  for (const auto& box : cadBoxes) {
    if (box.inBank)
      continue;
    const auto& a = box.axes;
    if (!box.envelope.empty()) {
      // the box in the frame of its envelope (both given in the bank frame): position E (c - c_env), axes A E^T
      const auto& [lvEnvelope, env] = envelopeLVs.at(box.envelope);
      const auto& e = env->axes;
      std::array<double,3> local;
      std::array<std::array<double,3>,3> axes;
      for (int i = 0; i < 3; i++) {
        local[i] = 0;
        for (int j = 0; j < 3; j++) {
          local[i] += e[i][j] * (box.centre[j] - env->centre[j]);
          axes[i][j] = a[i][0] * e[j][0] + a[i][1] * e[j][1] + a[i][2] * e[j][2];
        }
      }
      const G4RotationMatrix inEnvelope(G4ThreeVector(axes[0][0], axes[0][1], axes[0][2]), G4ThreeVector(axes[1][0], axes[1][1], axes[1][2]),
                                        G4ThreeVector(axes[2][0], axes[2][1], axes[2][2]));   // columns: the box axes
      place(createCadBoxLV(box), local[0], local[1], local[2], lvEnvelope, cadColour(box.material), -11, 0,
            new G4RotationMatrix(inEnvelope.inverse()));
      continue;
    }
    const BankTransform bank = banks->getBankTransform(box.bankId);
    const auto position = bank.toGlobal(box.centre);
    const G4RotationMatrix inBank(G4ThreeVector(a[0][0], a[0][1], a[0][2]), G4ThreeVector(a[1][0], a[1][1], a[1][2]),
                                  G4ThreeVector(a[2][0], a[2][1], a[2][2]));   // columns: the box axes
    auto rotation = new G4RotationMatrix((toG4Rotation(bank) * inBank).inverse());
    worldItems.push_back({createCadBoxLV(box), G4ThreeVector(position[0], position[1], position[2]), rotation,
                          cadColour(box.material), -11, box.bankId});
  }
  if (!cadBoxes.empty())
    printf("GeoBCSBanks: %d CAD component boxes, %d envelopes, %d materials (%s)\n", static_cast<int>(cadBoxes.size()),
           static_cast<int>(cadEnvelopes.size()), static_cast<int>(cadMaterials.size()), getParameterString("cad_components").c_str());

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
    worldItems.push_back({new G4LogicalVolume(new G4Box("AbsorberStrip", 0.5 * thickness, 0.5 * width, 0.5 * length),
                                              getParameterMaterial("absorber_strip_material"), "AbsorberStrip"),
                          G4ThreeVector(position[0], position[1], position[2]), rotation, BLACK, -6, bankId});
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
    worldItems.push_back({new G4LogicalVolume(new G4Box("AbsorberPlate", 0.5 * width, 0.5 * thickness, 0.5 * length),
                                              getParameterMaterial("absorber_strip_material"), "AbsorberPlate"),
                          G4ThreeVector(position[0], position[1], position[2]), rotation, BLACK, -7, bankId});
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
    worldItems.push_back({new G4LogicalVolume(new G4Tubs("SourceCylinder", 0, sourceRadius,
                                                         0.5 * getParameterDouble("show_source_cylinder_height_cm") * Units::cm, 0., 2 * M_PI),
                                              world_material, "SourceCylinder"),
                          G4ThreeVector(getParameterDouble("show_source_cylinder_x_meters") * Units::m,
                                        getParameterDouble("show_source_cylinder_y_meters") * Units::m,
                                        getParameterDouble("show_source_cylinder_z_meters") * Units::m),
                          rotation, G4Colour(1, 1, 1), -9, -1});
  }

  const double tankWall = getParameterDouble("vacuum_tank_wall_mm")*Units::mm;
  const double rIn = getParameterDouble("vacuum_tank_inner_radius_m")*Units::m;
  const double zFront = getParameterDouble("vacuum_tank_front_m")*Units::m;
  const double zRear = getParameterDouble("vacuum_tank_rear_m")*Units::m;

  // The mother volumes of the bank groups (with_bank_mothers): cylinders around the beam axis, from the bounding boxes
  // of their content in the world (z range, largest radius; the inner radius: the distance of the beam axis from the
  // content, from the convex hull of the bounding boxes projected along z), with a gap of 1 mm
  struct BankMother {
    std::string name;
    double zMin = 1e99, zMax = -1e99, rMin = 1e99, rMax = 0, zCentre = 0;
    G4LogicalVolume* lv = nullptr;
  };
  std::array<BankMother, 3> mothers;
  mothers[0].name = "front"; mothers[1].name = "middle"; mothers[2].name = "rear";
  auto groupOf = [](int bankId) { return bankId < 0 ? -1 : (bankId == 0 ? 2 : (bankId <= 4 ? 1 : 0)); };
  const bool withMothers = getParameterBoolean("with_bank_mothers");
  if (withMothers) {
    for (const auto& item : worldItems) {
      const int group = groupOf(item.bankId);
      if (group < 0)
        continue;
      auto& m = mothers[group];
      G4ThreeVector low, high;
      item.lv->GetSolid()->BoundingLimits(low, high);
      const G4RotationMatrix objectRotation = item.rotation ? item.rotation->inverse() : G4RotationMatrix();
      std::vector<std::array<double,2>> xy;
      for (int corner = 0; corner < 8; corner++) {
        const G4ThreeVector c = objectRotation * G4ThreeVector(corner & 1 ? high.x() : low.x(), corner & 2 ? high.y() : low.y(),
                                                               corner & 4 ? high.z() : low.z()) + item.position;
        m.zMin = std::min(m.zMin, c.z());
        m.zMax = std::max(m.zMax, c.z());
        m.rMax = std::max(m.rMax, std::hypot(c.x(), c.y()));
        xy.push_back({c.x(), c.y()});
      }
      m.rMin = std::min(m.rMin, distanceFromConvexHull(xy));
    }
    const double gap = 1.0 * Units::mm;
    for (auto& m : mothers) {
      m.zMin -= gap;
      m.zMax += gap;
      m.rMax += gap;
      m.rMin = (m.rMin - gap < 10.0 * Units::mm) ? 0.0 : m.rMin - gap;   // a ring only if it leaves a useful opening
      m.zCentre = 0.5 * (m.zMin + m.zMax);
    }
    if (mothers[0].zMax > mothers[1].zMin || mothers[1].zMax > mothers[2].zMin)
      throw std::runtime_error("GeoBCSBanks: with_bank_mothers: the mother volumes of the bank groups would overlap (front z "
                               + std::to_string(mothers[0].zMin) + " - " + std::to_string(mothers[0].zMax) + " mm, middle "
                               + std::to_string(mothers[1].zMin) + " - " + std::to_string(mothers[1].zMax) + " mm, rear "
                               + std::to_string(mothers[2].zMin) + " - " + std::to_string(mothers[2].zMax) + " mm)");
    if (tankWall > 0)
      for (const auto& m : mothers)
        if (m.rMax > rIn || m.zMin < zFront || m.zMax > zRear)
          throw std::runtime_error("GeoBCSBanks: with_bank_mothers: the " + m.name + " bank mother (radius "
                                   + std::to_string(m.rMax) + " mm, z " + std::to_string(m.zMin) + " - " + std::to_string(m.zMax)
                                   + " mm) does not fit in the vacuum tank");
    const G4Colour motherColour = getParameterBoolean("show_bank_mothers") ? G4Colour(0.3, 0.6, 1.0, 0.15) : INVISIBLE;
    for (auto& m : mothers) {
      const std::string name = "BankMother-" + m.name;
      m.lv = place(new G4Tubs(name, m.rMin, m.rMax, 0.5 * (m.zMax - m.zMin), 0., 2 * M_PI), world_material, 0, 0, m.zCentre,
                   lvWorld, motherColour, -14, 0, nullptr).logvol;
    }
  }
  for (const auto& item : worldItems) {
    G4LogicalVolume* mother = lvWorld;
    double zShift = 0;
    if (withMothers) {
      int group = groupOf(item.bankId);
      if (group < 0)   // the source cylinder: in the mother that contains its centre (if any)
        for (int k = 0; k < 3; k++)
          if (mothers[k].lv->GetSolid()->Inside(item.position - G4ThreeVector(0, 0, mothers[k].zCentre)) != kOutside)
            group = k;
      if (group >= 0) {
        mother = mothers[group].lv;
        zShift = mothers[group].zCentre;
      }
    }
    place(item.lv, item.position.x(), item.position.y(), item.position.z() - zShift, mother, item.colour, item.copyNumber, 0,
          item.rotation);
  }
  if (withMothers)
    for (const auto& m : mothers)
      printf("GeoBCSBanks: bank mother %-6s: z %.1f - %.1f mm, radius %.1f - %.1f mm, %d daughters\n", m.name.c_str(), m.zMin,
             m.zMax, m.rMin, m.rMax, static_cast<int>(m.lv->GetNoDaughters()));

  // Optional vacuum tank (outside all banks: checked in validateParameters)
  if (tankWall > 0) {
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


G4LogicalVolume *GeoBCS::createCadBoxLV(const CadComponents::Box& box) {
  // the volume name: CadComponent-<material key>-<set>:<CAD part>#<solid>.<box>
  const std::string name = "CadComponent-" + box.material + "-" + box.name;
  return new G4LogicalVolume(new G4Box(name, box.halfSize[0], box.halfSize[1], box.halfSize[2]),
                             getCadMaterial(box.material, box.fill), name);
}

G4Material *GeoBCS::getCadMaterial(const std::string& key, double fill) {
  // the material of a CAD component: parameter cad_<key>_material, with the density scaled by the fill of a homogenised
  // box. The factor is rounded to steps of 5% (relative), so that a set needs only a few materials (each material has
  // its own physics tables).
  const std::string base = getParameterString("cad_" + key + "_material");
  if (fill >= 0.975)
    return getParameterMaterial("cad_" + key + "_material");
  const double step = std::log(1.05);
  const double factor = std::exp(std::round(std::log(fill) / step) * step);
  char text[32];
  snprintf(text, sizeof(text), "%.4f", factor);
  // NCrystal cfg strings (without the NCrystal: prefix) take ;density=<f>x, the other material strings :scale_density=<f>
  const bool ncrystalCfg = base.find(".ncmat") != std::string::npos && base.rfind("NCrystal:", 0) != 0;
  const std::string name = ncrystalCfg ? base + ";density=" + text + "x" : base + ":scale_density=" + text;
  auto it = cadMaterials.find(name);
  if (it == cadMaterials.end())
    it = cadMaterials.emplace(name, NamedMaterialProvider::getMaterial(name)).first;
  return it->second;
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
  try {
    for (const auto& box : CadComponents::load(getParameterString("cad_components")))
      if (box.material != "B4C" && box.material != "Tufnol" && box.material != "Al")
        throw std::runtime_error("CAD component " + box.name + ": unknown material key " + box.material
                                 + " (known: B4C, Tufnol, Al: parameters cad_<key>_material)");
    for (const auto& box : CadComponents::loadPack(getParameterString("cad_pack_components")))
      if (box.material != "B4C" && box.material != "Tufnol" && box.material != "Al")
        throw std::runtime_error("CAD pack component " + box.name + ": unknown material key " + box.material);
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
