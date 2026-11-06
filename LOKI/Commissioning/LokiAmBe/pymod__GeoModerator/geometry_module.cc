#include "G4Interfaces/GeoConstructToPy.hh"

// The AmBe source in its HDPE moderator, for the moderator simulation (sb_lokiambe_moderatorsim) that makes the
// spectrum tables of the AmBe source (AmBeSourceGen): the neutrons leaving the moderator are recorded in RecordShell.

class GeoModerator final : public G4Interfaces::GeoConstructBase
{
public:
  GeoModerator();
  virtual ~GeoModerator(){}
  G4VPhysicalVolume* Construct() override;
  bool validateParameters() override;
};

SBLDPY_MODULE( _GeoModerator ) { G4Interfaces::geoToPy<GeoModerator>( m ); }

#include "G4Tubs.hh"
#include "G4Box.hh"
#include <iostream>

GeoModerator::GeoModerator()
  : GeoConstructBase("LokiAmBe/GeoModerator")
{
  // HDPE cylinder with a vertical axis (y), centred at the origin; the source is at the origin
  addParameterDouble("moderator_radius_cm", 13.0);
  addParameterDouble("moderator_height_cm", 26.0);
  addParameterString("moderator_material", "Polyethylene_CH2.ncmat;density=0.95gcm3");
  // stainless steel source capsule around the origin (0: none)
  addParameterDouble("capsule_radius_cm", 1.0);
  addParameterDouble("capsule_height_cm", 3.0);
  addParameterString("capsule_material", "G4_STAINLESS-STEEL");
}

bool GeoModerator::validateParameters()
{
  const double r = getParameterDouble("moderator_radius_cm"), h = getParameterDouble("moderator_height_cm");
  const double cr = getParameterDouble("capsule_radius_cm"), ch = getParameterDouble("capsule_height_cm");
  if (!(r > 0 && h > 0 && cr >= 0 && ch >= 0 && cr < r && ch < h)) {
    std::cout << "Error: invalid moderator / capsule dimensions" << std::endl;
    return false;
  }
  return true;
}

G4VPhysicalVolume* GeoModerator::Construct()
{
  const double r = getParameterDouble("moderator_radius_cm") * Units::cm;
  const double h = getParameterDouble("moderator_height_cm") * Units::cm;
  const double cr = getParameterDouble("capsule_radius_cm") * Units::cm;
  const double ch = getParameterDouble("capsule_height_cm") * Units::cm;
  auto vacuum = getMaterial("Vacuum");
  auto world = place(new G4Box("World", 2 * r + 10 * Units::cm, h + 10 * Units::cm, 2 * r + 10 * Units::cm),
                     vacuum, 0, 0, 0, 0, INVISIBLE);
  // G4Tubs has its axis along z: rotate it to the vertical (y) axis
  auto rot = new G4RotationMatrix();
  rot->rotateX(90 * Units::deg);
  // a 1 mm vacuum shell around the moderator: the neutrons leaving the moderator are recorded when they enter it
  auto shell = place(new G4Tubs("RecordShell", 0, r + 1 * Units::mm, 0.5 * h + 1 * Units::mm, 0, 360 * Units::deg),
                     vacuum, 0, 0, 0, world.logvol, INVISIBLE, 0, 0, rot);
  auto moderator = place(new G4Tubs("Moderator", 0, r, 0.5 * h, 0, 360 * Units::deg),
                         getParameterMaterial("moderator_material"), 0, 0, 0, shell.logvol,
                         G4Colour(0.9, 0.9, 0.9));
  if (cr > 0 && ch > 0)
    place(new G4Tubs("Capsule", 0, cr, 0.5 * ch, 0, 360 * Units::deg), getParameterMaterial("capsule_material"),
          0, 0, 0, moderator.logvol, G4Colour(0.5, 0.5, 0.5));
  return world.physvol;
}
