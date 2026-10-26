#ifndef G4GeoLoki_CadComponents_hh
#define G4GeoLoki_CadComponents_hh

#include <array>
#include <string>
#include <vector>

/// Optional components of the LOKI detector system taken from the CAD models, simplified to boxes.
///
/// A component set is a data file G4GeoLoki/data/cad_components_<name>.txt (a value containing a '/' is used as the
/// path of such a file directly). The files are generated from the STEP models (the cad_work scripts:
/// gen_components.py), one box per line:
///   <name> <material> <bankId> <mother> cx cy cz ax ay az bx by bz nx ny nz hx hy hz
/// with the box centre c [mm] and its three axes a, b, n (unit vectors) in the frame of the bank bankId (the frame of
/// the bank volume: x = depth, y = across the tubes, z = along the tubes), and the half lengths h [mm] along a, b, n.
/// mother "bank": the box is inside the bank volume (placed in it); "world": outside every bank volume, placed in the
/// world with the bank transform. So a component follows its bank (bank calibration, rear detector distance).
/// <material> is a key (e.g. B4C, Tufnol, Al), mapped to a Geant4 material by the geometry module.
/// Lines starting with '#' are comments.
///
/// Simplified sets (gen_simplified.py) use keyword lines, with homogenised boxes and envelopes:
///   box <name> <material> <fill> <bankId> <mother> cx cy cz ax ay az bx by bz nx ny nz hx hy hz
///   envelope <name> <bankId> cx cy cz ax ay az bx by bz nx ny nz hx hy hz
/// fill: the density factor of the box (the material of a CAD solid spread over its box, e.g. a plate with holes or a
/// hollow beam); mother: bank, world, or the name of an envelope of the same file. An envelope is a box of the world
/// material outside the bank volumes (placed in the world with the bank transform, like a world box) that holds boxes
/// (fewer volumes directly in the world). Both are given in the frame of their bank.
class CadComponents {
public:
  struct Box {
    std::string name, material;
    int bankId;
    bool inBank;                     // mother: the bank volume (true) or the world / an envelope (false)
    std::string envelope;            // the name of the envelope holding the box (empty: none)
    double fill = 1.0;               // density factor (homogenised box)
    std::array<double,3> centre;     // [Geant4 length units], bank frame
    std::array<std::array<double,3>,3> axes; // the box axes (rows) in the bank frame
    std::array<double,3> halfSize;   // [Geant4 length units]
  };
  struct Envelope {
    std::string name;                // <set>:<name>
    int bankId;
    std::array<double,3> centre;
    std::array<std::array<double,3>,3> axes;
    std::array<double,3> halfSize;
  };
  /// the boxes of a comma separated list of set names (or file paths); empty: none
  static std::vector<Box> load(const std::string& sets);
  /// the envelopes of the sets
  static std::vector<Envelope> loadEnvelopes(const std::string& sets);
  /// the boxes and envelopes of one data file
  static std::vector<Box> loadFile(const std::string& fileName, const std::string& setName,
                                   std::vector<Envelope>* envelopes = nullptr);
  /// the data file of a set name (empty if there is none); prefix "cad_components_" or "cad_pack_"
  static std::string findFile(const std::string& set, const std::string& prefix = "cad_components_");

  /// Pack components: boxes in the frame of every detector pack (the G4Para EmptyPackBox: x = depth, y = across the
  /// tubes, z = along the tubes), data files G4GeoLoki/data/cad_pack_<name>.txt, one box per line:
  ///   <name> <material> x0 x1 y0 y1 <zkind> a b      [mm]
  /// zkind "length": z in [-(L/2+a), L/2+a] (L: the straw length of the bank); "ends": z in [L/2+a, L/2+b] and its mirror
  /// An optional 10th value is the density factor (fill) of a homogenised box (default 1).
  struct PackBox {
    std::string name, material;
    double x0, x1, y0, y1;           // [Geant4 length units]
    bool alongLength;                // zkind length (true) or ends (false)
    double a, b;                     // [Geant4 length units]
    double fill = 1.0;
  };
  static std::vector<PackBox> loadPack(const std::string& sets);
};

#endif
