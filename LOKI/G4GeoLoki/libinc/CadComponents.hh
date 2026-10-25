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
class CadComponents {
public:
  struct Box {
    std::string name, material;
    int bankId;
    bool inBank;                     // mother: the bank volume (true) or the world (false)
    std::array<double,3> centre;     // [Geant4 length units], bank frame
    std::array<std::array<double,3>,3> axes; // the box axes (rows) in the bank frame
    std::array<double,3> halfSize;   // [Geant4 length units]
  };
  /// the boxes of a comma separated list of set names (or file paths); empty: none
  static std::vector<Box> load(const std::string& sets);
  /// the boxes of one data file
  static std::vector<Box> loadFile(const std::string& fileName, const std::string& setName);
  /// the data file of a set name (empty if there is none)
  static std::string findFile(const std::string& set);
};

#endif
