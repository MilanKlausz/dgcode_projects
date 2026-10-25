#include "G4GeoLoki/CadComponents.hh"
#include "Core/FindData.hh"
#include "Units/Units.hh"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::string CadComponents::findFile(const std::string& set) {
  if (set.find('/') != std::string::npos)
    return set;
  return Core::findData("G4GeoLoki", "cad_components_" + set + ".txt");
}

std::vector<CadComponents::Box> CadComponents::load(const std::string& sets) {
  std::vector<Box> boxes;
  std::stringstream items(sets);
  std::string set;
  while (std::getline(items, set, ',')) {
    set.erase(0, set.find_first_not_of(" \t"));
    set.erase(set.find_last_not_of(" \t") + 1);
    if (set.empty())
      continue;
    const std::string path = findFile(set);
    if (path.empty())
      throw std::runtime_error("CadComponents: unknown component set '" + set + "' (no data file G4GeoLoki/cad_components_"
                               + set + ".txt)");
    const auto more = loadFile(path, set);
    boxes.insert(boxes.end(), more.begin(), more.end());
  }
  return boxes;
}

std::vector<CadComponents::Box> CadComponents::loadFile(const std::string& fileName, const std::string& setName) {
  std::ifstream file(fileName);
  if (!file)
    throw std::runtime_error("CadComponents: can not open " + fileName);
  std::vector<Box> boxes;
  std::string line;
  int lineNumber = 0;
  while (std::getline(file, line)) {
    lineNumber++;
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream fields(line);
    Box box;
    std::string mother;
    fields >> box.name >> box.material >> box.bankId >> mother;
    double v[15];
    for (double& x : v)
      fields >> x;
    if (fields.fail() || box.bankId < 0 || box.bankId > 8 || (mother != "bank" && mother != "world"))
      throw std::runtime_error("CadComponents: invalid line " + std::to_string(lineNumber) + " of " + fileName);
    box.name = setName + ":" + box.name;
    box.inBank = (mother == "bank");
    for (int k = 0; k < 3; k++) {
      box.centre[k] = v[k] * Units::mm;
      box.halfSize[k] = v[12 + k] * Units::mm;
      for (int j = 0; j < 3; j++)
        box.axes[k][j] = v[3 + 3 * k + j];
    }
    // the axes must form a right-handed orthonormal frame (Geant4 rotations)
    const auto& a = box.axes;
    const double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                     + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if (std::fabs(std::fabs(det) - 1) > 1e-5)
      throw std::runtime_error("CadComponents: the axes of " + box.name + " are not orthonormal (" + fileName + ")");
    if (det < 0)                                       // a left-handed frame: flip the third axis (a box is symmetric)
      for (int j = 0; j < 3; j++)
        box.axes[2][j] = -box.axes[2][j];
    boxes.push_back(box);
  }
  return boxes;
}
