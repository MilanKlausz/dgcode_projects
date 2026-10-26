#include "G4GeoLoki/CadComponents.hh"
#include "Core/FindData.hh"
#include "Units/Units.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::string CadComponents::findFile(const std::string& set, const std::string& prefix) {
  if (set.find('/') != std::string::npos)
    return set;
  return Core::findData("G4GeoLoki", prefix + set + ".txt");
}

namespace {
  std::vector<std::string> setNames(const std::string& sets) {
    std::vector<std::string> names;
    std::stringstream items(sets);
    std::string set;
    while (std::getline(items, set, ',')) {
      set.erase(0, set.find_first_not_of(" \t"));
      set.erase(set.find_last_not_of(" \t") + 1);
      if (!set.empty())
        names.push_back(set);
    }
    return names;
  }
}

std::vector<CadComponents::PackBox> CadComponents::loadPack(const std::string& sets) {
  std::vector<PackBox> boxes;
  for (const auto& set : setNames(sets)) {
    const std::string path = findFile(set, "cad_pack_");
    if (path.empty())
      throw std::runtime_error("CadComponents: unknown pack component set '" + set + "' (no data file G4GeoLoki/cad_pack_" + set + ".txt)");
    std::ifstream file(path);
    if (!file)
      throw std::runtime_error("CadComponents: can not open " + path);
    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line)) {
      lineNumber++;
      if (line.empty() || line[0] == '#')
        continue;
      std::stringstream fields(line);
      PackBox box;
      std::string kind;
      fields >> box.name >> box.material >> box.x0 >> box.x1 >> box.y0 >> box.y1 >> kind >> box.a >> box.b;
      bool ok = !fields.fail();
      if (ok) {                                        // an optional 10th value: the fill
        double fill = 1.0;
        if (fields >> fill)
          box.fill = fill;
        else
          ok = fields.eof();
      }
      if (!ok || (kind != "length" && kind != "ends") || box.x1 <= box.x0 || box.y1 <= box.y0
          || (kind == "ends" && box.b <= box.a) || !(box.fill > 0 && box.fill <= 1))
        throw std::runtime_error("CadComponents: invalid line " + std::to_string(lineNumber) + " of " + path);
      box.name = set + ":" + box.name;
      box.alongLength = (kind == "length");
      box.x0 *= Units::mm; box.x1 *= Units::mm; box.y0 *= Units::mm; box.y1 *= Units::mm;
      box.a *= Units::mm; box.b *= Units::mm;
      boxes.push_back(box);
    }
  }
  return boxes;
}

namespace {
  std::vector<CadComponents::Box> loadSets(const std::string& sets, std::vector<CadComponents::Envelope>* envelopes) {
    std::vector<CadComponents::Box> boxes;
    for (const auto& set : setNames(sets)) {
      const std::string path = CadComponents::findFile(set);
      if (path.empty())
        throw std::runtime_error("CadComponents: unknown component set '" + set + "' (no data file G4GeoLoki/cad_components_"
                                 + set + ".txt)");
      const auto more = CadComponents::loadFile(path, set, envelopes);
      boxes.insert(boxes.end(), more.begin(), more.end());
    }
    return boxes;
  }
}

std::vector<CadComponents::Box> CadComponents::load(const std::string& sets) {
  return loadSets(sets, nullptr);
}

std::vector<CadComponents::Envelope> CadComponents::loadEnvelopes(const std::string& sets) {
  std::vector<Envelope> envelopes;
  loadSets(sets, &envelopes);
  return envelopes;
}

std::vector<CadComponents::Box> CadComponents::loadFile(const std::string& fileName, const std::string& setName,
                                                       std::vector<Envelope>* envelopes) {
  std::ifstream file(fileName);
  if (!file)
    throw std::runtime_error("CadComponents: can not open " + fileName);
  // the axes must form a right-handed orthonormal frame (Geant4 rotations); a left-handed frame: flip the third axis
  // (a box is symmetric)
  auto checkAxes = [&](std::array<std::array<double,3>,3>& a, const std::string& name) {
    const double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                     + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if (std::fabs(std::fabs(det) - 1) > 1e-5)
      throw std::runtime_error("CadComponents: the axes of " + name + " are not orthonormal (" + fileName + ")");
    if (det < 0)
      for (int j = 0; j < 3; j++)
        a[2][j] = -a[2][j];
  };
  std::vector<Box> boxes;
  std::vector<std::string> envelopeNames;
  std::string line;
  int lineNumber = 0;
  while (std::getline(file, line)) {
    lineNumber++;
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream fields(line);
    const auto invalid = [&]() {
      return std::runtime_error("CadComponents: invalid line " + std::to_string(lineNumber) + " of " + fileName);
    };
    std::string first;
    fields >> first;
    if (first == "envelope") {
      Envelope env;
      fields >> env.name >> env.bankId;
      double v[15];
      for (double& x : v)
        fields >> x;
      if (fields.fail() || env.bankId < 0 || env.bankId > 8)
        throw invalid();
      envelopeNames.push_back(env.name);
      env.name = setName + ":" + env.name;
      for (int k = 0; k < 3; k++) {
        env.centre[k] = v[k] * Units::mm;
        env.halfSize[k] = v[12 + k] * Units::mm;
        for (int j = 0; j < 3; j++)
          env.axes[k][j] = v[3 + 3 * k + j];
      }
      checkAxes(env.axes, env.name);
      if (envelopes)
        envelopes->push_back(env);
      continue;
    }
    Box box;
    std::string mother;
    if (first == "box")               // simplified sets: box <name> <material> <fill> <bankId> <mother> ...
      fields >> box.name >> box.material >> box.fill >> box.bankId >> mother;
    else {                            // the original format: <name> <material> <bankId> <mother> ...
      box.name = first;
      fields >> box.material >> box.bankId >> mother;
    }
    double v[15];
    for (double& x : v)
      fields >> x;
    if (fields.fail() || box.bankId < 0 || box.bankId > 8 || !(box.fill > 0 && box.fill <= 1))
      throw invalid();
    box.name = setName + ":" + box.name;
    box.inBank = (mother == "bank");
    if (mother != "bank" && mother != "world") {
      if (std::find(envelopeNames.begin(), envelopeNames.end(), mother) == envelopeNames.end())
        throw std::runtime_error("CadComponents: line " + std::to_string(lineNumber) + " of " + fileName + ": mother '" + mother
                                 + "' is neither bank, world nor an envelope defined before");
      box.envelope = setName + ":" + mother;
    }
    for (int k = 0; k < 3; k++) {
      box.centre[k] = v[k] * Units::mm;
      box.halfSize[k] = v[12 + k] * Units::mm;
      for (int j = 0; j < 3; j++)
        box.axes[k][j] = v[3 + 3 * k + j];
    }
    checkAxes(box.axes, box.name);
    boxes.push_back(box);
  }
  return boxes;
}
