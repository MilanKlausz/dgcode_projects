#include "G4GeoLoki/CalibMasks.hh"

#include <cmath>
#include <cassert>

// NOTE: All values are expressed in [mm] unit


CalibMasks::CalibMasksBase::CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, double elevationFromTubeFront, std::vector<double> pattern)
  : m_name(name),
    m_thickness(thickness),
    m_height(height),
    m_leftTubeEndDistance(leftTubeEndDistance),
    m_elevationFromTubeFront((assert(elevationFromTubeFront >= 0.05), elevationFromTubeFront)),
    m_pattern(pattern)
{
}

std::string CalibMasks::CalibMasksBase::getName() const{
  return this->m_name;
}
double CalibMasks::CalibMasksBase::getThickness() const{
  return this->m_thickness;
}
double CalibMasks::CalibMasksBase::getHeight() const{
  return this->m_height;
}
double CalibMasks::CalibMasksBase::getLeftTubeEndDistance() const{
  return this->m_leftTubeEndDistance;
}
double CalibMasks::CalibMasksBase::getElevationFromTubeFront() const{
  return this->m_elevationFromTubeFront;
}
std::vector<double> CalibMasks::CalibMasksBase::getPattern() const{
  return this->m_pattern;
}
double CalibMasks::CalibMasksBase::getWidth() const{
  return std::accumulate(this->m_pattern.begin(), this->m_pattern.end(), decltype(this->m_pattern)::value_type(0));
}

G4Material* CalibMasks::maskMaterial = NamedMaterialProvider::getMaterial("MAT_B4C:b10_enrichment=0.95");
G4Material* CalibMasks::maskBoxMaterial = NamedMaterialProvider::getMaterial("G4_Vacuum");
