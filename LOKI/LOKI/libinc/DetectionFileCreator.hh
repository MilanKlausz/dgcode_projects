#ifndef LOKI_DetectionFileCreator_hh
#define LOKI_DetectionFileCreator_hh

#include "mcpl.h"
#include "GriffDataRead/Setup.hh"
#include <iostream>

class DetectionFileCreator {
public:
  // weighted: store the weight of each detection event (simulations of weighted MCPL input), otherwise all weights are 1
  DetectionFileCreator(const char* fileName, GriffDataRead::StrMap& userData, bool weighted = false);
  DetectionFileCreator(const char* fileName);
  ~DetectionFileCreator();

  void addDetectionEvent(const int pixelId, const double tof, const double weight = 1.0);

private:
  std::string m_fileName;
  mcpl_outfile_t detMcpl;
  mcpl_particle_t *mcplParticle;

  void initiateMcplParticle();
};

#endif
