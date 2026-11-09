#ifndef GRAINITACALOSDACTION_H
#define GRAINITACALOSDACTION_H 1

#include "DDG4/Geant4SensDetAction.h"

#include <cmath>
#include <string>
#include <vector>

namespace dd4hep {
namespace sim {

  class GrainitaCaloSDData {
  public:
    GrainitaCaloSDData() = default;
    ~GrainitaCaloSDData() = default;

    Geant4Sensitive* sensitive{};
    int rawCollectionID = 0;
    std::string rawCollectionName = "GrainitaCalorimeterHitsRaw";
    int neighborCellSize = 5;
    bool useLightResponseFunction = true;

    //Fiber attenuation effect
    double fiberAttenuationLength = 529.6; // mm
    double fiberLength = 500.; // mm
    double innerRadius = 2200.; // mm
    double fiberReflection = 0.8;
    double longitudinalAttenuationNorm = 1. / (
      1. + fiberReflection * std::exp(-2*fiberLength / fiberAttenuationLength)
    ); // Normalize: f(fiberLength) = 1. 

    // Light response function parameters
    double responseX0 = 0.5;
    double responseFraction = 0.74;
    double responseAttLength1 = 2.45;
    double responseAttLength2 = 5.34;
    double norm = 1. / ( responseFraction*std::exp(-responseX0 / responseAttLength1) + (1-responseFraction) * std::exp(-responseX0 / responseAttLength2));
    
    double lightResponse(double distance) const {
      if(responseAttLength1<0 || responseAttLength2<0) return 1;

      if (distance < responseX0) {
        return 0.;
      }
      return norm * ( responseFraction*std::exp(-distance / responseAttLength1) + (1-responseFraction) * std::exp(-distance / responseAttLength2));
    }

    double longitudinalResponse(double globalR, double theta) const {
      if (fiberAttenuationLength <= 0. || fiberLength <= 0.) {
        return 1.;
      }
      const double sinTheta = std::sin(theta);
      if (std::abs(sinTheta) < 1.e-12) {
        return 1.;
      }
      const double zi = globalR - innerRadius / sinTheta;
      return longitudinalAttenuationNorm *
             (std::exp(-(fiberLength - zi) / fiberAttenuationLength) +
              fiberReflection * std::exp(-(fiberLength + zi) / fiberAttenuationLength));
    }

  };

  using GrainitaCaloSDAction = Geant4SensitiveAction<GrainitaCaloSDData>;

} // namespace sim
} // namespace dd4hep

#endif
