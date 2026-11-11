#ifndef GRAINITACALOSDACTION_H
#define GRAINITACALOSDACTION_H 1

#include "DDG4/Geant4SensDetAction.h"

#include <cmath>
#include <string>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWGridRhoPhiTheta_k4geo;
}
namespace sim {

  class GrainitaCaloSDData {
  public:
    GrainitaCaloSDData() = default;
    ~GrainitaCaloSDData() = default;

    Geant4Sensitive* sensitive{};
    int rawCollectionID = 0;
    std::string rawCollectionName = "GrainitaCalorimeterHitsRaw";
    // Arm length, in cells along phi and along theta, of the diamond-shaped
    // light-sharing window centred on the hit cell; 0 keeps only the hit cell.
    int neighborRadius = 3;
    const DDSegmentation::FCCSWGridRhoPhiTheta_k4geo* segmentation{};
    // Longitudinal attenuation along IP-pointing fibers, from the hit to the readout surface:
    //   readoutZ == 0 (barrel): cylinder R = outerRadius, fiber length = (outerRadius - R_hit) / sin(theta_hit)
    //   readoutZ  > 0 (endcap): plane |z| = readoutZ,     fiber length = (readoutZ - |z_hit|) / |cos(theta_hit)|
    double fiberAttenuationLength = 1000.; // mm
    double outerRadius = 2645.;            // mm
    double readoutZ = 0.;                  // mm

    /// Length of the IP-pointing fiber from the hit position (mm) to the readout surface (mm), >= 0.
    double fiberLength(double perp, double absZ, double mag) const {
      double length = (readoutZ > 0.) ? (readoutZ - absZ) * mag / absZ : (outerRadius - perp) * mag / perp;
      return length > 0. ? length : 0.;
    }

    // Light response function parameters
    bool useLightResponseFunction = true;
    // Light response function:
    //    y = slope*x + intersect (x<x0)
    //    y = norm*exp(-x/AttLength) (x>=x0)
    double x0 = 0.856;
    double AttLength = 3.4; // mm
    double slope = 0.93;
    double intersect = 0.206;
    double norm = 1. / std::exp(-1. * x0 / AttLength);

    double lightResponse(double distance) const {
      if (distance < x0) {
        return slope * distance + intersect;
      }
      return norm * std::exp(-distance / AttLength);
    }
  };

  using GrainitaCaloSDAction = Geant4SensitiveAction<GrainitaCaloSDData>;

} // namespace sim
} // namespace dd4hep

#endif
