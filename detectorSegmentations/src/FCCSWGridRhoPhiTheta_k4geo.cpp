#include "detectorSegmentations/FCCSWGridRhoPhiTheta_k4geo.h"
#include "DD4hep/Detector.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

namespace dd4hep {
namespace DDSegmentation {

  /// default constructor using an encoding string
  FCCSWGridRhoPhiTheta_k4geo::FCCSWGridRhoPhiTheta_k4geo(const std::string& cellEncoding)
      : FCCSWGridPhiTheta_k4geo(cellEncoding) {
    commonSetup();
  }

  FCCSWGridRhoPhiTheta_k4geo::FCCSWGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder)
      : FCCSWGridPhiTheta_k4geo(decoder) {
    commonSetup();
  }

  /// Initialization common to all ctors.
  void FCCSWGridRhoPhiTheta_k4geo::commonSetup() {
    // define type and description
    _type = "FCCSWGridRhoPhiTheta_k4geo";
    _description = "Rho-Phi-theta segmentation in the global coordinates";
    m_rhoBins.clear();

    // register all necessary parameters (additional to those registered in GridTheta_k4geo)
    registerParameter("grid_rho", "Grid size in rho", m_grid_rho, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("rho_bins", "Non-uniform rho bins", m_rhoBins, std::vector<double>{},
                      SegmentationParameter::NoUnit, true);
    registerParameter("rho_bin_count", "Number of bins for a uniform rho grid", m_rhoBinCount, 0,
                      SegmentationParameter::NoUnit, true);
    registerParameter("offset_R", "Offset in R", m_offsetR, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("offset_Z", "Endcap mode: |z| of the front plane (0: barrel mode, rho offset from offset_R)",
                      m_offsetZ, 0., SegmentationParameter::LengthUnit, true);
    registerParameter("min_theta", "Minimum physical theta", m_minTheta, 0., SegmentationParameter::AngleUnit, true);
    registerParameter("max_theta", "Maximum physical theta", m_maxTheta, M_PI, SegmentationParameter::AngleUnit, true);
    registerParameter("clear_extra_fields", "Clear CellID fields other than system, phi, theta and rho",
                      m_clearExtraFields, true);
    registerParameter("scale_phi_bins_with_sin_theta",
                      "Number of phi bins per theta ring = round(phi_bins * sin(theta)) (phi_bins: ring at 90 deg)",
                      m_scalePhiBinsWithSinTheta, false, SegmentationParameter::NoUnit, true);
    registerIdentifier("identifier_rho", "Cell ID identifier for rho", m_rhoID, "rho");

    m_thetaIndex = decoder()->index(fieldNameTheta());
    m_phiIndex = decoder()->index(fieldNamePhi());
    m_rhoIndex = decoder()->index(m_rhoID);
  }

  void FCCSWGridRhoPhiTheta_k4geo::clearExtraFields(CellID& cID) const {
    if (!m_clearExtraFields)
      return;

    for (const auto& field : decoder()->fields()) {
      const std::string& name = field.name();
      if (name == "system" || name == fieldNamePhi() || name == fieldNameTheta() || name == m_rhoID)
        continue;
      field.set(cID, 0);
    }
  }

  /// determine the local based on the cell ID
  Vector3D FCCSWGridRhoPhiTheta_k4geo::position(const CellID& cID) const {
    return positionFromRThetaPhi(rho(cID) * std::sin(theta(cID)), theta(cID), phi(cID));
  }

  /// determine the cell ID based on the position
  CellID FCCSWGridRhoPhiTheta_k4geo::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                                            const VolumeID& vID) const {
    CellID cID = vID;
    double lTheta = thetaFromXYZ(globalPosition);
    double lPhi = phiFromXYZ(globalPosition);
    double lrho =
        sqrt(radiusFromXYZ(globalPosition) * radiusFromXYZ(globalPosition) + globalPosition.Z * globalPosition.Z);

    // dd4hep::Segmentation::positionToBin: int(floor((position + 0.5 * cellSize - offset) / cellSize));
    // The phi grid depends on the theta ring (scale_phi_bins_with_sin_theta), hence theta is binned first.
    const int thetaBin = positionToBin(lTheta, gridSizeTheta(), offsetTheta());
    const int nPhi = phiBinsInThetaBin(thetaBin);
    decoder()->set(cID, m_thetaIndex, thetaBin);
    decoder()->set(cID, m_phiIndex,
                   wrapPhiBin(positionToBin(lPhi, 2 * M_PI / static_cast<double>(nPhi), offsetPhi()), nPhi));

    // std::cout<<" Position to cellID "<<std::endl;
    // std::cout<<" Print rho bins: "<<std::endl;
    // for(size_t i=0; i<m_rhoBins.size(); i++) std::cout<<m_rhoBins[i]<<'\t';
    // std::cout<<std::endl;
    // std::cout<<" Global position: ("<<globalPosition.X<<", "<<globalPosition.Y<<", "<<globalPosition.Z<<")
    // "<<std::endl; std::cout<<" Input: rho "<<lrho<<", grid rho: "<<m_grid_rho<<", offsetR "<<m_offsetR<<", theta
    // "<<lTheta<<", offset rho "<<m_offsetR/std::sin(lTheta)<<std::endl; std::cout<<" rho ID:  "<<positionToBin(lrho,
    // m_grid_rho, m_offsetR/std::sin(lTheta)+m_grid_rho/2. )<<std::endl;

    // Two ways for rho segmentation:
    if (m_rhoBins.size() < 2)
      // positionToBin(position, cellsize, offset): int(floor((position + 0.5 * cellSize - offset) / cellSize));
      decoder()->set(cID, m_rhoIndex, positionToBin(lrho, m_grid_rho, rhoOffset(lTheta) + m_grid_rho / 2.));
    else {
      // Note: a bug was found in dd4hep::Segmentation, and fixed in PR: https://github.com/AIDASoft/DD4hep/pull/1625
      // Following function can be used after updating dd4hep to latest ver.
      // int rhoID = Segmentation::positionToBin(lrho, m_rhoBins, m_offsetR/std::sin(lTheta));  //Fixed in DD4hep main
      // branch, no tag yet.

      std::vector<double>::const_iterator bin =
          std::upper_bound(m_rhoBins.begin(), m_rhoBins.end(), lrho - rhoOffset(lTheta));
      int rhoID = bin - m_rhoBins.begin() - 1;

      // Add a protection for cell ID, if the position is out of bin range
      // WARNING: may introduce bug. Better to well-define bin edges.
      if (rhoID < 0)
        rhoID = 0;
      if (rhoID > static_cast<int>(m_rhoBins.size()) - 2)
        rhoID = static_cast<int>(m_rhoBins.size()) - 2;
      decoder()->set(cID, m_rhoIndex, rhoID);
    }

    clearExtraFields(cID);
    return cID;
  }

  /// determine the distance to IP rho based on the cell ID
  double FCCSWGridRhoPhiTheta_k4geo::rho(const CellID cID) const {
    CellID rhoValue = decoder()->get(cID, m_rhoIndex);
    double m_theta = theta(cID);

    // dd4hep::Segmentation::binToPosition: bin * cellSize + offset;

    // std::cout<<" CellID to position"<<std::endl;
    // std::cout<<" Print rho bins: "<<std::endl;
    // for(size_t i=0; i<m_rhoBins.size(); i++) std::cout<<m_rhoBins[i]<<'\t';
    // std::cout<<std::endl;
    // std::cout<<" rho ID val: "<<rhoValue<<", theta val: "<<m_theta<<", offsetR "<<m_offsetR<<", rho =
    // "<<0.5*(m_rhoBins[rhoValue] + m_rhoBins[rhoValue+1])<<std::endl;

    if (m_rhoBins.size() < 2)
      return binToPosition(rhoValue, m_grid_rho, rhoOffset(m_theta) + m_grid_rho / 2.);
    else
      return rhoOffset(m_theta) + 0.5 * (m_rhoBins[rhoValue] + m_rhoBins[rhoValue + 1]);
  }

  double FCCSWGridRhoPhiTheta_k4geo::rhoOffset(double theta) const {
    if (m_offsetZ > 0.)
      return m_offsetZ / std::fabs(std::cos(theta)); // endcap: front plane |z| = offset_Z
    return m_offsetR / std::sin(theta);              // barrel: front cylinder R = offset_R
  }

  int FCCSWGridRhoPhiTheta_k4geo::wrapPhiBin(int bin, int nPhi) const {
    // phiFromXYZ() uses atan2(), hence phi bins are generally signed.  The first bin is
    // positionToBin(-pi), written as floor(0.5 - n/2 - offset/size) so that it is exact for zero offset.
    const double gridSize = 2 * M_PI / static_cast<double>(nPhi);
    const int firstPhiBin = static_cast<int>(std::floor(0.5 - 0.5 * nPhi - offsetPhi() / gridSize));
    return firstPhiBin + ((bin - firstPhiBin) % nPhi + nPhi) % nPhi;
  }

  int FCCSWGridRhoPhiTheta_k4geo::phiBinsInThetaBin(int thetaBin) const {
    if (!m_scalePhiBinsWithSinTheta)
      return phiBins();
    const double thetaCentre = binToPosition(thetaBin, gridSizeTheta(), offsetTheta());
    // sin(pi - theta) = sin(theta): the rings of the -z endcap mirror those of +z.
    const long nPhi = std::lround(phiBins() * std::sin(thetaCentre));
    return static_cast<int>(std::max(1L, nPhi));
  }

  double FCCSWGridRhoPhiTheta_k4geo::phi(const CellID cID) const {
    const int thetaBin = static_cast<int>(decoder()->get(cID, m_thetaIndex));
    const int phiBin = static_cast<int>(decoder()->get(cID, m_phiIndex));
    return binToPosition(phiBin, gridSizePhiInThetaBin(thetaBin), offsetPhi());
  }

  std::vector<double> FCCSWGridRhoPhiTheta_k4geo::cellDimensions(const CellID& cID) const {
    const int thetaBin = static_cast<int>(decoder()->get(cID, m_thetaIndex));
    return {gridSizePhiInThetaBin(thetaBin), gridSizeTheta()};
  }

  int FCCSWGridRhoPhiTheta_k4geo::firstThetaBin() const {
    const auto& thetaField = (*decoder())[m_thetaIndex];
    return std::max(positionToBin(m_minTheta, gridSizeTheta(), offsetTheta()), static_cast<int>(thetaField.minValue()));
  }

  int FCCSWGridRhoPhiTheta_k4geo::lastThetaBin() const {
    const auto& thetaField = (*decoder())[m_thetaIndex];
    const double lastTheta = (m_offsetZ > 0.) ? M_PI - m_minTheta : m_maxTheta;
    return std::min(positionToBin(lastTheta, gridSizeTheta(), offsetTheta()), static_cast<int>(thetaField.maxValue()));
  }

  bool FCCSWGridRhoPhiTheta_k4geo::thetaBinInCoverage(int bin) const {
    const auto& thetaField = (*decoder())[m_thetaIndex];
    if (bin < static_cast<int>(thetaField.minValue()) || bin > static_cast<int>(thetaField.maxValue()))
      return false;
    auto inRange = [&](double thetaLow, double thetaHigh) {
      return bin >= positionToBin(thetaLow, gridSizeTheta(), offsetTheta()) &&
             bin <= positionToBin(thetaHigh, gridSizeTheta(), offsetTheta());
    };
    if (m_offsetZ > 0.) // endcap: +z and mirrored -z coverage
      return inRange(m_minTheta, m_maxTheta) || inRange(M_PI - m_maxTheta, M_PI - m_minTheta);
    return inRange(m_minTheta, m_maxTheta);
  }

  // overrides the DDSegmentation::Segmentation::neighbours method
  void FCCSWGridRhoPhiTheta_k4geo::neighbours(const CellID& cID, std::set<CellID>& neighbours) const {
    neighboursInWindow(cID, neighbours, 1, 1, 1);
  }

  void FCCSWGridRhoPhiTheta_k4geo::neighboursInWindow(const CellID& cID, std::set<CellID>& neighbours, int rPhi,
                                                      int rTheta, int rRho) const {
    const int phiBin = static_cast<int>(decoder()->get(cID, m_phiIndex));
    const int thetaBin = static_cast<int>(decoder()->get(cID, m_thetaIndex));
    const int rhoBin = static_cast<int>(decoder()->get(cID, m_rhoIndex));
    // Clean the module and sector ID, so the neighbour finding only work by comparing the theta/phi/rho IDs.
    CellID base = cID;
    clearExtraFields(base);

    // Allowed (offset, bin) pairs along theta and rho, including the cell's own bin.
    // Phi candidates depend on the theta ring and are built in the loop below.
    std::vector<std::pair<int, int>> thetaCandidates{{0, thetaBin}}, rhoCandidates{{0, rhoBin}};

    // Theta is not periodic.  Use the bins reached by points at the detector
    // boundaries rather than testing bin centres: the first/last physical
    // cells can have centres just outside the boundary.
    for (int d = 1; d <= rTheta; ++d) {
      for (const int offset : {-d, d}) {
        if (thetaBinInCoverage(thetaBin + offset))
          thetaCandidates.push_back({offset, thetaBin + offset});
      }
    }

    // Non-uniform rho bins define their count directly.  A uniform grid needs
    // rho_bin_count from the compact file; the encoded field width is only a
    // storage limit and is not normally the number of physical layers.
    const auto& rhoField = (*decoder())[m_rhoIndex];
    const int minRhoBin = std::max(0, static_cast<int>(rhoField.minValue()));
    int maxRhoBin = static_cast<int>(rhoField.maxValue());
    if (m_rhoBins.size() >= 2)
      maxRhoBin = std::min(static_cast<int>(m_rhoBins.size()) - 2, maxRhoBin);
    else if (m_rhoBinCount > 0)
      maxRhoBin = std::min(m_rhoBinCount - 1, maxRhoBin);
    for (int d = 1; d <= rRho; ++d) {
      for (const int offset : {-d, d}) {
        if (rhoBin + offset >= minRhoBin && rhoBin + offset <= maxRhoBin)
          rhoCandidates.push_back({offset, rhoBin + offset});
      }
    }

    // Diamond (L1) window: |dPhi|/rPhi + |dTheta|/rTheta + |dRho|/rRho <= 1, evaluated in integers.
    // An axis with radius 0 only has offset 0, so its weight value is irrelevant.
    const long wPhi = std::max(rPhi, 1), wTheta = std::max(rTheta, 1), wRho = std::max(rRho, 1);
    // With scale_phi_bins_with_sin_theta the rings have different phi grids: in another ring, the phi offsets
    // are counted from the bin containing the phi of this cell's centre.  Rings with the same grid use phiBin
    // directly.
    const int nPhi = phiBinsInThetaBin(thetaBin);
    const double phiCentre = binToPosition(phiBin, 2 * M_PI / static_cast<double>(nPhi), offsetPhi());
    for (const auto& [dTheta, t] : thetaCandidates) {
      const int nPhiRing = phiBinsInThetaBin(t);
      const int phiBinRing =
          (nPhiRing == nPhi)
              ? phiBin
              : wrapPhiBin(positionToBin(phiCentre, 2 * M_PI / static_cast<double>(nPhiRing), offsetPhi()), nPhiRing);
      std::vector<std::pair<int, int>> phiCandidates{{0, phiBinRing}};
      if (nPhiRing > 1) {
        for (int d = 1; d <= rPhi; ++d) {
          phiCandidates.push_back({-d, wrapPhiBin(phiBinRing - d, nPhiRing)});
          phiCandidates.push_back({d, wrapPhiBin(phiBinRing + d, nPhiRing)});
        }
      }
      for (const auto& [dPhi, p] : phiCandidates) {
        for (const auto& [dRho, r] : rhoCandidates) {
          if (std::abs(dPhi) * wTheta * wRho + std::abs(dTheta) * wPhi * wRho + std::abs(dRho) * wPhi * wTheta >
              wPhi * wTheta * wRho)
            continue;
          CellID neighbour = base;
          decoder()->set(neighbour, m_phiIndex, p);
          decoder()->set(neighbour, m_thetaIndex, t);
          decoder()->set(neighbour, m_rhoIndex, r);
          // A phi window wider than the full circle wraps back onto the cell itself.
          if (neighbour != base)
            neighbours.insert(neighbour);
        }
      }
    }
  }

} // namespace DDSegmentation
} // namespace dd4hep
