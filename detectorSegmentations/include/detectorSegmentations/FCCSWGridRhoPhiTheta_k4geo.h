#ifndef DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H
#define DETECTORSEGMENTATIONS_GRIDRHOPHITHETA_K4GEO_H

// FCCSW
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

/** FCCSWGridRhoPhiTheta_k4geo Detector/detectorSegmentations/detectorSegmentations/FCCSWGridRhoPhiTheta_k4geo.h
 * FCCSWGridRhoPhiTheta_k4geo.h
 *
 *  Segmentation in rho, theta and phi.
 *  Based on FCCSWGridPhiTheta_k4geo, addition of Rho segmentation
 *
 */

namespace dd4hep {
namespace DDSegmentation {
  class FCCSWGridRhoPhiTheta_k4geo : public FCCSWGridPhiTheta_k4geo {
  public:
    /// default constructor using an arbitrary type
    FCCSWGridRhoPhiTheta_k4geo(const std::string& aCellEncoding);
    /// Default constructor used by derived classes passing an existing decoder
    FCCSWGridRhoPhiTheta_k4geo(const BitFieldCoder* decoder);

    /// destructor
    virtual ~FCCSWGridRhoPhiTheta_k4geo() = default;

    /**  Determine the global position based on the cell ID.
     *   @warning This segmentation has no knowledge of radius, so radius = 1 is taken into calculations.
     *   @param[in] aCellId ID of a cell.
     *   return Position (radius = 1).
     */
    virtual Vector3D position(const CellID& aCellID) const override;
    /**  Determine the cell ID based on the position.
     *   @param[in] aLocalPosition (not used).
     *   @param[in] aGlobalPosition position in the global coordinates.
     *   @param[in] aVolumeId ID of a volume.
     *   return Cell ID.
     */
    virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                          const VolumeID& aVolumeID) const override;

    /** Add the face-sharing neighbours in rho, phi and theta.
     *  Equivalent to neighboursInWindow(cellID, neighbours, 1, 1, 1).
     */
    virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

    /** Add the cells within a diamond-shaped window around cellID, i.e. bin
     *  offsets satisfying |dPhi|/rPhi + |dTheta|/rTheta + |dRho|/rRho <= 1
     *  (an axis with radius 0 is not varied).  With (1, 1, 1) this gives the
     *  face-sharing neighbours.  Phi neighbours wrap around at 2 pi; rho and
     *  theta do not.  The cell itself is never added.
     */
    void neighboursInWindow(const CellID& cellID, std::set<CellID>& neighbours, int rPhi, int rTheta, int rRho) const;

    /// First/last theta bin inside the physical coverage (and the theta field range).
    /// Barrel mode: [min_theta, max_theta]. Endcap mode (offset_Z > 0): the coverage is
    /// [min_theta, max_theta] U [pi - max_theta, pi - min_theta]; these return its outer limits,
    /// use thetaBinInCoverage() to test a bin.
    int firstThetaBin() const;
    int lastThetaBin() const;
    /// Whether a theta bin is inside the physical coverage (and the theta field range).
    bool thetaBinInCoverage(int bin) const;

    /** Number of phi bins in the ring of the given theta bin.
     *  Default: phi_bins for every ring.  With scale_phi_bins_with_sin_theta:
     *  max(1, round(phi_bins * sin(theta_c))), with theta_c the centre of the theta bin, i.e. phi_bins is the
     *  number of bins of the ring at theta = 90 deg and the phi cell width rho * sin(theta) * dPhi is about
     *  rho * 2 pi / phi_bins in every ring.
     *  offset_phi is the same absolute angle (centre of phi bin 0) in every ring.
     */
    int phiBinsInThetaBin(int thetaBin) const;
    /// Phi cell size in the ring of the given theta bin.
    inline double gridSizePhiInThetaBin(int thetaBin) const {
      return 2 * M_PI / static_cast<double>(phiBinsInThetaBin(thetaBin));
    }
    /// Azimuthal angle of the cell centre, using the phi grid of the cell's theta ring.
    /// Hides FCCSWGridPhiTheta_k4geo::phi(), which always uses phi_bins.
    double phi(const CellID cID) const;
    /// Cell dimensions (dPhi, dTheta), with dPhi of the cell's theta ring.
    virtual std::vector<double> cellDimensions(const CellID& cID) const override;

    /**  Find neighbours of the cell.
     *   Definition of neighbours is explained on slide 9:
     * https://indico.cern.ch/event/1475808/contributions/6219554/attachments/2966253/5218774/FCC_FullSim_HCal_slides.pdf
     *   @param[in] cID ID of a cell.
     *   @param[in] aDiagonal if true, will include neighbours from diagonal positions in the next and previous layers.
     *   return vector of neighbour cellIDs.
     */
    // std::vector<uint64_t> neighbours(const CellID cID, bool aDiagonal) const;

    /**  Find neighbours of the cell.
     *   Implement the signature from the Segmentation base class.
     */
    // virtual void neighbours(const CellID& cellID, std::set<CellID>& neighbours) const override;

  protected:
    /// Get rho from cellID
    double rho(const CellID cID) const;
    /// Distance from the IP, along the direction theta, to the front face of the detector:
    /// offset_R / sin(theta) for a cylinder (barrel), offset_Z / |cos(theta)| for a plane (endcap, offset_Z > 0).
    double rhoOffset(double theta) const;

  private:
    /// Clear all CellID fields except system, phi, theta and rho when enabled.
    void clearExtraFields(CellID& cellID) const;
    /// Map a phi bin onto the range [firstPhiBin, firstPhiBin + nPhi - 1] of a ring with nPhi phi bins.
    int wrapPhiBin(int bin, int nPhi) const;

    /// theta-dependent number of phi bins (see phiBinsInThetaBin()); false (default) keeps phi_bins everywhere
    bool m_scalePhiBinsWithSinTheta;

    /// the grid size in rho
    double m_grid_rho;
    std::vector<double> m_rhoBins;
    /// number of bins for a uniform rho grid (0 keeps the legacy field-range fallback)
    int m_rhoBinCount;

    /// physical theta coverage of the segmented detector
    double m_minTheta;
    double m_maxTheta;

    /// whether placement and other non-grid CellID fields are cleared
    bool m_clearExtraFields;

    /// the coordinate offset in rho / R
    double m_offsetR; // In cylinder case: offset depends on R not rho.
    /// endcap mode (> 0): |z| of the front plane from which rho is measured; 0 (default) = barrel mode with offset_R
    double m_offsetZ;

    /// the field name used for rho
    std::string m_rhoID;

    /// Initialization common to all ctors.
    void commonSetup();
    /// the field index used for theta
    int m_thetaIndex = -1;
    /// the field index used for phi
    int m_phiIndex = -1;
    /// the field index used for rho
    int m_rhoIndex = -1;
  };
} // namespace DDSegmentation
} // namespace dd4hep
#endif /* DETSEGMENTATION_GRIDRHOPHITHETA_H */
