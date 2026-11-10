//**************************************************************************
// \file    Grainita_ECAL_Endcaps_v02.cpp
// \brief:  Grainita ECAL endcaps geometry (development version).
//          Each endcap is a crystal layer divided into phi sectors and theta
//          modules (barrel-like construction, with optional tilt), with carbon
//          fiber plates on its front and back faces.
// \author: Fangyi Guo (fangyi.guo@cern.ch) start from Grainita_ECAL_Barrel_v02.cpp
// \date:   Nov. 2026
//**************************************************************************

// Includers from DD4hep
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"
#include <DD4hep/DetFactoryHelper.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

using namespace dd4hep;

// Build geometry
//
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  // A module side edge in the sector mid-plane for tilt in theta direction. 
  struct Edge {
    double x0;
    double slope;
    double at(double dz) const { return x0 + slope * dz; }
    double cosAngle() const { return 1. / std::sqrt(1. + slope * slope); }
    // Parallel edge whose face is displaced by the normal distance d (towards +x for d > 0)
    Edge shifted(double d) const { return Edge{x0 + d / cosAngle(), slope}; }
  };

  // A sector side wall in the transverse plane for tilt in phi direction. 
  // The wall is a plane parallel to the z axis.
  struct PhiLine {
    double k;
    double b;
    double y(double x) const { return k * x + b; }
    // Parallel line whose plane is displaced by the normal distance d (towards +y for d > 0)
    PhiLine shifted(double d) const { return PhiLine{k, b + d * std::sqrt(1. + k * k)}; }
  };


  std::cout << "--> Grainita_ECAL_Endcaps_v02::create_detector() start" << std::endl;

  // Get info from the xml file
  //
  // All geometry parameters are read from the <detector> element (no global constants):
  //   <dimensions zmin nphi ntheta/>
  //   <tilt angle/>   (optional, no tilt if absent)
  //   <crystal material thickness vis> <front rmin rmax/> <back rmin/> </crystal>
  //     front rmax: outer radius at zmin, defines the IP-pointing outer cone theta_out = atan(rmax / zmin)
  //   <backspace thickness/>   (optional, no back space if absent)
  //   <wall material thickness vis/>
  //   <plate material thickness vis/>
  xml_det_t x_det = e;
  std::string det_name = x_det.nameStr();
  int det_id = x_det.id();
  std::cout << "--> Going to create " << det_name << ", with ID: " << det_id << std::endl;

  xml_dim_t x_dim = x_det.dimensions();
  xml_comp_t x_crystal = x_det.child(_Unicode(crystal));
  xml_comp_t x_front = x_crystal.child(_Unicode(front));
  xml_comp_t x_back = x_crystal.child(_Unicode(back));
  xml_comp_t x_wall = x_det.child(_Unicode(wall));
  xml_comp_t x_plate = x_det.child(_Unicode(plate));

  // Define material
  Material air(description.material("Air"));
  Material MatCrystal(description.material(x_crystal.materialStr()));
  Material MatWall(description.material(x_wall.materialStr()));
  Material MatCarbonfiber(description.material(x_plate.materialStr()));

  // Detector global size. All radii are measured to the flat face (TGeoPgon convention).
  double zmin = x_dim.zmin();
  int nSec_phi = x_dim.attr<int>(_U(nphi));
  int nModule_theta = x_dim.attr<int>(_U(ntheta));
  double crystal_thick = x_crystal.thickness();
  double innerR_front = x_front.rmin();
  double outerR_front = x_front.rmax();
  double innerR_back = x_back.rmin();
  double Cseg_thick = x_wall.thickness();
  double Cframe_thick = x_plate.thickness();
  double module_tilt = x_det.hasChild(_Unicode(tilt)) ? xml_dim_t(x_det.child(_Unicode(tilt))).angle() : 0.;
  // Space behind the modules inside the sector, as the barrel BackSpace
  double back_space =
      x_det.hasChild(_Unicode(backspace)) ? xml_comp_t(x_det.child(_Unicode(backspace))).thickness() : 0.;
  if (back_space < 0)
    dd4hep::except("Grainita_ECAL_Endcaps_v02", "BackSpace (%g) must be >= 0", back_space);


  double dphi_sec = 2 * M_PI / nSec_phi;
  double phi_start = -dphi_sec / 2.; // flat face centred at phi = 0, as the barrel inner frame

  // |z| positions along the stack: front plate | sector (crystal layer + back space) | back plate
  // The "back" radii are those of the sector back face, adjacent to the back plate.
  double sector_depth = crystal_thick + back_space;
  double z_crystal_front = zmin + Cframe_thick;
  double z_sector_back = z_crystal_front + sector_depth;
  double zmax = z_sector_back + Cframe_thick;

  // Outer boundary: IP-pointing cone through the outer radius given at zmin (front face of the front plate).
  // Outer radii (flat faces) of the stack faces follow from it.
  double theta_out = atan(outerR_front / zmin);
  double outerR_sector_front = z_crystal_front * tan(theta_out);
  double outerR_sector_back = z_sector_back * tan(theta_out);

  if (innerR_front >= outerR_front || innerR_back >= outerR_sector_back)
    dd4hep::except("Grainita_ECAL_Endcaps_v02", "Inner radius must be smaller than outer radius (front %g/%g, back %g/%g)",
                   innerR_front, outerR_front, innerR_back, outerR_sector_back);

  std::cout << "  -- Dimension: |z| " << zmin << " - " << zmax << ", sector |z| " << z_crystal_front << " - "
            << z_sector_back << " (crystal " << crystal_thick << " + back space " << back_space << ")" << std::endl;
  std::cout << "     outer cone theta_out = " << theta_out << " rad; outer R at zmin " << outerR_front
            << ", sector front " << outerR_sector_front << ", sector back " << outerR_sector_back << std::endl;
  std::cout << "     inner R front " << innerR_front << ", back " << innerR_back << ", phi sectors " << nSec_phi
            << ", theta modules " << nModule_theta << std::endl;
  std::cout << "     tilt angle (deg): " << module_tilt * 180. / M_PI << std::endl;

  // ======== Phi tilt (same construction as the barrel) ======== //
  // In the transverse plane (sector local frame, x along the sector mid-line):
  //   inner vertices A(Rin, -dphi/2), D(Rin, +dphi/2) are kept on the regular polygon,
  //   outer vertices are rotated by +tilt: B(Rout_c, -dphi/2 + tilt), C(Rout_c, +dphi/2 + tilt),
  //   with Rout_c = Rout / cos(dphi/2) the polygon corner radius.
  // The side walls are the lines AB and DC (planes parallel to z), shared with the neighbour sectors.
  // The inner and outer faces stay at x = Rin and x = Rout (flat faces, as without tilt).
  // Rin / Rout: smallest inner and largest outer radius of the sector.
  double R_in_line = std::min(innerR_front, innerR_back);
  double R_out_line = outerR_sector_back;
  double R_out_corner = R_out_line / cos(dphi_sec / 2.);
  double y_A = R_in_line * tan(dphi_sec / 2.);
  double k_lo = -(R_out_corner * sin(dphi_sec / 2. - module_tilt) - y_A) /
                (R_out_corner * cos(dphi_sec / 2. - module_tilt) - R_in_line);
  double k_hi = (R_out_corner * sin(dphi_sec / 2. + module_tilt) - y_A) /
                (R_out_corner * cos(dphi_sec / 2. + module_tilt) - R_in_line);
  PhiLine wall_lo{k_lo, -y_A - R_in_line * k_lo}; // -y side
  PhiLine wall_hi{k_hi, y_A - R_in_line * k_hi};  // +y side
  double R_corner_max =
      std::max(std::hypot(R_out_line, wall_hi.y(R_out_line)), std::hypot(R_out_line, wall_lo.y(R_out_line)));
  std::cout << "  -- Phi walls: -y side k = " << wall_lo.k << ", b = " << wall_lo.b << "; +y side k = " << wall_hi.k
            << ", b = " << wall_hi.b << std::endl;
  std::cout << "     largest outer corner radius: " << R_corner_max << " (no tilt: " << R_out_corner << ")" << std::endl;

  // Crystal side walls in phi: displaced inwards by the normal wall thickness
  PhiLine wall_lo_c = wall_lo.shifted(Cseg_thick);
  PhiLine wall_hi_c = wall_hi.shifted(-Cseg_thick);

  // ======== Theta division ======== //
  // Theta range covering the whole sector, from its front and back faces.
  double theta_min = std::min(atan(innerR_front / z_crystal_front), atan(innerR_back / z_sector_back));
  double theta_max = theta_out;
  double theta_module = (theta_max - theta_min) / nModule_theta;
  std::cout << "  -- Theta range: " << theta_min << " - " << theta_max << " rad, step " << theta_module << std::endl;

  // Inner boundary follows the given radii; outer boundary is the IP-pointing cone theta_out (not tilted).
  // Internal theta walls: front end on the IP ray theta_i at the sector front face, back end on the
  // IP ray (theta_i - tilt) at the sector back face (as the barrel), i.e. tilted towards larger |cos(theta)|.
  Edge inner_edge{innerR_front, (innerR_back - innerR_front) / sector_depth};
  Edge outer_edge{outerR_sector_front, tan(theta_out)};
  auto theta_wall = [&](int i) {
    double theta = theta_min + i * theta_module;
    double xw_front = z_crystal_front * tan(theta);
    double xw_back = z_sector_back * tan(theta - module_tilt);
    return Edge{xw_front, (xw_back - xw_front) / sector_depth};
  };

  // Create the geometry
  DetElement ECAL(det_name, det_id);
  dd4hep::xml::setDetectorTypeFlag(e, ECAL);
  Volume worldVol = description.pickMotherVolume(ECAL);
  Assembly EcalEndcapsVol("EcalEndcaps");

  // ====== Carbon fiber plates on the front and back faces, shared by both endcaps ====== //
  // Plate outer radii: the outer cone at the plate face closer to the IP, so the plates stay inside the cone.
  PolyhedraRegular frontPlate(nSec_phi, phi_start, innerR_front, outerR_front, Cframe_thick);
  Volume frontPlate_vol("EcalEndcapFrontCarbonPlate", frontPlate, MatCarbonfiber);
  frontPlate_vol.setVisAttributes(description, x_plate.visStr());

  PolyhedraRegular backPlate(nSec_phi, phi_start, innerR_back, outerR_sector_back, Cframe_thick);
  Volume backPlate_vol("EcalEndcapBackCarbonPlate", backPlate, MatCarbonfiber);
  backPlate_vol.setVisAttributes(description, x_plate.visStr());

  // Sensitivity is optional: it is only set when the compact file provides a readout.
  if (sens.isValid())
    sens.setType("calorimeter");

  // Fill an EightPointSolid vertex array from the module (or sector) edges.
  // Local frame: x along the sector mid-line, y transverse, z along the global z axis.
  // The face closer to the IP is at local -dz on the +z side and at local +dz on the -z side.
  // The phi side faces are the planes y = lo.y(x) and y = hi.y(x).
  auto fill_vertices = [&](double* v, int side, double x_lo_front, double x_hi_front, double x_lo_back,
                           double x_hi_back, const PhiLine& lo, const PhiLine& hi) {
    auto face = [&](double* f, double x_lo, double x_hi) {
      double tmp[8] = {x_hi, hi.y(x_hi), x_hi, lo.y(x_hi), x_lo, lo.y(x_lo), x_lo, hi.y(x_lo)};
      std::copy(tmp, tmp + 8, f);
    };
    face(side > 0 ? v : v + 8, x_lo_front, x_hi_front);
    face(side > 0 ? v + 8 : v, x_lo_back, x_hi_back);
  };

  // Print the vertices of an EightPointSolid in the sector local frame (mm), IP-side face first.
  auto print_vertices = [&](const std::string& label, const double* v, int side, double z_front, double z_back) {
    std::cout << "  " << label << " vertices [x, y] (mm), R (mm), phi (deg):" << std::endl;
    for (int iface = 0; iface < 2; iface++) {
      const double* f = ((side > 0) == (iface == 0)) ? v : v + 8;
      std::cout << "    " << (iface == 0 ? "front" : "back") << " face, |z| = " << (iface == 0 ? z_front : z_back) / mm
                << std::endl;
      for (int i = 0; i < 8; i += 2)
        std::cout << "      [" << f[i] / mm << ", " << f[i + 1] / mm << "], R = " << std::hypot(f[i], f[i + 1]) / mm
                  << ", phi = " << atan2(f[i + 1], f[i]) * 180. / M_PI << std::endl;
    }
  };

  //========  Build and place the two endcaps (side = +1 for +z, -1 for -z)  ========  //
  for (int side : {+1, -1}) {
    std::string side_name = (side > 0) ? "pos" : "neg";

    // ------ Sector: one phi slice of the crystal layer (air), modules placed inside ------ //
    double vertices_sector[16];
    fill_vertices(vertices_sector, side, innerR_front, outerR_sector_front, innerR_back, outerR_sector_back, wall_lo,
                  wall_hi);
    EightPointSolid Sector("EcalEndcapSector_" + side_name, sector_depth / 2., vertices_sector);
    Volume Sector_vol("EcalEndcapSector_" + side_name, Sector, air);
    Sector_vol.setVisAttributes(description, x_det.visStr());
    if (side > 0) {
      std::cout << "Sector (side +z, sector 0):" << std::endl;
      print_vertices("Sector", vertices_sector, side, z_crystal_front, z_sector_back);
    }

    // ------ Modules along theta ------ //
    for (int it = 0; it < nModule_theta; it++) {
      Edge lo = (it == 0) ? inner_edge : theta_wall(it);
      Edge hi = (it == nModule_theta - 1) ? outer_edge : theta_wall(it + 1);

      // Module length along z:
      //  - innermost module (inner face parallel to the z axis): the crystal thickness;
      //  - other modules: every ray from the IP crossing the module has a crystal path >= crystal thickness
      //    (as the barrel). The smallest IP angle in the module is at the back end of its lower edge,
      //    hence h = thick * cos(atan(x_lo(h) / (z_front + h))), solved by fixed-point iteration.
      double height_module = crystal_thick;
      if (it > 0) {
        for (int iter = 0; iter < 100; iter++) {
          double h_new = crystal_thick * cos(atan(lo.at(height_module) / (z_crystal_front + height_module)));
          bool converged = std::fabs(h_new - height_module) < 1e-12 * crystal_thick;
          height_module = h_new;
          if (converged)
            break;
        }
      }

      double x_lo_front = lo.at(0.), x_hi_front = hi.at(0.);
      double x_lo_back = lo.at(height_module), x_hi_back = hi.at(height_module);

      if (side > 0) {
        std::cout << "Module #" << it << ": theta " << theta_min + it * theta_module << " - "
                  << theta_min + (it + 1) * theta_module << ", height " << height_module << ", x front "
                  << x_lo_front << " - " << x_hi_front << ", x back " << x_lo_back << " - " << x_hi_back
                  << std::endl;
      }

      // Crystal: all four side faces displaced inwards by the normal wall thickness Cseg_thick.
      // Front and back faces are not covered by the wall (carbon plate in front, as the barrel).
      Edge lo_c = lo.shifted(Cseg_thick);
      Edge hi_c = hi.shifted(-Cseg_thick);
      double xc_lo_front = lo_c.at(0.), xc_hi_front = hi_c.at(0.);
      double xc_lo_back = lo_c.at(height_module), xc_hi_back = hi_c.at(height_module);

      // Crystal must keep a positive size after the wall inset.
      if (xc_hi_front <= xc_lo_front || xc_hi_back <= xc_lo_back ||
          wall_hi_c.y(xc_lo_front) <= wall_lo_c.y(xc_lo_front) || wall_hi_c.y(xc_lo_back) <= wall_lo_c.y(xc_lo_back))
        dd4hep::except("Grainita_ECAL_Endcaps_v02",
                       "Module %d has no room for the crystal: x front %g - %g, x back %g - %g, wall %g", it,
                       x_lo_front, x_hi_front, x_lo_back, x_hi_back, Cseg_thick);

      double vertices_frame[16], vertices_crystal[16];
      fill_vertices(vertices_frame, side, x_lo_front, x_hi_front, x_lo_back, x_hi_back, wall_lo, wall_hi);
      fill_vertices(vertices_crystal, side, xc_lo_front, xc_hi_front, xc_lo_back, xc_hi_back, wall_lo_c, wall_hi_c);
      if (side > 0) {
        print_vertices("Frame", vertices_frame, side, z_crystal_front, z_crystal_front + height_module);
        print_vertices("Crystal", vertices_crystal, side, z_crystal_front, z_crystal_front + height_module);
      }

      std::string module_name = "EcalEndcapModule_" + side_name + "_" + std::to_string(it);
      EightPointSolid ModuleFrame(module_name, height_module / 2., vertices_frame);
      Volume ModuleFrame_vol(module_name, ModuleFrame, MatWall);
      ModuleFrame_vol.setVisAttributes(description, x_wall.visStr());

      std::string crystal_name = "EcalEndcapCrystal_" + side_name + "_" + std::to_string(it);
      EightPointSolid ModuleCrystal(crystal_name, height_module / 2., vertices_crystal);
      Volume ModuleCrystal_vol(crystal_name, ModuleCrystal, MatCrystal);
      ModuleCrystal_vol.setVisAttributes(description, x_crystal.visStr());
      if (sens.isValid())
        ModuleCrystal_vol.setSensitiveDetector(sens);

      PlacedVolume pv = ModuleFrame_vol.placeVolume(ModuleCrystal_vol);
      pv.addPhysVolID("stave", it);

      // Module front face aligned with the sector front face.
      Sector_vol.placeVolume(ModuleFrame_vol, Position(0., 0., -side * (sector_depth - height_module) / 2.));
    }

    // ------ Place plates and sectors in this endcap ------ //
    Assembly endcap_vol("EcalEndcap_" + side_name);
    endcap_vol.placeVolume(frontPlate_vol, Position(0., 0., side * (zmin + Cframe_thick / 2.)));
    endcap_vol.placeVolume(backPlate_vol, Position(0., 0., side * (z_sector_back + Cframe_thick / 2.)));

    for (int isec_phi = 0; isec_phi < nSec_phi; isec_phi++) {
      double phi = isec_phi * dphi_sec;
      Transform3D tr(RotationZ(phi), Translation3D(0., 0., side * (z_crystal_front + sector_depth / 2.)));
      PlacedVolume sector_pv = endcap_vol.placeVolume(Sector_vol, tr);
      sector_pv.addPhysVolID("sector", isec_phi);
    }

    PlacedVolume endcap_pv = EcalEndcapsVol.placeVolume(endcap_vol);
    endcap_pv.addPhysVolID("side", side);
    DetElement endcap_det(ECAL, "EcalEndcap_" + side_name, side);
    endcap_det.setPlacement(endcap_pv);
  }

  // Finalize geometry
  PlacedVolume EcalEndcaps_plv = worldVol.placeVolume(EcalEndcapsVol);
  EcalEndcaps_plv.addPhysVolID("system", det_id);
  ECAL.setPlacement(EcalEndcaps_plv);

  std::cout << "--> Grainita_ECAL_Endcaps_v02::create_detector() end" << std::endl;
  return ECAL;
}

DECLARE_DETELEMENT(Grainita_ECAL_Endcaps_v02, create_detector)

//**************************************************************************
