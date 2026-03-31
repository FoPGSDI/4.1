//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file duck_star_pgen.cpp
//  \brief Problem generator for duck-shaped TOV star.
//
//  Based on dyngr_tov.cpp. Maps a spherical TOV profile onto a duck-shaped
//  boundary defined by spherical harmonic deformation coefficients epsilon_lm.
//
//  New parameters (in <problem> block):
//    duck_lmax           - maximum ell for SH expansion (default 10)
//    duck_epsilon_file   - path to file with epsilon_lm coefficients
//    duck_amplitude      - overall scaling of deformation (default 1.0)
//
//  The epsilon file format:
//    # comment lines start with #
//    # R_0 and R_eq are read from header comments of the form:
//    #   # R_0 = <value>
//    # ell  m  epsilon_real  epsilon_imag
//    1  -1  -6.273e-04  -6.244e-02
//    ...
//
//  At each grid point (x,y,z), compute:
//    r     = sqrt(x^2 + y^2 + z^2)
//    theta = acos(z/r)
//    phi   = atan2(y, x)
//    R_duck(theta,phi) = R_0 * [1 + sum epsilon_lm * Y_lm(theta,phi)]
//    r_eff = r * R_0 / R_duck(theta,phi)
//  Then query the TOV solution at r_eff instead of r.

#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/adm.hpp"
#include "z4c/z4c.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "utils/tov/tov.hpp"
#include "utils/tov/tov_polytrope.hpp"
#include "utils/tov/tov_tabulated.hpp"
#include "utils/tov/tov_piecewise_poly.hpp"

#include <Kokkos_Random.hpp>

//----------------------------------------------------------------------------------------
// Real spherical harmonics Y_l^m(theta, phi)
//
// Convention: real-valued basis
//   m > 0:  Y_l^m = sqrt(2) * N_l^m * P_l^m(cos theta) * cos(m phi)
//   m = 0:  Y_l^0 = N_l^0 * P_l^0(cos theta)
//   m < 0:  Y_l^m = sqrt(2) * N_l^|m| * P_l^|m|(cos theta) * sin(|m| phi)
// where N_l^m = sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!)
//
// Associated Legendre polynomials P_l^m(x) computed via stable recurrence.
//----------------------------------------------------------------------------------------

// Maximum supported ell
constexpr int DUCK_LMAX_MAX = 10;

// Total number of (ell, m) coefficients for ell = 0..lmax
constexpr int n_sh_coeffs(int lmax) {
  return (lmax + 1) * (lmax + 1);
}

// Index mapping: (ell, m) -> flat index, with m in [-ell, ell]
// index = ell^2 + ell + m
KOKKOS_INLINE_FUNCTION
int sh_index(int ell, int m) {
  return ell * ell + ell + m;
}

//----------------------------------------------------------------------------------------
// Associated Legendre polynomial P_l^m(x) (without Condon-Shortley phase)
// Uses stable upward recurrence in l at fixed m.
//----------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
Real assoc_legendre(int ell, int m, Real x) {
  // Ensure m >= 0
  int am = (m >= 0) ? m : -m;

  if (am > ell) return 0.0;

  // Compute P_m^m(x) via the diagonal formula
  Real pmm = 1.0;
  if (am > 0) {
    Real somx2 = sqrt((1.0 - x) * (1.0 + x));
    Real fact = 1.0;
    for (int i = 1; i <= am; i++) {
      pmm *= fact * somx2;  // No Condon-Shortley phase
      fact += 2.0;
    }
  }

  if (ell == am) return pmm;

  // Compute P_{m+1}^m(x)
  Real pmmp1 = x * (2.0 * am + 1.0) * pmm;
  if (ell == am + 1) return pmmp1;

  // Upward recurrence: (l-m) P_l^m = (2l-1) x P_{l-1}^m - (l+m-1) P_{l-2}^m
  Real pll = 0.0;
  for (int ll = am + 2; ll <= ell; ll++) {
    pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + am - 1.0) * pmm) / (ll - am);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

//----------------------------------------------------------------------------------------
// Normalization factor for real spherical harmonics
// N_l^m = sqrt((2l+1)/(4pi) * (l-|m|)!/(l+|m|)!)
//----------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
Real sh_norm(int ell, int am) {
  // Compute (l-|m|)! / (l+|m|)! by cancellation
  Real ratio = 1.0;
  for (int k = ell - am + 1; k <= ell + am; k++) {
    ratio *= static_cast<Real>(k);
  }
  return sqrt((2.0 * ell + 1.0) / (4.0 * M_PI * ratio));
}

//----------------------------------------------------------------------------------------
// Evaluate real spherical harmonic Y_l^m(theta, phi)
//----------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
Real real_Ylm(int ell, int m, Real theta, Real phi) {
  int am = (m >= 0) ? m : -m;
  Real ct = cos(theta);
  Real plm = assoc_legendre(ell, am, ct);
  Real norm = sh_norm(ell, am);

  if (m > 0) {
    return sqrt(2.0) * norm * plm * cos(m * phi);
  } else if (m == 0) {
    return norm * plm;
  } else {
    return sqrt(2.0) * norm * plm * sin(am * phi);
  }
}

//----------------------------------------------------------------------------------------
// Duck deformation data structure
// Stored on device for use in Kokkos kernels
//----------------------------------------------------------------------------------------

// Maximum number of coefficients: (LMAX+1)^2 = 121
constexpr int MAX_NCOEFFS = (DUCK_LMAX_MAX + 1) * (DUCK_LMAX_MAX + 1);

struct DuckDeformation {
  Real R0;                          // Reference radius
  int lmax;                         // Maximum ell
  Real amplitude;                   // Overall scaling factor
  Real eps_real[MAX_NCOEFFS];       // Real part of epsilon_lm
  Real eps_imag[MAX_NCOEFFS];       // Imaginary part of epsilon_lm
  // NOTE: For real SH basis, we only use eps_real.
  // eps_imag is stored for completeness but not used in R(theta,phi).
  // If the input file uses complex SH, the user must convert to real SH
  // coefficients before running.

  // Evaluate R(theta, phi) / R0
  KOKKOS_INLINE_FUNCTION
  Real surface_ratio(Real theta, Real phi) const {
    Real sum = 0.0;
    for (int ell = 0; ell <= lmax; ell++) {
      for (int m = -ell; m <= ell; m++) {
        int idx = sh_index(ell, m);
        sum += eps_real[idx] * real_Ylm(ell, m, theta, phi);
      }
    }
    return 1.0 + amplitude * sum;
  }
};

//----------------------------------------------------------------------------------------
// Read epsilon_lm from file (host-side)
//----------------------------------------------------------------------------------------
DuckDeformation ReadDuckEpsilon(const std::string& filename, int lmax, Real amplitude) {
  DuckDeformation duck;
  duck.lmax = lmax;
  duck.amplitude = amplitude;
  duck.R0 = 1.0;  // default, overridden from file header

  // Zero out arrays
  for (int i = 0; i < MAX_NCOEFFS; i++) {
    duck.eps_real[i] = 0.0;
    duck.eps_imag[i] = 0.0;
  }

  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cout << "### FATAL ERROR in duck_star_pgen.cpp" << std::endl
              << "Cannot open duck epsilon file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  while (std::getline(infile, line)) {
    // Check for header comments with R_0
    if (line[0] == '#') {
      if (line.find("R_0") != std::string::npos && line.find("=") != std::string::npos) {
        size_t eq_pos = line.find("=");
        duck.R0 = std::stod(line.substr(eq_pos + 1));
      }
      continue;
    }

    // Parse data line: ell  m  eps_real  eps_imag
    std::istringstream iss(line);
    int ell, m;
    Real er, ei;
    if (!(iss >> ell >> m >> er >> ei)) continue;

    if (ell > lmax) continue;
    if (ell < 0 || abs(m) > ell) continue;

    int idx = sh_index(ell, m);
    duck.eps_real[idx] = er;
    duck.eps_imag[idx] = ei;
  }
  infile.close();

  std::cout << "Duck deformation loaded: lmax=" << duck.lmax
            << ", R0=" << duck.R0
            << ", amplitude=" << duck.amplitude << std::endl;

  // Print a few dominant coefficients for verification
  for (int ell = 0; ell <= std::min(lmax, 4); ell++) {
    for (int m = -ell; m <= ell; m++) {
      int idx = sh_index(ell, m);
      if (fabs(duck.eps_real[idx]) > 1e-6 || fabs(duck.eps_imag[idx]) > 1e-6) {
        std::cout << "  eps(" << ell << "," << m << ") = "
                  << duck.eps_real[idx] << " + " << duck.eps_imag[idx] << "i"
                  << std::endl;
      }
    }
  }

  return duck;
}

//----------------------------------------------------------------------------------------
// Prototypes for vector potential (same as standard TOV)
//----------------------------------------------------------------------------------------
template<class TOVEOS>
KOKKOS_INLINE_FUNCTION
static Real A1(const tov::TOVStar& tov_, const TOVEOS& eos, bool isotropic, Real pcut,
               Real magindex, Real x1, Real x2, Real x3);
template<class TOVEOS>
KOKKOS_INLINE_FUNCTION
static Real A2(const tov::TOVStar& tov_, const TOVEOS& eos, bool isotropic, Real pcut,
               Real magindex, Real x1, Real x2, Real x3);

// Prototype for user-defined history
void DuckTOVHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
// Main setup function
//----------------------------------------------------------------------------------------
template<class TOVEOS>
void SetupDuckTOV(ParameterInput *pin, Mesh* pmy_mesh_) {
  Real v_pert = pin->GetOrAddReal("problem", "v_pert", 0.0);
  bool isotropic = pin->GetOrAddBoolean("problem", "isotropic", false);
  bool minkowski = pin->GetOrAddBoolean("problem", "minkowski", false);

  // Duck-specific parameters
  int duck_lmax = pin->GetOrAddInteger("problem", "duck_lmax", 10);
  std::string duck_file = pin->GetOrAddString("problem", "duck_epsilon_file",
                                               "duck_epsilon.dat");
  Real duck_amplitude = pin->GetOrAddReal("problem", "duck_amplitude", 1.0);

  if (duck_lmax > DUCK_LMAX_MAX) {
    std::cout << "### WARNING: duck_lmax=" << duck_lmax
              << " exceeds maximum " << DUCK_LMAX_MAX
              << ". Clamping to " << DUCK_LMAX_MAX << std::endl;
    duck_lmax = DUCK_LMAX_MAX;
  }

  // Read duck deformation from file
  DuckDeformation duck_host = ReadDuckEpsilon(duck_file, duck_lmax, duck_amplitude);

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

  // Construct the TOV solution
  TOVEOS eos{pin};
  auto my_tov = tov::TOVStar::ConstructTOV(pin, eos);

  constexpr bool use_ye = tov::UsesYe<TOVEOS>;
  Real ye_atmo = pin->GetOrAddReal("mhd", "s0_atmosphere", 0.5);

  auto& w0_ = pmbp->pmhd->w0;
  int& nvars_ = pmbp->pmhd->nmhd;
  int& nscal_ = pmbp->pmhd->nscalars;

  // Capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  int &ie = indcs.ie;
  int &je = indcs.je;
  int &ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;

  auto &size = pmbp->pmb->mb_size;
  auto &adm = pmbp->padm->adm;
  auto &tov_ = my_tov;
  auto &eos_ = eos;
  auto duck_ = duck_host;  // Copied by value into lambda capture

  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);

  //-------------------------------------------------------------------------------------
  // Main initialization loop: set hydro primitives and ADM variables
  //-------------------------------------------------------------------------------------
  par_for("pgen_duck_tov", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n2-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Spherical coordinates from Cartesian
    Real r = sqrt(SQR(x1v) + SQR(x2v) + SQR(x3v));
    Real s = sqrt(SQR(x1v) + SQR(x2v));

    // Compute angular coordinates
    Real theta = (r > 0.0) ? acos(x3v / r) : 0.0;
    Real phi = atan2(x2v, x1v);

    // ---- Duck deformation: compute effective radius ----
    // R(theta, phi) / R_0 = 1 + sum eps_lm * Y_lm(theta, phi)
    Real R_ratio = duck_.surface_ratio(theta, phi);

    // Guard against non-positive surface function
    if (R_ratio <= 0.0) R_ratio = 1.0e-10;

    // Effective radius: r_eff = r * R_0 / R_duck = r / R_ratio
    // This maps the duck boundary to the sphere boundary
    Real r_eff = r / R_ratio;

    // Query TOV solution at effective radius
    Real rho, p, mass, alp, r_schw;
    Real vr = 0.0;
    Real p_pert_val = 0.0;
    Real ye = ye_atmo;
    auto &use_ye_ = use_ye;

    if (!isotropic) {
      tov_.GetPrimitivesAtPoint(eos_, r_eff, rho, p, mass, alp);
      if (r_eff <= tov_.R_edge) {
        Real x = r_eff / tov_.R_edge;
        vr = 0.5 * v_pert * (3.0 * x - x * x * x);
        if constexpr (use_ye) {
          ye = eos_.template GetYeFromRho<tov::LocationTag::Device>(rho);
        }
      }
    } else {
      tov_.GetPrimitivesAtIsoPoint(eos_, r_eff, rho, p, mass, alp);
      r_schw = tov_.FindSchwarzschildR(r_eff, mass);
      if (r_schw <= tov_.R_edge) {
        Real x = r_schw / tov_.R_edge;
        vr = 0.5 * v_pert * (3.0 * x - x * x * x);
        if constexpr (use_ye) {
          ye = eos_.template GetYeFromRho<tov::LocationTag::Device>(rho);
        }
      }
    }

    // Set hydrodynamic quantities
    w0_(m, IDN, k, j, i) = rho;
    w0_(m, IPR, k, j, i) = p;
    w0_(m, IVX, k, j, i) = vr * x1v / fmax(r, 1.0e-30);
    w0_(m, IVY, k, j, i) = vr * x2v / fmax(r, 1.0e-30);
    w0_(m, IVZ, k, j, i) = vr * x3v / fmax(r, 1.0e-30);
    auto &nvars = nvars_;
    auto &nscal = nscal_;
    if (use_ye && nscal >= 1) {
      w0_(m, nvars, k, j, i) = ye;
    }

    // ---- ADM variables: use SPHERICAL TOV metric at physical radius r ----
    // This is the key simplification: the metric is that of the unperturbed
    // spherical star. The duck-shaped density perturbation introduces an
    // O(epsilon^2) Hamiltonian constraint violation that Z4c will damp.
    Real rho_sph, p_sph, mass_sph, alp_sph;
    if (!isotropic) {
      tov_.GetPrimitivesAtPoint(eos_, r, rho_sph, p_sph, mass_sph, alp_sph);
    } else {
      Real r_schw_sph;
      tov_.GetPrimitivesAtIsoPoint(eos_, r, rho_sph, p_sph, mass_sph, alp_sph);
      r_schw_sph = tov_.FindSchwarzschildR(r, mass_sph);
    }

    adm.alpha(m, k, j, i) = alp_sph;
    if (minkowski) {
      adm.g_dd(m,0,0,k,j,i) = adm.g_dd(m,1,1,k,j,i) = adm.g_dd(m,2,2,k,j,i) = 1.0;
      adm.g_dd(m,0,1,k,j,i) = adm.g_dd(m,0,2,k,j,i) = adm.g_dd(m,1,2,k,j,i) = 0.0;
      adm.alpha(m,k,j,i) = 1.0;
    } else if (!isotropic) {
      Real fmet = 0.0;
      if (r > 0) {
        fmet = (1.0 / (1.0 - 2.0 * mass_sph / r) - 1.0) / (r * r);
      }
      adm.g_dd(m,0,0,k,j,i) = x1v*x1v*fmet + 1.0;
      adm.g_dd(m,0,1,k,j,i) = x1v*x2v*fmet;
      adm.g_dd(m,0,2,k,j,i) = x1v*x3v*fmet;
      adm.g_dd(m,1,1,k,j,i) = x2v*x2v*fmet + 1.0;
      adm.g_dd(m,1,2,k,j,i) = x2v*x3v*fmet;
      adm.g_dd(m,2,2,k,j,i) = x3v*x3v*fmet + 1.0;
      Real det = adm::SpatialDet(
              adm.g_dd(m,0,0,k,j,i), adm.g_dd(m,0,1,k,j,i),
              adm.g_dd(m,0,2,k,j,i), adm.g_dd(m,1,1,k,j,i),
              adm.g_dd(m,1,2,k,j,i), adm.g_dd(m,2,2,k,j,i));
      adm.psi4(m,k,j,i) = pow(det, 1.0/3.0);
    } else {
      Real fmet = 1.0;
      if (r > 0) {
        Real r_schw_sph2;
        Real rho_tmp, p_tmp, mass_tmp, alp_tmp;
        tov_.GetPrimitivesAtIsoPoint(eos_, r, rho_tmp, p_tmp, mass_tmp, alp_tmp);
        r_schw_sph2 = tov_.FindSchwarzschildR(r, mass_tmp);
        fmet = r_schw_sph2 / r;
      }
      Real psi4 = fmet * fmet;
      adm.g_dd(m,0,0,k,j,i) = adm.g_dd(m,1,1,k,j,i) = adm.g_dd(m,2,2,k,j,i) = psi4;
      adm.g_dd(m,0,1,k,j,i) = adm.g_dd(m,0,2,k,j,i) = adm.g_dd(m,1,2,k,j,i) = 0.0;
      adm.psi4(m,k,j,i) = psi4;
    }
    adm.beta_u(m,0,k,j,i) = adm.beta_u(m,1,k,j,i) = adm.beta_u(m,2,k,j,i) = 0.0;
    adm.vK_dd(m,0,0,k,j,i) = adm.vK_dd(m,0,1,k,j,i) = adm.vK_dd(m,0,2,k,j,i) = 0.0;
    adm.vK_dd(m,1,1,k,j,i) = adm.vK_dd(m,1,2,k,j,i) = adm.vK_dd(m,2,2,k,j,i) = 0.0;
  });

  //-------------------------------------------------------------------------------------
  // Magnetic field setup (identical to standard TOV)
  //-------------------------------------------------------------------------------------
  Real b_norm = pin->GetOrAddReal("problem", "b_norm", 0.0);
  Real pcut = pin->GetOrAddReal("problem", "pcut", 1e-6);
  Real magindex = pin->GetOrAddReal("problem", "magindex", 2);

  if (pin->GetOrAddBoolean("problem", "use_pcut_rel", false)) {
    Real pmax = eos_.template GetPFromRho<tov::LocationTag::Host>(tov_.rhoc);
    pcut = pcut * pmax;
  }

  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb = pmbp->nmb_thispack;
  DvceArray4D<Real> a1, a2, a3;
  Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

  auto &nghbr = pmbp->pmb->nghbr;
  auto &mblev = pmbp->pmb->mb_lev;

  par_for("pgen_potential", DevExeSpace(), 0,nmb-1,ks,ke+1,js,je+1,is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1f = LeftEdgeX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2f = LeftEdgeX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    Real x3f = LeftEdgeX(k-ks, nx3, x3min, x3max);

    a1(m,k,j,i) = A1(tov_, eos_, isotropic, pcut, magindex, x1v, x2f, x3f);
    a2(m,k,j,i) = A2(tov_, eos_, isotropic, pcut, magindex, x1f, x2v, x3f);
    a3(m,k,j,i) = 0.0;

    // AMR correction for A1 at x2-faces, x3-faces, and x2x3-edges
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;

    if ((nghbr.d_view(m,8 ).lev > mblev.d_view(m) && j==js) ||
        (nghbr.d_view(m,9 ).lev > mblev.d_view(m) && j==js) ||
        (nghbr.d_view(m,10).lev > mblev.d_view(m) && j==js) ||
        (nghbr.d_view(m,11).lev > mblev.d_view(m) && j==js) ||
        (nghbr.d_view(m,12).lev > mblev.d_view(m) && j==je+1) ||
        (nghbr.d_view(m,13).lev > mblev.d_view(m) && j==je+1) ||
        (nghbr.d_view(m,14).lev > mblev.d_view(m) && j==je+1) ||
        (nghbr.d_view(m,15).lev > mblev.d_view(m) && j==je+1) ||
        (nghbr.d_view(m,24).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,25).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,26).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,27).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,28).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,29).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,30).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,31).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,40).lev > mblev.d_view(m) && j==js && k==ks) ||
        (nghbr.d_view(m,41).lev > mblev.d_view(m) && j==js && k==ks) ||
        (nghbr.d_view(m,42).lev > mblev.d_view(m) && j==je+1 && k==ks) ||
        (nghbr.d_view(m,43).lev > mblev.d_view(m) && j==je+1 && k==ks) ||
        (nghbr.d_view(m,44).lev > mblev.d_view(m) && j==js && k==ke+1) ||
        (nghbr.d_view(m,45).lev > mblev.d_view(m) && j==js && k==ke+1) ||
        (nghbr.d_view(m,46).lev > mblev.d_view(m) && j==je+1 && k==ke+1) ||
        (nghbr.d_view(m,47).lev > mblev.d_view(m) && j==je+1 && k==ke+1)) {
      Real xl = x1v + 0.25*dx1;
      Real xr = x1v - 0.25*dx1;
      a1(m,k,j,i) = 0.5*(A1(tov_, eos_, isotropic, pcut, magindex, xl,x2f,x3f) +
                         A1(tov_, eos_, isotropic, pcut, magindex, xr,x2f,x3f));
    }

    // AMR correction for A2 at x1-faces, x3-faces, and x1x3-edges
    if ((nghbr.d_view(m,0 ).lev > mblev.d_view(m) && i==is) ||
        (nghbr.d_view(m,1 ).lev > mblev.d_view(m) && i==is) ||
        (nghbr.d_view(m,2 ).lev > mblev.d_view(m) && i==is) ||
        (nghbr.d_view(m,3 ).lev > mblev.d_view(m) && i==is) ||
        (nghbr.d_view(m,4 ).lev > mblev.d_view(m) && i==ie+1) ||
        (nghbr.d_view(m,5 ).lev > mblev.d_view(m) && i==ie+1) ||
        (nghbr.d_view(m,6 ).lev > mblev.d_view(m) && i==ie+1) ||
        (nghbr.d_view(m,7 ).lev > mblev.d_view(m) && i==ie+1) ||
        (nghbr.d_view(m,24).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,25).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,26).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,27).lev > mblev.d_view(m) && k==ks) ||
        (nghbr.d_view(m,28).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,29).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,30).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,31).lev > mblev.d_view(m) && k==ke+1) ||
        (nghbr.d_view(m,32).lev > mblev.d_view(m) && i==is && k==ks) ||
        (nghbr.d_view(m,33).lev > mblev.d_view(m) && i==is && k==ks) ||
        (nghbr.d_view(m,34).lev > mblev.d_view(m) && i==ie+1 && k==ks) ||
        (nghbr.d_view(m,35).lev > mblev.d_view(m) && i==ie+1 && k==ks) ||
        (nghbr.d_view(m,36).lev > mblev.d_view(m) && i==is && k==ke+1) ||
        (nghbr.d_view(m,37).lev > mblev.d_view(m) && i==is && k==ke+1) ||
        (nghbr.d_view(m,38).lev > mblev.d_view(m) && i==ie+1 && k==ke+1) ||
        (nghbr.d_view(m,39).lev > mblev.d_view(m) && i==ie+1 && k==ke+1)) {
      Real xl = x2v + 0.25*dx2;
      Real xr = x2v - 0.25*dx2;
      a2(m,k,j,i) = 0.5*(A2(tov_, eos_, isotropic, pcut, magindex, x1f,xl,x3f) +
                         A2(tov_, eos_, isotropic, pcut, magindex, x1f,xr,x3f));
    }
  });

  auto &b0 = pmbp->pmhd->b0;
  par_for("pgen_Bfc", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real dx3 = size.d_view(m).dx3;

    b0.x1f(m,k,j,i) = b_norm*((a3(m,k,j+1,i) - a3(m,k,j,i))/dx2 -
                       (a2(m,k+1,j,i) - a2(m,k,j,i))/dx3);
    b0.x2f(m,k,j,i) = b_norm*((a1(m,k+1,j,i) - a1(m,k,j,i))/dx3 -
                       (a3(m,k,j,i+1) - a3(m,k,j,i))/dx1);
    b0.x3f(m,k,j,i) = b_norm*((a2(m,k,j,i+1) - a2(m,k,j,i))/dx1 -
                       (a1(m,k,j+1,i) - a1(m,k,j,i))/dx2);

    if (i==ie) {
      b0.x1f(m,k,j,i+1) = b_norm*((a3(m,k,j+1,i+1) - a3(m,k,j,i+1))/dx2 -
                           (a2(m,k+1,j,i+1) - a2(m,k,j,i+1))/dx3);
    }
    if (j==je) {
      b0.x2f(m,k,j+1,i) = b_norm*((a1(m,k+1,j+1,i) - a1(m,k,j+1,i))/dx3 -
                           (a3(m,k,j+1,i+1) - a3(m,k,j+1,i))/dx1);
    }
    if (k==ke) {
      b0.x3f(m,k+1,j,i) = b_norm*((a2(m,k+1,j,i+1) - a2(m,k+1,j,i))/dx1 -
                           (a1(m,k+1,j+1,i) - a1(m,k+1,j,i))/dx2);
    }
  });

  auto &bcc_ = pmbp->pmhd->bcc0;
  par_for("pgen_Bcc", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real& w_bx = bcc_(m,IBX,k,j,i);
    Real& w_by = bcc_(m,IBY,k,j,i);
    Real& w_bz = bcc_(m,IBZ,k,j,i);
    w_bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    w_by = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    w_bz = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Sets initial conditions for duck-shaped TOV star in DynGRMHD
//  Compile with '-D PROBLEM=duck_star' to enroll as user-specific problem generator

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (!pmbp->pcoord->is_dynamical_relativistic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Duck TOV star problem can only be run when <adm> block is present"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  user_hist_func = &DuckTOVHistory;

  if (restart) {
    return;
  }

  // Select the right TOV template based on the EOS
  if (pmbp->pdyngr->eos_policy == DynGRMHD_EOS::eos_ideal) {
    SetupDuckTOV<tov::PolytropeEOS>(pin, pmy_mesh_);
  } else if (pmbp->pdyngr->eos_policy == DynGRMHD_EOS::eos_compose) {
    SetupDuckTOV<tov::TabulatedEOS>(pin, pmy_mesh_);
  } else if (pmbp->pdyngr->eos_policy == DynGRMHD_EOS::eos_hybrid) {
    SetupDuckTOV<tov::TabulatedEOS>(pin, pmy_mesh_);
  } else if (pmbp->pdyngr->eos_policy == DynGRMHD_EOS::eos_piecewise_poly) {
    SetupDuckTOV<tov::PiecewisePolytropeEOS>(pin, pmy_mesh_);
  } else {
    std::cout << "### WARNING in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Unknown EOS requested for duck TOV star problem" << std::endl
              << "Defaulting to fixed polytropic EOS" << std::endl;
    SetupDuckTOV<tov::PolytropeEOS>(pin, pmy_mesh_);
  }

  // Mesh block info for loop limits
  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;

  // Convert primitives to conserved
  pmbp->pdyngr->PrimToConInit(0, (n1-1), 0, (n2-1), 0, (n3-1));

  if (pmbp->pz4c != nullptr) {
    switch (indcs.ng) {
      case 2: pmbp->pz4c->ADMToZ4c<2>(pmbp, pin);
              pmbp->pz4c->ADMConstraints<2>(pmbp);
              break;
      case 3: pmbp->pz4c->ADMToZ4c<3>(pmbp, pin);
              pmbp->pz4c->ADMConstraints<3>(pmbp);
              break;
      case 4: pmbp->pz4c->ADMToZ4c<4>(pmbp, pin);
              pmbp->pz4c->ADMConstraints<4>(pmbp);
              break;
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// Vector potential functions (identical to standard TOV)
//----------------------------------------------------------------------------------------
template<class TOVEOS>
KOKKOS_INLINE_FUNCTION
static Real A1(const tov::TOVStar& tov_, const TOVEOS& eos, bool isotropic, Real pcut,
               Real magindex, Real x1, Real x2, Real x3) {
  Real r = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
  Real p, rho;
  if (!isotropic) {
    tov_.GetPandRho(eos, r, rho, p);
  } else {
    tov_.GetPandRhoIso(eos, r, rho, p);
  }
  return -x2*fmax(p - pcut, 0.0)*pow(1.0 - rho/tov_.rhoc, magindex);
}

template<class TOVEOS>
KOKKOS_INLINE_FUNCTION
static Real A2(const tov::TOVStar& tov_, const TOVEOS& eos, bool isotropic, Real pcut,
               Real magindex, Real x1, Real x2, Real x3) {
  Real r = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
  Real p, rho;
  if (!isotropic) {
    tov_.GetPandRho(eos, r, rho, p);
  } else {
    tov_.GetPandRhoIso(eos, r, rho, p);
  }
  return x1*fmax(p - pcut, 0.0)*pow(1.0 - rho/tov_.rhoc, magindex);
}

//----------------------------------------------------------------------------------------
// History function: track rho_max and alpha_min (same as standard TOV)
//----------------------------------------------------------------------------------------
void DuckTOVHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 2;
  pdata->label[0] = "rho-max";
  pdata->label[1] = "alpha-min";

  auto &w0_ = pm->pmb_pack->pmhd->w0;
  auto &adm = pm->pmb_pack->padm->adm;

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  Real rho_max = std::numeric_limits<Real>::max();
  Real alpha_min = -rho_max;
  Kokkos::parallel_reduce("DuckTOVHistSums",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &mb_max, Real &mb_alp_min) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    mb_max = fmax(mb_max, w0_(m,IDN,k,j,i));
    mb_alp_min = fmin(mb_alp_min, adm.alpha(m, k, j, i));
  }, Kokkos::Max<Real>(rho_max), Kokkos::Min<Real>(alpha_min));

#if MPI_PARALLEL_ENABLED
  if (global_variable::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &rho_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &alpha_min, 1, MPI_ATHENA_REAL, MPI_MIN, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&rho_max, &rho_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&alpha_min, &alpha_min, 1, MPI_ATHENA_REAL, MPI_MIN, 0, MPI_COMM_WORLD);
    rho_max = 0.;
    alpha_min = 0.;
  }
#endif

  pdata->hdata[0] = rho_max;
  pdata->hdata[1] = alpha_min;
}
