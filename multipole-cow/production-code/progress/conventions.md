# Conventions — Quacky Normal Modes Production Code

All conventions used in the production code for the "Quacky Normal Modes" paper.
This document is the single source of truth; any discrepancy with code comments should be resolved in favor of this file.

---

## Coordinate Systems

**Cartesian (x, y, z):** Used for mesh geometry, quadrupole tensors, inertia tensors, and all FEM computations.
Origin at mesh center of mass.

**Spherical (r, theta, phi):** Used for multipole expansions, spherical harmonic decomposition, and surface deformation parametrization.
Standard physics convention:
- r >= 0 (radial distance)
- theta in [0, pi] (polar angle from +z axis)
- phi in [0, 2pi) (azimuthal angle from +x axis in the x-y plane)

Conversion:
```
x = r sin(theta) cos(phi)
y = r sin(theta) sin(phi)
z = r cos(theta)
```

---

## Spherical Harmonics

**Definition:** Complex spherical harmonics Y_l^m(theta, phi) with Condon-Shortley phase included.

**Implementation:** `scipy.special.sph_harm_y(l, m, theta, phi)`

Note: `sph_harm_y` uses the physics convention (theta = polar, phi = azimuthal).
The older `sph_harm(m, l, phi, theta)` has swapped argument order and is avoided.

**Orthonormality:**
```
integral Y_l^m* Y_l'^m' dOmega = delta_{ll'} delta_{mm'}
```

**Conjugation:**
```
Y_l^{-m} = (-1)^m (Y_l^m)*
```

---

## Metric Signature

**Signature:** (-,+,+,+)

Used throughout all GR sections (TOV, tidal Love numbers, QNMs, post-Newtonian dynamics, waveform generation).

---

## Units

**Geometric units (G = c = 1):** Used in all general-relativistic calculations:
- TOV structure equations
- Tidal Love number ODE (h_2 master equation)
- Post-Newtonian orbital dynamics
- Waveform phasing
- QNM frequencies (dimensionless omega M)

**CGS units:** Used for physical quantities and dimensional estimates:
- Density: g cm^{-3}
- Pressure: dyn cm^{-2} = g cm^{-1} s^{-2}
- Mass: g (or M_sun = 1.989e33 g)
- Length: cm (or km)
- GW strain: dimensionless
- GW luminosity: erg s^{-1}

**Conversion factors:**
```
G  = 6.674e-8   cm^3 g^{-1} s^{-2}
c  = 2.998e10    cm s^{-1}
M_sun = 1.989e33 g
R_sun = 6.957e10 cm
```

---

## Multipole Moments

### Spherical Multipole Moments

**Definition:**
```
Q_l^m = integral rho(x) r^l Y_l^m(theta, phi) d^3x
```

Volume integral over the body interior.  The density rho is the mass density (g cm^{-3} in CGS, or dimensionless in code units where total mass M = 1).

For a discrete mesh with uniform density:
```
Q_l^m = (M / V) sum_tetra integral_{tetra} r^l Y_l^m dV
```

### Cartesian Quadrupole Tensor

**Definition (traceless):**
```
Q^C_{ij} = integral rho (3 x_i x_j - r^2 delta_{ij}) d^3x
```

Properties:
- Symmetric: Q^C_{ij} = Q^C_{ji}
- Traceless: Q^C_{ii} = 0
- 5 independent components (maps to l=2 spherical moments)

### Inertia Tensor

**Definition:**
```
I_{ij} = integral rho (r^2 delta_{ij} - x_i x_j) d^3x
```

Properties:
- Symmetric: I_{ij} = I_{ji}
- Positive definite
- Tr(I) = 2 integral rho r^2 d^3x

### Cross-Check Relation

The Cartesian quadrupole and inertia tensor are related by:
```
Q^C_{ij} = Tr(I) delta_{ij} - 3 I_{ij}
```

This identity serves as a validation check in all numerical computations.
Equivalently: I_{ij} = (1/3)[Tr(I) delta_{ij} - Q^C_{ij}], with Tr(I) = -(1/2) Tr(Q^C) ... but since Q^C is traceless, Tr(I) must be computed independently.

---

## Deformation Coefficients

**Surface parametrization:**
```
R(theta, phi) = R_0 [1 + sum_{l>=1, m} epsilon_{lm} Y_l^m(theta, phi)]
```

where:
- R_0: mean radius from spherical harmonic fit (monopole term)
- epsilon_{lm}: complex deformation coefficients (dimensionless)
- The l=0 term is absorbed into R_0

**Extraction method:** Regularized least-squares fit of the mesh surface radii to the SH basis, implemented in `build_sh_matrix()` and `fit_sh_coefficients()`.

**Reality condition:** Since R(theta, phi) is real:
```
epsilon_{l,-m} = (-1)^m (epsilon_{lm})*
```

**Normalization convention:** The epsilon_{lm} are defined so that |epsilon_{lm}| ~ 0.1--0.3 for the duck at low l, with the dominant power at l = 2 (body shape) and l = 4 (neck/bill structure).

---

## Love Numbers

### Scalar Love Number k_2

**Hinderer formula:** Starting from the metric perturbation h_2(r) satisfying the master ODE, define:
```
y = R_* h_2'(R_*) / h_2(R_*)
```
at the stellar surface R_*.  Then:
```
k_2 = (8/5) C^5 (1 - 2C)^2 [2 + 2C(y-1) - y]
      / {2C [6 - 3y + 3C(5y-8)]
         + 4C^3 [13 - 11y + C(3y-2) + 2C^2(1+y)]
         + 3(1-2C)^2 [2 - y + 2C(y-1)] ln(1-2C)}
```
where C = M/R_* is the compactness (geometric units).

### Dimensionless Tidal Deformability

```
Lambda = (2/3) k_2 C^{-5}
```

### Dimensionful Tidal Deformability

```
lambda = (2/3) k_2 R_*^5 / G
```

---

## Wigner 3j Symbols

**Implementation:** `sympy.physics.wigner.wigner_3j(j1, j2, j3, m1, m2, m3)`

Returns exact rational values (as sympy `Rational` objects).  Convert to float for numerical use.

**Key identity used (Gaunt integral):**
```
integral Y_{l1}^{m1*} Y_{l2}^{m2} Y_{l3}^{m3} dOmega
  = (-1)^{m1} sqrt((2l1+1)(2l2+1)(2l3+1) / (4pi))
    * wigner_3j(l1, l2, l3, 0, 0, 0)
    * wigner_3j(l1, l2, l3, -m1, m2, m3)
```

**Selection rules enforced by the 3j symbols:**
1. Triangle inequality: |l1 - l2| <= l3 <= l1 + l2
2. m-selection: -m1 + m2 + m3 = 0
3. Parity: l1 + l2 + l3 must be even (for the (0,0,0) symbol)

---

## GW Radiation

**Quadrupole formula (leading order):**
```
P_GW = (G / (45 c^5)) <I_ij^(3) I_ij^(3)>
```
where I_ij^(3) is the third time derivative of the (traceless) mass quadrupole moment, and angle brackets denote time averaging.

**Multipole formula:**
```
P_{l,m} = (G / c^{2l+1}) omega^{2l+2} |delta Q_l^m|^2
          * 4pi (l+1)(l+2) / [l(l-1) ((2l+1)!!)^2]
```

---

## Numerical Parameters

**Duck mesh:** DASSL rubber duck, 4732 vertices, 9460 faces.

**SH expansion:** l_max = 10 (default), adjustable.  Higher l_max captures finer features (bill tip, tail curl) but increases computational cost quadratically.

**FEM eigenvalue solver:** `scipy.sparse.linalg.eigsh` with shift-invert mode, k = 50 eigenvalues.

**TOV reference model:** Gamma = 2 polytrope, K = 100, rho_c = 1.28e-3 (geometric units).
