# Mathematical Conventions

Single source of truth for all mathematical conventions used in the reproduction
of "Higher multipoles of the cow" (Lehmann 2025, arXiv:2504.00506).

---

## 1. Coordinate System ("Cow Coordinates")

- **Handedness:** Right-handed Cartesian
- **Axes:** x = forward, y = up, z = right
- **Origin:** Center of mass (volumetric centroid) of the benchmark cow mesh
- **Benchmark units (BU):** Bounding box extents are (1.044, 0.6397, 0.3403) in (x, y, z)

> **IMPORTANT (mesh origin):** The raw mesh may NOT be centered at the center of
> mass. All coordinates must be translated to place the volumetric centroid at
> the origin before computing any moments.

---

## 2. Spherical Coordinates

| Symbol | Definition | Range |
|--------|-----------|-------|
| r | sqrt(x^2 + y^2 + z^2) | [0, inf) |
| theta | arccos(z / r) — polar angle (colatitude) from **z-axis** | [0, pi] |
| phi | arctan2(y, x) — azimuthal angle in **xy-plane** | [0, 2pi) |

For a Cartesian point (x, y, z):

    r     = sqrt(x**2 + y**2 + z**2)
    theta = arccos(z / r)
    phi   = arctan2(y, x)

> **WARNING:** The cow's vertical direction (y) is NOT the spherical polar axis.
> The polar axis is z (cow's right-hand direction).

> **[BLOCKING AMBIGUITY] Polar axis choice:**
> The paper does not explicitly state which Cartesian axis is the spherical
> polar axis. The convention above (polar axis = z) must be verified
> empirically by computing Q_l^m for all three axis choices (x, y, z as
> polar) and matching against the paper's tabulated values. Until this
> verification is complete, treat the polar axis assignment as provisional.

---

## 3. Spherical Harmonics Y_l^m

**Convention:** Condon-Shortley phase is INCLUDED.

**Normalization:** Orthonormal on S^2:

    integral |Y_l^m|^2 dOmega = 1

**scipy call (scipy >= 1.17):**

```python
from scipy.special import sph_harm_y
Y = sph_harm_y(n, m, theta, phi)   # n = l, theta = polar, phi = azimuthal
```

> **WARNING:** The old function `scipy.special.sph_harm(m, l, phi, theta)` has
> been REMOVED in scipy >= 1.17. Do NOT use it. Note the different argument
> order and convention in the old API.

**Conjugation symmetry:**

    Y_l^{-m} = (-1)^m conj(Y_l^m)

---

## 4. Multipole Moments Q_l^m

**Definition (interior, solid harmonics):**

    Q_l^m = integral_C d^3x rho(x) ||x||^l Y_l^m(x_hat)

Key points:
- Uses Y_l^m **directly** (NOT the complex conjugate)
- rho = 1 (uniform density)
- Computed about the center of mass, so the dipole (l = 1) vanishes

**Conjugation symmetry:**

    Q_l^{-m} = (-1)^m conj(Q_l^m)

---

## 5. Cartesian Quadrupole Q^C_ij

**Definition:**

    Q^C_ij = integral d^3x rho(x) (3 x_i x_j - r^2 delta_ij)

Key points:
- **No 1/2 factor** in front
- Traceless by construction: Tr(Q^C) = 0
- Symmetric: Q^C_ij = Q^C_ji

---

## 6. Inertia Tensor I_ij

**Definition:**

    I_ij = integral d^3x rho(x) (r^2 delta_ij - x_i x_j)

**Cross-check relation (verified algebraically from paper's numerical values):**

    Q^C_ij = Tr(I) delta_ij - 3 I_ij

---

## 7. Surface Map Decomposition

**Map F:** S^2 -> boundary(C) (sphere to cow boundary).

**Components:** f = (f_r, f_Delta_theta, f_Delta_phi)

| Component | Meaning |
|-----------|---------|
| f_r(Omega) | Radial distance of F(Omega) from origin |
| f_Delta_theta(theta, phi) | [F(theta, phi)]_theta - theta - const (monopole removed) |
| f_Delta_phi(theta, phi) | [F(theta, phi)]_phi - phi - const (monopole removed) |

**SH coefficients (uses CONJUGATE):**

    c_l^m = <f | Y_l^m> = integral f(Omega) conj(Y_l^m(Omega)) dOmega

> Note the distinction: multipole moments Q_l^m use Y_l^m directly, while
> surface map coefficients c_l^m use conj(Y_l^m).

---

## 8. Physical Units

| Quantity | Value | Notes |
|----------|-------|-------|
| Benchmark length unit L | 2.39 m | Cow length ~2.5 m; benchmark x-extent = 1.044 |
| Benchmark mass unit M | 13652 kg | From M/L^3 = 1 g/cm^3 |
| Cow density rho_physical | 1 g/cm^3 | Uniform throughout |

---

## 9. Numerical Benchmarks

These values from the paper serve as validation targets:

| Quantity | Value | Notes |
|----------|-------|-------|
| Q_0^0 | 0.0539 | Monopole moment |
| Volume | Q_0^0 * sqrt(4 pi) ~ 0.191 BU^3 | In benchmark units |
| Q_2^0 | -0.0029 | Quadrupole moment |
| <Q_triple_dot_ij Q_triple_dot^ij> | 0.00149 omega^6 | For y-axis rotation, benchmark units |

---

## 10. Known Ambiguities

### 10.1 [BLOCKING] Polar axis assignment

The paper does not state which Cartesian axis is the spherical polar axis.
**Resolution required:** Compute Q_l^m for each of the three possible polar
axis choices (x, y, z) and compare against the paper's tabulated values.
Mark all downstream results as provisional until resolved.

### 10.2 Mesh origin

The raw mesh coordinates are not guaranteed to be centered at the center of
mass. **Resolution:** Compute the volumetric centroid and translate all
vertices before computing any moments.

### 10.3 Table assignment (Section IV)

The paper's Section IV contains two tables of surface map coefficients:
- **First table** (f_r monopole = 1.0750): gradient flow map
- **Second table** (f_r monopole = 0.5825): harmonic map

Do not confuse the two when validating.
