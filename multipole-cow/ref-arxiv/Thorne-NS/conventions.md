# Conventions and Notation

## Paper: Thorne & Campolattaro (1967), ApJ 149, 591
### "Non-Radial Pulsation of General-Relativistic Stellar Models. I. Analytic Analysis for l >= 2"

---

## Units and Signature

- **Geometrized units**: c = G = 1
- **Metric signature**: (-,+,+,+) (timelike convention)
- **Follows**: Thorne (1966, 1967) "Relativistic Stellar Structure and Dynamics" (RSSD), but omits asterisks on geometrized quantities

## Coordinates

- Coordinate system: (t, r, theta, phi) — Schwarzschild-like coordinates
- Equilibrium line element (eq. 1):
  ds^2 = -e^nu dt^2 + e^lambda dr^2 + r^2(d theta^2 + sin^2 theta d phi^2)

## Metric Functions

| Symbol | Meaning |
|--------|---------|
| nu(r) | Gravitational potential: g_tt = -e^nu (replaces 2 Phi in RSSD) |
| lambda(r) | Radial metric function: e^{-lambda} = 1 - 2m/r (eq. 2) |
| m(r) | Mass inside radius r (eq. 3a) |
| rho(r) | Total mass-energy density |
| p(r) | Pressure |
| gamma | Adiabatic index = [(rho+p)/p](dp/drho)_{s,Z_i} (eq. 4) |

## Derivative Notation (defined on p. 596)

| Notation | Meaning |
|----------|---------|
| X' | Radial derivative: dX/dr |
| X'' | d^2 X/dr^2 |
| X_{,t} | Time derivative: dX/dt |
| X_{,tt} | d^2 X/dt^2 |
| X'_{,t} | Mixed derivative: d^2 X/dr dt |

## Spherical Harmonics

- Scalar: Y^l_M(theta, phi), even parity pi = (-1)^l
- Vector (even parity): Psi^l_{Mj} = d_j Y^l_M (eq. A1a)
- Vector (odd parity): Phi^l_{Mj} = epsilon_j^k d_k Y^l_M (eq. A1b)
- Tensor (even parity): Psi^l_{Mjk} = Y^l_{M|jk}, Phi^l_{Mjk} = gamma_{jk} Y^l_M (eq. A3a)
- Tensor (odd parity): chi^l_{Mjk} (eq. A3b)
- Antisymmetric tensor: epsilon_2^3 = -1/sin theta, epsilon_3^2 = sin theta (eq. A2)
- Unit sphere metric: gamma_{22} = 1, gamma_{33} = sin^2 theta (eq. A4)
- Covariant derivative slash notation: "|" denotes covariant differentiation w.r.t. gamma

## Perturbation Variables

### Gauge Choice
- **Regge-Wheeler gauge** throughout (Appendix A)
- Specialization to M = 0 (azimuthal index) and P_l(cos theta) (Legendre polynomials)

### Odd-Parity Perturbations (eq. 6)
| Variable | Role |
|----------|------|
| U(r,t) | Fluid displacement function |
| h_0(r,t) | Metric perturbation (t-phi component) |
| h_1(r,t) | Metric perturbation (r-phi component) |

### Even-Parity Perturbations (eq. 7)
| Variable | Role |
|----------|------|
| V(r,t) | Fluid displacement (angular, contravariant) |
| W(r,t) | Fluid displacement (radial, contravariant) |
| H_0(r,t) | Metric perturbation (tt component) |
| H_1(r,t) | Metric perturbation (tr component, "non-dynamical") |
| H_2(r,t) | Metric perturbation (rr component); H_2 = H_0 (eq. 8d) |
| K(r,t) | Metric perturbation (angular-angular component) |

## Eigenmode Decomposition

- Complex eigenfrequency (eq. 11): omega = sigma + i/tau
  - sigma: oscillation frequency (real part)
  - tau: damping time (positive = stable, radiating away energy)
- Convention (eq. 13): sigma >= 0
- Time dependence (eq. 12): H_0(r,t) = H(r) e^{i omega t}, etc.
  - H(r), K(r), W(r), V(r): complex radial eigenfunctions

## Key Equations

| Equation | Content |
|----------|---------|
| (1) | Equilibrium line element |
| (2) | e^{-lambda} = 1 - 2m/r |
| (3a-d) | TOV structure equations |
| (4) | Adiabatic index definition |
| (7b) | Even-parity metric perturbation matrix |
| (8a-d) | Initial-value equations for even-parity |
| (9a-c) | Propagation equations for even-parity |
| (14a-d) | Eigenequations (frequency-domain, main result) |
| (16) | Surface boundary condition: Delta p = 0 |
| (18) | Asymptotic wave form at r -> infinity |

## Boundary Conditions

### Center (r = 0, eq. 15a)
- K = A r^l + ..., H = A r^l + ..., W = B r^{l+1} + ..., V = -(B/l) r^l + ...
- Two free constants: A (metric) and B (fluid)
- Regularity: H_0 = K at r = 0 (eq. 10)

### Surface (r = R, eq. 15b)
- Regular series in (R - r)
- Four free constants: k_0, h_0, w_0, v_0
- Lagrangian pressure vanishes: Delta p = 0 (eq. 16)

### Infinity (r -> infinity, eq. 18)
- Incoming wave amplitude: C^{(I)}
- Outgoing wave amplitude: C^{(O)}
- Purely outgoing normal modes: C^{(I)} = 0 (QNMs)
- Standing waves: C^{(I)} = C^{(O)*} (eq. 31)

## Stability Criterion

- Stable QNM: omega_n^2 lies on positive real axis or upper half of complex plane
- Unstable QNM: omega_n^2 lies on negative real axis or lower half of complex plane
- Equivalently: tau_n > 0 means stable (damped); tau_n < 0 means unstable (growing)

## Physical Quantities

| Symbol | Meaning | Equation |
|--------|---------|----------|
| E_puls | Total pulsation energy | (29b) |
| E_kin | Kinetic energy | (25)-(26) |
| P | Radiated GW power | (30) |
| f(r) | Local observed frequency | (24) |
| M = m(R) | Total stellar mass | — |
| R | Stellar surface radius | — |

## Parity Classification

- **Even parity** (pi = (-1)^l, "electric-type"): causes density/pressure changes, produces pulsations and gravitational waves. This is the focus of the paper.
- **Odd parity** (pi = (-1)^{l+1}, "magnetic-type"): produces only stationary differential rotation and decoupled gravitational waves; no pulsation (Appendix B).
