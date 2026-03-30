# Quacky Normal Modes

## Paper Overview

**Premise (humorous):** What are the oscillation modes, tidal response, and gravitational wave signatures of a duck-shaped compact object? Can current GW observations rule out the duck hypothesis?

**Physics (rigorous):** We develop the general theory of oscillation modes and tidal deformability for non-spherical self-gravitating bodies, using a rubber duck as the worked example. The duck shape breaks rotational symmetry, lifting QNM degeneracies (gravitational Zeeman splitting), generating anisotropic Love numbers (the Love tensor), and modifying binary inspiral dynamics. We validate the formalism with numerical relativity simulations in AthenaK, culminating in the visual centerpiece: a duck-shaped TOV star melting to a sphere under its own gravity. We show that current PTA and LIGO observations cannot distinguish duck-shaped compact objects from spherical ones.

---

## Paper Structure

The paper follows a theory → numerics → observation arc:

### Part I: General Formalism (Sections 2–4)
- **Section 2:** Duck geometry, multipole characterization, deformation spectrum ε_ℓm
- **Section 3:** Quacky normal modes — perturbative splitting theory (Hadamard, Wigner 3j, Grotrian diagrams)
- **Section 4:** Tidal Love numbers — Love tensor, direction-dependent k_2(θ,φ), I-Love-Q breakdown, perturbative Λ correction
- **Section 5:** Binary duck PN formalism — quadrupole-monopole interaction, modified Kepler, system quadrupole, phase evolution

### Part II: Numerical Simulations (Sections 6–7)
- **Section 6:** QNM numerical results — FEM Helmholtz eigenvalues, cross-validation with perturbation theory
- **Section 7:** AthenaK numerical relativity — duck-shaped TOV star, evolution, duck melts to sphere, Ψ₄ extraction, waveform generation

### Part III: Experimental Verifications (Section 8)
- **Section 8:** Confrontation with GW data
  - 8.1: PTA band (supermassive ducks) — NANOGrav 15yr comparison, mock skymap, population degeneracy theorem
  - 8.2: LIGO band (stellar-mass ducks) — tidal deformability from inspiral, QNM spectroscopy from ringdown

### Section 1: Introduction (~2 pages)
### Section 9: Discussion and Conclusions (~1 page)

**Total: ~25 pages + figures/tables**

---

# Part I: General Formalism

---

## Section 2: Duck Geometry and Multipole Characterization

**Duck mesh:** DASSL rubber duck, 4732 vertices, 9460 faces, converted from OBJ via trimesh, normalized to unit scale.

**Multipole characterization** (existing validated pipeline in `multipole-cow/code/`):
- Q^C_duck: Cartesian quadrupole tensor (3×3, traceless, symmetric)
- I_duck: inertia tensor (3×3, symmetric), principal moments I_1, I_2, I_3
- Q_ℓ^m: spherical multipole moments up to ℓ=10
- R_eq = (3V/4π)^{1/3}: equivalent sphere radius
- κ_duck = Q^C_eigenvalue/(M R_eq²): shape-induced quadrupole parameter

**Surface deformation spectrum:** Parametrize the duck boundary as a deformed sphere:

$$R_{\rm duck}(\theta,\phi) = R_0\left[1 + \sum_{\ell \geq 1, m} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right]$$

The ε_ℓm are extracted via regularized least-squares SH fitting (existing code: `render_sphere_to_cow.py` → `build_sh_matrix()`, `fit_sh_coefficients()`).

**Verification:** Dipole vanishes at COM; Q^C traceless and symmetric; cross-check Q^C = Tr(I)δ − 3I.

**Figures/Tables:**
- Fig 1: Duck mesh and SH reconstruction at ℓ_max = 0, 2, 4, 8, full (existing renders)
- Fig 2: Deformation spectrum |ε_ℓm| vs ℓ (bar chart)
- Tab 1: Duck geometric and multipole parameters
- Tab 2: Deformation coefficients ε_ℓm up to ℓ=10

---

## Section 3: Quacky Normal Modes — Perturbative Splitting Theory

### 3.1 Problem Setup

The unperturbed problem: eigenvalue equation on a sphere of radius R_0. For the scalar Helmholtz model:

$$-\nabla^2 \psi_{n\ell m} = \omega_{n\ell}^2\,\psi_{n\ell m}, \quad \psi\big|_{r=R_0} = 0$$

Eigenvalues ω_{n,ℓ} = j_{n,ℓ}/R_0 (zeros of spherical Bessel j_ℓ). Each (n,ℓ) level has (2ℓ+1)-fold degeneracy.

For fluid modes (Cowling approximation): f-mode and p-mode frequencies from tabulated values (Kokkotas & Schmidt 1999).

### 3.2 Hadamard Formula for Domain Perturbation

The sensitivity of Dirichlet eigenvalues to boundary deformation δn(x):

$$\delta\omega_{n\ell m}^2 = -\int_{\partial D} \left|\frac{\partial \psi_{n\ell m}}{\partial n}\right|^2 \delta n\,dS$$

For the duck deformation δn = R_0 Σ ε_ℓ'm' Y_ℓ'^{m'}:

$$\delta\omega_{n\ell m}^2 = -R_0 \sum_{\ell',m'} \varepsilon_{\ell' m'} \int_{S^2} \left|\frac{\partial \psi_{n\ell m}}{\partial r}\bigg|_{R_0}\right|^2 Y_{\ell'}^{m'}(\theta,\phi)\,R_0^2\,d\Omega$$

### 3.3 Degenerate Perturbation Theory

Within each (n,ℓ) multiplet, diagonalize the (2ℓ+1) × (2ℓ+1) perturbation matrix:

$$V_{m_1, m_2}^{(n,\ell)} = -|R_{n\ell}'(R_0)|^2 R_0^3 \sum_{\ell',m'} \varepsilon_{\ell' m'}\int_{S^2} Y_\ell^{m_1*}\,Y_\ell^{m_2}\,Y_{\ell'}^{m'}\,d\Omega$$

The angular integral is the Gaunt integral, expressed via Wigner 3j symbols:

$$\int_{S^2} Y_\ell^{m_1*}\,Y_\ell^{m_2}\,Y_{\ell'}^{m'}\,d\Omega = (-1)^{m_1}\sqrt{\frac{(2\ell+1)^2(2\ell'+1)}{4\pi}} \begin{pmatrix} \ell & \ell & \ell' \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} \ell & \ell & \ell' \\ -m_1 & m_2 & m' \end{pmatrix}$$

### 3.4 Selection Rules

The 3j symbols enforce:
1. **Triangle inequality:** 0 ≤ ℓ' ≤ 2ℓ
2. **Parity:** ℓ' must be even
3. **m-selection:** m' = m_1 − m_2

**Consequence:** Only **even-ℓ' deformation harmonics** contribute:
- ε_{0,0}: overall shift (no splitting)
- ε_{2,m}: primary splitting (quadrupole deformation)
- ε_{4,m}: secondary splitting (hexadecapole)

### 3.5 Explicit ℓ=2 f-mode Splitting

The 5-fold degenerate ℓ=2 f-mode splits via the 5×5 perturbation matrix:

$$V_{m_1,m_2} = -|R_{n,2}'(R_0)|^2 R_0^3 \sum_{m'} \varepsilon_{2,m'}\,(-1)^{m_1}\sqrt{\frac{25 \cdot 5}{4\pi}} \begin{pmatrix} 2 & 2 & 2 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} 2 & 2 & 2 \\ -m_1 & m_2 & m' \end{pmatrix}$$

with (2 2 2 | 0 0 0) = √(2/35). Diagonalization gives 5 split frequencies.

### 3.6 Physical Interpretation

Which duck features drive the splitting:
- **ε_{2,0}:** axial elongation (body shape) → splits m=0 from m=±1, ±2
- **ε_{2,±1}:** bill-to-tail tilt asymmetry
- **ε_{2,±2}:** equatorial ellipticity (width vs depth)
- **ε_{4,m}:** finer structure (neck, bill tip, tail curl)

### 3.7 GW Emission from Duck QNMs

Each QNM mode (n,ℓ,m) radiates with power:

$$\dot{E}_{\ell,m} = \frac{G}{c^{2\ell+1}}\,\omega^{2\ell+2}\,|\delta Q_\ell^m|^2 \times \frac{4\pi(\ell+1)(\ell+2)}{\ell(\ell-1)[(2\ell+1)!!]^2}$$

For the duck: different m-substates have different ω, |δQ_ℓ^m|, and damping rates. The ringdown waveform:

$$h_+(t) + ih_\times(t) = \frac{1}{D_L}\sum_{n,\ell,m} A_{n\ell m}\,e^{-t/\tau_{n\ell m}}\,e^{i\omega_{n\ell m}t}\,{}_{-2}Y_{\ell m}(\iota,\varphi)$$

**Figures/Tables:**
- Fig 3: Perturbation matrix eigenvalues for ℓ=2,3,4 f-modes
- Fig 4: Grotrian diagram: sphere (degenerate) → duck (split)
- Tab 3: ℓ=2 f-mode splitting: perturbative prediction

---

## Section 4: Tidal Love Numbers of the Duck

### 4.1 Effective (Spherical) Love Number — Approach A

Replace the duck with its equivalent sphere, solve the h_2 master ODE:

$$h_2'' + \left\{\frac{2}{R} + \left[\frac{2M}{R^2} + 4\pi R(p-\rho)\right]e^\lambda\right\}h_2' - \left\{\frac{6e^\lambda}{R^2} - 4\pi\left[5\rho + 9p + (\rho+p)\frac{d\rho}{dp}\right]e^\lambda + \left(\frac{d\nu}{dR}\right)^2\right\}h_2 = 0$$

At surface, compute y = R_* h_2'(R_*)/h_2(R_*), then the Hinderer formula:

$$k_2 = \frac{8}{5}C^5(1-2C)^2[2+2C(y-1)-y]\left\{2C[6-3y+3C(5y-8)] + 4C^3[13-11y+C(3y-2)+2C^2(1+y)] + 3(1-2C)^2[2-y+2C(y-1)]\ln(1-2C)\right\}^{-1}$$

$$\Lambda = \frac{2}{3}k_2\,C^{-5}$$

**Newtonian k_2 for polytropes** (Brooker & Olle 1955):

| n | Γ | k_2^N |
|---|---|-------|
| 0 | ∞ | 0.750 |
| 1.0 | 2 | 0.260 |
| 1.5 | 5/3 | 0.143 |
| 3.0 | 4/3 | 0.0116 |

**Compact nuclear duck** (Γ=2, M ~ 1.4 M_☉, C ~ 0.15): k_2 ~ 0.08, Λ ~ 290–880. Comparable to GW170817: Λ = 190^{+390}_{-120}.

**Rubber duck** (elastic): k_2^elastic = (3/2)/[1 + 19μ_shear/(2ρgR)]. At terrestrial scale: k_2 ~ 2.4×10⁻⁸. Punchline: *a bathtub duck is 10⁸ times stiffer than a neutron star duck.*

### 4.2 The Love Tensor — Approach B

For non-spherical bodies, the tidal response is tensorial:

$$Q_{ij}^{\rm induced} = -\lambda_{ijkl}\,\mathcal{E}_{kl}$$

λ_ijkl is a symmetric 5×5 matrix (mapping traceless symmetric tensors to themselves) with up to **15 independent components**. Perturbative expansion:

$$\lambda_{ijkl} = \lambda_0\,\mathcal{P}_{ijkl}^{(0)} + \sum_{\ell',m'}\varepsilon_{\ell' m'}\,\delta\lambda_{ijkl}^{(\ell' m')} + \mathcal{O}(\varepsilon^2)$$

### 4.3 Direction-Dependent k_2(θ,φ) — Approach C

$$k_2(\theta,\phi) = \frac{3}{2R_{\rm eq}^5}\,\lambda_{ijkl}\,n_i n_j n_k n_l = k_2^{(0)} + \sum_m k_2^{(2,m)}Y_2^m + \sum_m k_2^{(4,m)}Y_4^m$$

**Visualization:** k_2 map on the duck surface. Bill (high local C) → suppressed k_2. Body (low C) → enhanced k_2.

### 4.4 Perturbative Correction to Λ — Approach E

The scalar Λ receives its leading shape correction at **second order** in ε:

$$\Lambda_{\rm duck} = \Lambda_{\rm sphere}\left[1 + \sum_{\ell\geq 1}\alpha_\ell\sum_m|\varepsilon_{\ell m}|^2 + \mathcal{O}(\varepsilon^3)\right]$$

Two sources: geometric (+15 Σ|ε_ℓm|² from R^{−5}) and physical (dk_2/dy × δy from Hadamard sensitivity). For Σ|ε_ℓm|² ~ 0.1: geometric correction alone is ~150%.

Mode-by-mode: ℓ=2 (body shape) dominant, ℓ=4 (neck/bill) secondary, ℓ≥6 negligible.

### 4.5 I-Love-Q Breakdown — Approach D

$$\ln\bar{I} = 1.47 + 0.0817\ln\bar{\Lambda} + 0.0149(\ln\bar{\Lambda})^2 + \cdots$$

Holds to O(1%) for spherical NSs. For the duck (ε ~ 0.3): expected O(100%) breakdown. The discrepancy is itself a "sphericity test" — the duck maximally violates I-Love-Q universality.

**Figures/Tables:**
- Fig 5: k_2(θ,φ) on duck surface
- Fig 6: I-Love-Q diagram: duck vs universal curve
- Tab 4: k_2 and Λ comparison: nuclear duck, rubber duck, NS (SLy/APR), BH
- Tab 5: Mode-by-mode ΔΛ/Λ from each ℓ

---

## Section 5: Binary Duck on Post-Newtonian Orbits

### 5.1 Quadrupole-Monopole Interaction

Two ducks (masses M_A, M_B) on quasi-circular orbits. The interaction (Poisson & Will, Eq. 9.41):

$$U_{QM} = -\frac{G M_B}{2 r^3} Q^{C,A}_{ij} n_i n_j$$

PN order: U_QM/E_N ~ (R/r)². For compact objects: 2PN [O(v⁴/c⁴)].

**Orientation cases:** tidally locked (Q_nn = const), freely spinning (Q_nn oscillates), precessing (general Euler angles).

### 5.2 Modified Orbital Dynamics

$$\omega^2 = \frac{G M_{\rm total}}{r^3}\left[1 + \frac{15\,Q_{nn}^{\rm eff}}{2 M_{\rm total}\, r^2}\right]$$

$$E(\omega) = E_{\rm Keplerian}(\omega)\left[1 + \delta_Q(\omega)\right]$$

### 5.3 System Quadrupole and GW Power

$$\mathcal{I}^{\rm sys}_{ij}(t) = \underbrace{\mu\left(x_i x_j - \tfrac{1}{3}r^2\delta_{ij}\right)}_{\rm orbital} + \underbrace{\sum_A \frac{1}{3}R_A(t)\,Q^{C,A}_{\rm body}\,R_A^T(t)}_{\rm body}$$

$$P_{\rm GW} = \frac{G}{45 c^5}\langle\dddot{\mathcal{I}}^{\rm sys}_{ij}\dddot{\mathcal{I}}^{{\rm sys},ij}\rangle = P_{\rm orb} + P_{\rm body} + P_{\rm cross}$$

### 5.4 Phase Evolution

$$\Psi(f) = \Psi_{\rm pp}(f) + \delta\Psi_Q(f), \quad \delta\Psi_Q \sim -\frac{75}{64\eta}\frac{\hat{Q}}{M r^2}\,v^{-1} \quad\text{[2PN]}$$

### 5.5 Waveform

$$h_+(t) = -\frac{2G\mu}{c^4 D_L}(\omega r)^2(1+\cos^2\iota)\cos\!\left(2\Phi(t) + \delta\Phi_Q\right)$$

**Figures/Tables:**
- Fig 7: GW power decomposition P_orb, P_body, P_cross vs f
- Fig 8: Phase dephasing ΔN(f) across LIGO and LISA bands
- Tab 6: PN coefficient comparison: standard vs duck quadrupole correction

---

# Part II: Numerical Simulations

---

## Section 6: QNM Numerical Results — FEM Eigenvalues

### 6.1 Helmholtz Eigenvalue Problem

Solve −∇²ψ = ω²ψ on the duck interior D with ψ|_{∂D} = 0.

**Method:**
1. Volumetric tetrahedralization via tetgen (`meshpy.tet`)
2. P1 FEM: stiffness K, mass M matrices
3. Generalized eigenvalue: K ψ = ω² M ψ
4. Solver: `scipy.sparse.linalg.eigsh` (shift-invert, k=50)
5. Mode classification: project onto Y_ℓ^m via `build_sh_matrix()`

### 6.2 Cross-Validation with Perturbation Theory

| Check | Perturbative (Sec 3) | FEM (Sec 6) |
|-------|---------------------|-------------|
| Sphere limit | Bessel zeros (exact) | Numerical eigenvalues |
| ℓ=2 splitting | 5×5 matrix | 5 split eigenvalues |
| ℓ=3 splitting | 7×7 matrix | 7 eigenvalues |
| Agreement | Linear in ε_ℓm | Exact (numerical) |

### 6.3 Convergence

Refine tet mesh at 2–3 resolutions, verify eigenvalue convergence.

**Figures/Tables:**
- Fig 9: FEM eigenvalue spectrum: duck vs sphere (horizontal bars showing splitting)
- Fig 10: 3D mode shapes: ℓ=2, m=−2..+2 on duck mesh
- Tab 7: First 20 eigenfrequencies: duck vs sphere, with (ℓ,m) labels
- Tab 8: Cross-validation: perturbative vs FEM splitting for ℓ=2,3,4

---

## Section 7: AthenaK Numerical Relativity Simulations

### 7.1 Duck-Shaped TOV Star Initial Data

**Responsibility: Hengrui Zhu**

**What we provide:**
1. Duck deformation coefficients ε_ℓm up to ℓ_max = 10 as text file (`duck_epsilon.dat`)
2. Surface function R(θ,φ)/R_0 = 1 + Σ ε_ℓm Y_ℓ^m(θ,φ)
3. EoS parameters: Γ = 2, K = 100, ρ_c = 1.28×10⁻³

**What Hengrui implements:**
1. New `pgen/duck_star.cpp` based on existing `pgen/dyngr_tov.cpp`
2. Key modification — map TOV profile onto duck shape:
   $$\rho(r, \theta, \phi) = \rho_{\rm TOV}\!\left(\frac{r \cdot R_0}{R(\theta,\phi)}\right), \quad P(r, \theta, \phi) = P_{\rm TOV}\!\left(\frac{r \cdot R_0}{R(\theta,\phi)}\right)$$
3. Constraint handling: accept O(ε²) Hamiltonian constraint violation, let Z4c damp it. Monitor constraint norms.
4. Perturbation: v_pert ~ −0.01 (radial cubic profile)

**AthenaK configuration (GR):** Based on `whisky_tov.athinput`:
```
<coord>
general_rel = true, m = 0.0, a = 0.0
<mhd>
gamma = 2.0, dyn_eos = ideal
<problem>
rhoc = 1.28e-3, kappa = 100.0, npoints = 10000.0, dr = 1e-3
v_pert = -0.01, b_norm = 0.0, user_hist = true
duck_lmax = 10, duck_epsilon_file = duck_epsilon.dat
<adm>
```

**AthenaK configuration (Newtonian):** Based on `polytropic_star.athinput`:
```
<problem>
rho_c = 1.0, R_star = 3.654, v_pert = -0.01
rho_floor = 1.0e-8, njeans = 16
duck_lmax = 10, duck_epsilon_file = duck_epsilon.dat
```

### 7.2 Evolution and Mode Extraction

**GR evolution:** Z4c gauge (1+log lapse, Gamma-driver shift), GR hydro (HLLC, PPM), Cartesian AMR 128³ base + 3–4 levels, ~100–200 T_dyn.

**Post-processing (`code/duck_qnm_extract.py`):**
1. Load .hst → ρ_c(t) time series
2. Hanning window, FFT → P(f)
3. Peak identification → ω_R; Lorentzian fit → ω_I
4. GR: Ψ₄ extraction at multiple radii, _{−2}Y_{ℓm} decomposition, damped sinusoid fit

### 7.3 The Duck Melts to a Sphere (CENTERPIECE)

With isotropic pressure only, the duck is NOT in equilibrium. It oscillates and relaxes:

1. **t = 0:** Duck-shaped density contours
2. **t ~ T_dyn:** Large-scale breathing oscillations
3. **t ~ few T_dyn:** High-ℓ features (bill tip, tail) smooth first
4. **t ~ 10–50 T_dyn:** Low-ℓ features (elongation) persist
5. **t → ∞:** Spherical equilibrium

**Timescales:**
- Newtonian: oscillation at QNM frequencies, no physical damping (numerical viscosity only)
- GR: GW damping τ_GW ~ c⁵/(Gω^{2ℓ+2}R^{2ℓ}) ~ 0.1–1 s for ℓ=2 f-mode. Higher ℓ damps faster → bill melts first.

**Visualization (Fig 11, CENTERPIECE):** Time-lapse panels of 3D density isosurface at t = 0, T/4, T/2, T, 2T, 5T, 10T. Color by δρ/ρ₀.

### 7.4 Waveform Generation

Extract h_+(t), h_×(t) from Ψ₄ integration. The duck ringdown waveform has (2ℓ+1)× more distinct (ω, τ) pairs than a sphere. Generate template waveforms for injection into LIGO noise curves (Section 8.2).

**Figures/Tables:**
- **Fig 11: Duck melts to sphere (CENTERPIECE)** — time-lapse density isosurface
- Fig 12: FFT power spectrum ρ_c(t): duck vs sphere, showing mode splitting
- Fig 13: Ψ₄ ringdown: (2,0), (2,±1), (2,±2) channels with distinct frequencies
- Fig 14: Complex frequency plane: ω_R vs ω_I for first 15 modes
- Tab 9: QNM frequencies: AthenaK vs FEM vs perturbative (cross-tier comparison)
- Tab 10: GW damping times per mode: duck vs sphere

---

# Part III: Experimental Verifications

---

## Section 8: Confrontation with Gravitational Wave Data

### 8.1 PTA Band — Supermassive Binary Ducks (nHz)

A population of supermassive binary ducks (M ~ 10⁶–10⁹ M_☉) produces a stochastic GWB in the nanohertz band.

**Population-level degeneracy theorem:** For any SMBHB model fitting NANOGrav 15yr data, a binary duck population with identical masses and merger rate predicts h_c(f) agreeing to better than 10⁻⁵.

**Proof:** The duck quadrupole correction is:

$$\frac{\delta(\dot{E})}{\dot{E}} \sim \kappa_{\rm duck}\left(\frac{\pi M_c f}{c^3}\right)^{4/3} \sim 10^{-5}\text{--}10^{-9}$$

at PTA frequencies (f ~ 1–100 nHz), which is 4–9 orders of magnitude below the ~30% measurement uncertainty on A_GWB = (2.4 ± 0.7) × 10⁻¹⁵.

**NANOGrav comparison:**
- Duck prediction: h_c = A_duck × (f/f_yr)^{−2/3}, A_duck = 2.4 × 10⁻¹⁵ (identical to SMBHB)
- Spectral shape: γ = 13/3 (identical — depends on inspiral dynamics, not body shape)
- Hellings-Downs: identical (depends on quadrupolar GW radiation, same for all sources)

**Electromagnetic indistinguishability:** At M ~ 10⁹ M_☉, duck angular size at z = 0.1 is θ ~ 6×10⁻⁸ arcsec — 5 orders of magnitude below EHT resolution.

**Mock PTA skymap:** Generate via `hasasia` + `healpy`. Mollweide projection with stochastic background + loudest 5 duck binary sources (duck icons) + PTA sensitivity contours.

**Punchline:** *The nanohertz GW background observed by NANOGrav, EPTA, PPTA, and CPTA is equally consistent with a cosmological population of inspiraling ducks.*

### 8.2 LIGO Band — Stellar-Mass Ducks (10–1000 Hz)

Two observational channels for stellar-mass ducks:

**(a) Inspiral — tidal deformability:**
- The duck's Λ enters the GW phase at 5PN (tidal effects): δΨ_tidal ∝ Λ̃ v¹⁰
- For a 1.4+1.4 M_☉ duck binary, Λ̃_duck differs from Λ̃_sphere by O(10–100%) due to shape (Section 4.4)
- Compare with GW170817 constraint: Λ̃ = 300^{+420}_{-230}
- The duck's anisotropic Love tensor (Section 4.2) means Λ̃ depends on the duck's orientation relative to the orbital plane — an effect absent for spherical NS

**(b) Ringdown — QNM spectroscopy:**
- Post-merger ringdown encodes the duck's QNM spectrum (Section 3)
- The split frequencies from gravitational Zeeman splitting (Section 3.5) produce a richer ringdown than a spherical remnant
- At LIGO sensitivity, the ℓ=2 f-mode is detectable for events at d < 100 Mpc with SNR > 8
- The m-dependent splitting Δω/ω ~ O(ε) ~ 0.3 is large enough to resolve with next-generation detectors (Cosmic Explorer, Einstein Telescope)
- **Key question:** Can QNM spectroscopy distinguish a duck from a sphere? Yes, if Δω/ω > 1/Q_mode (where Q = ω_R/(2ω_I) is the quality factor)

**Detectability analysis:**
- Inject duck ringdown waveforms (from AthenaK, Section 7.4) into simulated LIGO noise
- Compute SNR for various distances and orientations
- Determine the minimum SNR to resolve the m-splitting
- Compare with current (O4) and future (CE, ET) detector sensitivities

**Speculative:** *If a GW event is observed with anomalous QNM splitting inconsistent with Kerr or standard NS models, the duck hypothesis cannot be excluded.*

### 8.3 Breaking the Degeneracy (Speculative)

Methods that could in principle distinguish ducks from spherical objects:
1. **GW memory:** nonlinear memory from duck mergers encodes full multipole structure (requires SNR ~ 100+)
2. **Multi-band:** combining PTA + LISA + LIGO constrains the inspiral-to-merger transition, where duck effects are largest
3. **Tidal disruption events:** duck-shaped tidal field differs, but monopole dominates at disruption radius
4. **Direct imaging:** EHT angular resolution insufficient by 5 orders of magnitude

**Figures/Tables:**
- Fig 15: Duck GWB overlaid on NANOGrav 15yr violin plot
- Fig 16: Mock PTA skymap (Mollweide projection)
- Fig 17: Injected duck ringdown waveform in LIGO noise, showing m-splitting
- Tab 11: Binary duck population parameters reproducing NANOGrav signal
- Tab 12: Detectability of duck QNM splitting: SNR vs distance for O4, CE, ET

---

## Specification for Hengrui Zhu (AthenaK Initial Data)

### We provide:
1. **ε_ℓm file** (`duck_epsilon.dat`): `# ell  m  epsilon_real  epsilon_imag` format
2. **Surface function** R(θ,φ)/R_0 on a fine grid
3. **EoS:** Γ = 2, K = 100, ρ_c = 1.28×10⁻³

### Hengrui implements:
1. `pgen/duck_star.cpp` based on `dyngr_tov.cpp`: read ε_ℓm, compute R(θ,φ), map TOV via `GetPrimitivesAtPoint(eos, r·R_0/R(θ,φ), ...)`
2. Constraint: Option A (accept O(ε²) violation, Z4c damps it). Monitor norms.
3. Output: standard binary + history (.hst with ρ_max, α_min)

### Hengrui does NOT handle:
- Post-processing, FFT, Ψ₄ decomposition, figures (our Python scripts)
- Newtonian runs (we modify polytropic_star pgen)

---

## Complete Figures and Tables

| # | Content | Section |
|---|---------|---------|
| **Part I: Formalism** | | |
| Fig 1 | Duck mesh + SH reconstruction at ℓ_max = 0, 2, 4, 8, full | 2 |
| Fig 2 | Deformation spectrum |ε_ℓm| vs ℓ | 2 |
| Fig 3 | Perturbation matrix eigenvalues for ℓ=2,3,4 f-modes | 3 |
| Fig 4 | Grotrian diagram: sphere → duck level splitting | 3 |
| Fig 5 | Direction-dependent k_2(θ,φ) on duck surface | 4 |
| Fig 6 | I-Love-Q diagram: duck vs universal curve | 4 |
| Fig 7 | GW power decomposition P_orb, P_body, P_cross vs f | 5 |
| Fig 8 | Phase dephasing ΔN(f): LIGO + LISA bands | 5 |
| **Part II: Numerics** | | |
| Fig 9 | FEM eigenvalue spectrum: duck vs sphere | 6 |
| Fig 10 | 3D mode shapes on duck mesh | 6 |
| **Fig 11** | **Duck melts to sphere: time-lapse (CENTERPIECE)** | **7** |
| Fig 12 | FFT of ρ_c(t): duck vs sphere, mode splitting | 7 |
| Fig 13 | Ψ₄ ringdown channels | 7 |
| Fig 14 | Complex frequency plane ω_R vs ω_I | 7 |
| **Part III: Observations** | | |
| Fig 15 | Duck GWB on NANOGrav 15yr violin plot | 8.1 |
| Fig 16 | Mock PTA skymap (Mollweide) | 8.1 |
| Fig 17 | Duck ringdown in LIGO noise | 8.2 |
| | | |
| Tab 1 | Duck geometric/multipole parameters | 2 |
| Tab 2 | Deformation coefficients ε_ℓm | 2 |
| Tab 3 | Perturbative QNM splitting predictions | 3 |
| Tab 4 | Love numbers: nuclear duck, rubber duck, NS, BH | 4 |
| Tab 5 | Mode-by-mode ΔΛ/Λ from each ℓ | 4 |
| Tab 6 | PN coefficients: standard vs duck correction | 5 |
| Tab 7 | FEM eigenfrequencies: duck vs sphere (first 20) | 6 |
| Tab 8 | Cross-tier QNM comparison: perturbative vs FEM vs AthenaK | 6–7 |
| Tab 9 | QNM frequencies from AthenaK | 7 |
| Tab 10 | GW damping times per mode | 7 |
| Tab 11 | Binary duck population for NANOGrav | 8.1 |
| Tab 12 | Duck QNM detectability: SNR vs distance | 8.2 |

---

## Key References

### Stellar Oscillations and QNMs
1. Kokkotas & Schmidt, Living Rev. Rel. 2, 2 (1999)
2. Andersson & Kokkotas, MNRAS 299, 1059 (1998)
3. Cowling, MNRAS 101, 367 (1941)
4. Chandrasekhar, ApJ 139, 664 (1964)
5. Chandrasekhar, "Ellipsoidal Figures of Equilibrium" (Yale, 1969)

### Mode Splitting
6. Yoshida & Lee, MNRAS 317, L1 (2000)
7. Gaertig & Kokkotas, PRD 78, 064063 (2008)
8. Passamonti et al., MNRAS 394, 730 (2009)
9. Doneva et al., PRD 88, 044052 (2013)
10. Lamb, Proc. London Math. Soc. 13, 189 (1881)

### Domain Perturbation Theory
11. Hadamard, "Leçons sur le Calcul des Variations" (1910)
12. Henry & Douanla, J. Math. Anal. Appl. 430, 995 (2015)

### AthenaK
13. Stone et al. (2024)
14. Zhu et al. (2024)
15. Fields et al. (2024)

### Tidal Love Numbers and I-Love-Q
16. Love, Proc. R. Soc. A 82, 73 (1909)
17. Hinderer, ApJ 677, 1216 (2008)
18. Flanagan & Hinderer, PRD 77, 021502 (2008)
19. Damour & Nagar, PRD 80, 084035 (2009)
20. Binnington & Poisson, PRD 80, 084018 (2009)
21. Yagi & Yunes, Science 341, 365 (2013); PRD 88, 023009 (2013)
22. Landry & Poisson, PRD 91, 104018 (2015)
23. Pani et al., PRD 92, 024010 (2015)
24. Brooker & Olle, MNRAS 115, 101 (1955)

### Binary PN Dynamics
25. Poisson & Will, "Gravity" (Cambridge, 2014) — Ch. 9–11
26. Blanchet, Living Rev. Rel. 17, 2 (2014)
27. Poisson, PRD 57, 5287 (1998)
28. Barker & O'Connell, PRD 12, 329 (1975)

### PTA and GWB
29. Agazie et al. (NANOGrav), ApJL 951, L8 (2023)
30. Antoniadis et al. (EPTA), A&A 678, A50 (2023)
31. Reardon et al. (PPTA), ApJL 951, L6 (2023)
32. Phinney, arXiv:astro-ph/0108028 (2001)
33. Sesana et al., MNRAS 390, 192 (2008)
34. Mingarelli et al., Nature Astronomy 1, 886 (2017)
35. Hellings & Downs, ApJ 265, L39 (1983)
36. Taylor et al., ApJ 819, L6 (2016)
37. Hazboun et al., PRD 100, 104028 (2019)

### LIGO and Tidal Effects
38. Abbott et al. (LIGO/Virgo), PRL 119, 161101 (2017) — GW170817
39. Abbott et al. (LIGO/Virgo), PRL 121, 161101 (2018) — GW170817 tidal constraints
40. Hinderer et al., PRD 81, 123016 (2010) — tidal effects in GW signals

### Original Works
41. Lehmann, arXiv:2504.00506 (2025) — "Higher multipoles of the cow"
42. Tolman, Phys. Rev. 55, 364 (1939); Oppenheimer & Volkoff, Phys. Rev. 55, 374 (1939)
