# Quasinormal Modes of the Duck: Gravitational Zeeman Splitting and Shape Relaxation in Numerical Relativity

## Paper Overview

**Premise (humorous):** What happens when a duck-shaped self-gravitating body oscillates and relaxes under its own gravity?

**Physics (rigorous):** The non-spherical geometry of the duck breaks the rotational symmetry of the eigenvalue problem for stellar oscillations, lifting the (2ℓ+1)-fold degeneracy of each QNM multiplet — the gravitational analog of Zeeman splitting. We compute this splitting analytically (perturbation theory), numerically (FEM eigenvalues), and via full numerical relativity (AthenaK time-domain evolution). The NR simulation captures the "duck melting to a sphere" — the visual centerpiece of the paper.

**Additional section:** A population of supermassive binary ducks produces a nanohertz GW signal indistinguishable from the PTA-observed stochastic background, demonstrating that current observations cannot rule out the duck hypothesis.

---

## Paper Structure

### Section 1: Introduction
- The spherical cow approximation and its discontents (Lehmann 2025)
- Extend to the duck: non-spherical geometry as a symmetry-breaking perturbation
- Preview: Zeeman splitting of QNMs, duck melting, PTA indistinguishability
- ~2 pages

### Section 2: Duck Geometry and Multipole Characterization
- Duck mesh (DASSL rubber duck, 4732 vertices, 9460 faces)
- Multipole moments: Q^C, I, Q_ℓ^m up to ℓ=10 (existing validated pipeline)
- Surface deformation coefficients ε_ℓm via SH fitting
- Equivalent sphere radius R_eq, principal axes, ellipticity
- Shape-induced quadrupole parameter κ_duck = Q^C/(M R_eq²)
- ~3 pages

### Section 3: Perturbative QNM Splitting (Analytical — Tier C)
- Degenerate perturbation theory for eigenvalue shifts
- Hadamard-type formula for domain perturbation
- Selection rules from angular momentum coupling
- First-order frequency shifts δω_{n,ℓ,m} as function of ε_ℓm
- Which duck features (bill, tail, body) dominate the splitting
- ~4 pages

### Section 4: FEM Eigenvalue Problem (Numerical — Tier B)
- Helmholtz equation on duck interior domain
- Volumetric tetrahedralization, P1 FEM
- Eigenvalue spectrum: duck vs sphere
- Mode classification via SH projection
- ~3 pages

### Section 5: AthenaK Numerical Relativity Simulation (Tier A)
- **5.1** Duck-shaped TOV star initial data (Hengrui Zhu)
- **5.2** Evolution and mode extraction
- **5.3** The duck melts to a sphere — shape relaxation dynamics
- **5.4** GW emission from duck QNMs
- ~5 pages

### Section 6: The Duck Hypothesis and Pulsar Timing Arrays
- Binary duck population model in the nanohertz band
- Characteristic strain spectrum h_c(f) matching NANOGrav 15yr
- Mock PTA skymap
- "Current observations cannot distinguish binary ducks from SMBHBs"
- ~3 pages

### Section 7: Discussion and Conclusions
- ~1 page

**Total: ~21 pages + figures/tables**

---

## Section 3: Perturbative QNM Splitting — Detailed Methodology

### 3.1 Problem Setup

The duck boundary is parametrized as a deformed sphere:

$$R_{\rm duck}(\theta,\phi) = R_0\left[1 + \sum_{\ell \geq 1, m} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right]$$

where R_0 = R_eq is the equivalent sphere radius and ε_ℓm are the deformation coefficients extracted from the SH fitting pipeline (existing code: `render_sphere_to_cow.py` → `build_sh_matrix()`, `fit_sh_coefficients()`).

The unperturbed problem is the eigenvalue equation on a sphere of radius R_0. For the scalar Helmholtz problem (simplest model):

$$-\nabla^2 \psi_{n\ell m} = \omega_{n\ell}^2\,\psi_{n\ell m}, \quad \psi\big|_{r=R_0} = 0$$

The eigenvalues are ω_{n,ℓ} = j_{n,ℓ}/R_0 where j_{n,ℓ} is the n-th zero of the spherical Bessel function j_ℓ. Each (n,ℓ) level has (2ℓ+1)-fold degeneracy.

For fluid modes (Cowling approximation), the unperturbed eigenvalues are the f-mode and p-mode frequencies of a uniform-density or polytropic sphere, which are tabulated (Kokkotas & Schmidt, Living Rev. Rel. 2, 1999).

### 3.2 Hadamard Formula for Domain Perturbation

The Hadamard formula gives the sensitivity of Dirichlet eigenvalues to boundary deformation. For a domain D with boundary ∂D deformed by δn(x) in the outward normal direction:

$$\delta\omega_{n\ell m}^2 = -\int_{\partial D} \left|\frac{\partial \psi_{n\ell m}}{\partial n}\right|^2 \delta n\,dS$$

For our spherical-to-duck deformation with δn = R_0 Σ ε_ℓ'm' Y_ℓ'^{m'}:

$$\delta\omega_{n\ell m}^2 = -R_0 \sum_{\ell',m'} \varepsilon_{\ell' m'} \int_{S^2} \left|\frac{\partial \psi_{n\ell m}}{\partial r}\bigg|_{R_0}\right|^2 Y_{\ell'}^{m'}(\theta,\phi)\,R_0^2\,d\Omega$$

### 3.3 Degenerate Perturbation Theory

Within each degenerate (n,ℓ) multiplet, the first-order correction requires diagonalizing the (2ℓ+1) × (2ℓ+1) perturbation matrix:

$$V_{m_1, m_2}^{(n,\ell)} = -R_0 \sum_{\ell',m'} \varepsilon_{\ell' m'} \int_{S^2} \frac{\partial \psi_{n\ell m_1}^*}{\partial r}\bigg|_{R_0} \frac{\partial \psi_{n\ell m_2}}{\partial r}\bigg|_{R_0} Y_{\ell'}^{m'}(\theta,\phi)\,R_0^2\,d\Omega$$

**Key simplification:** On the sphere, ψ_{n,ℓ,m}(r,θ,φ) = R_{n,ℓ}(r) Y_ℓ^m(θ,φ), so:

$$\frac{\partial \psi_{n\ell m}}{\partial r}\bigg|_{R_0} = R_{n\ell}'(R_0)\,Y_\ell^m(\theta,\phi)$$

The angular integral becomes:

$$\int_{S^2} Y_\ell^{m_1*}(\Omega)\,Y_\ell^{m_2}(\Omega)\,Y_{\ell'}^{m'}(\Omega)\,d\Omega = (-1)^{m_1}\sqrt{\frac{(2\ell+1)^2(2\ell'+1)}{4\pi}} \begin{pmatrix} \ell & \ell & \ell' \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} \ell & \ell & \ell' \\ -m_1 & m_2 & m' \end{pmatrix}$$

This uses the Gaunt integral expressed via Wigner 3j symbols.

### 3.4 Selection Rules

The 3j symbols enforce:
1. **Triangle inequality:** |ℓ − ℓ| ≤ ℓ' ≤ ℓ + ℓ, i.e., 0 ≤ ℓ' ≤ 2ℓ
2. **Parity:** ℓ + ℓ + ℓ' must be even → ℓ' must be even
3. **m-selection:** m' = m_1 − m_2

**Consequence:** Only the **even-ℓ' deformation harmonics** (ℓ' = 0, 2, 4, ...) contribute to first-order splitting within an (n,ℓ) multiplet:
- ε_{0,0} shifts all levels equally (no splitting, just overall shift)
- ε_{2,m} (quadrupole deformation) produces the primary splitting
- ε_{4,m} (hexadecapole) contributes secondary splitting
- Higher even-ℓ' contribute progressively smaller corrections

### 3.5 Explicit Formula for ℓ=2 f-mode Splitting

The ℓ=2 f-mode (5-fold degenerate on sphere) splits under the duck's quadrupole deformation (ε_{2,m}). The perturbation matrix is 5×5:

$$V_{m_1,m_2} = -|R_{n,2}'(R_0)|^2 R_0^3 \sum_{m'} \varepsilon_{2,m'}\,(-1)^{m_1}\sqrt{\frac{25 \cdot 5}{4\pi}} \begin{pmatrix} 2 & 2 & 2 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} 2 & 2 & 2 \\ -m_1 & m_2 & m' \end{pmatrix}$$

The 3j symbol (2 2 2 | 0 0 0) = √(2/35).

Diagonalizing this 5×5 matrix gives the 5 split frequencies. The eigenvalues are linear combinations of ε_{2,m'} weighted by the 3j symbols.

### 3.6 Physical Interpretation

The duck's deformation spectrum ε_ℓm encodes which geometric features drive the splitting:
- **ε_{2,0}:** axial elongation (body shape) → splits m=0 from m=±1, ±2
- **ε_{2,±1}:** tilt asymmetry (e.g., bill-to-tail asymmetry in one plane)
- **ε_{2,±2}:** equatorial ellipticity (width vs depth)
- **ε_{4,m}:** finer structure (neck, bill tip, tail curl)

The dominant contribution to splitting comes from the largest ε_ℓm coefficients. For the rubber duck, the dominant deformation is expected to be ε_{2,0} (prolate along the bill-tail axis) and ε_{2,2} (equatorial asymmetry from the round body).

### 3.7 Implementation

**Inputs:** ε_ℓm from Phase 0 (SH fitting of duck mesh), sphere eigenvalues (Bessel zeros for Helmholtz, or tabulated f/p-mode values for fluid).

**Computation:**
1. Build perturbation matrix V for each (n,ℓ) multiplet using Wigner 3j (from `sympy.physics.wigner.wigner_3j` or `py3nj`)
2. Diagonalize V → eigenvalues give frequency shifts δω²
3. The split frequencies are ω_{n,ℓ,m} = √(ω_{n,ℓ}² + δω²_m)

**Output:** Table of split frequencies for ℓ=2, 3, 4 f-modes; Grotrian diagram.

---

## Section 4: FEM Eigenvalue Problem — Detailed Methodology

### 4.1 Problem

Solve -∇²ψ = ω²ψ on the duck interior D with ψ|_{∂D} = 0 (Dirichlet).

### 4.2 Method

1. **Volume mesh:** Tetrahedralize duck interior using tetgen (via `meshpy.tet`) from the surface mesh `duck.off`
2. **P1 FEM:** Linear (hat function) basis on tetrahedra
   - Stiffness matrix: K_ij = ∫_D ∇φ_i · ∇φ_j dV
   - Mass matrix: M_ij = ∫_D φ_i φ_j dV
   - Both assembled element-by-element, standard formulas for linear tets
3. **Generalized eigenvalue problem:** K ψ = ω² M ψ
4. **Solver:** `scipy.sparse.linalg.eigsh(K, k=50, M=M, sigma=0.0)` (shift-invert for smallest eigenvalues)
5. **Mode classification:** Project each eigenvector onto Y_ℓ^m basis at vertex positions (reusing `build_sh_matrix()` from `render_sphere_to_cow.py`):
   $$a_{\ell m}^{(n)} = \sum_i \psi_n(x_i)\,Y_\ell^{m*}(\theta_i,\phi_i)\,\Delta V_i$$
   Angular power spectrum: P_ℓ = Σ_m |a_{ℓm}|²

### 4.3 Comparison with Sphere

Run the same FEM on a sphere of radius R_eq (or use analytical Bessel zeros). The sphere eigenvalues have (2ℓ+1)-fold degeneracy. The duck eigenvalues show the splitting. Differences validate the perturbative prediction from Section 3.

### 4.4 Convergence

Refine the tet mesh (increase vertex count) and verify eigenvalue convergence. Report the first 20 eigenvalues at 2-3 mesh resolutions.

---

## Section 5: AthenaK NR Simulation — Detailed Methodology

### 5.1 Duck-Shaped TOV Star Initial Data

**Responsibility: Hengrui Zhu**

**What Hengrui needs from us:**
1. Duck deformation coefficients ε_ℓm up to ℓ_max = 10 (output of Phase 0)
2. The surface function R(θ,φ) = R_0 [1 + Σ ε_ℓm Y_ℓ^m(θ,φ)]
3. Specification of the EoS and central density

**What Hengrui provides:**
1. Modified `pgen/duck_star.cpp` based on existing `pgen/dyngr_tov.cpp`
2. The key modification: instead of querying the TOV solution at radius r = |x|, query it at the rescaled radius r_eff = r / R(θ,φ) * R_0, where R(θ,φ) is the duck surface function
3. This maps the spherical TOV profile onto the duck shape:
   $$\rho(r, \theta, \phi) = \rho_{\rm TOV}\!\left(\frac{r \cdot R_0}{R(\theta,\phi)}\right)$$
   $$P(r, \theta, \phi) = P_{\rm TOV}\!\left(\frac{r \cdot R_0}{R(\theta,\phi)}\right)$$
4. Constraint solver: after setting the fluid data, run ADMToZ4c + ADMConstraints to initialize the spacetime consistently
5. Perturbation: small velocity perturbation v_pert ~ -0.01 (radial, cubic profile as in existing whisky_tov.athinput) to excite modes

**EoS specification:**
- **Stellar mass duck (full GR):** Polytrope with Γ = 2, K = 100, ρ_c = 1.28×10⁻³ (standard TOV test star, M ≈ 1.4 M_☉, R ≈ 9.6 km in code units, C ≈ 0.15)
- **Supermassive duck (Newtonian):** Lane-Emden n=1.5 (Γ = 5/3), K ≈ 0.4001, ρ_c = 1.0 (existing validated setup)

**AthenaK configuration (GR run):**

Based on existing `whisky_tov.athinput`:
```
<coord>
general_rel = true
m = 0.0
a = 0.0

<mhd>
gamma = 2.0
dyn_eos = ideal

<problem>
rhoc = 1.28e-3
kappa = 100.0
npoints = 10000.0
dr = 1e-3
v_pert = -0.01
b_norm = 0.0
user_hist = true
# New parameters for duck shape:
duck_lmax = 10
duck_epsilon_file = duck_epsilon.dat

<adm>
```

**AthenaK configuration (Newtonian run):**

Based on existing `polytropic_star.athinput`:
```
<problem>
rho_c = 1.0
R_star = 3.654
v_pert = -0.01
rho_floor = 1.0e-8
njeans = 16
duck_lmax = 10
duck_epsilon_file = duck_epsilon.dat
```

### 5.2 Evolution and Mode Extraction

**GR evolution:**
- Z4c gauge: 1+log lapse, Gamma-driver shift
- GR hydro: HLLC Riemann solver, PPM reconstruction
- Grid: Cartesian AMR, typical resolution 128³ base with 3-4 AMR levels
- Boundary: outflow
- Duration: ~100-200 dynamical times (T_dyn = √(R³/GM) ≈ 0.5 ms for NS)

**Diagnostics extracted at each timestep:**
- ρ_max (maximum density — oscillation tracer)
- α_min (minimum lapse — GR diagnostic)
- Central density ρ_c(t)
- Multipole moments of density perturbation at several radii: decompose δρ(r,θ,φ,t) into Y_ℓ^m

**Post-processing (Python, `code/duck_qnm_extract.py`):**
1. Load history file (.hst) → ρ_c(t) time series
2. Remove secular drift, apply Hanning window
3. FFT → power spectrum P(f)
4. Identify peaks → QNM frequencies ω_R
5. Fit each peak to Lorentzian → extract damping rate ω_I (width)
6. For GR runs: also extract Ψ₄ at extraction radii, decompose into (ℓ,m), fit damped sinusoids

**Newtonian evolution:**
- Same procedure but with Newtonian hydro + multigrid Poisson
- No Ψ₄ extraction (no GWs in Newtonian gravity)
- GW damping rates computed analytically from the mode shapes (Section 5.4)

### 5.3 The Duck Melts to a Sphere

**The key visual result.** With only isotropic pressure (no anisotropic stress), the duck-shaped star is NOT in equilibrium — it will oscillate and relax toward spherical symmetry.

**Dynamics:**
1. **t = 0:** Duck-shaped TOV star with density contours following the duck geometry
2. **t ~ T_dyn:** Large-scale oscillations begin; the duck "breathes" and deforms
3. **t ~ few × T_dyn:** High-ℓ features (bill tip, tail) smooth out first (higher ℓ → higher frequency → faster damping)
4. **t ~ 10-50 × T_dyn:** Low-ℓ features (overall elongation) persist longer
5. **t → ∞:** Approaches spherical equilibrium

**What determines the melting timescale?**
- **Newtonian:** No GW damping. Relaxation is via numerical viscosity (resolution-dependent) or explicit viscosity if added. The oscillation frequencies are the QNMs. In practice, the duck oscillates indefinitely in the absence of dissipation — the "melting" is the envelope of the oscillation, not a monotonic relaxation.
- **GR:** GW emission provides physical damping. The damping time τ_GW for the ℓ=2 f-mode is:
  $$\tau_{\rm GW} \sim \frac{c^5}{G\omega^{2\ell+2}R^{2\ell}} \sim 0.1\text{-}1 \text{ s}$$
  (for a 1.4 M_☉ NS f-mode at ~2 kHz)

  Higher ℓ modes damp faster: τ_GW ∝ ω^{-(2ℓ+2)}. So the bill (high-ℓ features) melts first, consistent with the visual.

**Visualization plan:**
- **Fig 7 (centerpiece):** Time-lapse panels of 3D density isosurface at t = 0, T/4, T/2, T, 2T, 5T, 10T showing the duck gradually losing its features
- Render using matplotlib Poly3DCollection (existing `render_mesh()` from `render_sphere_to_cow.py`) or VTK for higher quality
- Color by local density perturbation δρ/ρ₀

### 5.4 GW Emission from Duck QNMs

Each QNM mode (n,ℓ,m), when excited, radiates gravitational waves. The radiated power is:

$$\dot{E}_{\ell,m} = \frac{G}{c^{2\ell+1}}\,\omega^{2\ell+2}\,|\delta Q_\ell^m|^2 \times \frac{4\pi(\ell+1)(\ell+2)}{\ell(\ell-1)[(2\ell+1)!!]^2}$$

where δQ_ℓ^m is the time-varying multipole moment associated with mode (n,ℓ,m).

**For a sphere:** all m-substates have the same |δQ_ℓ^m| → same damping rate.

**For the duck:** the broken symmetry means:
1. Different m-substates have different ω_{n,ℓ,m} → different radiation frequencies
2. Different mode shapes → different |δQ_ℓ^m| → different amplitudes and damping rates
3. Mode mixing → a nominal (ℓ,m) mode radiates in multiple multipole channels

**Ringdown waveform from duck QNM excitation:**

$$h_+(t) + ih_\times(t) = \frac{1}{D_L}\sum_{n,\ell,m} A_{n\ell m}\,e^{-t/\tau_{n\ell m}}\,e^{i\omega_{n\ell m}t}\,{}_{-2}Y_{\ell m}(\iota,\varphi)$$

For the duck, the number of distinct (ω, τ) pairs is (2ℓ+1) times larger than for a sphere at each (n,ℓ). The Grotrian diagram (Fig 10) illustrates this spectral richness.

**From AthenaK GR run:** Extract Ψ₄ at multiple extraction radii → decompose into _{-2}Y_{ℓm} → fit each (ℓ,m) channel to damped exponential → measure (ω_R, ω_I = 1/τ) per mode.

---

## Section 6: PTA Indistinguishability — Methodology Specification

*[Detailed methodology to be provided by the PTA research subagent — see separate report]*

**Core argument:**
1. A population of supermassive binary ducks (M ~ 10⁸⁻⁹ M_☉) inspiraling in the nanohertz band produces a stochastic GW background
2. The characteristic strain spectrum follows h_c(f) = A × (f/f_yr)^{-2/3} (same spectral shape as circular SMBHB inspiral)
3. The duck's permanent quadrupole modifies individual source power by δP/P ~ (v/c)⁴ ~ 10⁻⁴ to 10⁻⁸ in the PTA band — completely undetectable
4. Therefore: the NANOGrav 15yr GWB signal (A_GWB ~ 2.4×10⁻¹⁵, γ ~ 13/3) is equally consistent with a population of binary ducks as with standard SMBHBs
5. **Punchline:** "We cannot distinguish binary ducks from binary black holes with current PTA observations"

**Deliverables:**
- Mock PTA skymap (Mollweide projection) showing strain from individual loud binary duck sources
- Overlay duck GWB spectrum on NANOGrav 15yr violin plot
- Brief parameter estimation: what binary duck population (mass function, merger rate) reproduces the observed A_GWB

---

## Cross-Tier Validation Plan

| Check | Tier C (Perturbative) | Tier B (FEM) | Tier A (AthenaK) |
|-------|----------------------|--------------|-------------------|
| Sphere limit | Analytical eigenvalues (Bessel zeros) | FEM on sphere mesh | Spherical TOV evolution |
| ℓ=2 splitting | 5×5 matrix diagonalization | 5 split eigenvalues | 5 peaks in FFT |
| ℓ=3 splitting | 7×7 matrix | 7 eigenvalues | 7 peaks |
| Frequency values | δω linear in ε_ℓm | Exact (numerical) | Exact (NR) |
| Agreement | Should match Tier B for small ε | Reference standard | Should match Tier B to resolution |
| Damping rates | Analytical formula | Not available | Measured from Ψ₄ (GR) |

---

## Figures and Tables Plan

| # | Content | Section |
|---|---------|---------|
| Fig 1 | Duck mesh and SH reconstruction at ℓ_max = 0, 2, 4, 8, full (existing) | 2 |
| Fig 2 | Deformation spectrum |ε_ℓm| vs ℓ (bar chart) | 2 |
| Fig 3 | Perturbation matrix eigenvalues for ℓ=2,3,4 f-modes | 3 |
| Fig 4 | Grotrian diagram: sphere (degenerate) → duck (split) | 3 |
| Fig 5 | FEM eigenvalue spectrum: duck vs sphere (horizontal bars) | 4 |
| Fig 6 | 3D mode shapes: ℓ=2, m=-2..+2 on duck mesh (colored by eigenvector) | 4 |
| **Fig 7** | **Duck melts to sphere: time-lapse density isosurface (CENTERPIECE)** | **5** |
| Fig 8 | FFT power spectrum of ρ_c(t): duck vs sphere, showing mode splitting | 5 |
| Fig 9 | Ψ₄ ringdown: (2,0), (2,±1), (2,±2) channels (GR run) | 5 |
| Fig 10 | Complex frequency plane: ω_R vs ω_I for first 15 modes | 5 |
| Fig 11 | Duck GWB spectrum overlaid on NANOGrav 15yr violin plot | 6 |
| Fig 12 | Mock PTA skymap: Mollweide projection of strain from binary duck population | 6 |
| Tab 1 | Duck geometric and multipole parameters | 2 |
| Tab 2 | Deformation coefficients ε_ℓm up to ℓ=10 | 2 |
| Tab 3 | ℓ=2 f-mode splitting: perturbative vs FEM vs AthenaK | 3-5 |
| Tab 4 | First 20 eigenfrequencies: duck vs sphere | 4 |
| Tab 5 | GW damping times per mode: duck vs sphere | 5 |
| Tab 6 | Binary duck population parameters reproducing NANOGrav signal | 6 |

---

## Specification for Hengrui Zhu (AthenaK Initial Data)

### What we provide:
1. **Duck deformation coefficients** ε_ℓm as a text file (`duck_epsilon.dat`):
   ```
   # ell  m  epsilon_real  epsilon_imag
   2  0  -0.268  0.000
   2  1   0.001  0.003
   2  2   0.145 -0.012
   ...
   ```
   (To be computed from Phase 0 SH fitting)

2. **The surface function** R(θ,φ)/R_0 = 1 + Σ ε_ℓm Y_ℓ^m(θ,φ), evaluated on a fine (θ,φ) grid if needed

3. **EoS parameters:** Γ = 2, K = 100, ρ_c = 1.28×10⁻³ (standard test star)

### What Hengrui implements:
1. **New problem generator** `pgen/duck_star.cpp` based on `dyngr_tov.cpp`
   - Read ε_ℓm from file
   - Compute R(θ,φ) using Y_ℓ^m (need SH evaluation on the grid — can use a lookup table or inline computation)
   - Map TOV profile: query `GetPrimitivesAtPoint(eos, r * R_0 / R(θ,φ), ...)` instead of `GetPrimitivesAtPoint(eos, r, ...)`
   - This gives ρ(r,θ,φ) = ρ_TOV(r·R_0/R(θ,φ)), P similarly
   - Apply velocity perturbation as in standard setup

2. **Constraint handling:** The deformed initial data violates the Hamiltonian constraint (by O(ε²)). Two options:
   - **Option A (simpler):** Accept O(ε²) constraint violation — it will be damped by Z4c constraint damping during evolution. Suitable if ε is small enough.
   - **Option B (rigorous):** Solve the constraint equations on the deformed data via the elliptic solver (conformal thin-sandwich or XCTS). This requires more infrastructure.
   - **Recommendation:** Start with Option A. For ε ~ 0.1-0.3 (duck deformation), the constraint violation is O(0.01-0.1). Monitor constraint norms during evolution.

3. **Output:** Standard AthenaK binary output + history file with ρ_max, α_min at each timestep

### What Hengrui does NOT need to worry about:
- Post-processing (FFT, mode extraction, Ψ₄ decomposition) — handled by our Python scripts
- Newtonian runs — we can modify the existing polytropic star pgen ourselves
- Analysis and figures — all done post-simulation

---

## Key References

### Stellar Oscillations and QNMs
1. Kokkotas & Schmidt, Living Rev. Rel. 2, 2 (1999) — "Quasi-Normal Modes of Stars and Black Holes" [comprehensive review]
2. Andersson & Kokkotas, MNRAS 299, 1059 (1998) — "Towards gravitational wave asteroseismology"
3. Cowling, MNRAS 101, 367 (1941) — the Cowling approximation
4. Chandrasekhar, ApJ 139, 664 (1964) — "General variational principle governing oscillations"
5. Chandrasekhar, "Ellipsoidal Figures of Equilibrium" (Yale, 1969) — oscillations of non-spherical bodies

### Mode Splitting from Deformation/Rotation
6. Yoshida & Lee, MNRAS 317, L1 (2000) — inertial modes of slowly rotating relativistic stars
7. Gaertig & Kokkotas, PRD 78, 064063 (2008) — oscillations of rapidly rotating NS
8. Passamonti et al., MNRAS 394, 730 (2009) — coupling of polar and axial modes in rotating NS
9. Doneva et al., PRD 88, 044052 (2013) — universal relations for NS oscillations
10. Lamb, Proc. London Math. Soc. 13, 189 (1881) — oscillations of ellipsoidal figures [analytical baseline]

### Domain Perturbation Theory
11. Hadamard, "Leçons sur le Calcul des Variations" (1910) — eigenvalue sensitivity to boundary shape
12. Henry & Douanla, J. Math. Anal. Appl. 430, 995 (2015) — domain perturbation for eigenvalues [modern treatment]
13. Kozlov, Maz'ya & Movchan, "Asymptotic Analysis of Fields in Multi-Structures" (2000) — spectral theory on perturbed domains

### AthenaK and Numerical Methods
14. Stone et al. (2024) — AthenaK framework
15. Zhu et al. (2024) — AthenaK NR solver (Z4c formalism)
16. Fields et al. (2024) — GR hydro/MHD in dynamical spacetimes

### I-Love-Q and Tidal Deformability
17. Yagi & Yunes, Science 341, 365 (2013) — I-Love-Q relations
18. Flanagan & Hinderer, PRD 77, 021502 (2008) — tidal Love numbers in GW
19. Hinderer, ApJ 677, 1216 (2008) — Love number computation

### PTA and GWB
20. Agazie et al. (NANOGrav), ApJL 951, L8 (2023) — NANOGrav 15yr GWB evidence
21. Antoniadis et al. (EPTA), A&A 678, A50 (2023) — EPTA DR2
22. Phinney, arXiv:astro-ph/0108028 (2001) — GWB from compact binary population
23. Sesana et al., MNRAS 390, 192 (2008) — SMBHB population models
24. Burke-Spolaor et al., A&AR 27, 5 (2019) — PTA review
25. Mingarelli et al., Nature Astronomy 1, 886 (2017) — individual SMBHB sources in PTA
26. Hellings & Downs, ApJ 265, L39 (1983) — Hellings-Downs correlation
27. Taylor et al., ApJ 819, L6 (2016) — PTA sky sensitivity
28. Hazboun et al., PRD 100, 104028 (2019) — hasasia PTA sensitivity curves
29. Sesana, Vecchio & Volonteri, MNRAS 394, 2255 (2009) — resolvable SMBHB sources
30. Kelley et al., MNRAS 471, 4508 (2017) — GWB from Illustris SMBHB population
31. Reardon et al. (PPTA), ApJL 951, L6 (2023) — PPTA GWB search

### Original Works
32. Lehmann, arXiv:2504.00506 (2025) — "Higher multipoles of the cow"
33. Tolman, Phys. Rev. 55, 364 (1939); Oppenheimer & Volkoff, Phys. Rev. 55, 374 (1939) — TOV equation

---

## Section 6: The Nanohertz Duck — PTA Indistinguishability (Full Specification)

### 6.1 Core Argument

A population of supermassive binary ducks (M ~ 10⁶–10⁹ M_☉) inspiraling in the nanohertz band produces a stochastic GW background **rigorously indistinguishable** from the signal observed by PTAs. The duck's permanent quadrupole Q^C modifies individual binary GW emission at 2PN order, where the correction is:

$$\frac{\delta(\dot{E})}{\dot{E}} \sim \kappa_{\rm duck}\left(\frac{\pi M_c f}{c^3}\right)^{4/3}$$

**Numerical estimates in the PTA band:**

| f (nHz) | v/c | (v/c)⁴ | δ_duck (κ=1) |
|---------|-----|---------|--------------|
| 1 | 0.005 | 6×10⁻¹⁰ | ~10⁻⁹ |
| 10 | 0.011 | 1.5×10⁻⁸ | ~10⁻⁸ |
| 100 | 0.024 | 3.3×10⁻⁷ | ~10⁻⁷ |

These corrections are 4–9 orders of magnitude below the ~30% measurement uncertainty on A_GWB.

### 6.2 Population-Level Degeneracy Theorem

**Theorem (informal):** For any SMBHB population model (mass function Φ(M), merger rate dn/dz, eccentricity e(f), mass ratio distribution p(q)) that fits the NANOGrav 15yr data, there exists a binary duck population with identical mass function and merger rate whose predicted h_c(f) agrees to better than 10⁻⁵.

**Proof sketch:** The GWB characteristic strain is:

$$h_c^2(f) = \frac{4G^{5/3}}{3\pi^{1/3}c^2}f^{-4/3}\int_0^{z_{\max}} dz\,\frac{dn}{dz}\,\frac{[\mathcal{M}_c(1+z)]^{5/3}}{(1+z)^{1/3}\,d_L^2(z)}\,[1 + \delta_{\rm duck}(f,z)]$$

Since δ_duck < 10⁻⁵ for all f in the PTA band, the integral with and without the duck correction differ by less than 10⁻⁵. The population parameters (Φ, dn/dz, p(q)) enter identically for ducks and black holes of the same mass. QED.

### 6.3 NANOGrav 15yr Comparison

The observed signal:
- A_GWB = (2.4 ± 0.7) × 10⁻¹⁵ at f_ref = 1/yr = 31.7 nHz
- Spectral index γ = 13/3 (power law h_c ∝ f^{−2/3})
- Hellings-Downs inter-pulsar correlations at 3.5–4σ

The duck prediction:
- h_c(f) = A_duck × (f/f_yr)^{−2/3} with A_duck = A_SMBHB = 2.4 × 10⁻¹⁵
- Spectral shape identical (circular GW-driven inspiral gives γ = 13/3 regardless of body shape)
- Hellings-Downs correlation identical (depends only on quadrupolar GW radiation pattern, same for all sources)

### 6.4 Electromagnetic Indistinguishability

At supermassive scales (M ~ 10⁹ M_☉):
- Duck physical size R ~ R_Schwarzschild ~ 20 AU
- At z = 0.1 (d_L ~ 450 Mpc): angular size θ ~ 6×10⁻⁸ arcsec
- EHT resolution: ~20 μarcsec → 5 orders of magnitude too coarse
- Gravitational lensing: multipole signature falls off as r^{−(ℓ+1)}, unresolvable beyond ~100 R

**Conclusion:** No electromagnetic observation can resolve the duck shape at cosmological distances.

### 6.5 Mock PTA Skymap Methodology

**Tools:** `hasasia` (PTA sensitivity), `healpy` (HEALPix sky maps), `holodeck` or custom Monte Carlo (population)

**Procedure:**
1. Generate duck binary population: draw N ~ 10⁶ sources from Φ(M) × p(q) × dn/dz
2. Assign sky positions uniformly (isotropic GWB)
3. Compute h₀ for each source: h₀ = (4/d_L)(G M_c/c²)^{5/3}(πf)^{2/3}/c
4. Identify N_res ~ 0–5 individually resolvable sources (h₀ > h_c/√N_f)
5. Create HEALPix strain map and render Mollweide projection
6. Overlay PTA sensitivity contours from `hasasia`

### 6.6 Figures for This Section

**Fig 11: Violin plot overlay**
- NANOGrav 15yr free-spectrum posteriors (gray violins, publicly available)
- Duck GWB power law (orange line): h_c = 2.4×10⁻¹⁵ × (f/f_yr)^{−2/3}
- "Duck-corrected" power law (dashed, visually identical to solid)
- Label: "The duck model is indistinguishable from the SMBHB model within NANOGrav credible intervals"

**Fig 12: Mock duck skymap**
- Mollweide projection (equatorial coordinates)
- Color: log(h_c) stochastic background, scale 10⁻¹⁶ to 10⁻¹⁴
- Markers: loudest 5 binary duck sources (duck icons)
- Contours: PTA 3σ sky sensitivity
- Inset: zoom on loudest source

### 6.7 The Punchline

> *We cannot distinguish a supermassive binary duck from a supermassive binary black hole based on PTA observations, electromagnetic observations, or any currently feasible measurement. The nanohertz gravitational wave background observed by NANOGrav, EPTA, PPTA, and CPTA is equally consistent with a cosmological population of inspiraling ducks.*

### 6.8 How Could One Break the Degeneracy? (Speculative)

For completeness and humor:
1. **GW memory from duck mergers:** nonlinear memory encodes full multipole structure, but requires SNR ~ 100+ (not yet achievable)
2. **QNM spectroscopy:** ringdown encodes duck spectrum (Topic 2), but PTA frequency resolution (Δf ~ 2 nHz) is too coarse
3. **Tidal disruption events:** duck-shaped tidal field differs from spherical, but monopole dominates at disruption radius
4. **Direct EHT imaging:** angular resolution 5 orders of magnitude insufficient
