# Higher Multipoles of the Duck: Gravitational Wave Physics of Anatidae

## Context

April 1st paper extending Lehmann (2025, arXiv:2504.00506) "Higher Multipoles of the Cow" to a rubber duck. Three rigorous research threads: (1) binary duck inspiral with higher-multipole PN corrections across all mass scales, (2) quasinormal modes of a duck-shaped self-gravitating body — full GR for stellar mass, Newtonian for supermassive, (3) equation of state requirements for duck matter. The premise is humorous; the methodology must be theorem-level rigorous or numerically sound.

**This document is a research specification — topic hunting and method details. No code implementation in this round.**

---

## Available Infrastructure

| Tool | Location | Capability |
|------|----------|------------|
| Multipole pipeline | `multipole-cow/code/` | mesh → Q_lm → Q^C → I → GW power (validated ~2%) |
| Duck mesh | `data/duck.off` | 4732 verts, 9460 faces, watertight, unit scale |
| AthenaK (main) | `/data/haiyangw/athenak/` | GR hydro (Z4c NR), Newtonian hydro, multigrid Poisson self-gravity |
| AthenaK (oscillation) | `/data/haiyangw/claude/oscillation/single-star/` | Polytropic star oscillations, Lane-Emden IC, validated FMG convergence |

---

## Phase 0: Duck Baseline (shared prerequisite)

Run the existing multipole pipeline on the duck mesh to establish all baseline quantities.

**Required outputs:**
- Q^C_duck (3×3 Cartesian quadrupole), I_duck (3×3 inertia tensor), principal moments
- Q_lm up to ℓ=10 (spherical multipoles)
- R_eq = (3V/4π)^{1/3} (equivalent sphere radius)
- ε_lm (SH surface deformation coefficients — reuse `render_sphere_to_cow.py` fitting pipeline)
- κ_duck = Q^C_eigenvalues / (M R_eq²) (shape-induced quadrupole parameter)
- GW power for single spinning duck (reuse `gw_radiation.py`)

**Verification:** Dipole vanishes at COM; Q^C traceless, symmetric; cross-check Q^C = Tr(I)δ − 3I.

---

## Topic 1: Binary Duck on Post-Newtonian Orbits

### 1.1 Physical Setup

Two ducks (masses M_A, M_B) on quasi-circular orbits. Each duck has a permanent body-frame quadrupole Q^C_body computed from the mesh. Results shown across three mass regimes:

| Regime | Mass scale | Detector band | Key physics |
|--------|-----------|---------------|-------------|
| Stellar | 1–3 M_☉ | LIGO/Virgo/KAGRA (10–1000 Hz) | NS-like, tidal deformability, nuclear EoS |
| Intermediate | 10²–10⁴ M_☉ | LISA/DECIGO gap | Transition regime |
| Supermassive | 10⁶–10⁹ M_☉ | LISA (0.1–100 mHz) | Long inspiral, large accumulated phase |

### 1.2 Quadrupole-Monopole Interaction

The Newtonian interaction between body A's quadrupole and body B's monopole (Poisson & Will, "Gravity", Eq. 9.41; Barker & O'Connell 1975):

$$U_{QM} = -\frac{G M_B}{2 r^3} Q^{C,A}_{ij} n_i n_j$$

where Q^C_A is duck A's quadrupole rotated to the inertial frame, **n** is the unit separation vector.

**PN order counting:** U_QM/E_N ~ (R/r)². For compact objects where R ~ GM/c², this is 2PN [O(v⁴/c⁴)]. For extended ducks with physical radius R_phys >> R_Schwarzschild, the effect is parametrically larger — it enters at an effective order (R_phys/r)² which can be O(1) for close separations.

**Orientation cases:**
- *Tidally locked:* duck body x-axis along **n** → Q_nn = const (dominant configuration for close binaries)
- *Freely spinning:* duck spins about principal axis at Ω_spin ≠ ω_orb → Q_nn oscillates → additional GW frequency components
- *Precessing:* general Euler angle evolution → rich phenomenology

### 1.3 Modified Orbital Dynamics

**Modified Kepler's law** (circular orbits):

$$\omega^2 = \frac{G M_{\rm total}}{r^3}\left[1 + \frac{15\,Q_{nn}^{\rm eff}}{2 M_{\rm total}\, r^2}\right]$$

**Binding energy** with quadrupole correction:

$$E(\omega) = E_{\rm Keplerian}(\omega)\left[1 + \delta_Q(\omega)\right]$$

where δ_Q encodes the fractional correction from the permanent quadrupoles. This is computed by inverting the modified Kepler relation to get r(ω), then substituting into E(r) = E_Kep(r) + U_QM(r).

**Energy balance and inspiral:** dE/dt = −P_GW governs the frequency evolution. The quadrupole correction modifies both E(ω) and P_GW(ω), giving a modified chirp rate dω/dt.

### 1.4 GW Emission: System Quadrupole

The total mass quadrupole of the binary system:

$$\mathcal{I}^{\rm sys}_{ij}(t) = \underbrace{\mu\left(x_i x_j - \tfrac{1}{3}r^2\delta_{ij}\right)}_{\rm orbital} + \underbrace{\sum_A \frac{1}{3}R_A(t)\,Q^{C,A}_{\rm body}\,R_A^T(t)}_{\rm body}$$

For tidally locked ducks, the body term rotates at ω_orb → contributes at frequency 2ω (same as orbital quadrupole) → coherent addition/subtraction.

**Power decomposition:**

$$P_{\rm GW} = \frac{G}{45 c^5}\langle\dddot{\mathcal{I}}^{\rm sys}_{ij}\dddot{\mathcal{I}}^{{\rm sys},ij}\rangle = P_{\rm orb} + P_{\rm body} + P_{\rm cross}$$

- P_orb: standard Peters-Mathews (leading order)
- P_body: from duck's permanent quadrupole alone
- P_cross: interference between orbital and body quadrupoles (can be positive or negative!)

### 1.5 Higher Multipole Contributions

**(a) Current quadrupole** S_ij (from spin angular momentum distribution):
- Enters at 0.5PN relative to mass quadrupole (Blanchet, Living Rev. Rel. 17, 2014)
- For a rigid duck spinning at Ω about principal axis: S_ij involves the inertia tensor
- Contributes GW power: P_current = G/(45c⁷)⟨S⃛_ij S⃛^ij⟩

**(b) Mass octupole** I_ijk (from Q_3^m):
- Also enters at 0.5PN relative
- For equal-mass binary: vanishes by symmetry (parity)
- For unequal-mass: extract from existing compute_Qlm(ell=3)
- Contributes: P_octupole = G/(189c⁷)⟨I⃛_ijk I⃛^ijk⟩

### 1.6 Phase Evolution

In the stationary phase approximation (matched filtering observable):

$$\Psi(f) = \Psi_{\rm pp}(f) + \delta\Psi_Q(f)$$

The quadrupole correction to the phase (Poisson 1998, PRD 57 5287):

$$\delta\Psi_Q \sim -\frac{75}{64\eta}\frac{\hat{Q}}{M r^2}\,v^{-1}$$

where Q̂ is a mass-weighted combination of both ducks' quadrupoles, η = μ/M is the symmetric mass ratio, and v = (πMf)^{1/3}.

**Key observable:** the number of dephasing cycles ΔN = δΨ/(2π) accumulated over the detector band. For stellar-mass ducks in LIGO band (10–1000 Hz), this determines detectability. For LISA-band supermassive ducks, the enormous number of cycles (~10⁵) makes even tiny fractional corrections measurable.

### 1.7 Duck Tidal Deformability

**(a) Shape-induced quadrupole parameter** (analogous to spin-induced κ for NS):

$$\kappa_{\rm duck} = \frac{Q^C_{\rm eigenvalue}}{M R_{\rm eq}^2}$$

Compare: Kerr BH has κ = 1; neutron stars have κ ~ 2–14 (Laarakkers & Poisson, ApJ 512, 1999).

**(b) Elastic tidal Love number** k_2 (for deformable rubber duck):

$$k_2 = \frac{3/2}{1 + 19\mu_{\rm shear}/(2\rho g R)}$$

For rubber: E ~ 0.01 GPa, ν ~ 0.5 → μ_shear ~ 3.3 MPa.

**(c) Dimensionless tidal deformability** (LIGO observable):

$$\Lambda = \frac{2}{3}k_2\left(\frac{c^2 R}{G M}\right)^5$$

This is directly comparable to the Λ measured from GW170817 (Λ ~ 100–800 for NS).

### 1.8 Method Summary

| Step | Input | Method | Output |
|------|-------|--------|--------|
| Body moments | duck.off | Existing pipeline | Q^C, I, Q_lm, κ_duck |
| U_QM(r) | Q^C, M, r, orientation | Analytical (Poisson & Will) | Interaction energy |
| ω(r), E(r) | U_QM | Modified Kepler inversion | Orbital dynamics |
| I^sys(t) | Orbital params + Q^C | Analytical rotation | Time-dependent quadrupole |
| P_GW | I^sys(t) | Numerical 3rd derivative (existing method) | Power decomposition |
| Ψ(f) | E(ω), P_GW(ω) | Energy balance integration | Phase evolution |
| h_+, h_× | Ψ(f), amplitudes | SPA or time-domain ODE | Waveform |
| k_2, Λ | Rubber material properties | Elastic sphere formula | Tidal deformability |

### 1.9 Figures & Tables

| # | Content |
|---|---------|
| Fig 1 | Binary duck schematic: two ducks on circular orbit with labeled axes and quadrupole ellipsoids |
| Fig 2 | |δE/E_N| vs r/R_eff (log-log): tidally locked, freely precessing, point-particle |
| Fig 3 | GW power decomposition vs f: P_orb, P_body, P_cross, P_current (log-log), each labeled with PN order |
| Fig 4 | Dephasing ΔN(f) in cycles: stellar-mass (LIGO), intermediate, supermassive (LISA) — three panels |
| Fig 5 | Time-domain h_+(t): duck binary vs point-particle, with inset showing phase difference |
| Fig 6 | Λ_duck vs rubber Young's modulus, overlaid with Λ_NS bands from GW170817 |
| Tab 1 | Duck body parameters: M, V, R_eq, I_1, I_2, I_3, Q^C eigenvalues, κ_duck per axis |
| Tab 2 | GW power breakdown at r = 5, 10, 20 R_eff |
| Tab 3 | PN coefficient comparison: standard terms vs duck quadrupole correction |
| Tab 4 | k_2 and Λ for rubber duck vs NS (SLy, APR) vs BH (k_2 = 0) |

### 1.10 Verification

- V1: Q^C → 0 recovers Peters-Mathews P = (32/5)(G⁴/c⁵)(μ²M³)/r⁵
- V2: Spheroidal Q^C = diag(Q,−Q/2,−Q/2) matches Poisson & Will Eq. 9.41 analytically
- V3: Equal-mass + parity → mass octupole vanishes
- V4: Energy conservation dE/dt = −P_GW along inspiral (numerical check)
- V5: P_body/P_orb scales as (R/r)² (dimensional analysis)
- V6: r → ∞ limit with single spinning duck recovers existing validated GW power
- V7: SPA phase agrees with time-domain integration to O(1/N_cycles)
- V8: All physical units consistent (CGS throughout, as in existing gw_radiation.py)

---

## Topic 2: Quasinormal Modes of a Self-Gravitating Duck

### 2.1 Physical Picture

A duck-shaped self-gravitating body has natural oscillation modes (QNMs). Unlike black hole QNMs (vacuum perturbations of Kerr, solved via Teukolsky), these are *stellar* QNMs — coupled fluid + gravitational perturbations, closer to neutron star oscillation theory.

**Two regimes** (per user specification):
- **Stellar mass (1–3 M_☉): Full GR** — Z4c formalism in AthenaK with GR hydrodynamics
- **Supermassive (10⁶+ M_☉): Newtonian** — self-gravitating gas cloud with multigrid Poisson, as in the existing polytropic star setup

**Central question:** How does the duck's non-spherical geometry lift the (2ℓ+1)-fold degeneracy of each QNM multiplet? This is the gravitational analog of Zeeman splitting in atomic physics — the duck shape plays the role of the magnetic field.

### 2.2 Tier A: AthenaK Time-Domain Evolution (Primary Method)

This is the most rigorous approach: evolve the full nonlinear equations and extract modes from the ringdown.

#### 2.2a Stellar Mass Duck (Full GR)

**Framework:** AthenaK Z4c + GR hydro (Zhu et al. 2024, Fields et al. 2024)

**Initial data construction:**
1. Start from a TOV solution for a polytropic star (Γ = 2, K chosen for M ~ 1.4 M_☉, R ~ 10 km)
2. Deform the density profile using duck SH coefficients:
   $$\rho(r, \theta, \phi) = \rho_{\rm TOV}\!\left(\frac{r}{R(\theta,\phi)}\right), \quad R(\theta,\phi) = R_0\!\left[1 + \sum_{\ell,m} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right]$$
3. Solve the constraint equations (Hamiltonian + momentum) on the deformed data using the elliptic solver, or use a conformal thin-sandwich approach
4. Set pressure from the EoS: P = Kρ^Γ

**Problem generator:** New `pgen/duck_star_gr.cpp` based on existing `pgen/dyngr_tov.cpp`
- Read ε_lm from file (output of Phase 0)
- Construct R(θ,φ) on the AthenaK grid
- Map TOV profile to duck shape
- Apply small velocity perturbation to excite modes

**Evolution:**
- Z4c gauge: 1+log lapse, Gamma-driver shift
- GR hydro with HLLC Riemann solver
- AMR: refine around the duck surface (steep density gradient)
- Evolve for ~100 dynamical times to resolve mode frequencies

**Wave extraction:**
- Extract Ψ₄ at fixed extraction radii using spin-weighted spherical harmonic decomposition (infrastructure exists in JAX-NR: `wave_extraction/weyl_scalar_extraction.py` — adapt to AthenaK output format)
- Decompose into (ℓ,m) components
- Each (ℓ,m) channel: fit ω_R + i ω_I (real frequency + damping rate)

**Alternative: SpECTRE ringdown pipeline** (`spectre/Evolution/Ringdown/`) could be adapted for post-processing AthenaK output if the data format is converted.

#### 2.2b Supermassive Duck (Newtonian)

**Framework:** AthenaK Newtonian hydro + multigrid Poisson self-gravity (existing validated setup)

**Initial data:** Lane-Emden polytrope (n=1.5, Γ=5/3) deformed by duck ε_lm coefficients, exactly as in the existing `polytropic_star.athinput` but with non-spherical initial profile.

**Problem generator:** New `pgen/duck_star_newt.cpp` based on existing `pgen/polytropic_star.cpp` in `/data/haiyangw/claude/oscillation/single-star/`

**Mode extraction:** FFT of density/velocity time series at multiple points → frequency peaks → angular decomposition for (ℓ,m) classification.

**Key physics difference:** No gravitational wave damping in Newtonian gravity. The modes are purely oscillatory (real frequencies). Damping comes only from numerical viscosity or explicit dissipation terms. The GW damping rate must be computed separately (see §2.4).

### 2.3 Tier B: FEM Eigenvalue Problem (Complementary)

Solve the scalar Helmholtz equation on the duck interior as a membrane model:

$$-\nabla^2 \psi_n = \omega_n^2\,\psi_n, \quad \psi\big|_{\partial\mathcal{D}} = 0$$

**Method:**
- Volumetric tetrahedralization of duck interior (tetgen via meshpy.tet or pygalmesh)
- P1 FEM: assemble stiffness K and mass M matrices
- Generalized eigenvalue problem: K ψ = ω² M ψ
- Solve with scipy.sparse.linalg.eigsh (shift-invert for smallest eigenvalues)
- Classify modes by projecting eigenvectors onto Y_lm (reusing `build_sh_matrix()`)

**Value:** Fast, mesh-resolution controlled, gives the mode splitting pattern without any fluid/GR physics. Isolates the purely geometric effect of the duck shape on the eigenspectrum.

### 2.4 Tier C: Perturbative Mode Splitting (Analytical)

Treat the duck as a perturbed sphere and compute first-order frequency shifts.

**Duck boundary parametrization:**

$$R_{\rm duck}(\theta,\phi) = R_0\left[1 + \sum_{\ell\geq 2,m} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right]$$

where ε_lm are extracted from the existing SH fitting pipeline.

**First-order degenerate perturbation theory:** Within each unperturbed (n,ℓ) multiplet (2ℓ+1 degenerate states), construct and diagonalize the perturbation matrix:

$$V_{m,m'} = \sum_{\ell'',m''} \varepsilon_{\ell'',m''}\,\langle n,\ell,m | \hat{V}_{\ell'',m''} | n,\ell,m'\rangle$$

**Selection rules** (from angular momentum coupling):
- V is nonzero only when |ℓ − ℓ''| ≤ ℓ ≤ ℓ + ℓ'' and m' = m + m''
- The ℓ'' = 2 (quadrupole) deformation couples each ℓ-multiplet to itself and to ℓ ± 2
- Higher ℓ'' deformations give higher-order couplings

**Matrix elements:** Computed via Wigner 3j symbols (from sympy.physics.wigner or py3nj). The radial part involves overlap integrals of unperturbed eigenfunctions with the boundary perturbation.

**Output:** Predicted frequency shifts δω_{n,ℓ,m} for each mode, expressed as linear combinations of ε_lm. This gives analytical insight into which duck features (bill, tail, body) drive the splitting.

### 2.5 Emission to Infinity

Each QNM, when excited, radiates gravitational waves. The key question: how does the duck's broken symmetry make the radiation anisotropic?

**For a sphere:** all m-substates of a given (n,ℓ) have the same damping rate τ_GW. The radiated power in mode (ℓ,m) is:

$$\dot{E}_{\ell,m} = \frac{G}{c^{2\ell+1}}\,\omega^{2\ell+2}\,|\delta Q_\ell^m|^2 \times C_\ell$$

where δQ_ℓ^m is the time-varying multipole from the oscillation, and C_ℓ = 4π(ℓ+1)(ℓ+2)/[ℓ(ℓ−1)(2ℓ+1)!!²].

**For the duck:** the broken symmetry means:
1. Frequencies split: ω_{n,ℓ,m} ≠ ω_{n,ℓ,m'} → different radiation frequencies
2. Mode shapes change: |δQ_ℓ^m| differs between m-substates → different amplitudes
3. Mode mixing: a duck mode nominally labeled (ℓ,m) has admixtures of other (ℓ',m') → radiation in multiple multipole channels

**Ringdown waveform:** After excitation (e.g., from a perturbation or merger), the GW signal is:

$$h(t) = \sum_{n,\ell,m} A_{n\ell m}\,e^{-t/\tau_{n\ell m}}\,\cos(\omega_{n\ell m}\,t + \phi_{n\ell m})\,{}_{-2}Y_{\ell m}(\iota,\varphi)$$

For the duck, the number of distinct (ω, τ) pairs is larger than for a sphere, and the mode amplitudes A_{nℓm} depend on the excitation geometry relative to the duck's shape.

**Extraction from AthenaK:**
- GR runs: extract Ψ₄(t) at multiple (ℓ,m) → fit each channel to damped sinusoids → measure (ω_R, ω_I) per mode
- Newtonian runs: extract density perturbation δρ(t) → project onto Y_lm → FFT for frequencies (damping rates computed separately from GW formula)

### 2.6 Connection to Topic 3

The relaxation of a duck-shaped star toward spherical symmetry (when only isotropic pressure is present) IS the excitation of QNMs. The characteristic timescale for losing duck features at each ℓ is set by:
- The f-mode frequency at that ℓ (oscillation timescale ~ 1/ω_f)
- The GW damping time τ_GW (for compact/stellar-mass ducks)
- The viscous damping time τ_visc (for supermassive gas-cloud ducks)

This provides a natural bridge: Topic 2 gives the frequencies, Topic 3 asks what internal physics prevents the relaxation.

### 2.7 Figures & Tables

| # | Content |
|---|---------|
| Fig 7 | AthenaK density snapshots: duck-shaped star at t = 0, t = T/4, t = T/2 (GR and Newtonian side by side) |
| Fig 8 | FFT power spectrum: duck vs sphere, showing mode splitting (Newtonian) |
| Fig 9 | Ψ₄ ringdown: (ℓ,m) = (2,0), (2,±1), (2,±2) channels showing distinct frequencies (GR) |
| Fig 10 | Grotrian diagram: sphere (degenerate levels) → duck (split), lines colored by ℓ |
| Fig 11 | 3D mode shape visualization: f-mode ℓ=2 with m = −2,−1,0,+1,+2 on duck mesh |
| Fig 12 | Complex frequency plane: ω_R vs ω_I for first 15 modes, duck vs sphere |
| Tab 5 | Duck deformation coefficients ε_lm up to ℓ = 10 |
| Tab 6 | First 20 eigenfrequencies: duck vs sphere (Helmholtz, FEM) |
| Tab 7 | f-mode splitting comparison: Tier A (AthenaK GR) vs Tier A (AthenaK Newt.) vs Tier B (FEM) vs Tier C (perturbative) |
| Tab 8 | GW damping times τ_{n,ℓ,m} for the duck, compared with spherical star |

### 2.8 Verification

- **Sphere limit:** Set ε_lm = 0 → recover (2ℓ+1)-fold degeneracy and known eigenvalues (Bessel zeros for Helmholtz; tabulated f/p-mode values for fluid modes)
- **Ellipsoid:** Known mode splitting from Lamb (1881) → verify perturbative and numerical results
- **Cross-tier:** Tier C (perturbative) ≈ Tier A/B (numerical) for near-spherical deformations
- **AthenaK convergence:** Resolution study (refine grid by 2×, verify ω converges)
- **GR constraints:** Monitor Hamiltonian and momentum constraints throughout evolution
- **Energy conservation:** Total energy (kinetic + potential + internal) conserved to truncation error

---

## Topic 3: Equation of State of Duck Matter

### 3.1 The Central Question

What internal physics is required for a self-gravitating body to maintain a duck shape against gravitational collapse toward spherical symmetry?

The answer depends dramatically on mass scale:

| Mass | ρ_mean | Support mechanism | Key physics |
|------|--------|-------------------|-------------|
| 1 M_☉ | ~10¹⁴ g/cm³ | Nuclear degeneracy | Needs anisotropic stress ~10⁵× nuclear |
| 10³ M_☉ | ~10⁸ g/cm³ | Electron degeneracy + radiation | White dwarf regime |
| 10⁶ M_☉ | ~18 g/cm³ | Radiation pressure | Supermassive star |
| 10⁸ M_☉ | ~10⁻⁴ g/cm³ | Radiation pressure | Gas cloud, L ~ L_Edd |
| 10⁹ M_☉ | ~10⁻⁶ g/cm³ | Radiation pressure | Near-vacuum, GR unstable |

### 3.2 Compactness Analysis

**Global compactness:** C = GM/(R_eq c²)

**Bounds:**
- Buchdahl (1959): C < 4/9 ≈ 0.444 for isotropic perfect fluid
- Schwarzschild: C < 1/2 (black hole forms)
- Anisotropic generalization (Bowers & Liang 1974): bound modified by anisotropic stress

**Duck-specific: directional compactness.** Define C(θ,φ) = GM/(R(θ,φ) c²). The duck's bill has the smallest radial extent → highest local compactness. As M increases:

1. The bill region violates Buchdahl first → "bill collapse"
2. The body remains subcritical longer
3. At sufficiently high M, the entire duck collapses

**Method:** Compute R(θ,φ) from the SH surface expansion, then map C(θ,φ) over the surface for each mass scale. Identify the critical mass M_crit(θ,φ) where each direction first violates Buchdahl:

$$M_{\rm crit}(\theta,\phi) = \frac{4}{9}\frac{R(\theta,\phi)\,c^2}{G}$$

The global critical mass is M_crit = min over all angles.

### 3.3 TOV Analysis: Spherical Baseline

For the equivalent sphere (monopole approximation): what EoS supports mass M at radius R_eq?

**TOV equation** (full GR):

$$\frac{dP}{dr} = -\frac{(\rho + P/c^2)(m(r) + 4\pi r^3 P/c^2)\,G}{r^2(1 - 2Gm(r)/(rc^2))}$$

**Newtonian limit** (supermassive regime where C ≪ 1):

$$\frac{dP}{dr} = -\frac{Gm(r)\rho(r)}{r^2}$$

**Uniform density analytical solution:**
- P_c = (2π/3) G ρ² R² (Newtonian)
- P_c = ρc² [√(1−2C) − √(1−2C(r/R)²)] / [√(1−2C(r/R)²) − √(1−2C)] (GR, Schwarzschild interior)

**EoS families to test:**

| EoS | Formula | Physical context | Regime |
|-----|---------|-----------------|--------|
| Polytrope (stiff) | P = Kρ^{5/3} | Non-relativistic degenerate | Stellar |
| Polytrope (soft) | P = Kρ^{4/3} | Relativistic degenerate | Stellar |
| SLy (tabulated) | Numerical table | Nuclear matter, soft | NS |
| APR (tabulated) | Numerical table | Nuclear matter, stiff | NS |
| MIT bag | P = (ρ−4B)c²/3 | Strange quark matter | Exotic NS |
| Radiation | P = aT⁴/3 | Photon gas | Supermassive |
| Ideal gas | P = nkT | Thermal support | Supermassive |

**Method:** Integrate TOV for each EoS family, producing M(R) curves. Overlay the duck's R_eq(M) locus (assuming self-similar scaling of duck shape). Identify which EoS families can support a duck at each mass scale.

### 3.4 Anisotropic Pressure Support — The "Duck Matter Problem"

**The key physics:** A spherical star needs only isotropic pressure P. A duck-shaped star requires anisotropic stress to resist gravitational relaxation toward a sphere. The stress tensor must have:

$$T^{ij} = P\,\delta^{ij} + \Pi^{ij}$$

where Π^{ij} is the anisotropic (deviatoric) stress that maintains the duck shape.

**Required anisotropic stress** for deformation mode ε_lm:

At leading order (perturbative, Newtonian), the gravitational potential perturbation from the shape deformation is:

$$\delta\Phi_\ell(r) = -\frac{4\pi G\rho}{2\ell+1}\left[\frac{r^{\ell+2}}{\ell+3} - \frac{R^{2\ell+3}}{(\ell+3)r^{\ell+1}}\right]\varepsilon_{\ell m}$$

The required pressure perturbation to maintain hydrostatic equilibrium:

$$\delta P_\ell \sim \rho\,\frac{\partial\delta\Phi_\ell}{\partial r}\,\varepsilon_{\ell m}\,R$$

The shear stress required scales as:

$$\sigma_{\rm required} \sim G\rho^2 R^2\,\varepsilon_{\ell m}$$

**Connection to neutron star mountains** (Ushomirsky, Cutler & Bildsten, MNRAS 319, 2000; Johnson-McDaniel & Owen, PRD 88, 2013):

Maximum sustainable ellipticity of a neutron star crust:

$$\varepsilon_{\rm max} \sim \frac{\sigma_{\rm break}}{G\rho^2 R^2} \sim 10^{-7}\text{ to }10^{-5}$$

where σ_break is the breaking strain of nuclear pasta/crystal.

**The duck has ε ~ 0.3** (from aspect ratio and SH decomposition). Therefore:

$$\frac{\sigma_{\rm required}}{\sigma_{\rm nuclear}} \sim \frac{\varepsilon_{\rm duck}}{\varepsilon_{\rm NS,max}} \sim \frac{0.3}{10^{-6}} \sim 3 \times 10^5$$

**This is the central result:** *Duck matter must be approximately 10⁵ times stronger than the strongest known material in the universe (nuclear pasta).* This "duck matter problem" is the comedic-but-rigorous punchline of Topic 3.

**Possible exotic resolutions** (to be discussed tongue-in-cheek):
- Quark matter with enhanced shear modulus (σ ~ 10× nuclear) — insufficient by ~10⁴
- Crystalline color-superconducting quark matter — possibly sufficient (Mannarelli et al. 2007)
- Cosmic strings woven into duck shape — beyond known physics
- Magnetic field support (magnetar fields B ~ 10¹⁵ G) — B²/(8π) ~ 4×10²⁹ dyne/cm², comparable to nuclear but still insufficient
- Rapid rotation (centrifugal support) — oblate, not duck-shaped

### 3.5 Supermassive Duck

At M ~ 10⁶–10⁹ M_☉, the physics shifts entirely from nuclear to radiation/thermal.

**Mean density** (at equivalent sphere radius, assuming duck scales self-similarly with mass):

$$\bar{\rho} = \frac{3M}{4\pi R_{\rm eq}^3}$$

For R_eq scaling with the Schwarzschild radius (maximally compact): ρ ∝ M⁻². Actual supermassive objects have R >> R_s, so densities are even lower.

**Radiation-supported duck (supermassive star analog):**

Central temperature from radiation pressure balance:

$$P_c = \frac{2\pi}{3}G\rho^2 R^2 = \frac{1}{3}aT_c^4$$

$$T_c = \left(2\pi G\rho^2 R^2 / a\right)^{1/4}$$

Luminosity ~ L_Edd = 4πGMm_pc/σ_T = 1.26 × 10³⁸ (M/M_☉) erg/s.

**Directional Eddington limit:** For the duck, the effective gravity g_eff(θ,φ) and the projected area vary with direction:

$$L_{\rm Edd}(\theta,\phi) \propto R(\theta,\phi)^2\,g_{\rm eff}(\theta,\phi)$$

The bill has smaller R → lower L_Edd → **the bill blows off first** under radiation pressure.

**GR instability (Chandrasekhar 1964):** A radiation-dominated body has adiabatic index Γ = 4/3, which is marginally stable in Newtonian gravity but unstable in GR:

$$\Gamma_{\rm crit} = \frac{4}{3} + K\,\frac{GM}{Rc^2}$$

where K ≈ 1.12 for uniform density. Any GR correction destabilizes the supermassive duck. The duck shape (non-radial perturbation modes) provides additional instability channels beyond the radial Chandrasekhar mode.

**AthenaK validation:** Evolve a supermassive duck with Newtonian hydro + self-gravity + simple cooling prescription:
- Does it maintain shape with anisotropic stress? (verify EoS requirements)
- How fast does it relax to sphere with isotropic pressure? (connects to Topic 2 QNMs)
- At what point does it fragment? (Jeans instability in bill/tail?)

### 3.6 Figures & Tables

| # | Content |
|---|---------|
| Fig 13 | Local compactness map: 3D duck surface colored by C(θ,φ) at M = 1.4 M_☉ |
| Fig 14 | Mass-radius diagram: TOV curves (SLy, APR, MIT bag) + duck locus + Buchdahl/Schwarzschild lines |
| Fig 15 | Required shear stress vs mass (log-log): σ_nuclear, σ_quark, σ_exotic lines + duck requirement |
| Fig 16 | Duck matter phase diagram: (ρ, T) plane with regions for nuclear, quark, radiation, exotic |
| Fig 17 | Supermassive duck fate: M vs ρ with stability boundaries (Chandrasekhar, Jeans, Eddington) |
| Fig 18 | Directional Eddington map: Mollweide projection of L_Edd(θ,φ)/L_Edd(sphere) — bill region highlighted |
| Tab 9 | Compactness at key mass scales: 1.4, 2.0, 10, 10³, 10⁶, 10⁹ M_☉ + critical mass for bill collapse |
| Tab 10 | Required P_c for each EoS at M = 1.4 M_☉ duck |
| Tab 11 | Duck ε_lm vs NS mountain ε_max for nuclear/quark/exotic matter (the "duck matter gap") |
| Tab 12 | Supermassive duck: T_c, L/L_Edd, τ_KH, Γ_crit for M = 10⁶, 10⁷, 10⁸, 10⁹ M_☉ |

### 3.7 Verification

- Sphere limit: recover Buchdahl bound, known TOV M_max for each EoS
- Known NS benchmarks: SLy → M_max ~ 2.05 M_☉, R_1.4 ~ 11.7 km
- Newtonian scaling: P_c = (2π/3)Gρ²R² (verified analytically)
- Supermassive: T_c and L agree with Hoyle & Fowler (1963) for supermassive stars
- AthenaK: duck relaxation timescale matches f-mode period from Topic 2

---

## Cross-Topic Connections

The three topics are not independent — they form a coherent narrative:

```
Topic 3 (EoS)           Topic 2 (QNMs)            Topic 1 (Binary)
     │                       │                          │
     │  "What holds the      │  "How does it            │  "What GW signal
     │   duck together?"     │   vibrate?"              │   does it produce?"
     │                       │                          │
     ├── Required σ for  ────┤── Relaxation to      ────┤── Permanent Q^C
     │   duck shape          │   sphere = f-mode        │   modifies inspiral
     │                       │   excitation             │
     ├── Compactness     ────┤── GR vs Newt. QNMs   ────┤── Tidal deformability
     │   determines regime   │   depend on C            │   Λ depends on EoS
     │                       │                          │
     └── Supermassive    ────┘── Duck-to-sphere      ────┘── LISA-band ducks
         duck fate               relaxation GWs           with long inspirals
```

**Key unifying result:** The GW signature of a binary duck system encodes information about the duck's internal structure (EoS, via Λ) AND its shape (via κ_duck and QNM spectrum). A sufficiently sensitive GW detector could, in principle, distinguish a duck from a cow.

---

## Implementation Roadmap

### Stage 1: Duck baseline (Python, builds on existing code)
- Run existing pipeline on duck.off → all moments, geometry
- Extract ε_lm deformation spectrum
- Single spinning duck GW power

### Stage 2: Binary duck (Python, new modules)
- Orbital dynamics with Q-M interaction
- System quadrupole and GW power decomposition
- Phase evolution and waveforms across mass scales
- Love number and tidal deformability

### Stage 3: EoS analysis (Python, new modules)
- Compactness analysis and local compactness map
- TOV solver with EoS library
- Anisotropic stress requirements ("duck matter problem")
- Supermassive duck: radiation support, Eddington, GR stability

### Stage 4: QNMs — Newtonian (AthenaK + Python)
- New pgen: duck-shaped polytrope with Newtonian self-gravity
- Based on existing `/data/haiyangw/claude/oscillation/single-star/polytropic_star.athinput`
- Evolve, extract time series, FFT for mode frequencies
- Compare with FEM Helmholtz eigenvalues and perturbative predictions

### Stage 5: QNMs — Full GR (AthenaK + Python)
- New pgen: duck-shaped TOV star with Z4c evolution
- Based on existing `pgen/dyngr_tov.cpp`
- Ψ₄ extraction and (ℓ,m) decomposition
- GW damping rates per mode

### Stage 6: Integration
- Unified driver and figure generation
- Cross-tier QNM comparison table
- Update RESEARCH_NOTE.md and CLAUDE.md
- Draft paper outline

---

## Key References

**Binary PN dynamics:**
1. Poisson & Will, "Gravity: Newtonian, Post-Newtonian, Relativistic" (Cambridge, 2014) — Ch. 9 (multipole interactions), Ch. 10 (PN EOM), Ch. 11 (GW emission)
2. Blanchet, Living Rev. Rel. 17, 2 (2014) — PN waveform coefficients, multipole formalism
3. Poisson, PRD 57, 5287 (1998) — spin-induced quadrupole in GW phase
4. Mikoczi, Vasuth & Gergely, PRD 71, 124043 (2005) — generalized spin orientations
5. Barker & O'Connell, PRD 12, 329 (1975) — quadrupole-monopole interaction in GR
6. Flanagan & Hinderer, PRD 77, 021502 (2008) — tidal Love numbers in GW
7. Hinderer, ApJ 677, 1216 (2008) — Love number computation
8. Laarakkers & Poisson, ApJ 512, 282 (1999) — spin-induced quadrupole of NS

**QNMs and stellar oscillations:**
9. Kokkotas & Schmidt, Living Rev. Rel. 2, 2 (1999) — QNM review (stars and BH)
10. Andersson & Kokkotas, MNRAS 299, 1059 (1998) — f-mode and w-mode universality
11. Cowling, MNRAS 101, 367 (1941) — the Cowling approximation
12. Lamb, Proc. London Math. Soc. 13, 189 (1881) — oscillations of ellipsoidal figures
13. Chandrasekhar, "Ellipsoidal Figures of Equilibrium" (Yale, 1969)

**AthenaK and numerical methods:**
14. Stone et al. (2024) — AthenaK framework
15. Zhu et al. (2024) — AthenaK NR solver
16. Fields et al. (2024) — GR hydro/MHD in dynamical spacetimes

**EoS, structure, and NS mountains:**
17. Tolman, Phys. Rev. 55, 364 (1939); Oppenheimer & Volkoff, Phys. Rev. 55, 374 (1939) — TOV
18. Buchdahl, Phys. Rev. 116, 1027 (1959) — compactness bound
19. Ushomirsky, Cutler & Bildsten, MNRAS 319, 902 (2000) — NS mountains
20. Johnson-McDaniel & Owen, PRD 88, 044004 (2013) — maximum ellipticity
21. Bowers & Liang, ApJ 188, 657 (1974) — anisotropic stars
22. Mannarelli, Rajagopal & Sharma, PRD 76, 074026 (2007) — crystalline color superconductor
23. Read et al., PRD 79, 124032 (2009) — parametrized EoS

**Supermassive stars:**
24. Hoyle & Fowler, MNRAS 125, 169 (1963) — supermassive stars
25. Chandrasekhar, ApJ 140, 417 (1964) — GR instability (Γ = 4/3)

**Original work:**
26. Lehmann, arXiv:2504.00506 (2025) — "Higher multipoles of the cow"

---

## Fun Ratings

### Topic 1: Binary Duck on PN Orbits

| Sub-topic | Fun | Notes |
|-----------|-----|-------|
| Quadrupole-monopole interaction | 6/10 | Solid physics, but essentially textbook PN with Q^C plugged in. The duck adds flavor but the math is Poisson & Will chapter 9 verbatim. |
| Tidally locked vs precessing orientations | 7/10 | Imagining two ducks tidally locked bill-to-bill on a circular orbit is delightful. The precessing case ("tumbling ducks") is even better. |
| GW power decomposition (P_orb + P_body + P_cross) | 7/10 | The cross-term P_cross being *negative* (destructive interference between the duck's body and its orbit) is a genuinely fun result to compute. |
| Phase dephasing across LIGO/LISA bands | 5/10 | Important for the paper's rigor, but it's a standard SPA integral. Less visually amusing. |
| Duck tidal deformability / Love number | 9/10 | Comparing a rubber duck's Lambda to GW170817 NS constraints is peak April 1st energy. "The rubber duck is 10^12 times more deformable than a neutron star" is a line that writes itself. |
| **Topic 1 overall** | **7/10** | Good backbone for the paper. Rigorous and publishable, moderate fun. The Love number comparison is the highlight. |

### Topic 2: Quasinormal Modes of a Self-Gravitating Duck

| Sub-topic | Fun | Notes |
|-----------|-----|-------|
| "Gravitational Zeeman splitting" framing | 10/10 | This analogy is *beautiful*. The duck shape as a symmetry-breaking field lifting QNM degeneracy --- this alone could be a talk title. |
| AthenaK full GR evolution of a duck-shaped star | 9/10 | Running a production NR code on a duck is inherently absurd and wonderful. The density snapshots (Fig 7) will be iconic. |
| Watching the duck relax to a sphere | 10/10 | The AthenaK simulation of a duck losing its duckness --- bill smoothing out, tail disappearing --- is the most visually compelling result in the entire paper. *The duck melts under its own gravity.* |
| FEM Helmholtz eigenvalues (Tier B) | 4/10 | Technically clean but not funny. It's just a Laplacian eigenvalue problem on a funny domain. |
| Perturbative splitting via Wigner 3j (Tier C) | 6/10 | Elegant physics, but the fun is abstract. The "which duck features drive the splitting" angle (bill vs tail vs body) saves it. |
| Grotrian diagram (Fig 10) | 8/10 | A Grotrian energy-level diagram for a *duck* is an instant classic figure for talks. |
| m-dependent damping rates | 7/10 | "The duck rings down at different rates depending on which direction you poke it" is a nice punchline. |
| **Topic 2 overall** | **8.5/10** | The strongest topic for visual impact, talk material, and genuine novelty. The Zeeman analogy and the melting duck simulation are the stars. |

### Topic 3: Equation of State of Duck Matter

| Sub-topic | Fun | Notes |
|-----------|-----|-------|
| "Bill collapse" --- the bill violates Buchdahl first | 9/10 | The local compactness map with the bill glowing red as the first region to collapse is hilarious and rigorous. |
| The "duck matter problem" (sigma ~ 10^5 x nuclear) | 10/10 | This is THE punchline of the paper. "Maintaining a duck shape at nuclear density requires matter 100,000 times stronger than anything in the known universe." Instant classic. |
| Exotic resolutions (cosmic strings, crystalline color superconductor) | 9/10 | Deadpan enumeration of increasingly exotic matter models that *still can't hold a duck together* is comedic gold. |
| TOV analysis with EoS families | 4/10 | Necessary scaffolding, but it's standard NS structure calculation with a duck label. |
| Mass-radius diagram with duck locus | 6/10 | Nice figure, but M-R diagrams are so ubiquitous in NS physics that even a duck on one feels routine. |
| Supermassive duck: "the bill blows off first" | 10/10 | The directional Eddington limit causing the bill to ablate under radiation pressure is *perfect*. The Mollweide projection (Fig 18) showing the bill highlighted as the weak point is a top-tier figure. |
| GR instability of a radiation-supported duck | 7/10 | Chandrasekhar meets Anatidae. Solid but the Gamma = 4/3 + O(C) story is well-trodden. |
| **Topic 3 overall** | **8/10** | The duck matter problem and bill-blows-off-first are both potential paper-title-level results. Strong comedic payoff. |

### Overall Rankings

| Rank | Topic | Fun | Why |
|------|-------|-----|-----|
| 1 | **Topic 2: QNMs** | 8.5/10 | Gravitational Zeeman splitting + watching the duck melt in AthenaK = unbeatable combo of novelty and visual impact |
| 2 | **Topic 3: EoS** | 8/10 | "Duck matter problem" is the best single punchline; bill collapse and bill-blows-off are runner-ups |
| 3 | **Topic 1: Binary** | 7/10 | Most rigorous and publishable, but less inherently funny. The Love number comparison is the standout moment |

**Best single result across all topics:** The duck matter problem (Topic 3, section 3.4) --- sigma_required / sigma_nuclear ~ 3 x 10^5. It's quantitative, surprising, and the deadpan exotic-matter discussion elevates it.

**Best visual:** AthenaK duck melting to a sphere (Topic 2, Fig 7). A time-lapse of a self-gravitating duck losing its features under its own gravity is the figure that gets shared on social media.

