# Multipole Expansion of the Cow --- Reproduction Study

## Preamble: Guidelines for These Notes

*Read this section whenever updating.*

### What These Notes Should Be

Treat this as a publishable research paper, except:

- **Uncertainty is explicit** --- use markers (`[HYPOTHESIS]`, `[PRELIMINARY]`, `[SOLID]`), not hedging prose
- **Gaps are visible** --- mark them (`[BLOCKING]`, `[FUTURE]`), don't smooth over
- **Sections can be unbalanced** --- developed where we have results, skeletal where we don't
- **Abandoned paths are documented** --- in appendix, not omitted
- **Structure is provisional** --- may need revision as understanding evolves

**A good research note = a publishable paper with explicit, fully-enumerated holes.**

### Bidirectional Criterion

- **Forward:** Every marker, if resolved, should advance the paper
- **Backward:** Every loose end should be captured by a marker

If you see an unmarked gap, add a marker. If a marker wouldn't help when resolved, remove it.

### Update Guidelines

- **Extend existing sections** by default; new sections fragment the narrative
- **Add new section** only when content is genuinely a new thread
- **Revise in place** when information changes; don't append "UPDATE: actually..."
- **Prune to appendix** when abandoning an approach; don't delete, move
- **Restructure** when the narrative no longer fits; flag it, don't do silently

### Anti-patterns

- Smoothing gaps with hedging language instead of marking
- Hiding uncertainty to make it "read better"
- Creating new section for every thought
- Orphaned content that doesn't connect to narrative

---

## Notes

### Thesis / Research Question

Reproduce the multipole expansion of the cow (Lehmann 2025, arXiv:2504.00506) from first principles: derive the mathematics, implement the numerics, and validate against the published results. `[HYPOTHESIS]`

### Motivation

The paper introduces a framework for systematically improving the spherical cow approximation (SCA) via spherical multipole expansion. Reproducing this work serves as:
1. A validation of the published results
2. An exercise in computational physics (mesh processing, spherical harmonics, rigid body mechanics)
3. A foundation for potential extensions

### Scope of Reproduction

The paper has three main computational components:

**A. Multipole expansion of bovine potentials (Sec III)**
- Compute spherical multipole moments $Q_\ell^m$ of a 3D cow mesh
- Compute the Cartesian quadrupole tensor $Q_{ij}^C$
- Derive gravitational wave emission power via the quadrupole formula
- Compute the inertia tensor and spindown timescale
- `[FUTURE: verify numerical values against Table I]`

**B. Multipole expansion of bovine surfaces (Sec IV)**
- Distance gradient flow method: inflate the cow to a star-shaped domain, project to sphere
- Harmonic map method: compute conformal map from cow surface to sphere via libigl
- Expand resulting maps $(f_r, f_{\Delta\theta}, f_{\Delta\phi})$ in spherical harmonics
- `[FUTURE: verify against Tables II and III]`

**C. Bovine rigid body mechanics (Sec V)**
- Cow tipping analysis: compute minimum tipping force given the geometry
- Requires inertia tensor + surface geometry from parts A and B
- `[FUTURE: reproduce force estimates]`

### Mathematical Conventions

`[BLOCKING: agent must converge on conventions before implementation]`

Coordinate system ("cow coordinates"):
- $x$: forward direction of the cow
- $y$: orthogonal to ground (up)
- $z$: cow's right

Benchmark units:
- Bounding box: $(1.044, 0.6397, 0.3403)$ in $(x, y, z)$
- Density: $\rho = 1$ (constant over the cow volume)

Spherical multipole moment:
$$Q_\ell^m \equiv \int_{\mathcal{C}} \mathrm{d}^3\mathbf{x}\, \rho(\mathbf{x})\, \|\mathbf{x}\|^\ell\, Y_\ell^m(\hat{\mathbf{x}})$$

Cartesian quadrupole:
$$Q_{ij}^C = \int \mathrm{d}^3\mathbf{x}\, \rho(\mathbf{x})\, (3 x_i x_j - r^2 \delta_{ij})$$

### Possible Approaches

- **Mesh source:** libigl tutorial data `cow.off` (same as paper)
- **Spherical harmonics:** Use `scipy.special.sph_harm` or `pyshtools`
- **Mesh processing:** `trimesh`, `libigl` Python bindings, or `igl`
- **Gradient flow:** Custom implementation using signed distance function
- **Harmonic map:** Via libigl's harmonic map routines

### Numerical Results

`[SOLID]` Cartesian quadrupole Q^C validated within 1.2% of paper values.
`[SOLID]` Inertia tensor I validated within 0.7% of paper values.
`[SOLID]` GW power coefficient ⟨Q⃛_ij Q⃛^ij⟩/ω⁶ ≈ 0.00144, paper gives 0.00149 (~3% error).
`[SOLID]` Physical GW power Ė ≈ 5.37×10⁻⁴¹ erg/s, paper gives 5.5×10⁻⁴¹ (~2%).
`[SOLID]` Spindown timescale ≈ 1.98×10⁴⁹ s, paper gives 1.9×10⁴⁹ (~4%).
`[SOLID]` All 17 validation tests pass.

`[PRELIMINARY]` Spherical Q_ℓ^m: surface-integral method has normalization offset vs paper convention. The `√(4π/(2ℓ+1))` factor in the potential formula suggests the paper absorbs this into Q_ℓ^m. Cartesian Q^C (convention-independent) validates correctly.

`[SOLID]` 3D sphere-to-mesh transition rendered for cow, bunny, and duck at ℓ_max = 0, 1, 2, 4, 8, 16, full. Progression clearly shows increasing detail from sphere to recognizable shape. See `results/sphere_to_cow_transition.png`, `results/sphere_to_bunny_transition.png`, and `results/sphere_to_duck_transition.png`. The rendering script (`code/render_sphere_to_cow.py`) is mesh-agnostic and accepts any .off file. The duck mesh (DASSL rubber duck, 4732 vertices, 9460 faces) was obtained from GitHub and converted from OBJ to OFF via trimesh.

`[FUTURE: surface map methods]` Gradient flow and harmonic map methods (Sec IV Tables II/III) are stubbed but not yet producing validated coefficients. These require more complex mesh processing (SDF computation, cotangent Laplacian, stereographic projection).

### April 1st Paper: Gravitational Wave Physics of the Duck

`[HYPOTHESIS]` Three research threads extending the multipole framework to a rubber duck mesh (DASSL duck, 4732 vertices, 9460 faces, converted from OBJ via trimesh). Full specification in `progress/duck_gw_specification.md`.

**Topic 1: Binary Duck on Post-Newtonian Orbits**
- Quadrupole-monopole interaction U_QM = −(GM_B/2r³) Q^C_ij n_i n_j modifies orbital dynamics at effective 2PN order (Poisson & Will Ch. 9)
- System quadrupole includes orbital + body permanent quadrupole → power decomposition P_orb + P_body + P_cross
- Phase evolution δΨ_Q enters at 2PN in SPA; measurable dephasing across LIGO (stellar) and LISA (supermassive) bands
- Duck tidal deformability Λ and shape-induced quadrupole κ_duck — compare with GW170817 constraints
- `[FUTURE: implement binary_duck.py, binary_gw.py, duck_waveform.py]`

**Topic 2: Quasinormal Modes of a Self-Gravitating Duck**
- Stellar mass (1–3 M_☉): full GR via AthenaK Z4c + GR hydro, based on existing `pgen/dyngr_tov.cpp`
- Supermassive (10⁶+ M_☉): Newtonian via AthenaK multigrid Poisson, based on existing polytropic star setup
- Duck shape lifts (2ℓ+1)-fold QNM degeneracy — gravitational Zeeman splitting
- Three-tier verification: AthenaK time-domain (Tier A), FEM Helmholtz eigenvalues (Tier B), perturbative splitting via Wigner 3j (Tier C)
- Anisotropic GW emission: broken symmetry gives m-dependent damping rates
- `[FUTURE: implement duck_star pgen in AthenaK, duck_helmholtz.py, duck_mode_splitting.py]`

**Topic 3: Equation of State of Duck Matter**
- Local compactness C(θ,φ) = GM/R(θ,φ)c² varies over duck surface; bill collapses first (highest C)
- TOV analysis with multiple EoS families (SLy, APR, MIT bag, polytrope, radiation)
- The "duck matter problem": maintaining duck shape at NS density requires anisotropic stress ~10⁵× nuclear σ_break (ε_duck ~ 0.3 vs ε_NS_max ~ 10⁻⁶)
- Supermassive duck: radiation-supported, directional Eddington limit (bill blows off first), GR instability at Γ_crit = 4/3 + O(C)
- `[FUTURE: implement duck_compactness.py, duck_tov.py, duck_anisotropic.py, duck_supermassive.py]`

**Cross-topic connections:** Duck relaxation to sphere = QNM excitation (Topic 2↔3); GW signal encodes both shape κ_duck and EoS Λ (Topics 1↔3); a sufficiently sensitive detector could distinguish a duck from a cow.

---

## Appendix

### Abandoned Approaches

[None yet.]
