# Duck-Shaped TOV Star: Instructions for Hengrui

## Overview

We need AthenaK simulations of a duck-shaped neutron star oscillating and
relaxing to a sphere. The purpose is to extract quasi-normal mode (QNM)
frequencies from a non-spherical compact object and compare them against
perturbative predictions (gravitational Zeeman splitting of QNM degeneracies).

The attached `duck_star_pgen.cpp` is a draft problem generator based on
`dyngr_tov.cpp`. The only physics modification is the **radial rescaling** that
maps the spherical TOV profile onto the duck shape. Everything else (metric
setup, magnetic fields, history output, EOS dispatch) is identical to the
standard TOV pgen.

## What the pgen does

1. Solves the standard TOV equations for a Gamma=2 polytrope (same as
   `dyngr_tov.cpp`).
2. At each grid point `(x, y, z)`, computes spherical coordinates
   `(r, theta, phi)`.
3. Evaluates the duck surface function using real spherical harmonics:

       R_duck(theta, phi) = R_0 * [1 + amplitude * sum_{l,m} eps_{lm} Y_lm(theta, phi)]

4. Computes the effective radius `r_eff = r / R_ratio` where
   `R_ratio = R_duck / R_0`.
5. Queries the TOV solution at `r_eff` instead of `r` to get density and
   pressure. This stretches/compresses the radial profile directionally,
   creating a duck-shaped density distribution.
6. Sets the **metric** from the spherical TOV at the physical radius `r` (not
   `r_eff`). This introduces an O(epsilon^2) Hamiltonian constraint violation
   that Z4c damps.

## Input files needed

### duck_epsilon.dat

Format:
```
# Duck SH deformation coefficients epsilon_lm
# R(theta,phi) = R_0 [1 + sum epsilon_lm Y_lm(theta,phi)]
# R_0 = 0.2937620957
# R_eq = 0.4359705304
# ell  m  epsilon_real  epsilon_imag
   1   -1  -6.273383998814e-04  -6.243833123852e-02
   1    0  +2.027182372167e+00  +0.000000000000e+00
   ...
```

This file is already provided in the `production-code/` directory. The `R_0`
value is read from the header comment. The coefficients use the real spherical
harmonic convention (the imaginary column is present for completeness but is
not used by the pgen -- we convert to real SH basis internally).

**Important:** The epsilon values in this file are the raw fit to the duck mesh
at unit scale. For a physically reasonable perturbation, use the
`duck_amplitude` parameter to scale them down. For example:
- `duck_amplitude = 1.0`: full duck shape (large deformation, nonlinear)
- `duck_amplitude = 0.1`: 10% duck (linear regime, good for perturbation
  theory comparison)
- `duck_amplitude = 0.01`: 1% duck (firmly linear, for convergence tests)

## Suggested athinput parameters

### Full GR run (primary)

```
<comment>
problem  = Duck-shaped TOV star

<job>
basename = duck_tov

<mesh>
nghost = 4
nx1    = 128
x1min  = -51.2
x1max  = 51.2
ix1_bc = outflow
ox1_bc = outflow

nx2    = 128
x2min  = -51.2
x2max  = 51.2
ix2_bc = outflow
ox2_bc = outflow

nx3    = 128
x3min  = -51.2
x3max  = 51.2
ix3_bc = outflow
ox3_bc = outflow

<mesh_refinement>
refinement = adaptive
num_levels = 4
refine_interval = 1

<meshblock>
nx1  = 16
nx2  = 16
nx3  = 16

<time>
evolution  = dynamic
integrator = rk3
cfl_number = 0.4
nlim       = -1
tlim       = 10000
ndiag      = 1

<coord>
general_rel = true
m           = 0.0
a           = 0.0
excise      = false

<z4c>
diss        = 0.1
chi_div_floor = 0.1

<mhd>
eos         = ideal
dyn_eos     = ideal
dyn_error   = reset_floor
reconstruct = ppmx
rsolver     = hlle
dfloor      = 1.0e-10
tfloor      = 1.0e-8
dthreshold  = 1.02
gamma       = 2.0
dyn_scratch = 1
fofc        = true

<adm>

<problem>
rhoc              = 1.28e-3
kappa             = 100.0
npoints           = 10000
dr                = 1e-3
b_norm            = 0.0
pcut              = 1e-6
magindex          = 1
user_hist         = true
v_pert            = -0.01
duck_lmax         = 10
duck_epsilon_file = duck_epsilon.dat
duck_amplitude    = 0.1

<output1>
file_type   = rst
dt          = 1000.0

<output2>
file_type   = hst
dt          = 0.1
data_format = %20.15e

<output3>
file_type   = bin
variable    = mhd_w
dt          = 10.0

<output4>
file_type   = bin
variable    = z4c
dt          = 100.0
```

### Key parameter choices

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `duck_amplitude` | 0.1 | Linear regime for comparison with perturbation theory. Try 0.01 for convergence, 1.0 for full nonlinear duck. |
| `v_pert` | -0.01 | Gentle radial kick to excite f-modes cleanly. |
| `tlim` | 10000 M | ~50 ms, ~100 dynamical times. Enough to resolve QNM frequencies. |
| `hst dt` | 0.1 | High-cadence history for FFT (Nyquist freq ~5 kHz). Reduce to 0.01 if disk allows. |
| `bin dt` | 10.0 | 3D snapshots every ~0.05 ms for time-lapse visualization. |
| AMR levels | 4 | Effective resolution ~0.2 M across the star. |

### Newtonian comparison run (optional, lower priority)

If a Newtonian self-gravity solver is available (e.g., `polytropic_star` pgen),
the same radial rescaling approach applies. Use:
```
<problem>
rho_c              = 1.0
R_star             = 3.654
v_pert             = -0.01
rho_floor          = 1.0e-8
njeans             = 16
duck_lmax          = 10
duck_epsilon_file  = duck_epsilon.dat
duck_amplitude     = 0.1
```
The Newtonian run serves as a cross-check (no GW damping, only numerical
viscosity). QNM frequencies should match the Cowling approximation.

## What output we need

### Essential
1. **History file** (`.hst`): columns for `rho-max` and `alpha-min` at high
   cadence. This is the primary observable for FFT-based mode extraction.
2. **3D binary snapshots** (`mhd_w`): density, pressure, velocity fields at
   regular intervals for visualization (the "duck melts to sphere" figure).

### Highly desired
3. **Psi4 extraction**: If the Weyl scalar extraction module is available,
   extract Psi4 on coordinate spheres at r = 40, 60, 80, 100 M. Decompose
   into spin-weighted spherical harmonics _{-2}Y_{lm} for (l,m) up to l=4.
   This gives us the GW signal directly.
4. **Constraint norms**: L2 norm of the Hamiltonian constraint as a function
   of time, to verify that Z4c damps the initial violation.

### Nice to have
5. **z4c snapshots** at coarser intervals (dt=100) for checking gauge dynamics.
6. **Restart files** at regular intervals in case we need to extend runs.

## Constraint handling recommendations

The duck-shaped initial data is **not** a self-consistent solution to the
Einstein constraint equations. We are perturbing only the matter sector
(density, pressure) while keeping the metric spherical. This is a standard
approach in NR for studying perturbations of equilibrium configurations.

**Expected behavior:**
- Initial Hamiltonian constraint violation: ||H||_2 ~ O(epsilon^2). For
  `duck_amplitude = 0.1`, this is ~1% level.
- Z4c constraint damping should reduce ||H|| by ~1 order of magnitude within
  t ~ 50 M.
- If the constraint violation grows rather than decays, this indicates a
  resolution or gauge problem.

**If constraints blow up:**
1. Increase resolution (add AMR levels or increase base grid).
2. Reduce `duck_amplitude` to stay in the linear regime.
3. Try `diss = 0.2` in the `<z4c>` block for stronger Kreiss-Oliger dissipation.

## Resolution and AMR suggestions

### Minimum viable run
- Base: 64^3, 3 AMR levels, MeshBlock 16^3
- Effective resolution: ~0.4 M (R_*/20)
- Cost: ~few hundred GPU-hours on A100

### Production run
- Base: 128^3, 4 AMR levels, MeshBlock 16^3
- Effective resolution: ~0.2 M (R_*/40)
- Cost: ~few thousand GPU-hours on A100

### Convergence test
Run at 3 resolutions (0.4 M, 0.2 M, 0.1 M) with `duck_amplitude = 0.01` to
verify that QNM frequencies converge. Expected: 2nd-order convergence with PPMx.

### AMR refinement criteria
Refine on density gradient. The duck surface is a strong density discontinuity
that needs to be well-resolved. The AMR should automatically place the finest
levels around the stellar surface.

## Code integration notes

1. The pgen file `duck_star_pgen.cpp` should be placed in
   `src/pgen/duck_star.cpp` and compiled with `-D PROBLEM=duck_star`.
2. The `DuckDeformation` struct uses fixed-size arrays (121 entries for
   lmax=10) to avoid dynamic memory allocation inside Kokkos lambdas.
3. The real spherical harmonic evaluation is done inline via recurrence
   relations -- no external library needed.
4. The `duck_epsilon.dat` file must be in the run directory (or provide
   the full path in the athinput file).

## Contact

For questions about the physics setup, deformation coefficients, or
post-processing pipeline, contact Haiyang. The Python post-processing scripts
(FFT, Psi4 decomposition, waveform generation) will be provided separately
once the first simulation data is available.
