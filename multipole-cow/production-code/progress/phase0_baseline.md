# Phase 0: Duck Baseline — Quacky Normal Modes

**Date:** 2026-03-30

## Mesh Properties

- **File:** `data/duck.off`
- **Vertices:** 4732
- **Faces:** 9460
- **Volume:** V = 3.4710431392e-01
- **Center of mass:** (0.000318, 0.010456, -0.029563)
- **Bounding box size:** (0.892083, 1.000000, 0.963750)
- **Equivalent radius:** R_eq = (3V/4pi)^(1/3) = 0.4359705304
- **Mean radius (SH fit):** R_0 = 0.2937620957

## Cartesian Quadrupole Tensor Q^C

```
  [-1.1770676965e-02  +1.0285007841e-05  -9.7492724397e-05]
  [+1.0285007841e-05  +4.3183292223e-03  -9.0605795369e-03]
  [-9.7492724397e-05  -9.0605795369e-03  +7.4523477430e-03]
```

Eigenvalues: +1.5080670775e-02, -3.3093917489e-03, -1.1771279026e-02

## Inertia Tensor I

```
  [+3.5889958210e-02  -3.4283359471e-06  +3.2497574799e-05]
  [-3.4283359471e-06  +3.0526956147e-02  +3.0201931790e-03]
  [+3.2497574799e-05  +3.0201931790e-03  +2.9482283307e-02]
```

Principal moments (eigenvalues): 3.5890158897e-02, 3.3069529804e-02, 2.6939508963e-02

## Dimensionless Shape Parameters

- kappa_duck = Q^C_eigenvalues / (M * R_eq^2) = (+0.228584, -0.050162, -0.178422)

## Spherical Multipole Moments Q_lm

- Q_00 = 9.7916319160e-02+0.0000000000e+00j (monopole, = V * Y_00 normalization)
- Q_10 = 5.20e-18+0.00e+00j (should vanish at COM)

Full Q_lm up to l=10 stored in `duck_baseline.json`.

## SH Deformation Spectrum

Radial deformation: R(theta,phi) = R_0 [1 + sum epsilon_lm Y_lm]

| l | P_l = sum_m |eps_lm|^2 | sqrt(P_l) |
|---|---|---|
| 1 | 4.11726625e+00 | 2.029105 |
| 2 | 8.11923385e+00 | 2.849427 |
| 3 | 1.25740847e+01 | 3.545996 |
| 4 | 5.93612371e+00 | 2.436416 |
| 5 | 5.05430546e+00 | 2.248178 |
| 6 | 2.53066032e+00 | 1.590805 |
| 7 | 1.58139814e+00 | 1.257537 |
| 8 | 4.37583750e-01 | 0.661501 |
| 9 | 1.80995857e-01 | 0.425436 |
| 10 | 5.76905974e-02 | 0.240189 |

## Output Files

- `production-code/duck_baseline.json` -- all numerical results
- `production-code/duck_epsilon.dat` -- epsilon_lm in tabular format
- `results/duck_deformation_spectrum.png` -- bar chart of |epsilon_lm|
- `results/duck_power_spectrum.png` -- P_l vs l

## Status

Phase 0 baseline computation complete. All moments, eigenvalues, and
deformation coefficients computed and saved. Ready for Phase 1 (normal mode analysis).
