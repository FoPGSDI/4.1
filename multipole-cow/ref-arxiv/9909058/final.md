# Kokkotas & Schmidt (1999) -- Key Results for Quacky Normal Modes

**Paper:** "Quasi-Normal Modes of Stars and Black Holes," Living Rev. Rel. 2, 2 (1999)
**arXiv:** gr-qc/9909058
**Markdown output:** `md_output/`

## Relevance to Our Work

This is the canonical review of QNMs for compact objects. It provides the theoretical framework we adapt for the duck-shaped star: mode classification, perturbation equations, and frequency/damping scalings.

## Key Results Extracted

### 1. Mode Classification (Section 5)

The paper establishes the mode taxonomy we use:

- **f-mode (fundamental):** No nodes inside star. Frequency proportional to mean density. For a uniform-density Newtonian star:
  $$\omega^2 = \frac{2\ell(\ell-1)}{2\ell+1} \frac{M}{R^3}$$
  Typical NS values: f ~ 2.87 kHz, damping time ~ 0.11 s.

- **p-modes (pressure/acoustic):** Infinitely many, pressure is restoring force. Frequencies > 4--7 kHz for p_1. Damping times ~ seconds.

- **g-modes (gravity):** Buoyancy-driven, require non-zero Schwarzschild discriminant. Frequencies < 100 Hz, damping times ~ years.

- **w-modes (spacetime):** Unique to GR, no significant fluid motion. High frequency (5--12 kHz), very rapid damping (~ 0.01--0.1 ms). Three subfamilies: curvature, trapped, and interface modes.

- **r-modes (rotational):** Appear in rotating stars, frequency ~ 2m*Omega / [l(l+1)]. Generically CFS-unstable for l >= 2.

### 2. Perturbation Equations (Section 5.1)

Polar perturbations of a relativistic star decompose into three coupled wave equations for spacetime variables (S, F) and density perturbation (H):
- S and F propagate at light speed
- H propagates at sound speed c_s
- Outside the star: reduces to Zerilli equation

Key approximations:
- **Cowling approximation:** Freeze spacetime perturbations, keep only fluid equation for H
- **Inverse Cowling approximation (ICA):** Freeze fluid, keep spacetime equations for S, F

### 3. Table of QNM Properties (Table 2)

| Mode  | Frequency   | Damping time |
|-------|-------------|--------------|
| f     | 2.87 kHz    | 0.11 sec     |
| p_1   | 6.57 kHz    | 0.61 sec     |
| g_1   | 19.85 Hz    | years        |
| w_1   | 12.84 kHz   | 0.024 ms     |
| w_II  | 8.79 kHz    | 0.016 ms     |

### 4. Stability Criterion

Non-radial stability governed by Schwarzschild discriminant:
$$S(r) = \frac{dp}{dr} - \frac{\Gamma_1 p}{\rho + p}\frac{d\rho}{dr}$$
S > 0 everywhere => stable QNMs. Related Brunt-Vaisala frequency N^2 = g*S(r).

### 5. Axial Perturbation Equation

Single wave equation for toroidal modes:
$$-\frac{1}{c^2}\frac{\partial^2 X}{\partial t^2} + \frac{\partial^2 X}{\partial r_*^2} + \frac{e^\nu}{r^3}\left[\ell(\ell+1)r + r^3(\rho - p) - 6M\right] = 0$$
Reduces to Regge-Wheeler equation outside the star.

## How We Use This

- The f-mode frequency formula omega^2 ~ M/R^3 is the basis for our "quack frequency" scaling
- Mode classification (f, p, g, w) provides the framework for identifying duck QNMs
- Cowling approximation justifies our Newtonian treatment of supermassive duck oscillations
- The w-mode spectrum provides the GR correction hierarchy
