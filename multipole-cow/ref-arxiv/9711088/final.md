# Andersson & Kokkotas (1998) -- Key Results for Quacky Normal Modes

**Paper:** "Towards gravitational-wave asteroseismology," MNRAS 299, 1059 (1998)
**arXiv:** gr-qc/9711088
**Markdown output:** `md_output/`

## Relevance to Our Work

This paper establishes the "universal relations" between QNM frequencies and stellar parameters (M, R). These empirical fitting formulas are what we adapt for duck-shaped stars to show how shape deformation modifies the coefficients.

## Key Results Extracted

### 1. f-mode Frequency -- Universal Relation (Eq. 1)

$$\omega_f \text{ (kHz)} \approx 0.78 + 1.635 \left(\frac{\bar{M}}{\bar{R}^3}\right)^{1/2}$$

where $\bar{M} = M / (1.4 M_\odot)$ and $\bar{R} = R / (10\text{ km})$.

This confirms omega_f ~ sqrt(mean density), with a linear fit that holds across 12 realistic EOS to within a few percent.

### 2. f-mode Damping Time (Eq. 2)

$$\frac{1}{\tau_f \text{ (s)}} \approx \frac{\bar{M}^3}{\bar{R}^4}\left[22.85 - 14.65\left(\frac{\bar{M}}{\bar{R}}\right)\right]$$

Derived from quadrupole formula scaling: tau_f ~ R(R/M)^3.

### 3. p-mode Frequency (Eq. 3)

$$\omega_p \text{ (kHz)} \approx \frac{1}{\bar{M}}\left(1.75 + 5.59\frac{\bar{M}}{\bar{R}}\right)$$

Less robust than f-mode relation; p-modes are sensitive to interior matter distribution.

### 4. w-mode Frequency and Damping (Eqs. 4--5)

$$\omega_w \text{ (kHz)} \approx \frac{1}{\bar{R}}\left[20.92 - 9.14\left(\frac{\bar{M}}{\bar{R}}\right)\right]$$

$$\frac{1}{\tau_w \text{ (ms)}} \approx \frac{1}{\bar{M}}\left[5.74 + 103\left(\frac{\bar{M}}{\bar{R}}\right) - 67.45\left(\frac{\bar{M}}{\bar{R}}\right)^2\right]$$

w-mode frequency inversely proportional to stellar size; damping depends on compactness.

### 5. Parameter Estimation Accuracy (Table 2)

Combining f-mode and w-mode frequencies: errors (5%, 2%) in (R, M).
Combining f-mode frequency and damping: errors (6.5%, 17.6%).
Best combination (w-mode freq + f-mode damping): errors (3.2%, 1.9%).

### 6. Black Hole vs Neutron Star Discrimination

BH fundamental QNM: f ~ 12 kHz * (M_sun/M), tau ~ 0.05 ms * (M/M_sun).
NS f-mode at similar frequency but damping ~1000x longer => easy to distinguish.

## How We Use This

- The universal relations (Eqs. 1--5) are the baseline: for a spherical star, these predict QNM frequencies from (M, R)
- For the duck, the shape breaks spherical symmetry => mode splitting. The deviation from these universal relations quantifies the "duckiness"
- The fitting formulas give us the zero-order frequencies that get split by the duck multipole moments
- The parameter estimation framework shows what precision is needed to detect duck vs sphere differences
