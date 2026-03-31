# Hinderer (2008) -- Key Results for Quacky Normal Modes

**Paper:** "Tidal Love numbers of neutron stars," ApJ 677, 1216 (2008)
**arXiv:** 0711.2420
**Markdown output:** `md_output/`

## Relevance to Our Work

This paper provides the complete relativistic calculation of the tidal Love number k_2 for neutron stars. The Love number quantifies the quadrupolar deformability of a star in response to an external tidal field -- directly relevant to how a duck-shaped star responds to tidal forces from a companion.

## Key Results Extracted

### 1. Definition of Love Number (Eq. 3--4)

The tidal deformability parameter lambda relates the induced quadrupole moment Q_ij to the external tidal field E_ij:
$$Q_{ij} = -\lambda \mathcal{E}_{ij}$$

The dimensionless Love number k_2 is:
$$k_2 = \frac{3}{2} G \lambda R^{-5}$$

### 2. The h_2 ODE (Eq. 14)

The key differential equation for the metric perturbation H(r) inside the star:
$$H'' + H'\left[\frac{2}{r} + e^\lambda\left(\frac{2m(r)}{r^2} + 4\pi r(p - \rho)\right)\right] + H\left[-\frac{6e^\lambda}{r^2} + 4\pi e^\lambda\left(5\rho + 9p + \frac{\rho+p}{dp/d\rho}\right) - \nu'^2\right] = 0$$

Boundary condition at r=0: H(r) = a_0 r^2 [1 + O(r^2)].

### 3. The k_2 Formula (Eq. 20)

The master formula for k_2 in terms of compactness C = M/R and y = R*H'(R)/H(R):

$$k_2 = \frac{8C^5}{5}(1-2C)^2[2+2C(y-1)-y] \times \left\{2C[6-3y+3C(5y-8)] + 4C^3[13-11y+C(3y-2)+2C^2(1+y)] + 3(1-2C)^2[2-y+2C(y-1)]\log(1-2C)\right\}^{-1}$$

### 4. Newtonian Limit (Eq. 22)

$$k_2^N = \frac{1}{2}\left(\frac{2-y}{y+3}\right)$$

For n=1 polytrope: k_2^N(n=1) = -1/2 + 15/(2*pi^2) ~ 0.260.

### 5. Fitting Formula (Eq. 25)

For 0.5 <= n <= 1.0 and 0.1 <= M/R <= 0.24:
$$k_2 \approx \frac{3}{2}\left(-0.41 + \frac{0.56}{n^{0.33}}\right)\left(\frac{M}{R}\right)^{-0.003}$$

k_2 depends more strongly on polytropic index n (central condensation) than on compactness M/R.

### 6. Table of Love Numbers

| n   | M/R      | k_2    |
|-----|----------|--------|
| 0.5 | 1e-5     | 0.4491 |
| 0.5 | 0.1      | 0.251  |
| 0.5 | 0.2      | 0.095  |
| 1.0 | 1e-5     | 0.2599 |
| 1.0 | 0.1      | 0.122  |
| 1.0 | 0.2      | 0.0459 |

Relativistic values are lower than Newtonian, especially for higher n (more centrally condensed stars).

### 7. LIGO Measurability

LIGO II can constrain lambda <= 2.01e37 g cm^2 s^2 (90% confidence) for 1.4 M_sun NS-NS inspiral at 50 Mpc.

## How We Use This

- The h_2 ODE and k_2 formula provide the framework for computing Love numbers of duck-shaped stars
- For a non-spherical star, the Love number becomes a tensor (not scalar) -- the duck shape introduces anisotropic tidal response
- The Newtonian limit k_2^N gives us the baseline for Newtonian duck calculations
- The compactness dependence quantifies relativistic corrections to the duck's tidal deformability
- The fitting formula provides quick estimates for comparing duck vs sphere Love numbers
