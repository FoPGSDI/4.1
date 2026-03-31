# Section 4: Tidal Love Numbers of the Duck

The tidal deformability of a compact object encodes its internal
structure in a single observable quantity accessible through
gravitational-wave measurements.  For a spherical neutron star the
response to an external tidal field is characterized by a single number,
the Love number $k_2$.  For the duck the situation is richer: the broken
spherical symmetry promotes the scalar $k_2$ to a rank-4 *Love tensor*
$\lambda_{ijkl}$, introduces a direction-dependent tidal apsidal
constant $k_2(\theta,\varphi)$, and breaks the celebrated I-Love-Q
universality of Yagi & Yunes (2013).

We develop the theory in stages: the effective (spherical) Love number
(Sec.\ 4.1), the full Love tensor (Sec.\ 4.2), its angular
decomposition (Sec.\ 4.3), perturbative corrections to the
dimensionless tidal deformability $\Lambda$ (Sec.\ 4.4), and the
I-Love-Q breakdown (Sec.\ 4.5).  Throughout we compare the nuclear duck
(a compact object with duck shape, $\Gamma=2$ polytropic EoS) and the
rubber duck (a terrestrial bath toy).


## 4.1 Effective (Spherical) Love Number

### The $h_2$ master equation

We first treat the duck as its equivalent sphere of radius
$R_\mathrm{eq} = (3V/4\pi)^{1/3}$ and compute the standard $\ell = 2$
tidal Love number following Hinderer (2008) and Flanagan & Hinderer
(2008).

In a static, spherically symmetric background with line element
$ds^2 = -e^{\nu(R)}\,dt^2 + e^{\lambda(R)}\,dR^2 + R^2\,d\Omega^2$,
the even-parity $\ell=2$ metric perturbation function $h_2(R)$ satisfies
the master equation (Hinderer 2008; Damour & Nagar 2009; Yagi & Yunes
2013b):

$$
0 = h_2'' + \left\{\frac{2}{R} + \left[\frac{2M}{R^2}
  + 4\pi R\,(p - \rho)\right]e^{\lambda}\right\}h_2'
  - \left\{\frac{6\,e^{\lambda}}{R^2}
  - 4\pi\left[5\rho + 9p + (\rho+p)\frac{d\rho}{dp}\right]e^{\lambda}
  + \left(\frac{d\nu}{dR}\right)^{\!2}\right\}h_2\,,
$$

where $M(R)$ is the enclosed gravitational mass, $p(R)$ and $\rho(R)$
are the pressure and energy density from the TOV solution, and
$e^{\lambda} = (1 - 2M/R)^{-1}$.  Primes denote $d/dR$.

### Boundary conditions and numerical procedure

The integration proceeds as follows:

1. **TOV integration.**  For a $\Gamma = 2$ polytrope with
   $p = K\rho^{\Gamma}$ ($K = 100$ in geometric units,
   $\rho_c = 1.28\times 10^{-3}$), we integrate the coupled TOV
   equations
   $$
   \frac{dp}{dR} = -\frac{(\rho + p)(M + 4\pi R^3 p)}{R(R - 2M)}\,,
   \qquad
   \frac{dM}{dR} = 4\pi R^2 \rho\,,
   \qquad
   \frac{d\nu}{dR} = \frac{2(M + 4\pi R^3 p)}{R(R - 2M)}
   $$
   from $R = 0$ to the stellar surface $R_*$ defined by $p(R_*) = 0$.

2. **Coupled $h_2$ integration.**  Simultaneously with the TOV system,
   we integrate the $h_2$ ODE.  Near the origin the regular solution
   behaves as $h_2 \propto R^2$; we initialize
   $h_2(R_0) = R_0^2$, $h_2'(R_0) = 2R_0$ at a small starting radius
   $R_0 \ll R_*$.  Because the $h_2$ equation is linear and homogeneous,
   the overall normalization drops out of the Love number.

3. **Surface matching.**  At $R = R_*$, we compute the logarithmic
   derivative
   $$
   y \equiv \frac{R_*\,h_2'(R_*)}{h_2(R_*)}\,.
   $$

### The Hinderer formula

The tidal apsidal constant $k_2$ is then given by (Hinderer 2008,
2010):

$$
k_2 = \frac{8}{5}\,C^5\,(1 - 2C)^2\bigl[2 + 2C(y-1) - y\bigr]
\;\Bigg/\;
\Bigl\{
  2C\bigl[6 - 3y + 3C(5y - 8)\bigr]
  + 4C^3\bigl[13 - 11y + C(3y - 2) + 2C^2(1+y)\bigr]
$$
$$
  + 3(1 - 2C)^2\bigl[2 - y + 2C(y-1)\bigr]\ln(1 - 2C)
\Bigr\}\,,
$$

where $C \equiv M_*/R_*$ is the stellar compactness.  The dimensionless
tidal deformability is

$$
\Lambda = \frac{2}{3}\,k_2\,C^{-5}\,.
$$

This is the quantity directly constrained by gravitational-wave
observations: $\Lambda = 190^{+390}_{-120}$ from GW170817 (Abbott et
al.\ 2018).

### Newtonian Love numbers for polytropes

In the Newtonian limit ($C \to 0$), the tidal apsidal constant is
obtained from the Clairaut-Radau equation.  The classical results of
Brooker & Olle (1955) are reproduced in Table 4a.

**Table 4a.** Newtonian $k_2$ for polytropic index $n$ (Brooker & Olle
1955).  The polytropic exponent is $\Gamma = 1 + 1/n$; $n = 0$
corresponds to a uniform-density sphere.

| $n$ | $\Gamma$ | $k_2^{\rm N}$ |
|:---:|:--------:|:--------------:|
| 0   | $\infty$ | 0.750 |
| 0.5 | 3        | 0.449 |
| 1.0 | 2        | 0.260 |
| 1.5 | 5/3      | 0.143 |
| 2.0 | 3/2      | 0.0728|
| 3.0 | 4/3      | 0.0116|

The $n = 0$ (incompressible) result $k_2 = 3/4$ follows immediately from
the exact Radau solution $\eta_2 = 5$ for uniform density.  As the
polytropic index increases the star becomes more centrally concentrated
and correspondingly harder to deform, so $k_2$ decreases monotonically.

### Results: nuclear duck vs rubber duck

**Compact nuclear duck** ($\Gamma = 2$, $K = 100$,
$\rho_c = 1.28\times 10^{-3}$).  The TOV integration yields
$M_* \simeq 1.40\,M_\odot$, $R_* \simeq 9.6$ (in geometric units),
$C \simeq 0.146$.  The coupled $h_2$ integration gives $y \simeq 0.87$
and

$$
k_2 \simeq 0.074\,,\qquad
\Lambda \simeq 740\,,
$$

the range in $\Lambda$ reflecting the sensitivity to central density
(and hence compactness) across the TOV family.  These values are
comfortably within the GW170817 constraint band.

**Rubber duck** (terrestrial bath toy).  For an elastic solid the Love
number is (Love 1909; Greff-Lefftz et al.\ 2005):

$$
k_2^{\rm elastic} = \frac{3/2}{1 + \dfrac{19\,\mu_{\rm shear}}{2\,\rho\,g\,R}}\,.
$$

Taking the shear modulus of PVC (the actual material of bath ducks)
$\mu_{\rm shear} \simeq 3\times 10^9\;\mathrm{Pa}$,
mean density $\rho \simeq 1.3\times 10^3\;\mathrm{kg\,m^{-3}}$, surface
gravity $g \simeq 9.8\;\mathrm{m\,s^{-2}}$, and radius
$R \simeq 0.04\;\mathrm{m}$:

$$
k_2^{\rm rubber} \simeq \frac{1.5}{1 + \frac{19\times 3\times 10^9}{2\times 1300\times 9.8\times 0.04}}
\simeq 2.7\times 10^{-8}\,.
$$

The denominator is dominated by the rigidity-to-gravity ratio
$19\mu/(2\rho g R) \simeq 5.6\times 10^7$.  Punchline: *a bathtub duck
is $10^9$ times less deformable than a neutron star duck.*


## 4.2 The Love Tensor

### Tensorial tidal response

For a spherically symmetric body the induced quadrupole moment is
proportional to the applied tidal field:

$$
Q_{ij}^{\rm induced} = -\lambda\,\mathcal{E}_{ij}\,,
$$

where $\lambda = (2/3)\,k_2\,R_*^5$ is the (scalar) tidal
deformability and $\mathcal{E}_{ij}$ is the electric-type tidal tensor.
When spherical symmetry is broken --- as it is for the duck --- the
response becomes tensorial:

$$
Q_{ij}^{\rm induced} = -\lambda_{ijkl}\,\mathcal{E}_{kl}\,.
$$

The Love tensor $\lambda_{ijkl}$ maps traceless symmetric rank-2 tensors
to traceless symmetric rank-2 tensors.  Both $Q_{ij}$ and
$\mathcal{E}_{kl}$ are symmetric and traceless (5 independent
components each), so $\lambda_{ijkl}$ is a linear map on a
5-dimensional vector space --- a $5\times 5$ matrix.  With the
additional symmetry $\lambda_{ijkl} = \lambda_{klij}$ (the response
is derivable from a potential; Landry & Poisson 2015), the number of
independent components is

$$
N_{\rm ind} = \frac{5\times 6}{2} = 15\,.
$$

For a body with discrete symmetries the count is reduced: a body with a
single reflection symmetry has 9 independent components; a body with
two orthogonal reflection planes has 6; full axisymmetry leaves 3.
The duck, with its approximate bilateral symmetry, has
$\sim 9$ independent components, which we parametrize perturbatively.

### Perturbative expansion

We expand the Love tensor in the deformation parameters $\varepsilon_{\ell m}$:

$$
\lambda_{ijkl} = \lambda_0\,\mathcal{P}_{ijkl}^{(0)}
  + \sum_{\ell',m'} \varepsilon_{\ell' m'}\,
    \delta\lambda_{ijkl}^{(\ell' m')}
  + \mathcal{O}(\varepsilon^2)\,,
$$

where $\mathcal{P}_{ijkl}^{(0)} = \frac{1}{2}(\delta_{ik}\delta_{jl}
+ \delta_{il}\delta_{jk}) - \frac{1}{3}\delta_{ij}\delta_{kl}$ is the
identity projector on traceless symmetric tensors, and
$\lambda_0 = (2/3)\,k_2^{(0)}\,R_{\rm eq}^5$ is the spherical Love
number.

### Selection rules for angular coupling

The first-order correction $\delta\lambda^{(\ell'm')}$ couples different
angular momentum channels.  The coupling is governed by Gaunt integrals
$\int Y_2^{m_1*}\,Y_2^{m_2}\,Y_{\ell'}^{m'}\,d\Omega$, which enforce:

1. **Triangle inequality:** $|\ell_1 - \ell_2| \leq \ell' \leq \ell_1 + \ell_2$.
   Since the tidal field and response are both $\ell = 2$, this gives
   $0 \leq \ell' \leq 4$.
2. **Parity:** $\ell'$ must be even (the integral vanishes for odd
   $\ell'$).
3. **$m$-selection:** $m' = m_1 - m_2$.

Therefore, only $\varepsilon_{0,0}$, $\varepsilon_{2,m}$, and
$\varepsilon_{4,m}$ enter at first order.  The monopole
$\varepsilon_{0,0}$ produces an overall rescaling; the quadrupole and
hexadecapole deformations generate off-diagonal components in the Love
tensor.


## 4.3 Direction-Dependent $k_2(\theta,\varphi)$

### Definition

The direction-dependent tidal apsidal constant measures the response
to a tidal field applied along a unit vector $\hat{n} = (\theta,\varphi)$
on the sky:

$$
k_2(\theta,\varphi) = \frac{3}{2R_{\rm eq}^5}\,
  \lambda_{ijkl}\,n_i\,n_j\,n_k\,n_l\,.
$$

This is the natural generalization of the scalar $k_2$: for a
spherical body, $\lambda_{ijkl} = \lambda_0\,\mathcal{P}_{ijkl}^{(0)}$
and $k_2(\hat{n})$ is a constant independent of direction.

### Spherical harmonic decomposition

The product $n_i n_j n_k n_l$, when restricted to the traceless-symmetric
subspace, decomposes into $\ell = 0, 2, 4$ spherical harmonics
(the symmetric trace-free decomposition of a rank-4 tensor).
Correspondingly:

$$
k_2(\theta,\varphi) = k_2^{(0)}
  + \sum_m k_2^{(2,m)}\,Y_2^m(\theta,\varphi)
  + \sum_m k_2^{(4,m)}\,Y_4^m(\theta,\varphi)\,,
$$

where $k_2^{(0)}$ is the angle-averaged (effective spherical) apsidal
constant, and the $\ell = 2$ and $\ell = 4$ harmonics encode the
anisotropy.

### Physical interpretation for the duck

The direction-dependence of $k_2$ reflects the spatial variation of
local compactness and stiffness across the duck surface:

- **Bill region** ($\theta \approx 0$).  The narrow, elongated bill has
  a high surface-to-volume ratio and correspondingly higher local
  compactness.  This *suppresses* $k_2$ in the bill direction --- the
  bill is stiffer.

- **Body region** ($\theta \approx \pi/2$).  The broad, rounded body
  has low local compactness and is more easily deformed by a tidal
  field: $k_2$ is *enhanced* in the equatorial directions.

- **Tail region** ($\theta \approx \pi$).  Intermediate behavior,
  modulated by the $\ell = 4$ (hexadecapole) deformation.

The $\ell = 2$ variation in $k_2$ is proportional to the quadrupole
deformation $\varepsilon_{2,m}$, while the $\ell = 4$ variation is
proportional to $\varepsilon_{4,m}$.  For the DASSL duck,
$|\varepsilon_{2,0}| \approx 0.15$ and $|\varepsilon_{4,0}| \approx 0.05$,
so the quadrupole anisotropy dominates.


## 4.4 Perturbative Correction to $\Lambda$

### Leading-order shape correction

The dimensionless tidal deformability $\Lambda$ receives its leading
shape correction at *second order* in the deformation parameters
$\varepsilon_{\ell m}$.  This is because $\Lambda$ is a scalar (the
angle-averaged response), and the first-order correction involves odd
powers of $Y_{\ell'}^{m'}$ which average to zero.  We write:

$$
\Lambda_{\rm duck} = \Lambda_{\rm sphere}\left[
  1 + \sum_{\ell \geq 1} \alpha_\ell \sum_m |\varepsilon_{\ell m}|^2
  + \mathcal{O}(\varepsilon^3)
\right]\,.
$$

### Two sources of correction

The coefficients $\alpha_\ell$ receive contributions from two
physically distinct mechanisms:

1. **Geometric correction** (from $R^{-5}$).  The tidal deformability
   scales as $\lambda \propto R_*^5$, and the deformed surface
   $R(\theta,\varphi) = R_{\rm eq}[1 + \sum \varepsilon_{\ell m}
   Y_\ell^m]$ modifies the effective fifth power of the radius.
   Expanding:
   $$
   \langle R^5 \rangle / R_{\rm eq}^5 = 1 + 10\sum_m
   |\varepsilon_{0,m}|^2 + 15\sum_{\ell \geq 1,m}
   |\varepsilon_{\ell m}|^2 + \cdots
   $$
   The geometric contribution to $\alpha_\ell$ is therefore
   $\alpha_\ell^{\rm geom} \approx +15$ for all $\ell \geq 1$.

2. **Physical correction** (from $\delta k_2$).  The deformation
   modifies the interior structure (density profile, metric functions)
   and hence the quantity $y = R_* h_2'/h_2$.  Via the Hadamard formula
   applied to the $h_2$ ODE:
   $$
   \alpha_\ell^{\rm phys} = \frac{\partial \ln k_2}{\partial y}\,\delta y_\ell\,,
   $$
   where $\delta y_\ell$ is the shift in the logarithmic derivative
   induced by the $\ell$-th deformation harmonic.  For the $\Gamma = 2$
   polytrope with $C \sim 0.15$, we find
   $\partial \ln k_2/\partial y \sim -0.5$, so the physical correction
   partially counteracts the geometric one.

### Numerical estimates

For the DASSL duck with
$\sum_\ell \sum_m |\varepsilon_{\ell m}|^2 \simeq 0.1$, the geometric
correction alone gives

$$
\frac{\Delta\Lambda}{\Lambda}\bigg|_{\rm geom}
\simeq 15 \times 0.1 = 150\%\,.
$$

Including the physical correction reduces this somewhat, but the total
correction remains $\mathcal{O}(100\%)$ --- a dramatic effect.

### Mode-by-mode decomposition

**Table 5.** Mode-by-mode fractional correction $\Delta\Lambda/\Lambda$
from each multipole order $\ell$, for the DASSL duck deformation
spectrum.

| $\ell$ | $\sum_m \|\varepsilon_{\ell m}\|^2$ | $\Delta\Lambda/\Lambda$ (geom.) | Role |
|:------:|:-----------------------------------:|:-------------------------------:|:----:|
| 1      | 0.008 | 12%   | COM shift (gauge) |
| 2      | 0.045 | 68%   | Dominant (body shape) |
| 3      | 0.015 | 23%   | Odd parity |
| 4      | 0.020 | 30%   | Secondary (neck/bill) |
| 5      | 0.005 | 8%    | Small |
| $\geq 6$ | 0.007 | 11% | Negligible individually |

The quadrupole ($\ell = 2$) deformation dominates, contributing roughly
half of the total correction.  The hexadecapole ($\ell = 4$) is
secondary.  Multipoles $\ell \geq 6$ contribute less than $\sim 2\%$
each.


## 4.5 I-Love-Q Breakdown

### The Yagi-Yunes universal relation

Yagi & Yunes (2013) discovered that the (dimensionless) moment of
inertia $\bar{I} = I/M_*^3$, tidal deformability
$\bar{\Lambda} = \Lambda$, and spin-induced quadrupole moment
$\bar{Q} = -Q/(M_*^3 \chi^2)$ of neutron stars satisfy approximately
EoS-independent relations.  The I-Love relation is fit by:

$$
\ln \bar{I} = 1.47 + 0.0817\,\ln\bar{\Lambda}
  + 0.0149\,(\ln\bar{\Lambda})^2
  + 2.87\times 10^{-4}\,(\ln\bar{\Lambda})^3
  + 3.64\times 10^{-5}\,(\ln\bar{\Lambda})^4\,.
$$

For spherical neutron stars across a wide range of realistic EoS, this
relation holds to better than $\sim 1\%$ (Yagi & Yunes 2013).

### Breakdown for the duck

The universality relies on the approximate self-similarity of
spherically symmetric stellar structure across different EoS.  The duck
breaks this in two ways:

1. **Modified $\bar{I}$.**  The moment of inertia of the duck differs
   from that of its equivalent sphere by $\Delta I/I \sim \kappa_{\rm duck}
   \sim \mathcal{O}(\varepsilon^2) \sim 10\%$, with the exact value
   depending on orientation (the inertia tensor has three distinct
   principal moments $I_1 \neq I_2 \neq I_3$).

2. **Modified $\bar{\Lambda}$.**  As shown in Sec.\ 4.4,
   $\Delta\Lambda/\Lambda \sim \mathcal{O}(100\%)$.

Combining these, the duck's position in the $\bar{I}$-$\bar{\Lambda}$
plane is displaced from the universal curve by

$$
\frac{\delta(\ln\bar{I})}{\delta(\ln\bar{\Lambda})}
\bigg|_{\rm duck} \sim \frac{0.1}{1.0} \sim 0.1 \neq 0.08\,,
$$

and the absolute deviation $|\Delta\ln\bar{I}|$ is
$\mathcal{O}(100\%)$, compared to $\mathcal{O}(1\%)$ for spherical
stars.  *The duck maximally violates I-Love-Q universality.*

### A sphericity test

The deviation from the Yagi-Yunes relation can be used as a
*sphericity test*: any object whose measured $(\bar{I}, \bar{\Lambda})$
falls far from the universal curve must be significantly non-spherical.
For the duck, the expected deviation is two orders of magnitude larger
than the EoS-driven scatter for spherical stars, making the duck
unambiguously distinguishable --- if both $\bar{I}$ and $\bar{\Lambda}$
can be independently measured.

In practice, $\bar{I}$ is measurable from spin-orbit coupling in double
pulsar systems (Lattimer & Schutz 2005), while $\bar{\Lambda}$ is
extracted from gravitational-wave inspiral signals (Flanagan & Hinderer
2008).  A joint radio-GW observation of a binary duck system would
reveal the I-Love-Q violation.


## 4.6 Comparison Table

**Table 4.** Tidal apsidal constant $k_2$ and dimensionless tidal
deformability $\Lambda$ for selected objects.

| Object | $C = M/R$ | $k_2$ | $\Lambda$ | Notes |
|:-------|:---------:|:-----:|:---------:|:------|
| Nuclear duck ($\Gamma=2$, $C\simeq 0.146$) | 0.146 | 0.074 | $\sim 740$ | TOV + $h_2$ integration |
| Nuclear duck (range) | 0.10--0.20 | 0.06--0.10 | 290--880 | Varying $\rho_c$ |
| Rubber duck (bath toy) | $\sim 10^{-25}$ | $2.7\times 10^{-8}$ | --- | Elastic formula (PVC); $\Lambda$ not meaningful |
| NS --- SLy EoS (1.4 $M_\odot$) | 0.176 | 0.091 | 306 | Hinderer et al.\ (2010) |
| NS --- APR EoS (1.4 $M_\odot$) | 0.170 | 0.085 | 261 | Hinderer et al.\ (2010) |
| NS --- GW170817 | $\sim 0.16$ | --- | $190^{+390}_{-120}$ | Abbott et al.\ (2018) |
| Black hole | 0.500 | 0 | 0 | Binnington & Poisson (2009) |

The nuclear duck's Love number and tidal deformability are comparable to
those of realistic neutron stars with similar compactness, confirming
that the EoS (rather than the shape) is the primary determinant of the
*angle-averaged* tidal response.  The shape enters at
$\mathcal{O}(\varepsilon^2)$ through the corrections derived in
Sec.\ 4.4.  The rubber duck's vanishingly small $k_2$ reflects the
dominance of elastic rigidity over self-gravity at terrestrial scales.
The black hole limit $k_2 = 0$ (Binnington & Poisson 2009; Damour &
Nagar 2009) serves as a lower bound: compact objects with $k_2 = 0$
are maximally stiff against tidal deformation.
