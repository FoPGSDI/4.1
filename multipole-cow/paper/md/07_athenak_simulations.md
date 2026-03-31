# Section 7: AthenaK Numerical Relativity Simulations

We validate the perturbative predictions of Sections 3--6 with fully nonlinear
numerical relativity simulations using the AthenaK code
(Stone et al.\ 2024; Zhu et al.\ 2024).  The central result of this section is the
visual centerpiece of the paper: a duck-shaped TOV star, released from rest,
melts to a sphere under its own gravity while radiating gravitational waves
at the split QNM frequencies predicted by Section 3.

## 7.1 Duck-Shaped TOV Star Initial Data

### Spherical baseline

We begin with a standard Tolman--Oppenheimer--Volkoff solution for a
$\Gamma = 2$ polytrope with adiabatic constant $K = 100$ (geometric units,
$G = c = M_\odot = 1$) and central density $\rho_c = 1.28 \times 10^{-3}$.
This produces a stable equilibrium star with gravitational mass
$M \approx 1.4\,M_\odot$, circumferential radius $R_* \approx 8.13\,M$
(compactness $C \equiv M/R_* \approx 0.17$), and dynamical timescale
$T_{\rm dyn} = 2\pi/\omega_0 \approx 0.5\,$ms.

The TOV solver integrates the structure equations outward from the center
using $N_{\rm pts} = 10{,}000$ radial points with step $\delta r = 10^{-3}$,
producing profiles $\rho_{\rm TOV}(r)$, $P_{\rm TOV}(r)$, $m(r)$, and the
lapse $\alpha(r)$.

### Mapping onto the duck shape

The duck-shaped initial data is constructed by a radial rescaling of the
spherical TOV profile.  Given the duck surface function from Section 2,
$$
R_{\rm duck}(\theta,\phi) \;=\; R_0\!\left[1 + \sum_{\ell=1}^{\ell_{\max}}
\sum_{m=-\ell}^{\ell} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right],
$$
we define the **effective radial coordinate** at each grid point
$(r, \theta, \phi)$:
$$
r_{\rm eff}(r,\theta,\phi) \;=\; r \cdot \frac{R_0}{R_{\rm duck}(\theta,\phi)}.
$$
The hydrodynamic primitives are then assigned by querying the TOV solution at
the rescaled radius:
$$
\rho(r,\theta,\phi) = \rho_{\rm TOV}(r_{\rm eff}), \qquad
P(r,\theta,\phi) = P_{\rm TOV}(r_{\rm eff}).
$$
The mapping preserves the EoS relation $P = K\rho^\Gamma$ at every point.
Points with $r_{\rm eff} > R_*$ (outside the rescaled stellar surface) are
assigned atmosphere values $\rho = \rho_{\rm floor}$, $P = P_{\rm floor}$.

The deformation coefficients $\varepsilon_{\ell m}$ are read from a data file
(`duck_epsilon.dat`) containing real spherical harmonic coefficients up to
$\ell_{\max} = 10$.  This file is produced by the SH fitting pipeline of
Section 2.  For the duck, the dominant deformation amplitudes are
$|\varepsilon_{1,0}| \approx 2.0$ (dipole/COM shift),
$|\varepsilon_{2,0}| \approx 2.7$ (prolate elongation), and
$|\varepsilon_{3,0}| \approx 3.5$ (bill--tail asymmetry).  An overall
amplitude parameter $\epsilon$ can be used to uniformly scale all coefficients
to study the linear regime ($\epsilon \ll 1$) or the fully nonlinear duck
($\epsilon = 1$).

### Metric initial data

The spatial metric and extrinsic curvature are initialized from the
*spherical* TOV solution evaluated at the *physical* radius $r$ (not $r_{\rm eff}$):
$$
\gamma_{ij} = \delta_{ij} + \frac{x_i x_j}{r^2}\left(\frac{1}{1 - 2m(r)/r} - 1\right),
\qquad \alpha = \alpha_{\rm TOV}(r), \qquad \beta^i = 0, \qquad K_{ij} = 0.
$$
This is the metric of a spherical star with the same total mass.  The
duck-shaped density perturbation introduces a Hamiltonian constraint violation
of order $\mathcal{O}(\varepsilon^2)$, since the constraint equations are
quadratic in the matter source terms and we have perturbed only the matter.
For the scaled runs with $\epsilon \lesssim 0.3$, this violation is
$\lesssim 10^{-2}$ in the $L^2$ norm and is rapidly damped by the Z4c
constraint-damping terms (typical e-folding time $\sim 5M$).

### Velocity perturbation

A small radial velocity perturbation excites the fundamental modes:
$$
v^r(r) = \frac{v_{\rm pert}}{2}\left(3\frac{r}{R_*} - \frac{r^3}{R_*^3}\right),
\qquad v_{\rm pert} = -0.01,
$$
applied only inside $r \leq R_*$.  The cubic profile vanishes at both the
center and the surface, avoiding spurious reflections.

---

## 7.2 Evolution Setup

### Spacetime evolution

We evolve the Einstein equations using the Z4c formulation with standard
gauge conditions:
- **Lapse:** $1{+}\log$ slicing, $\partial_t \alpha = -2\alpha K$.
- **Shift:** Gamma-driver condition,
  $\partial_t \beta^i = \frac{3}{4}\tilde{\Gamma}^i - \eta\,\beta^i$,
  with damping parameter $\eta = 2/M$.

Z4c provides constraint damping through its auxiliary variables
$\Theta$ and $\hat{\Gamma}^i$, which exponentially suppress the
$\mathcal{O}(\varepsilon^2)$ constraint violation from the non-self-consistent
initial data within a few crossing times.

### GR hydrodynamics

The general-relativistic hydrodynamic equations are solved using:
- **Reconstruction:** PPMx (piecewise parabolic method, extrema-preserving).
- **Riemann solver:** HLLE.
- **Time integration:** RK3 with CFL number 0.4.
- **EoS:** Ideal gas, $\Gamma = 2$.
- **Floors:** $\rho_{\rm floor} = 10^{-10}$, $T_{\rm floor} = 10^{-8}$,
  with density threshold 1.02 for floor activation.
- **Error policy:** `reset_floor` (reset to atmosphere on conserved-to-primitive
  failure).
- **Magnetic field:** $B = 0$ (unmagnetized).

### Grid and AMR

The computational domain is a Cartesian box
$[-L, L]^3$ with $L = 102.4\,M$.  Boundary conditions are
reflecting on inner faces and outflow (diode) on outer faces.
The base resolution is $128^3$ with MeshBlock size $64^3$, giving a base
grid spacing $\Delta x = 1.6\,M$.  Three levels of adaptive mesh refinement
(AMR) refine on the density gradient, achieving an effective resolution of
$\Delta x_{\rm fine} = 0.2\,M \approx R_*/40$ across the star.

### Duration

The simulation is evolved for $t_{\rm final} = 10{,}000\,M \approx 50\,$ms
$\approx 100\,T_{\rm dyn}$, sufficient to observe multiple oscillation
cycles and to measure damping rates for the lowest-order modes.  History
data ($\rho_{\max}$, $\alpha_{\min}$) are recorded at every time step
(dt = $10^{-5}$).  3D binary snapshots of hydrodynamic primitives are
written at intervals $\Delta t = 10\,M$ for post-processing.

---

## 7.3 The Duck Melts to a Sphere

The duck-shaped TOV star is not in hydrostatic equilibrium: the isotropic
pressure gradient $-\nabla P$ does not balance the gravitational acceleration
$-\rho\,\nabla\Phi$ on the non-spherical boundary.  The star therefore
oscillates and relaxes toward its unique spherical equilibrium configuration.
This process --- the duck melting to a sphere --- is the visual centerpiece
of the paper (Fig.\ 11).

### Hierarchy of relaxation timescales

The relaxation is not uniform: high-$\ell$ surface features decay faster
than low-$\ell$ ones.  This is a direct consequence of the mode-dependent
GW damping timescale:
$$
\tau_{\rm GW}(\ell) \;\sim\; \frac{c^5}{G\,\omega^{2\ell+2}\,R^{2\ell}}.
$$
For the $\ell = 2$ f-mode of a $1.4\,M_\odot$ star,
$\tau_{\rm GW} \sim 0.1$--$1\,$s, while for $\ell = 4$,
$\tau_{\rm GW}$ is shorter by a factor $\sim (\omega R/c)^4 \sim 10^{-4}$.
In practice, the dominant damping mechanism at early times is the
nonlinear coupling and redistribution of energy among modes, which operates
on the dynamical timescale $T_{\rm dyn}$.

### Morphological evolution

The evolution proceeds through a characteristic sequence:

1. **$t = 0$:** The initial duck-shaped density contours show the bill,
   body, and tail clearly delineated in the $\rho = 0.5\rho_c$ isosurface.

2. **$t \sim T_{\rm dyn}/4$:** Large-scale breathing oscillation begins.
   The star contracts globally as the excess gravitational potential energy
   from the non-spherical shape converts to kinetic energy.

3. **$t \sim T_{\rm dyn}$:** The highest-$\ell$ features --- the bill tip
   ($\ell \gtrsim 6$) and tail curl ($\ell \sim 5$--$8$) --- have already
   smoothed.  The star resembles a prolate spheroid with a residual bump
   where the head was.

4. **$t \sim 2$--$5\,T_{\rm dyn}$:** The quadrupolar ($\ell = 2$) elongation
   oscillates coherently.  The density isosurface alternates between
   prolate and oblate configurations as the f-mode rings.

5. **$t \sim 10$--$50\,T_{\rm dyn}$:** Only the lowest-order modes
   ($\ell = 2$, possibly $\ell = 3$) remain excited.  The oscillation
   amplitude decays through GW emission and numerical viscosity.

6. **$t \to \infty$:** The star settles to the unique spherical TOV
   equilibrium, confirming the stability of the background solution.

### Physical interpretation

The melting process illustrates the spectral decomposition theorem of
Section 3 in the time domain.  The initial duck shape, expanded as
$\delta\rho/\rho_0 = \sum_{\ell,m} a_{\ell m}(t)\,Y_\ell^m$,
evolves as
$$
a_{\ell m}(t) = a_{\ell m}(0)\,e^{-t/\tau_\ell}\cos(\omega_\ell t + \phi_\ell),
$$
where each mode decays independently at its own frequency and damping rate.
The bill melts first because it is dominated by $\ell \geq 6$ harmonics,
which have the shortest damping times.  The body elongation persists longest
because $\ell = 2$ has the longest $\tau_{\rm GW}$.

---

## 7.4 Mode Extraction

### Central density oscillations

The maximum rest-mass density $\rho_{\max}(t)$, recorded in the history
file at high cadence ($\Delta t = 10^{-5}\,M$), traces the radial ($\ell = 0$)
oscillation mode.  We extract the QNM spectrum via:

1. **Windowing:** Apply a Hanning window to the time series
   $\delta\rho(t) = \rho_{\max}(t) - \langle\rho_{\max}\rangle$ to suppress
   spectral leakage.
2. **FFT:** Compute the power spectral density $P(f) = |\tilde{\delta\rho}(f)|^2$.
3. **Peak identification:** Locate peaks in $P(f)$ corresponding to the
   fundamental ($F$) and pressure ($p_1$, $p_2$, ...) modes.
4. **Damping rate:** Fit each peak to a Lorentzian profile
   $P(f) \propto [(f - f_R)^2 + f_I^2]^{-1}$ to extract the imaginary
   frequency $\omega_I = 2\pi f_I$.

For the duck, the key observable is the **splitting** of the $\ell = 2$
f-mode peak.  In the spherical limit, a single peak appears at
$f_{F} \approx 1.58\,$kHz.  For the duck, this peak splits into up to
$2\ell + 1 = 5$ sub-peaks, with frequency separations
$\Delta f/f \sim \mathcal{O}(\varepsilon)$ as predicted by Section 3.5.

### Weyl scalar $\Psi_4$ decomposition

The Newman--Penrose scalar $\Psi_4$ is extracted on coordinate spheres at
radii $r_{\rm ext} = 40, 60, 80, 100\,M$ and decomposed into spin-weighted
spherical harmonics:
$$
\Psi_4(t, r_{\rm ext}, \theta, \phi) = \sum_{\ell \geq 2}\sum_{m=-\ell}^{\ell}
\Psi_4^{\ell m}(t, r_{\rm ext})\;{}_{-2}Y_{\ell m}(\theta, \phi).
$$
Each $(\ell, m)$ channel is fit to a damped sinusoid:
$$
\Psi_4^{\ell m}(t) = A_{\ell m}\,e^{-t/\tau_{\ell m}}\,
\sin(\omega_{\ell m}\,t + \phi_{\ell m}),
$$
yielding the complex QNM frequency
$\tilde{\omega}_{\ell m} = \omega_{\ell m} + i/\tau_{\ell m}$ for each
angular channel independently.

For the spherical star, modes with different $m$ but the same $\ell$ have
identical complex frequencies.  For the duck, the $(2\ell+1)$-fold degeneracy
is lifted: the five $\ell = 2$ channels yield five distinct
$(\omega_{2m}, \tau_{2m})$ pairs, constituting direct numerical evidence for
the gravitational Zeeman effect.

### Cross-validation

The extracted frequencies are compared against:
- **Perturbative prediction** (Section 3): $\omega_{n\ell m}^{(\rm pert)}$
  from degenerate perturbation theory.
- **FEM eigenvalues** (Section 6): $\omega_{n\ell m}^{(\rm FEM)}$ from
  Helmholtz eigenvalue analysis.
- **Linear perturbation theory** (Cowling approximation): $\omega_{n\ell}^{(\rm Cowling)}$
  for the unperturbed sphere.

Agreement between the three tiers validates the perturbative framework and
confirms that the AthenaK simulation resolves the QNM spectrum correctly.
Results are collected in Table 9 (AthenaK frequencies) and Table 8
(cross-tier comparison).

---

## 7.5 Waveform Generation

### Gravitational wave strain

The GW strain is obtained by double time-integration of $\Psi_4$:
$$
h_+(t) - i\,h_\times(t) = -\int_0^t dt'\int_0^{t'} dt''\;\Psi_4(t''),
$$
with fixed-frequency integration (FFI) to suppress low-frequency drift.
The resulting ringdown waveform decomposes as:
$$
h_+(t) + i\,h_\times(t) = \frac{1}{D_L}\sum_{n,\ell,m}
\mathcal{A}_{n\ell m}\,e^{-t/\tau_{n\ell m}}\,
e^{i\omega_{n\ell m}\,t}\;{}_{-2}Y_{\ell m}(\iota, \varphi),
$$
where $D_L$ is the luminosity distance, $(\iota, \varphi)$ the sky
orientation of the observer, and the sum runs over all excited QNM
overtones $(n)$, angular degrees $(\ell)$, and azimuthal orders $(m)$.

### The duck ringdown signature

For a spherical star, the ringdown waveform has a simple structure:
each $\ell$ contributes a single damped sinusoid (all $m$ substates
share the same $\omega$ and $\tau$).  For the duck, the broken spherical
symmetry produces a richer signal:

- The $\ell = 2$ sector alone contributes **five** distinct damped
  sinusoids with frequencies $\omega_{2,m}$ ($m = -2, \ldots, +2$),
  each with its own damping time $\tau_{2,m}$.
- The interference between these closely spaced frequencies produces
  **beating** in the waveform envelope, with beat period
  $T_{\rm beat} \sim 2\pi/\Delta\omega \sim T_{\rm dyn}/\varepsilon$.
- The relative amplitudes $\mathcal{A}_{n\ell m}$ depend on the
  orientation of the duck's symmetry axes relative to the observer,
  introducing an anisotropic radiation pattern absent in the spherical case.

### Template waveforms

We generate template waveforms at representative distances and orientations
for injection into simulated detector noise (Section 8.2).  For a
$1.4\,M_\odot$ duck at $D_L = 50\,$Mpc observed along the duck's
bill axis ($\iota = 0$), only $m = \pm 2$ modes are excited (by selection
rules of ${}_{-2}Y_{\ell m}$), producing a clean two-frequency ringdown.
Edge-on observation ($\iota = \pi/2$) excites all five $m$-substates,
maximizing the spectroscopic information content.

The resulting waveforms, along with their Fourier transforms overlaid on
detector noise curves, are presented in Figures 13 and 17.
