# Quasi-Normal Modes of Relativistic Stars

Pulsating stars are important sources of information for
astrophysics. Nearly every star undergoes some kind of pulsation
during its evolution from the early stages of formation until the very
late stages, usually the catastrophic creation of a compact object
(white dwarf, neutron star or black hole). Pulsations of supercompact
objects are of great importance for relativistic astrophysics since
these pulsations are accompanied by the emission of gravitational
radiation. Neutron star oscillations were also proposed to explain the
quasi-periodic variability found in radio-pulsar and X-ray burster
signals [vanhorn, MVH88]. In this chapter we shall discuss
various features of neutron star non-radial pulsations i.e.\ the
various modes of pulsation,  mode excitation, detection probability
and the possibility to  extract information (to estimate, for example,
the radius, mass and stellar equation of state) from the detection of
the associated gravitational waves. It is not in our plans to discuss
rotating relativistic stars; the interested reader should refer to
another review in this journal [Nick98]. Radial oscillations are
also not discussed since they are not interesting for gravitational
wave research.

## Stellar Pulsations: The Theoretical Minimum

For the study of stellar oscillations we shall consider a spherically
symmetric and static spacetime which can be described by the
Schwarzschild solution outside the star, see
equation ([ref:metric1]). Inside the star, assuming that the stellar
material is behaving like an ideal fluid, we define the energy momentum
tensor

$$
  T_{\mu\nu}=(\rho +p)u_\mu u_\nu + p g_{\mu \nu},
$$

where $p(r)$ is the pressure, $\rho(r)$ is the total energy
density. Then from the conservation of the energy-momentum and the
condition for hydrostatic equilibrium we can derive the
Tolman-Oppenheimer-Volkov (TOV) equations for the interior of a
spherically symmetric star in equilibrium. Specifically,

$$
  e^{-\lambda} = 1 - {2 m(r) \over r},
$$

and the “mass inside radius $r$” is represented by

$$
  m(r) = 4\pi \int_0^r \rho r^2 dr.
$$

This means that the total mass of the star is $M=m(R)$, with $R$ being
the star’s radius. To determine a stellar model we must solve

$$
  {dp \over dr} = - {\rho + p \over 2} {d\nu \over dr},
$$

where

$$
  {d\nu \over dr} = {2 e^\lambda(m + 4\pi p r^3) \over r^2}.
$$

These equations should of course, be supplemented with  an equation of
state $p = p (\rho, \dots)$ as input. Usually is sufficient to use a
one-parameter equation of state to model neutron stars, since the
typical thermal energies are much smaller than the Fermi energy. The
polytropic equation of state $p=K \rho^{1+1/N}$ where $K$ is the
polytropic constant and $N$ the polytropic exponent, is used in most
of the studies. The existence of a unique global solution of the
Einstein equations for a given equation of state and a given value of
the central density has been proven by Rendall and
Schmidt [RS91].

If we assume a small variation in the fluid or/and in the spacetime
we must deal with the perturbed Einstein equations

$$
  \delta \left( G^{\mu}_{\nu} -  {{8 \pi G}\over {c^4}}
    T^{\mu}_{\nu}\right) = 0,
$$

and the variation of the fluid equations of motion

$$
  \delta \left( T^{\mu}_{\nu ; \mu} \right) = 0,
$$

while the perturbed metric will be given by equation ([ref:scw_pert]).

Following the procedure of the previous section one can decompose the
perturbation equations into spherical harmonics. This decomposition
leads to two classes of oscillations according to the parity of the
harmonics (exactly as for the black hole case). The first ones called
{\em even} (or spheroidal, or polar) produce spheroidal deformations
on the fluid, while the second are the {\em odd} (or toroidal, or
axial) which produce toroidal deformations.

For the {\em  polar} case one can use certain combinations of the
metric perturbations as unknowns, and the linearized field equations
inside the star will be equivalent to the following system  of three
wave equations for unknowns $S, F, H$:

$$
  -{1\over c^2}{\partial^2S\over \partial^2t}
  +{\partial^2S\over \partial^2r_*} + L_1(S,F,\ell) = 0,

$$

$$
  -{1\over c^2}{\partial^2F\over \partial^2t}
  +{\partial^2F\over \partial^2r_*} + L_2(S,F,H,\ell) = 0,

$$

$$
  -{1\over (c_s)^2}{\partial^2H\over \partial^2t}
  +{\partial^2H\over \partial^2r_*} + L_3(H,H’,S,S’,F,F’,\ell) = 0,

$$

and the constraint

$$
  {\partial^2F\over \partial^2r_*}+L_4(F,F’,S,S’,H,\ell) = 0.

$$

The linear functions $L_i$, ($i=1, 2, 3, 4$) depend on the background
model and their explicit form can be found in [KES93, aaks]. The
functions $S$ and $F$ correspond to the perturbations of the spacetime
while the  function $H$ is proportional to the density perturbation
and is only defined on the background star. With $c_s$ we define the
speed of sound and with a prime we denote differentiation  with
respect to $r_*$:

$$
  {\partial \over {\partial r_*}} = e^{(v-\lambda)/2}
  {\partial \over {\partial r}}.
$$

Outside the star there are only perturbations of the spacetime. These
are described by a single wave equation, the Zerilli equation
mentioned in the previous section, see equations ([ref:qnmwave])
and ([ref:zerpot]). In [KES93] it was shown that (for background
stars whose boundary density is positive) the above system -- together
with the geometrical transition conditions at the boundary of the star
and regularity conditions at the center -- admits a well posed Cauchy
problem. The constraint is preserved under the evolution. We see that
two variables propagate along light characteristics and the density
$H$ propagates with the sound velocity of the background star.

It is possible to eliminate the constraint -- first done by
Moncrief [Monc74c] -- if one solves the
constraint ([ref:starconstr]) for $H$ and puts the corresponding
expression into $L_2$. (The characteristics for $F$ change then to
sound characteristics inside the star and light characteristics
outside.) This way one has just to solve two coupled wave equations
for $S$ and $F$ with unconstrained data, and to calculate $H$ using the
constraint from the solution of the two wave equations. Again the
explicit form of the equation can be found in [aaks].

Turning next to quasi-normal modes in the spirit of
section [ref:section_2], we can Laplace transform  the two wave
equations and obtain a system of ordinary differential equations which
is of fourth order. The Green function can be constructed from
solutions of the homogeneous equations (having the appropriate
behavior at the center and infinity)  and its analytic continuation
may have poles  defining the quasi-normal mode frequencies.

From the form of the above equations one can easily see two limiting
cases. Let us first assume that the gravitational field is very
weak. Then equation ([ref:starper1]) and ([ref:starper2]) can be
omitted (actually $S \to 0$ in the weak field
limit [Thorne69b, aaks]) and we find that one equation is enough
to describe (with acceptable accuracy) the oscillations of the
fluid. This approach is known as the Cowling
approximation [cowl41]. Inversely, we can assume that the
coupling between the two equations ([ref:starper1])
and ([ref:starper2]) describing the spacetime perturbations with the
equation ([ref:starper3]) is weak and consequently derive  all the
features of the spacetime perturbations from only the two of
them. This is what is called the “inverse Cowling approximation”
(ICA) [AKS96].

For the {\em axial} case the perturbations reduce to a single wave
equation for the spacetime perturbations which describes toroidal
deformations

$$
  -{1\over c^2}{\partial^2 X\over \partial^2t}
  +{\partial^2X\over \partial^2r_*}
  + {e^v\over r^3}\left[\ell(\ell+1)r+r^3(\rho-p)-6M\right] = 0,

$$

where $X\sim h_{r\phi}$. Outside the star, pressure and density are
zero and this equation is reduced to the Regge-Wheeler equation, see
equations ([ref:qnmwave]) and ([ref:zerpot]). In Newtonian theory, if
the star is non-rotating and the static model is a perfect fluid
(i.e.\ shear stresses are absent), the {\em axial} oscillations are a
trivial solution of zero frequency to the perturbation equations and
the variations of pressure and density are zero. Nevertheless, the
variation of the velocity field is not zero and produces
non-oscillatory eddy motions. This means that there are no oscillatory
velocity fields. In the relativistic case the picture is
identical [TC67] nevertheless; in this case there are still QNMs,
the ones that we will describe later as “{\em spacetime or
$w$-modes}” [kokkotas94].

When the star is set in slow rotation then the axial modes are no longer
degenerate, but instead a new family of modes emerges, the so-called
$r$-modes. An interesting property of these modes that has been
pointed out by Andersson [nils98, FM98] is that these modes are
generically unstable due the Chandrasekhar-Friedman-Schutz
instability [Chandra70, FS78] and furthermore it has been
shown [AKS98, LOM98] that these modes can potentially restrict
the rotation period of  newly formed neutron stars and also that they
can radiate away detectable amounts of gravitational
radiation [Owen98]. The equations describing the perturbations of
slowly rotating relativistic stars have been derived by
Kojima [kojima92, kojima93], and Chandrasekhar and
Ferrari [chandra91b].

## Mode Analysis

The study of stellar oscillations in a general relativistic context
already has a history of 30 years. Nevertheless, recent results have
shown remarkable features which had previously been overlooked.

Until recently most studies treated the stellar oscillations in a
nearly Newtonian manner, thus practically ignoring the dynamical
properties of the spacetime [TC67, Thorne68, PT69, Thorne69a,
  Thorne69b, lindblom83, MVS83, DL85, MVH88]. The spacetime was used
as the medium upon which the gravitational waves, produced by the
oscillating star, propagate. In this way all the families of modes
known from Newtonian theory were found for relativistic stars while in
addition the damping times due to gravitational radiation were
calculated.

Inspired by a simple but instructive model  [kokkotas86], Kokkotas
and Schutz showed the existence of a new family of
modes: the $w$-modes [kokkotas92]. These are {\em spacetime modes} and their
properties, although different, are closer to the black hole QNMs than
to the standard fluid stellar modes. The main characteristics of the
$w$-modes are high frequencies accompanied with very rapid
damping. Furthermore, these modes hardly excite any fluid motion. The
existence of these modes has been verified by subsequent
work [leins93, AKS95]; a part of the spectrum was found earlier
by Kojima [kojima88] and it has been shown that they exist also
for odd parity (axial) oscillations [kokkotas94]. Moreover,
sub-families of $w$-modes have been found for both the polar and axial
oscillations i.e.\ the {\em interface} modes found by Leins et
al. [leins93] (see also [AKK96]), and the {\em trapped}
modes found by Chandrasekhar and Ferrari [chandra91c] (see
also [kokkotas94, kojima95, AKK96]). Recently, it has been proven
that one can reveal all the properties of the $w$-modes even if one
“freezes” the fluid oscillations (Inverse Cowling
Approximation) [AKS96]. In the rest of this section we shall
describe the features of both  families of oscillation modes, fluid and spacetime, for the
case $\ell=2$.

### Families of Fluid Modes

For non-rotating stars the fluid modes exist {\em only for polar}
oscillations. Here we will describe the properties of the most
important modes for gravitational wave emission. These are the
fundamental, the pressure and the gravity modes; this division has
been done in a phenomenological way by Cowling [cowl41]. For an
extensive discussion of other families of fluid modes we refer the
reader to [gau95, gau96] and [MVS83, MVH88]. Tables of
frequencies and damping times of neutron star oscillations for twelve
equations of state can be found in a recent work [AK98] which
verifies and extends earlier work [lindblom83]. In Table 2 we
show characteristic frequencies and damping times of various QNM
modes for a typical neutron star.

| mode | frequency | damping time |
|---|---|---|
| $f$ | 2.87 kHz | 0.11 sec |
| $p_1$ | 6.57 kHz | 0.61 sec |
| $g_1$ | 19.85 Hz | years |
| $w_1$ | 12.84 kHz | 0.024 ms |
| $w_{II}$ | 8.79 kHz | 0.016 ms |

- The $f$-mode ({\em f}undamental) is a {\em stable} mode which
  exists only for non-radial oscillations. The frequency is
  proportional to the mean density of the star and it is nearly
  independent of the details of the stellar structure. An exact
  formula for the frequency can be derived for Newtonian uniform
  density stars
  $$
  \omega^2 = {{2\ell(\ell-1)}\over {2\ell+1}} {M\over R^3} \ .
  $$
  This relation is approximately correct also for the relativistic
  case [AKK96] (see also the discussion in
  section [ref:section_5_4]). The $f$-mode eigenfunctions have {\em no
  nodes} inside the star, and they grow towards the surface. A
  typical neutron star has an $f$-mode with a frequency of $1.5-3$ kHz
  and the damping time of this oscillation is less than a second
  ($0.1-0.5$ sec). Detailed data for the frequencies and damping times
  (due to gravitational radiation) of the $f$-mode for various
  equations of state can be found in [lindblom83,
  AK98]. Estimates for the damping times due to viscosity can be
  found in [CL87, CLS90].
- The $p$-modes ({\em p}ressure or acoustic) exist for both radial
  and non-radial oscillations. There are infinitely many of them. The
  pressure is the restoring force and it experiences substantial
  fluctuations when these modes are excited. Usually, the radial
  component of the fluid displacement vector is significantly larger
  than the tangential component. The oscillations are thus nearly
  radial. The frequencies depend on the travel time of an acoustic
  wave across the star. For a neutron star the frequencies are
  typically higher than $4-7$ kHz ($p_1$-mode) and the damping times
  for the first few $p$-modes are of the order of a few seconds. Their
  frequencies and damping times  increase with the order of the
  mode. Detailed data for the frequencies and damping times (due to
  gravitational radiation) of the $p_1$-mode for various equations of
  state can be found in [AK98].
- The {$g$-modes} ({\em g}ravity) arise because gravity tends to
  smooth out material inhomogeneities along equipotential
  level-surfaces and buoyancy is the restoring force. The changes in
  the  pressure are very small along the star. Usually, the tangential
  components of the fluid displacement vector are dominant in the
  fluid motion. The $g$-modes require a {\em non-zero Schwarzschild
  discriminant} in order to have non-zero frequency, and if they
  exist there are infinitely many of them. If the perturbation is
  stable to convection, the $g$-modes will be stable ($\omega^2>0$);
  if unstable to convection the $g$-modes are unstable ($\omega^2<0$);
  and if marginally stable to convection, the $g$-mode frequency
  vanishes. For typical neutron stars they have frequencies smaller
  than a hundred Hz (the frequency decreases with the order of the
  mode), and they usually damp out in time much longer than a few days
  or even years. For an extensive discussion about $g$-modes in
  relativistic stars refer to [Finn86, Finn87]; and for a study
  of the instability of the $g$-modes of rotating stars to
  gravitational radiation reaction refer to [Lai98].
- The {$r$-modes} ({\em r}otational) in a non-rotating star are
  purely toroidal (axial) modes with vanishing frequency. In a
  rotating star, the displacement vector acquires spheroidal
  components and the frequency in the rotating frame, to first order
  in the rotational frequency $\Omega$ of the star, becomes
  $$
  \omega_r = {{2 m \Omega}\over {\ell (\ell+1)}}.

  $$
  An inertial observer measures a frequency of
  $$
  \omega_i = \omega_r - m\Omega.

  $$
  From ([ref:omega_r]) and ([ref:omega_i]) it can be deduced that a
  counter-rotating (with respect to the star, as defined in the
  co-rotating frame) $r$-mode appears as co-rotating with the star to
  a distant inertial observer. Thus, all $r$-modes with $\ell \geq 2$
  are generically unstable to the emission of gravitational radiation,
  due to the Chandrasekhar-Friedman-Schutz (CFS)
  mechanism [Chandra70, FS78]. The instability is active as long
  as its growth-time is shorter than the damping-time due to the
  viscosity of neutron star matter. Its effect is to slow down, within
  a year, a rapidly rotating neutron star to slow rotation rates and
  this explains why only slowly rotating pulsars are associated with
  supernova remnants [AKS98, LOM98, KS98]. This suggests that the
  $r$-mode instability might not allow millisecond pulsars to be
  formed after an accretion induced collapse of a white
  dwarf [AKS98]. It seems that millisecond pulsars can only be
  formed by the accretion induced spin-up of old, cold neutron
  stars. It is also possible that the gravitational radiation emitted
  due to this instability  by a newly formed neutron star could be
  detectable by the advanced versions of the gravitational wave
  detectors presently under construction [Owen98]. Recently,
  Andersson, Kokkotas and Stergioulas [AKSt98] have suggested
  that the $r$-instability might be responsible for stalling the
  neutron star spin-up in strongly accreting Low Mass X-ray Binaries
  (LMXBs). Additionally, they suggested that the gravitational waves
  from the neutron stars, in such LMXBs, rotating at the instability
  limit may well be detectable. This idea was also suggested  by
  Bildsten [Bild98] and studied in detail by
  Levin [Levin98].

<!-- TODO: figure content -->

### Families of Spacetime or $w$-Modes

The spectra of the three known families of $w$-modes are different but
the spectrum of each family is similar {\em both} for polar and axial
stellar oscillations. As we have mentioned earlier they are clearly
modes of the spacetime and from numerical calculations appear to be stable.

- The \mathbf{curvature modes} are the standard
  $w$-modes [kokkotas92]. They are the most  important for
  astrophysical applications. They are clearly related to the
  spacetime curvature and exist for all relativistic stars. Their main
  characteristic is the rapid damping of the oscillations. The damping
  rate increases as the compactness of the star decreases: For nearly
  Newtonian stars (e.g.\ white dwarfs) these modes have not been
  calculated due to numerical instabilities in the various codes, but
  this case is of marginal importance due to the very fast damping
  that these modes will undergo. One of their main characteristics is
  the absence of significant fluid motion (this is a common feature
  for all families of $w$-modes). Numerical studies have indicated the
  existence of an infinite number of modes; model problems suggest
  this too [kokkotas86, BS93, amodel]. For a typical neutron star
  the frequency of the first $w$-mode is around $5-12$ kHz and
  increases with the order of the mode. Meanwhile, the typical damping
  time is of the order of a few tenths of a millisecond and decreases
  slowly with the order of the mode.
- The \mathbf{trapped modes} exist only for supercompact stars ($R
  \le 3M$) i.e.\ when the surface of the star is inside the peak of
  the gravitational field’s potential barrier [chandra91c,
  kokkotas94]. Practically, the first few curvature modes become
  trapped as the star becomes more and more compact, and even the
  $f$-mode shows similar behavior [kojima95, AKK96]. The trapped
  modes, as with all the spacetime modes, do not induce any
  significant fluid motions and there are only a finite number of them
  (usually less than seven or so). The number of trapped modes
  increases as  the potential well becomes deeper, i.e.\ with
  increasing compactness of the star. Their damping is quite slow
  since the gravitational waves have to penetrate the potential
  barrier. Their frequencies can be of the order of a few hundred Hz
  to a few kHz, while their damping times can be of the order of a few
  tenths of a second. In general no realistic equations of state are
  known that would allow the formation of a sufficiently compact star
  for the trapped modes to be relevant.
- The \mathbf{interface modes} [leins93] are extremely rapidly
  damped modes. It seems that there is only a finite number of such
  modes ($2-3$ modes only) [AKK96], and they are in some ways
  similar to the modes for acoustic waves scattered off a hard
  sphere. They do not induce any significant fluid motion and their
  frequencies can be from 2 to 15 kHz for typical neutron stars while
  their damping times are of the order of less than a tenth of a
  millisecond.

## Stability

The stability of radial oscillations for non-rotating stars in general
relativity is well understood. Especially, the stability of static
spherically symmetric stars can be determined by examining the
mass-radius relation for a sequence of equilibrium stellar models, see
for example Chapter 24 in [mtw]. The radial perturbations are
described by a Sturm-Liouville second order equation with the
frequency of the mode being the eigenvalue $\omega^2$, then for real
$\omega$ the modes will be stable while for imaginary $\omega$ they
will be unstable [Chandra64], see also Chapter 17.2
in [ST83].

The stability of the non-radially pulsating stars (Newtonian or
relativistic) is determined by the Schwarzschild discriminant

$$
  S(r) = {{d p}\over {d r}} - {{\Gamma_1 p}\over {\rho+p}}
  {{d \rho}\over {d r}},
$$

where $\Gamma_1$ is the star’s adiabatic index. This can be understood
if we define the local buoyancy force $f$ per unit volume acting on a
fluid element displaced a small radial distance $\delta r$ to be

$$
  f \sim -g(r) S(r) \delta r,
$$

where $g$ is the local acceleration of gravity. When $S$ is negative
in some region the buoyancy force is positive and the star is unstable
against convection, while when $S$ is positive the buoyancy force is
restoring and the star is stable against convection. Another way of
discussing the stability is through the so-called Brunt-Väisälä
frequency $N^2=g S(r)$ which is the characteristic frequency of the
local fluid oscillations. Following earlier discussions when $N^2$ is
positive, the fluid element undergoes oscillations, while when $N^2$
is negative the fluid is locally unstable. In other words, in
Newtonian theory stability to non-radial oscillations  can be
guaranteed only if $S>0$ everywhere within the star [Cox]. In
general relativity [DI73], this is a sufficient condition, and so
if $S>0$ the quasi-normal modes are stable. For an extensive
discussion of stellar instabilities for both non-rotating and rotating
stars  (which are actually more interesting for the gravitational wave
astronomy) refer to [BFS87, lindblom97, Nick98].

For completeness the same applies as outlined at the end of
section [ref:section_3_3]. A model calculation of Price and
Husain [PH92], however indicated that the nearly Newtonian
quasi-normal modes might be a basis for the fluid
perturbations. Further mathematical investigation is needed to clarify
this issue.