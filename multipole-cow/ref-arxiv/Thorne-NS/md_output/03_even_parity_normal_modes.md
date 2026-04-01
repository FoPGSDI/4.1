# EVEN-PARITY NORMAL MODES OF PULSATION

The study of the even-parity pulsations of the equilibrium configuration is simplified
considerably by analyzing the pulsations into normal modes. For pulsations the only pulsations
one expects to occur in nature are pulsations which generate outgoing gravitational
waves; and these are damped primarily (but not entirely) in models in which the
gravitational radiation is purely outgoing. Because of radiation damping, such normal
modes will not have real frequencies and eigenvalues; rather, their eigenfrequencies
and frequencies will be complex. However, assuming completeness of the set of
complex eigenfunctions with outgoing waves, any real pulsation with purely outgoing
radiation can be expressed as a linear combination of the complex normal modes.

In this section we shall discuss the eigenvalue equations and the boundary conditions
which determine the normal modes of even parity—including real normal modes (standing
gravitational waves), as well as complex normal modes with various mixtures of
outgoing and incoming gravitational radiation. All discussion of the relationship between the normal modes and the real pulsations of an equilibrium configuration will
be delayed until § IV.

## a) Eigenequations

We concentrate our attention on normal modes which belong to a particular even-parity
spherical harmonic $(l,\, M = 0,\, \varpi = [-1]^{l})$ and which have a particular complex
frequency,

$$\omega = \sigma + i/\tau \;. \tag{11}$$

These normal modes can be characterized by their complex radial eigenfunctions $H(r)$,
$K(r)$, $W(r)$, $V(r)$

$$\begin{aligned}
  H_{0}(r,t) &= H(r)e^{i\omega t} \;, & K(r,t) &= K(r)e^{i\omega t} \;, \\
  W(r,t)     &= W(r)e^{i\omega t} \;, & V(r,t) &= V(r)e^{i\omega t} \;.
\end{aligned} \tag{12}$$

As a matter of convention we restrict ourselves to modes with

$$\sigma \geq 0 \;. \tag{13}$$

The remaining modes can be obtained from these by complex conjugation.

Here and throughout the remainder of this paper we ignore the auxiliary,
"non-dynamical" functions $H_{1}$ and $H_{2}$ of the perturbed metric (7b). They can be
computed at any time from the initial-value equations (8).

The differential equations which the eigenfunctions $H$, $K$, $W$, $V$ satisfy are obtained
by substituting expressions (12) into the unsolved initial-value equation (8a) and into
the propagation equations (9):

$$\begin{aligned}
  H' + r^{-1}e^{\lambda}\bigl[l(l+1)/2 + 1 + 4\pi r^{2}(p - \rho)\bigr]H
  &= r K'' \\
  &\quad + e^{\lambda}(3 - 5m/r - 4\pi r^{2}\rho)K
     - r^{-1}e^{\lambda}\bigl[l(l+1)/2 - 1 \\
  &\quad + 8\pi r^{2}(\rho + p)\bigr]K
     + 8\pi r^{-1}(\rho + p)e^{\lambda/2}W'
     + 8\pi r^{-1}\rho' e^{\lambda/2}W \\
  &\quad + 8\pi l(l+1)r^{-1}(\rho + p)e^{\lambda}V \;,
\end{aligned} \tag{14a}$$

$$\begin{aligned}
  {}-e^{\nu-\lambda}K''
  + 2r^{-1}e^{\nu}\bigl[-1 + m/r + 2\pi r^{2}(\rho - p)\bigr]K'
  - \omega^{2}K \\
  + r^{-2}e^{\nu}\bigl[l(l+1) - 2 + 8\pi r^{2}(\rho + p - \gamma p)\bigr]K \\
  + r^{-2}e^{\nu}\bigl[2e^{-\lambda} - 4\pi r^{2}(\rho + p + \gamma p)\bigr]H \\
  - 8\pi r^{-2}e^{\nu-\lambda/2}(\rho + p - \gamma p)W'
  - 8\pi r^{-2}e^{\nu-\lambda/2}(\rho' - p')W \\
  - 8\pi l(l+1)r^{-2}e^{\nu}(\rho + p - \gamma p)V
  &= 0 \;,
\end{aligned} \tag{14b}$$

$$\begin{aligned}
  {}-\bigl\{\gamma p e^{\nu/2}\bigl[r^{-2}e^{-\lambda/2}W'
    + l(l+1)r^{-2}V - H/2 - K\bigr]\bigr\}' \\
  - \omega^{2}r^{-2}e^{(\lambda-\nu)/2}(\rho + p)W
  + \tfrac{1}{2}(\rho + p)e^{\nu/2}(r^{-2}e^{-\lambda/2}\nu')'W \\
  - l(l+1)r^{-2}(\rho + p)(e^{\nu/2})'V
  - \tfrac{1}{2}(\rho + p)e^{-\nu}(He^{3\nu/2})' \\
  + (\rho + p)(Ke^{\nu/2})'
  &= 0 \;,
\end{aligned} \tag{14c}$$

$$\begin{aligned}
  {}-\omega^{2}e^{-\nu}(\rho + p)V
  + l(l+1)r^{-2}\gamma p V
  + r^{-2}\gamma p e^{-\lambda/2}W'
  + r^{-2}p' e^{-\lambda/2}W \\
  - \tfrac{1}{2}(\rho + p + \gamma p)H - \gamma p K
  &= 0 \;.
\end{aligned} \tag{14d}$$

Of these eigenequations, only the three propagation equations (14b)--(14d) involve the
eigenfrequency, $\omega$. Note that the real and imaginary parts of the eigenfunctions are
coupled in the eigenequations (14) only by means of the squared eigenfrequency $\omega^2$.
This permits the existence of purely real, standing-wave eigenfunctions.

The eigenequations (14) comprise a fifth-order differential system; for any given
choice of the complex frequency, $\omega$, there are five linearly independent, complex solutions
to equations (14). Physical boundary conditions outlined in the next two sections
make four of the five solutions unacceptable. The remaining solution is physically acceptable for all choices of $\omega$, providing one permits arbitrary admixtures of ingoing and
outgoing radiation. The restriction to purely outgoing radiation restricts $\omega$ to take on a
discrete set of values, as we shall see in § IIIc.

## b) Boundary Conditions at the Star's Center and Surface

By expanding the eigenequations (14) about $r = 0$ one finds that only two of the
five independent solutions are physically acceptable at the center of the star. (See
Appendix D for analysis.) The two acceptable solutions have the form:

$$K = Ar^l + \cdots, \qquad H = Ar^l + \cdots, \qquad W = Br^{l+1} + \cdots, \tag{15a}$$

$$V = -(B/l)r^l + \cdots, \qquad \text{near } r = 0\;.$$

Here $A$ and $B$ are independent, freely specifiable constants. These acceptable solutions
are characterized by small fluid motions and by small perturbations in the density, the
pressure, and the geometry of spacetime. The unacceptable solutions—which are discussed in Appendix D—are characterized by perturbations in the density, the pressure,
and/or the geometry of spacetime which diverge as $r = 0$.

An expansion of the eigenequations (14) about $r = R$ (surface of star) reveals that
four of the five independent solutions are physically acceptable there. (See Appendix D
for analysis.) The acceptable solutions have the following expansion upon $(R - r)$:

$$\begin{aligned}
  K &= k_0 + k_1(R-r) + \cdots, &
  H &= h_0 + h_1(R-r) + \cdots, \\
  W &= w_0 + w_1(R-r) + \cdots, &
  V &= v_0 + v_1(R-r) + \cdots,
\end{aligned} \tag{15b}$$

near $r = R$.

Of the constants which enter into this expansion only four are freely specifiable. The
remaining constants are fixed by the asymptotic form of the eigenequations and by the
demand that the Lagrangian perturbation in the pressure vanish at the star's surface:

$$\Delta p \equiv \lim_{r \to R^-} \gamma p \left[
  \frac{e^{-\lambda/2}}{r^2} W' - \frac{l(l+1)}{r} V
  + \frac{1}{2} H_1 + K
\right] P_l(\cos\theta) = 0\;. \tag{16}$$

For example, if the pressure-density relation near the star's surface is polytropic and the
adiabatic index behaves smoothly,

$$p = \rho^{1+1/n} + \cdots, \qquad \mu = \bar{\mu}(R-r)^n + \cdots, \tag{17a}$$

$$\gamma = \gamma_1 + \gamma_2(R-r) + \cdots,$$

then the four constants $k_0$, $h_0$, $w_0$ and $w_0$ are freely specifiable in the acceptable solutions
(15b). The eigenequation (14d) fixes $v_0$ in terms of these four arbitrary constants

$$\omega^2 R(R - 2M)^{-1} v_0 = R^{-1}(1 - 2M/R)^{1/2} (N+1) w_0\;, \tag{17b}$$

and the eigenequations (14a)--(14d) all conspire to fix the remaining constants in the
ascending power series. Here $M \equiv m(r = R)$ is the total mass of the star. (In this paper
$M$ is used both as the star's total mass and as the spherical harmonic projection index.)

In general, the one physically unacceptable solution at the star's surface is that
solution for which the Lagrangian change in the pressure does not vanish. For the
example of a polytropic pressure-density relation (eq. [17a]) this unacceptable solution has
fluid displacements which diverge at the star's surface ($W \sim V \sim [R - r]^{-8}$); see
Appendix D for details.

For any arbitrary choice of the complex pulsation frequency, $\omega$, there will be just one
solution to the eigenequations (14) which is well-behaved both at the center of the star
and at the surface. That solution is the unique linear mixture of the two solutions
well behaved at the center, which contains none of the solution unacceptable at the
surface. Henceforth we confine our attention to that unique solution, for any given
choice of $\omega$, which is acceptable both at the center and surface.

## c) Behavior of the Physically Acceptable Eigensolution at $r = \infty$

For any given choice of $\omega$ the physically acceptable eigensolution will be
characterized by the behavior it takes on outside the star by some particular mixture of
ingoing and outgoing gravitational radiation. As is shown in Appendix E, the gravitational
waves have the asymptotic form (solution to eigenqs. [14a] for large $r$):

$$\begin{aligned}
K &= C^{(\mathrm{I})} g_l^{(+)}(r)\, e^{-i(\omega t + \text{SM}\ln r)} + C^{(\mathrm{O})} g_l^{(-)}(r)\, e^{-i(\omega t - \text{SM}\ln r)} \\
H &= i a C^{(\mathrm{I})} g_l^{(+)}(r)\, e^{-i(\omega t + \text{SM}\ln r)} - i a C^{(\mathrm{O})} g_l^{(-)}(r)\, e^{-i(\omega t - \text{SM}\ln r)} \\
&\quad \text{for } r \gg M \text{ and } \{R \to r\}^{-1}
\end{aligned} \tag{18}$$

Here $C^{(\mathrm{I})}$ and $C^{(\mathrm{O})}$ are complex constants which are fixed by the demand that the
external perturbation in the geometry of spacetime join smoothly at the star's surface
to the internal perturbation

$$K,\; K',\; H \text{ and } H \text{ continuous at } r = R \;. \tag{19}$$

Since the time dependence of the normal modes is $e^{i\omega t}$ (cf. eq. [12]), *the constant $C^{(\mathrm{I})}$
is the amplitude for incoming gravitational waves, while $C^{(\mathrm{O})}$ is the amplitude for outgoing
gravitational waves*. For any given choice of the complex eigenfrequency, $\omega$, and of the
spherical harmonic index, $l$, one has a uniquely determined ratio of outgoing amplitude
to incoming amplitude:

$$C^{(\mathrm{O})}/C^{(\mathrm{I})} \text{ is a unique function of } \omega \text{ and } l \;. \tag{20}$$

From previous experience with one-dimensional analyses one expects this—that, for
any given choice of $l$, there should be a discrete spectrum of complex eigenfrequencies
for which $C^{(\mathrm{O})}/C^{(\mathrm{I})}$ is infinite. These eigenfrequencies and the corresponding eigenfunctions make up normal modes with purely outgoing radiation. Similarly, those discrete
eigenfrequencies and functions for which $C^{(\mathrm{O})}/C^{(\mathrm{I})}$ is zero are normal modes with
purely incoming radiation. Finally, those normal modes with real eigenfrequencies
($r = \omega$ in eq. [11]) have real eigenfunctions (cf. eqs. [14]) and represent physically
acceptable, standing waves with equal inward and outward fluxes.

Most of the eigenfunctions (18) far from the star might appear at first sight to be
physically unacceptable. In general, they contain a divergence in the metric perturbation
tensor (eqs. [5], [7b], and [12]) which diverges at large $r$ rather than having the $1/r$
behavior that one expects of a radiation field. The divergence of the metric perturbation
at large $r$ arises from two sources:

1. There is an exponential divergence associated with the imaginary part, $i/r$, of $\omega$.
This exponential divergence is physically acceptable. For $r > 0$ the divergence occurs in the outgoing-wave term and is a result of the exponential damping of stable pulsations of the star.
For $r < 0$ the divergence occurs in the incoming-wave term and is associated with an
input of energy from infinity which increases exponentially with time.

2. There is also a divergence of order $r$ in the amplitude of the metric perturbation
(18). This divergence, which is present even for standing waves ($r = \infty$), has its origin
in the Regge-Wheeler choice of gauge. It is a coordinate-dependent divergence which
can be removed by an infinitesimal coordinate transformation (change of gauge) and,
consequently, it is *a divergence with no physical reality*. We have verified that it has no
physical reality by calculating—using a computer program of Thorne and Zimmerman
(1967)—the Riemann curvature invariants for the metric perturbation (18). We find
that, aside from exponential divergences associated with the damping time, $\tau$, the
perturbation in the curvature invariants has the acceptable form

$$\delta R_I \sim 1/r \;.$$

Edelstein (1967), in studying non-radial perturbations of the Schwarzschild geometry,
has independently verified that the divergence of order $r$ in $\delta g_{\mu\nu}/g_{\mu\nu}$ is non-physical. His
approach was to exhibit explicitly a coordinate transformation (change of gauge)
which made the amplitude of the exterior eigenfunctions go to zero for large $r$. We
shall not have need of this change of gauge in the present paper, but it is comforting
to know of its existence.

Summarizing the results of § III: To each choice of the complex eigenfrequency,
$\omega_n$, there is a unique eigenfunction which is physically acceptable at both $r = 0$ and
the surface of the star. That eigenfunction has the form (18) far from the star, which
is always a physically acceptable form if one allows for the possibility of both incoming
and outgoing gravitational waves.
