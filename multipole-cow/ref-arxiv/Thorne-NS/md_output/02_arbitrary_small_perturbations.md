# ARBITRARY, SMALL PERTURBATIONS

## (a) The Equilibrium Configuration

Throughout this paper we shall use the conventions and notation of Thorne (1966, 1967), "Relativistic Stellar Structure and Dynamics"; cited henceforth as RSSD) except that we shall adopt geometrized units ($c = G = 1$). We shall omit the asterisks ($*$) used in RSSD for geometrized quantities, and we shall use the symbol $\nu$ in place of $2\Phi$ for $\ln g_{tt}$.

As is discussed in RSSD, a relativistic equilibrium configuration can be described by a coordinate system $(t, r, \theta, \phi)$, with respect to which the geometry of spacetime is given by the line element

$$ds^{2} = -(e^{\nu})\,(dt)^{2} + e^{\lambda}\,(dr)^{2} + r^{2}(d\theta^{2} + \sin^{2}\theta\,d\phi^{2})\;.$$ (1)

Here $\nu$ and $\lambda$ are functions of the radius $r$, and $m(r)$ is defined as the "interior mass" $m(r)$, by

$$e^{-\lambda} = 1 - 2m/r\;.$$ (2)

In order to specify the hydrodynamic structure of an equilibrium configuration, one generally gives—in addition to the "gravitational potential," $\nu(r)$, and the mass inside radius $r$, $m(r)$—also the total density of mass-energy, $\rho$, and the pressure, $p$, as functions of $r$. The quantities $r$, $m$, $\rho$, and $p$ are related to each other by the mass equation,

$$m = \int_0^r 4\pi r^2 \rho \, dr \;,$$ (3a)

the Tolman-Oppenheimer-Volkoff equation of hydrostatic equilibrium,

$$\frac{dp}{dr} = -\frac{(\rho + p)(m + 4\pi r^3 p)}{r(r - 2m)} \;,$$ (3b)

the source equation for $\nu$,

$$d\nu/dr = -2(\rho + p)^{-1}(dp/dr) \;,$$ (3c)

and an equation of state,

$$p = p(\rho, Z_1, \ldots, Z_N) \;.$$ (3d)

Here $s$ is the entropy per baryon and $Z_1, \ldots, Z_N$ are the fractional abundances of the various nuclear species. The adiabatic index, $\gamma \equiv \Gamma_1$ in the notation of RSSD—which governs the response of the stellar material to adiabatic compressions is related to the equation of state by

$$\gamma = \bigl[(\rho + p)/p\bigr](\partial p/\partial \rho)_{s,\, Z_1,\, \ldots,\, Z_N} \;.$$ (4)

In order to construct a relativistic stellar model, one solves equations (3a–3d), together with certain equations of thermal structure and with gas characteristic equations outlined in RSSD. (We do not write down here the remaining structure equations because they play no role in the theory of pulsation.) Throughout the remainder of this paper we assume that somebody has given us an equilibrium configuration in which the structure parameters $r$, $m$, $\rho$, and $p$ satisfy equations (3a–3d), and we examine the behavior of that configuration under non-radial perturbations.

## (b) Fluid Displacement and Perturbation in the Geometry of Spacetime

The small-amplitude motion of our perturbed configuration is described by the 3-vector displacement, $\xi(t, r, \theta, \phi)$, of the fluid with respect to the equilibrium coordinate system $(t, r, \theta, \phi)$. As a result of the fluid motion, the geometry of spacetime around and inside the equilibrium configuration is no longer described by the line element (1). Rather, the geometry fluctuates in a manner described by ten functions, $h_{\mu\nu} \equiv h_{\mu\nu}$, of $(t, r, \theta, \phi) = (x^0, x^1, x^2, x^3)$:

$$ds^2 = (g_{\mu\nu}^{(0)} + h_{\mu\nu})\,dx^\mu\,dx^\nu \;.$$ (5)

The entire theory of non-radial pulsations consists of the study of the "equations of motion" which govern the thirteen functions, $\xi_i$ and $h_{\mu\nu}$, of $(t, r, \theta, \phi)$.

In order to make the equations of motion as simple as possible, we analyze $\xi_i$ and $h_{\mu\nu}$ into vectorial and tensorial spherical harmonics in a manner first employed by Regge and Wheeler (1957). (See Appendix A for details.) Each spherical harmonic is characterized by the usual indices $l$ and $M$, and by a parity which can be $(-1)^l$ or $(-1)^{l+1}$. Group-theoretic considerations—or, alternatively, an analysis of the equations of motion—reveal that for small-amplitude motions there is no coupling between the various harmonics $(l, m)$. Consequently, we can confine our attention to the various harmonics individually.

For small-amplitude motion in a given spherical harmonic, the equations of motion are simplified considerably by making a particular choice of "gauge"—the Regge-Wheeler choice—for the gravitational field. Making such a choice of gauge corresponds to removing the coordinate arbitrariness in $(t, r, \theta, \phi)$, which arbitrariness is of the same size as the perturbation, $h_{\mu\nu}$, in the geometry of spacetime.

The equations of motion can be simplified still further by restricting one's attention to spherical harmonics with $M = 0$. Once these are understood, the spherical-harmonic motions ($l$, $M \neq 0$, $r$) can be obtained from ($l$, $M = 0$, $r$) by suitable rotations about the center of the star.

Having specialized to a particular spherical harmonic ($l$, $M$, $r$), having made the Regge-Wheeler choice of gauge for that spherical harmonic, and having restricted attention to $M = 0$, we obtain the following vastly simplified forms for the fluid displacement, $\xi$, and the metric perturbation, $h_{\mu\nu}$ (see Appendix A for derivation):

If $r = (-1)^{l+1}$—"odd parity"; "magnetic-type" perturbations—then the covariant components of the fluid-displacement vector take the form

$$\xi_{t} = \xi_{r} = 0 \;, \qquad \xi_{\theta} = U(r,t)\,\sin\theta\,\partial_{\phi} P_{l}(\cos\theta) \;; \tag{6a}$$

and the only non-vanishing components of the metric perturbation are

$$\begin{aligned}
h_{t\phi} &= h_{\phi t} = h_{0}(r,t)\,\sin\theta\,\partial_{\theta} P_{l}(\cos\theta) \;,\\
h_{r\phi} &= h_{\phi r} = h_{1}(r,t)\,\sin\theta\,\partial_{\theta} P_{l}(\cos\theta) \;.
\end{aligned} \tag{6b}$$

If $r = (-1)^{l}$—"even parity"; "electric-type" perturbations—then the contravariant components of the fluid displacement vector take the form

$$\xi^{r} = r^{-2} e^{-\lambda/2} V(r,t)\, P_{l}(\cos\theta) \;, \qquad \xi^{\theta} = -r^{-2} W(r,t)\,\partial_{\theta} P_{l}(\cos\theta) \;, \qquad \xi^{\phi} = 0 \;; \tag{7a}$$

and the metric perturbation takes the form

$$h_{\mu\nu} =
\left(
\begin{array}{c|cccc}
   & t & r & \theta & \phi \\
  \hline
  t & H_{0} e^{\nu} & H_{1} & 0 & 0 \\
  r & H_{1} & H_{2} e^{\lambda} & 0 & 0 \\
  \theta & 0 & 0 & r^{2} K & 0 \\
  \phi & 0 & 0 & 0 & r^{2} \sin^{2}\!\theta\, K
\end{array}
\right)
P_{l}(\cos\theta) \;. \tag{7b}$$

Here $V$, $W$, $H_{0}$, $H_{1}$, $H_{2}$, and $K$ are functions of $r$ and $t$.

For odd-parity perturbations the equations of motion are a set of coupled equations for the fluid-displacement function $U(r,t)$ and the metric perturbation functions $h_{0}(r,t)$, $h_{1}(r,t)$ of equations (6). For even-parity perturbations the equations of motion are a set of coupled equations for the fluid-displacement functions $V(r,t)$, $W(r,t)$ and the metric perturbation functions $H_{0}(r,t)$, $H_{1}(r,t)$, $H_{2}(r,t)$, $K(r,t)$ of equations (7). The odd and even cases are considered separately in the next two sections.

## (c) Odd-Parity Motions

In Appendix B the Einstein field equations are calculated and discussed for the odd-parity case. From that discussion one learns that the odd-parity motions of a star are not characterized by pulsations which emit gravitational waves; rather, they are characterized by a stationary, differential rotation of the fluid inside the star and by gravitational waves which do not couple to the star at all.

This result should not be surprising. Pulsations can occur only if the perturbation causes a change in the star's internal density and pressure, $\rho$ and $p$. However, $\rho$ and $p$ are scalar fields; and all *scalar* spherical harmonics are of even parity. (Only vector and tensor spherical harmonics can have odd parity.) Therefore, a perturbation of odd tensor spherical harmonics cannot change the star's density or pressure distributions and therefore cannot cause the star to pulsate.

This explanation for the non-existence of odd-parity pulsations can be restated in more physical terms as follows: Odd-parity gravitational waves are "transverse" waves in the sense that they only couple to stars which can support anisotropic stresses. (Anisotropic stresses are characterized by tensors, whose spherical harmonics can have odd parity.) Consequently, by idealizing our equilibrium configurations as made of a perfect fluid (purely isotropic stress, $-p$), we rule out, ab initio, the possibility of odd-parity pulsations and of the generation of odd-parity gravitational waves.

Because odd-parity pulsations are non-existent, we shall confine our analysis to even-parity motions throughout the remainder of this paper.

## (d) Even-Parity Motions

The equations of motion which govern the time evolution of even-parity motions are derived from Einstein's field equations in Appendix C. For a given harmonic $(l,\,M=0,\,s=(-1)^{l})$ these equations of motion are a set of coupled equations in the fluid displacement functions $V(r,t)$, $W(r,t)$ and for the metric perturbation functions $H_0(r,t)$, $H_1(r,t)$, $H_2(r,t)$, $K(r,t)$ (cf. eqs. [3]). Initial-value equations, which determine the time evolution of $K$, $W$, $V$ from initial data, have been specified; and propagation equations, which determine the time evolution of $K$, $W$, $V$. For $l\geq 2$ the *initial-value equations* are (cf. Appendix C)

$$\begin{aligned}
H_0' + r^{-1}e^{\lambda}\bigl[l(l+1)/2 + 1 + 4\pi r^2(\rho-p)\bigr]H_0 \\
\quad = rK'' + e^{\lambda}(3-5m/r-4\pi r^2 p)K'
  - r^{-1}e^{\lambda}\bigl[l(l+1)/2 - 1 \\
\qquad + 8\pi r^2(\rho+p)\bigr]K
  + 8\pi r^{-1}(\rho+p)e^{\lambda/2}W' \\
\qquad + 8\pi r^{-1}\rho'\,e^{\lambda/2}W
  + 8\pi l(l+1)r^{-1}(\rho+p)e^{\nu}V \;,
\end{aligned} \tag{8a}$$

$$\begin{aligned}
l(l+1)H_1 = 2r^2 K'_{,t} - 2re^{\lambda}(-1+3m/r+4\pi r^2 p)K_{,t}
  - 2rH_{0,t} \\
  + 16\pi(\rho+p)e^{\lambda/2}W_{,t} \;,
\end{aligned} \tag{8b}$$

$$H_{1,t} = e^{\nu}H_0' + r^{-1}e^{\lambda+\nu}(2m/r+8\pi r^2 p)H_0 - e^{\nu}K' \;, \tag{8c}$$

$$H_2 = H_0 \;, \tag{8d}$$

and the time derivatives of equations (8a) and (8d). For $l\geq 2$ the *propagation equations* are

$$\begin{aligned}
K_{,tt} - e^{-\lambda}K'' + 2r^{-1}e^{-\lambda}\bigl[-1+m/r+2\pi r^2(\rho-p)\bigr]K' \\
+ r^{-2}e^{\nu}\bigl[l(l+1)-2+8\pi r^{-2}(\rho+p-\gamma p)\bigr]K
  + r^{-2}e^{\nu}\bigl[2e^{-\lambda}-4\pi r^2(\rho+p+\gamma p)\bigr]H_0 \\
- 8\pi r^{-2}(\rho+p-\gamma p)e^{\nu-\lambda/2}W'
  - 8\pi r^{-2}(\rho'-p')e^{\nu-\lambda/2}W \\
- 8\pi l(l+1)r^{-2}(\rho+p-\gamma p)e^{\nu}V = 0 \;,
\end{aligned} \tag{9a}$$

$$\begin{aligned}
r^{-2}(\rho+p)e^{(\lambda-\nu)/2}W_{,tt}
  = \bigl\{\gamma p e^{\nu/2}\bigl[r^{-2}e^{-\lambda/2}W'
    + l(l+1)r^{-2}V - \tfrac{1}{2}H_0 - K\bigr]\bigr\}' \\
+ \tfrac{1}{2}(\rho+p)e^{\nu/2}(r^{-2}e^{-\lambda/2}\nu')'W
  - l(l+1)r^{-2}(\rho+p)(e^{\nu/2})'V \\
- \tfrac{1}{2}(\rho+p)e^{-\nu}(H_0 e^{3\nu/2})'
  + (\rho+p)(Ke^{\nu/2})' = 0 \;,
\end{aligned} \tag{9b}$$

$$\begin{aligned}
e^{-\nu}(\rho+p)V_{,tt} + l(l+1)r^{-2}\gamma p V
  + r^{-2}\gamma p\,e^{-\lambda/2}W' + r^{-2}p'\,e^{-\lambda/2}W \\
- \tfrac{1}{2}(p+\rho+\gamma p)H_0 - \gamma p K = 0 \;.
\end{aligned} \tag{9c}$$

Here and throughout this paper primes denote radial derivatives, and commas followed by $t$'s denote time derivatives: $X' \equiv \partial X/\partial r$, $X'' \equiv \partial^{2}X/\partial r^{2}$, $X_{,t} \equiv \partial X/\partial t$, and $X_{,tt} \equiv \partial^{2}X/\partial t^{2}$.

In deriving equations (8) and (9) we have assumed at several points that $l \neq 0$ and $l \neq 1$; and *we shall maintain this assumption throughout the remainder of the paper*. This assumption is not a serious limitation on our analysis, since we are interested primarily in the emission of gravitational waves by pulsating stars, and gravitational waves must have $l \geq 2$.

The structure of the equations of motion (8) and (9) reflects the fact that even-parity motions have three degrees of freedom—two, $W$ and $V$, associated with the motion of the fluid; and one, $K$, associated with the gravitational waves.[^1]

[^1]: If one specifies at an initial moment of coordinate time the radial distributions of $W$, $V$, $K$ and $W_{,t}$, $V_{,t}$, $K_{,t}$, then one can use the three propagation equations (9)—plus the initial value equation (8a) with the boundary condition

    $$H_0 = K \quad \text{at } r = 0 \tag{10}$$

    (cf. eq. [15a])—to propagate $W$, $V$, $K$ and $K_{,t}$. At any moment of coordinate time the functions $H_0$, $H_1$, and $H_2$, which are needed to fix completely the perturbed geometry of spacetime, are determined from $W$, $V$, $K$, $V_{,t}$, $K_{,t}$ by the initial-value equations (8).

    The dynamical evolution problem (eqs. [9]) would be much simpler if one could solve the initial-value equations (8a) for $H_2$ as a linear combination of $K'$, $K$, $K_{,t}$, $W'$, $W$, $V''$, $V'$, $V$, and then use that solution, plus equation (8a) itself, to eliminate $H_0$ and $H_2$ from the dynamical equations (9). Unfortunately, equations (8a) cannot be so inverted and, therefore, equation (8a) is needed along with equations (9) to fix the dynamical propagation of the three degrees of freedom. The only way one can incorporate equation (8a) into the propagation equations (9) directly is by increasing equations (9) from second-order partial differential equations to third-order partial differential equations. We prefer not to do this.
