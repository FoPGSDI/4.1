# Slowly Rotating, Isolated Neutron Stars: $\mathcal{O}(\chi^2)$

Let us now look at slowly-rotating NS solutions at quadratic order in spin. Following Sec. [ref:sec:linear], we first discuss the differential equations that describe the solution and then we solve them in the exterior region. We then discuss the asymptotic behaviors of the solutions at the NS center, obtain the interior solutions numerically, and match it to the exterior solution at the NS surface.

## Einstein Equations and Exterior Solutions

At quadratic order in spin, the $\theta$-component of the equation of motion $\nabla^{\mu} T_{\mu\theta}^\text{mat}=0$, valid only inside the star, yields

$$
\xi_2 = -\frac{R^2 e^{-\lambda} (3 h_2+e^{-\nu} R^2 \omega_1^2)}{3(M+4\pi p R^3)}\,.

$$

The $(\theta,\theta)-(\phi,\phi)$, $(R,\theta)$ and $(R,R)$ components of the Einstein Equations give respectively,

$$\begin{aligned}

m_2 &=& -R e^{-\lambda} h_2  +\frac{1}{6} R^4 e^{-(\nu+\lambda)} \left[ R e^{-\lambda} \left(\frac{d \omega_1}{dR} \right)^2
+ 16 \pi R \omega_1^2 (\rho +p) \right]\,, \\

\frac{dK_2}{dR} &=& -\frac{dh_2}{dR} + \frac{R-3M-4\pi pR^3}{R^2} e^{\lambda} h_2
+ \frac{R-M+4\pi p R^3}{R^3} e^{2\lambda} m_2\,, \\

\frac{dh_2}{dR} &=& -\frac{R-M+4 \pi p R^3}{R}e^{\lambda} \frac{dK_2}{dR} +\frac{3-4\pi (\rho +p) R^2}{R} e^{\lambda} h_2 +\frac{2}{R} e^{\lambda} K_2 +\frac{1+8\pi p R^2}{R^2} e^{2\lambda} m_2 +\frac{R^3}{12}e^{-\nu} \left( \frac{d\omega_1}{dR} \right)^2  \\
& & -\frac{4 \pi (\rho +p) R^4 \omega_1^2}{3R}e^{-\nu+\lambda}\,.
\end{aligned}$$

By imposing asymptotic flatness at spatial infinity, one finds the exterior solutions  [hartle1967]

$$\begin{aligned}

h_2^\mathrm{ext} &=& \frac{1}{M_* R^3} \left( 1 + \frac{M_*}{R} \right) S^2  + AQ_2^2 \left( \frac{R}{M_*}-1 \right)  \\
             &=&  \frac{1}{M_* R^3} \left( 1 + \frac{M_*}{R} \right) S^2 -\frac{3A R^2}{M_* (R-2 M_*)} \left[ 1- 3 \frac{M_*}{R} + \frac{4}{3} \frac{M_*^2}{R^2} + \frac{2}{3} \frac{M_*^3}{R^3} + \frac{R}{2 M_*} f(R)^2 \ln f(R) \right]\,, \\

K_2^\mathrm{ext} &=& - \frac{1}{M_* R^3} \left( 1+\frac{2M_*}{R} \right) S^2 + \frac{2AM_*}{\sqrt{R (R-2 M_*)}}Q_2^1\left( \frac{R}{M_*}-1 \right) - AQ_2^2 \left( \frac{R}{M_*}-1 \right)  \\
              &=& - \frac{1}{M_* R^3} \left( 1+\frac{2M_*}{R} \right) S^2 + \frac{3AR}{M_*} \left[1+ \frac{M_*}{R} - \frac{2}{3} \frac{M_*^2}{R^2} + \frac{R}{2M_*} \left( 1- \frac{2 M_*^2}{R^2} \right) \ln f(R)  \right]\,, \\

m_2^\mathrm{ext} &=&  - \frac{1}{M_* R^2} \left( 1 -7\frac{M_*}{R} + 10 \frac{M_*^2}{R^2} \right)S^2 +\frac{3 A R^2}{M_*} \left[ 1 - 3 \frac{M_*}{R} + \frac{4}{3} \frac{M_*^2}{R^2} + \frac{2}{3} \frac{M_*^3}{R^3} + \frac{R}{2 M_*} f(R)^2 \ln f(R)  \right]\,,
\end{aligned}$$

with $f(R) \equiv 1- 2M_*/R$, $Q_2^2$ and $Q_2^1$ the associated Legendre functions of the second kind and $A$ an integration constant that is to be determined by matching with the interior solution at the NS surface.

The spin-induced quadrupole moment $Q^\mathrm{(rot)}$ can be read off from the coefficient of the $P_2(\cos\theta)/R^3$ term in the Newtonian potential [hartlethorne]:

$$
Q^\mathrm{(rot)} = - \frac{S^2}{M_*} - \frac{8}{5} A M_*^3\,.

$$

Notice that the quadrupole moment depends both on the magnitude of the spin angular momentum $S$ and the integration constant $A$, determined after matching the interior and exterior linear- and quadratic-order in spin solutions at the NS surface. The quadrupole moment represents the quadrupolar deformation of a body away from sphericity, with $Q^\mathrm{(rot)}<0$ corresponding to an oblate deformation. Notice also that the first term of Eq. (Eq. quadrupole) is identical to the relation one obtains for BH, which means that $A\to0$ in the GR test-particle limit.

![](QbarM-PRD.eps)

## Interior Solutions

Let us begin by Taylor-expanding Eqs. (Eq. xi2)--(Eq. h2R) about the NS center and solving the expanded equations to obtain

\allowdisplaybreaks
$$\begin{aligned}

h_2 (R) &=&  B R^2 + \mathcal{O}(R^4)\,, \\

K_2(R) &=& - B R^2 + \mathcal{O}(R^4)\,, \\

m_2(R) &=& - B R^3 + \mathcal{O}(R^5)\,, \\

\xi_2(R) &=& - \frac{3 B + e^{-\nu_c} \omega_c^2}{4\pi (\rho_c +3p_c)}  R + \mathcal{O}(R^3), \quad (R \to 0^+)\,,  \\
\end{aligned}$$

where $B$ is a constant that determines the NS quadrupole moment. As before, the constant $\nu_c$ is defined as $\nu_c \equiv \nu(r_\epsilon)$.

We numerically solve the evolution Eqs. (Eq. k2R) and (Eq. h2R) with the initial conditions of Eqs. (Eq. h20) and (Eq. k20), using an adaptive 4th-order Runge-Kutta algorithm [gsl]. As before, when solving these equations we must impose the following boundary conditions, such that $h_2$ and $K_2$ are continuous at the NS surface:

$$
h_2^\mathrm{int} (\mathcal{R}_*) = h_2^\mathrm{ext} (\mathcal{R}_*), \quad K_2^\mathrm{int} (\mathcal{R}_*) = K_2^\mathrm{ext} (\mathcal{R}_*)\,.
$$

These matching conditions determine the constants $A$ in Eqs. (Eq. h2-ext) and (Eq. k2-ext) and $B$ in Eqs. (Eq. h20) and (Eq. k20).

In practice, we follow [hartle1967,kent-CSNS] and first solve the interior solution as a sum of a particular solution, with some test-value for $B$, and the product of an undetermined constant and the homogeneous solution. We then fix this undetermined constant, together with $A$, by requiring that the interior and exterior solutions match at the NS surface. We have checked the results obtained through this method by solving the equations using the Riccati method [dieci1,dieci2,takata].

Figure [ref:fig:QbarMC] shows the dimensionless rotationally-induced quadrupole moment $\bar{Q}$ as functions of $M_*$ and $C$, where $\bar{Q}$ is defined by

$$

\bar{Q} \equiv -\frac{Q^\mathrm{(rot)}}{M_*^3 \chi^2}\,,
$$

where we recall that the dimensionless spin parameter $\chi$ is defined by $\chi \equiv S/M_*^2$. This $\bar{Q}$ is the same as the dimensionless quadrupole moment $a$ in [poisson-quadrupole]. As in the $\bar{I}$ case, the $\bar{Q}$ curves for realistic EoSs approach each other as $C$ increases. Moreover, these curves also approach the $\bar{Q}$ value for a BH as the compactness approaches $0.5$. As before, however, the NS sequence does not go to a BH solution for any finite choice of central density.

## Rotational Love Number

With the quadratic isolated NS solutions at hand, we can now introduce the rotational Love number [mora-will]. In general, Love numbers represent the deformability of a NS away from sphericity. The rotational Love number, in particular, refers to the deformability of a NS due to its spin.

Love numbers are defined in a {*buffer zone*}, the region ${\cal{R}} \gg R \gg \mathcal{R}_*$, where ${\cal{R}}$ is the radius of curvature of the source of the perturbation. For example, the $(t,t)$ component of the metric can be expanded in the buffer zone as [mora-will,Yunes:2005nn,Yunes:2006iw,JohnsonMcDaniel:2009dq,hinderer-love,Chatziioannou:2012gq]

$$\begin{aligned}

\frac{1-g_{tt}}{2} &=& - \frac{M_*}{R} - \frac{4\pi}{5} \frac{Q^\mathrm{(rot)}}{R^3} \sum_m Y_{2m} ( \hat{\Omega} ) Y_{2m}^* ( \hat{n} ) + \mathcal{O} \left( \frac{\mathcal{R}_*^{4}}{R^4} \right)  \\
&& + \frac{4\pi}{15} \mathcal{E}^\mathrm{(rot)} R^2 \sum_m Y_{2m} ( \hat{\Omega} ) Y_{2m}^* (\hat{n} )  + \mathcal{O}\left(\frac{R^{3}}{{\cal{R}}^3}\right)  \\
&=& - \frac{M_*}{R} - \frac{Q^\mathrm{(rot)}}{R^3} P_2 (\hat{\Omega} \cdot \hat{n}) + \mathcal{O} \left( \frac{\mathcal{R}_*^{4}}{R^4} \right)  \\
&& + \frac{1}{3} \mathcal{E}^\mathrm{(rot)} R^2 P_2 (\hat{\Omega} \cdot \hat{n})  + \mathcal{O}\left(\frac{R^{3}}{{\cal{R}}^3}\right)\,.
\end{aligned}$$

The quantity $\mathcal{E}^\mathrm{(rot)}$ is related to the trace of the rotationally-induced, electric, quadrupole tidal tensor, i.e. the quadrupolar contribution of the centrifugal potential. In the Newtonian limit, this quantity reduces to $\mathcal{E}^\mathrm{(rot)} = \Omega_*^2$ [mora-will]. As usual, $Y_{2m} ( \hat{\Omega} )$ are the $\ell =2$ spherical harmonics in the $\hat{\Omega}$ direction, where $\hat{n}$ is the principal axis of the perturbation, which in this case corresponds to the unit vector of the spin angular momentum $\hat{S}$.

The $\ell=2$ rotational Love number $\lambda^\mathrm{(rot)}$ is defined by [mora-will,berti-iyer-will]

$$
\lambda^\mathrm{(rot)} \equiv - \frac{Q^\mathrm{(rot)}}{\mathcal{E}^\mathrm{(rot)}} = -\frac{Q^\mathrm{(rot)}}{\Omega_{*}^2}\,,
$$

where the second equality uses the Newtonian expression for ${\cal{E}}^{\mathrm{(rot)}}$. As defined here, $\lambda^\mathrm{(rot)}$ has unit of (mass)$^5$ or (length)$^5$ (recall that we use geometric units throughout this paper, where $c=1=G$), and thus, there are 2 natural ways of normalizing it [mora-will,berti-iyer-will];

$$\begin{aligned}
k_2^\mathrm{(rot)} & \equiv & \frac{3}{2} \frac{\lambda^\mathrm{(rot)}}{\mathcal{R}_*^5}\,,  \\
\bar{\lambda}^\mathrm{(rot)} &\equiv& \frac{\lambda^\mathrm{(rot)}}{M_*^5} = \frac{2}{3} k_2^\mathrm{(rot)} C^{-5}\,.
\end{aligned}$$

By using Eqs. (Eq. I), (Eq. bar-I) and (Eq. bar-Q), one can rewrite $\bar{\lambda}^\mathrm{(rot)}$ as

$$
\bar{\lambda}^\mathrm{(rot)} = \bar{I}^2 \bar{Q}\,.
$$

In this paper, we refer to $k_2^\mathrm{(rot)}$ as the $\ell=2$ rotational apsidal constant, while we refer to $\bar{\lambda}^\mathrm{(rot)}$ as the $\ell=2$ dimensionless rotational Love number.