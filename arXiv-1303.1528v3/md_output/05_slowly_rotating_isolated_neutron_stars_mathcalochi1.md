# Slowly Rotating, Isolated Neutron Stars: $\mathcal{O}(\chi^1)$

Let us now focus on constructing slowly-rotating, isolated NS solutions. In this section, we only consider axisymmetric perturbations and construct NSs to linear order in spin. We will first discuss the differential equation that needs to be solved, and then we will solve them in the exterior region modulo an integration constant. After this, we discuss the asymptotic behavior of the solution at the NS center, which can then be used as an initial condition to solve the equations in the interior region. Finally, we determine the integration constant by matching the interior and exterior solutions at the NS surface.

## Einstein Equations and Exterior Solutions

At linear order in $\epsilon$, the only non-vanishing component of the Einstein Equations is the $(t,\phi)$ one:

$$
\frac{d^2 \omega_1}{d R^2} + 4\frac{1- \pi R^2 (\rho +p) e^{\lambda}}{R}\frac{d \omega_1}{d R}  -16 \pi (\rho +p) e^{\lambda} \omega_1=0\,.

$$

Solving this equation in the exterior (i.e. setting $p = 0 = \rho$), one finds [hartle1967]

$$
\omega_1^\mathrm{ext} = \Omega_* - \frac{2S}{R^3} = \Omega_* \left( 1- \frac{2I}{R^3} \right)\,,

$$

where we have defined the moment of inertia by

$$
I \equiv \frac{S}{\Omega_*}\,.

$$

This quantity characterizes how fast a body can spin given a fixed spin angular momentum $S$. Notice that the exterior solution depends on two constants $\Omega_{*}$ and $S$. The former must be specified {*a priori*}, just like $\rho_{c}$, and it describes how fast the NS is rotating. The latter is determined by matching this exterior solution to an interior solution at the NS surface.

## Interior Solutions

Before we can solve for the interior solution, we first need initial conditions at the NS center. Taylor-expanding Eq. (Eq. omega1RR) about the NS center, we find that the interior solution must asymptotically behave as

$$

\omega_1(R) = \omega_c + \frac{8\pi}{5} (\rho_c + p_c) \omega_c R^2 + \mathcal{O}(R^3) \quad (R \to 0^+)\,.
$$

This solution contains a single constant, $\omega_c$, because we have eliminated another constant by requiring regularity of the solution at the NS center. The constant $\omega_{c}$ determines the NS spin angular momentum $S$, or equivalently, the NS moment of inertia $I$; in particular, $I$ increases as $\omega_c$ increases.

![](IbarM-PRD.eps)

We numerically solve Eq. (Eq. omega1RR) with the initial condition in Eq. (Eq. omega10) via an adaptive 4th-order Runge-Kutta method [gsl]. In solving this equation, one can take advantage of its homogeneity, its scale-invariance, as done e.g. in [yunespsaltis,kent-CSNS]. Once the interior solution has been found, we match it to the exterior one in Eq. (Eq. omega-ext) at the NS surface $R=\mathcal{R}_*$. The matching ensures that the solution is continuous and differentiable at the NS surface:

$$
\omega_1^{\mathrm{int}}(\mathcal{R}_*) = \omega_1^{\mathrm{ext}}(\mathcal{R}_*), \quad  \omega_1’{}^{\mathrm{int}}(\mathcal{R}_*) = \omega_1’{}^{\mathrm{ext}}(\mathcal{R}_*)\,,

$$

where we use the superscript “int” to refer to interior quantities. Through these conditions, we determine $S$ (or equivalently $I$)  and $\omega_c$ as a function of $\Omega_{*}$. In practice, due to the scale invariance of Eq. (Eq. omega1RR), the exterior solution can be divided by $\Omega_{*}$ and thus it only depends on the single constant $I$. Similarly, the interior solution can be obtained for $\omega_{1}^{\mathrm{int}}/\Omega_{*}$ as a function of a single constant $\bar{\omega}_{c} = \omega_{c}/\Omega_{*}$. Therefore, the conditions in Eq. (Eq. BC-omega) uniquely determine $I$ and $\omega_{c}$. This then determines the full solution, and thus also $S$, up to the overall constant of proportionality $\Omega_{*}$.

The moment of inertia can be expressed entirely as a function of the interior solution. From Eqs. (Eq. tt-zeroth)--(Eq. omega1RR), (Eq. omega-ext) and (Eq. BC-omega), $I$ takes the form [hartle1967,kalogera-psaltis]

$$
I = \frac{8\pi}{3} \frac{1}{\Omega_*} \int_0^{\mathcal{R}_*} \frac{e^{-(\nu^\mathrm{int} + \lambda^\mathrm{int} )/2} R^5 (\rho + p) \omega_1^{\mathrm{int}}}{R-2M(R)} dR\,.

$$

In the Newtonian limit (superscript “N”), Eq. (Eq. I) reduces to [hartle1967]

$$
I^\text{N} = \frac{8\pi}{3} \int^{\mathcal{R}_*}_0 R^4 \rho (R)dR\,,

$$

For later convenience, we define the dimensionless moment of inertia $\bar{I}$

$$

\bar{I} \equiv \frac{I}{M_*^3}\,.
$$

Figure [ref:fig:IbarMC] shows $\bar{I}$ as a function of the NS mass $M_*$ and compactness $C$. We have verified that the moment of inertia $I$ obtained here agrees exactly with previous results in the literature [pani-NS-EDGB]. Observe that the different $\bar{I}$ curves for realistic EoSs approach each other as $C$ increases. Moreover, observe that all these curves approach the value of $\bar{I}$ for a BH as $C \to 0.5$, shown with a solid cross in Fig. [ref:fig:IbarMC]. Of course, none of the NS sequences considered here will ever lead to a BH solution for any finite choice of central density.