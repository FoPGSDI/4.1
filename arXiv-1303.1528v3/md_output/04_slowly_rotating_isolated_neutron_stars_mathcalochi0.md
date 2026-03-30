# Slowly Rotating, Isolated Neutron Stars: $\mathcal{O}(\chi^0)$

In this section, we construct non-rotating, isolated NS solutions, which will later be used as background solutions to construct slowly-rotating, tidally-deformed NS solutions in Secs. [ref:sec:linear], [ref:sec:quadratic] and [ref:sec:tidal].

## Einstein Equations and Exterior Solutions

The $(t,t)$ and $(R,R)$ components of the Einstein Equations yield

$$\begin{aligned}

\frac{d M}{dR} &=& 4 \pi R^2 \rho\,, \\

\frac{d \nu}{dR} &=& 2\frac{4 \pi R^3 p + M}{R(R-2M)}\,,
\end{aligned}$$

respectively. Combining the $R$-component of the equation of motion $\nabla^\mu T_{\mu R}^\text{mat}=0$ and Eq. (Eq. rr-zeroth), one obtains the Tolman-Oppenheimer-Volkoff (TOV) equation:

$$
\frac{dp}{dR} = -\frac{(4\pi R^3 p + M) (\rho + p)}{R(R-2M)}\,.

$$

Equations (Eq. tt-zeroth), (Eq. rr-zeroth) and (Eq. TOV-zeroth) together with the equation of state $p=p(\rho)$ close the system of differential equations.

The exterior solutions to the above equations can be obtained by setting $\rho=0=p$. One finds [hartle1967]

$$
\nu^\mathrm{ext} (R) = -\lambda^\mathrm{ext} (R) = \ln \left(1 - \frac{2 M_*}{R} \right)\,.

$$

We use the superscripts “ext” to refer to exterior quantities.

## Interior Solutions

First, we solve Eqs. (Eq. tt-zeroth) and (Eq. TOV-zeroth) together with the equation of state with initial conditions

\allowdisplaybreaks
$$\begin{aligned}
\rho (r_\epsilon) &=& \rho_c + \mathcal{O}(r_\epsilon^2)\,, \\
p (r_\epsilon) &=& p_c + \mathcal{O}(r_\epsilon^2)\,, \\
M(r_\epsilon) &=& \frac{4\pi}{3} \rho_c r_\epsilon^3 + \mathcal{O}(r_\epsilon^5)\,,
\end{aligned}$$

where $\rho_c$ and $p_c$ are the central density and pressure respectively, and $r_\epsilon$ corresponds to the core radius which we take to be $r_\epsilon = 100 \; {\rm{cm}} \ll \mathcal{R}_*$. We have checked that all of our results are independent of the choice of $r_{\epsilon}$ provided this is a very small number relative to the NS radius. We solve Eqs. (Eq. tt-zeroth) and (Eq. TOV-zeroth)  outwards from $r=r_\epsilon$ until $p$ vanishes. The NS radius $\mathcal{R}_*$ and the NS mass $M_*$ are then defined by $p(\mathcal{R}_*)=0$ and $M_* = M(\mathcal{R}_*)$ respectively. For later convenience, we introduce the NS compactness

$$
C \equiv \frac{M_*}{\mathcal{R}_*}\,.
$$

Notice that The central pressure $p_c=p(\rho_{c})$ is determined from the EoS, once $\rho_c$ is chosen. The central density $\rho_c$ is then a free parameter that effectively determines the mass and radius of the NS.

With this solutions, we can then solve Eq. (Eq. rr-zeroth). One approach is to use the boundary condition [see Eq. (Eq. ext-nu-lambda)]

$$
e^{\nu(\mathcal{R}_*)} = 1- \frac{2 M_*}{\mathcal{R}_*}
$$

at the NS surface as an initial condition and then integrate inwards toward the core. Another approach is to use the fact that Eq. (Eq. rr-zeroth) is shift invariant, as done e.g. in [kent-CSNS]. All throughout this paper, numerical solutions to the initial value problem are obtained with an adaptive 4th-order Runge-Kutta method [gsl].

![](MR-new-PRD.eps)

![](profile-PRD.eps)

Figure [ref:fig:MR] shows the mass-radius relation for various EoSs. We have checked that the mass-radius relation for the SLy EoS agrees with that shown in [shibata-fitting]. As anticipated, all EoSs lead to NSs with maximum mass larger than $1.93 \; M_\odot$ (the black dashed horizontal line), which is the lower bound for the mass of J1614-2230 [1.97NS]. We do not show the mass-radius relation for the polytropic EoSs because the I-Love-Q relations that we present in Sec. [ref:sec:I-Love-Q] only depend on the NS compactness and do not depend on the mass-radius relation. Figure [ref:fig:profile] shows the interior profile of the NS density (top) and pressure (bottom) as functions of the radial coordinate for a compactness of $C=0.17$, which corresponds to a NS with $M_* = 1.4 M_\odot$ and $\mathcal{R}_* \approx 12.1 \; {\rm{km}}$ for the APR EoS.