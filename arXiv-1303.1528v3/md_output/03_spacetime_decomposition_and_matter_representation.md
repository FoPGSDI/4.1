# Spacetime Decomposition and Matter Representation

In this paper, we consider uniformly rotating NSs that are slightly deformed either due to rotation or tidal fields. Such solutions can be numerically constructed perturbatively in a slow-rotation and tidal-deformation expansion, taking the non-rotating, isolated solution as a background. In this section, we explain the metric decomposition employed here and the stress-energy tensor we will use to describe NSs.

## Metric Decomposition

We choose Boyer-Lindquist type coordinates $(t,r,\theta,\phi)$ and decompose the metric as

$$\begin{aligned}
ds^2 &=& -e^{\bar{\nu}(r)} \left[1+2 \epsilon^2 \bar{h}_2 (r) \alpha Y_{2 m}(\theta,\phi) \right] dt^2  \\
& & + e^{\bar{\lambda}(r)} \left[ 1+\frac{2 \epsilon^2 \bar{m}_2 (r) \alpha Y_{2 m}(\theta,\phi)}{r-2\bar{M}(r)} \right] dr^2  \\
& & + r^2 \left[ 1+2 \epsilon^2 \bar{K}_2(r) \alpha Y_{2 m}(\theta,\phi) \right]  \\
& & \times \left\{ d \theta^2 + \sin^2 \theta \left[ d\phi - \epsilon [\Omega_* - \bar{\omega}_1 (r) P_1’ (\cos \theta) ] dt \right]^2 \right\}  \\
& & + \mathcal{O}(\epsilon^3)\,,

\end{aligned}$$

where $\bar{M}(r)$ is defined by

$$
\bar{M}(r) \equiv \frac{\left[ 1-e^{-\bar{\lambda}(r)} \right]r}{2}\,,
$$

$P_\ell (\cos \theta)$ is the $\ell$-th order Legendre polynomial, $P_{1}’ = d P_{1}/d(\cos{\theta})$ and $Y_{\ell m}(\theta, \phi)$ is the spherical harmonic function. The quantity $\epsilon$ here is a book-keeping parameter that we will later set to unity and we only introduce to remind ourselves of the order of the approximation. Terms linear in $\epsilon$ are induced only by linear-order in rotation effects, while tidal-deformation effects enter at ${\cal{O}}(\epsilon^{2})$. We will work here to quadratic order in $\epsilon$.

A slow-rotation expansion is quite appropriate to model old NSs. Recycled millisecond pulsars, the fastest NSs observed to date, have angular velocities in the kHz, but this number is small relative to the NS mass, i.e. $M_{*} \Omega_{*} \lesssim 0.01$, where $\Omega_*$ is the NS angular velocity. For the fastest millisecond pulsar J1939+2134 [fastest-pulsar], with period $1.5 \; {\rm{ms}}$, the dimensionless spin parameter, defined via $\chi \equiv S/M_{*}^{2} = I \Omega_{*}/M_{*}^{2}$, is still small $\chi \lesssim 0.3$, using a Newtonian expression for the moment of inertia. Thus, a slow-rotation expansion is well-justified, especially when carried out to second order. This approximation, however, would break down if considering newly-born NSs, which are likely to be differentially rotating, much hotter and with much larger magnetic fields. Notice also that the NSs that will source GWs in the band of ground-based detectors are expected to have significantly smaller spins than that. This is because NSs spin-down [bildsten-cutler] as they inspiral and ground-based detectors will only be sensitive to the last 17 minutes of the orbit before coalescence.

The free functions in our metric decomposition are $\bar{\nu}$ and $\bar{\lambda}$ at $\mathcal{O}(\epsilon^0)$, $\bar{\omega}_1$ at $\mathcal{O}(\epsilon)$ and $\bar{h}_2$, $\bar{K}_2$ and $\bar{m}_2$ at $\mathcal{O}(\epsilon^2)$. The leading-order correction due to slow rotation enters at $\mathcal{O}(\epsilon)$, while that due to tidal deformations enters at $\mathcal{O}(\epsilon^2)$. For the former, we restrict ourselves to axisymmetric perturbations; at $\mathcal{O}(\epsilon)$ only the $(\ell,m)=(1,0)$ mode survives, while at $\mathcal{O}(\epsilon^2)$ only the $(\ell,m)=(0,0)$ and $(\ell,m)=(2,0)$ modes survive. For the latter, we are only interested in the spin and tidal, {*quadrupolar*} deformations, and thus we only keep $\ell=2$ modes in Eq. (Eq. metric-ansatz-rth), but allow for all $m$ modes. Henceforth, we set the constant $\alpha=2 \sqrt{\pi/5}$ so that $\alpha Y_{20}(\theta, \phi) = P_2 (\cos \theta)$.

As pointed out by Hartle [hartle1967], one needs to be careful about choosing coordinates when deriving and solving perturbed equations. A perturbative analysis is valid only if perturbed quantities are much smaller than the unperturbed one. If one were to carry out calculations in $(t,r,\theta,\phi)$ coordinates, such conditions would be violated in certain situations. For example, in the region of spacetime outside the unperturbed star but inside the perturbed star, the ratio of the perturbed pressure (or density) to that of the unperturbed pressure (or density) diverges, which violates our perturbative treatment.

In order to overcome this problem, we transform the radial coordinate via [hartle1967]

$$
r(R,\theta) = R + \epsilon^2 \xi_2(R) \alpha Y_{2 m}(\theta, \phi) + {\cal{O}}(\epsilon^{3})\,,
$$

where $\xi_{2}(R)$ is such that

$$
\rho [r(R, \theta, \phi)]=\rho(R) = \rho^{(0)}(R)\,.
$$

In other word, the new radial coordinate $R$ is chosen such that $\rho [r(R, \theta,\phi)]$ is identical to the unperturbed density $\rho^{(0)}(r)$. By construction, the density and pressure in these new coordinates contain only the unperturbed contributions. Notice that $\xi_{2} Y_{2m}$ is well-defined only inside the star and we take it to be constant outside. This means that the exterior metric in $(t,r,\theta,\phi)$ coordinates can be obtained simply by replacing $R \to r$ in the exterior metric in $(t,R,\theta,\phi)$ coordinates.

The transformed metric in $(t,R,\theta,\phi)$ coordinates can be found in [kent-CSNS] for the axisymmetric case.  Henceforth, we will relabel the metric coefficients via

\allowdisplaybreaks
$$\begin{aligned}
\nu(R) &\equiv& \bar{\nu}(r)  = \bar{\nu}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
\lambda(R) &\equiv& \bar{\lambda}(r)  = \bar{\lambda}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
\omega_{1}(R) &\equiv& \bar{\omega}_{1}(r)  = \bar{\omega}_{1}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
h_{2}(R) &\equiv& \bar{h}_{2}(r)  = \bar{h}_{2}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
m_{2}(R) &\equiv& \bar{m}_{2}(r)  = \bar{m}_{2}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
K_{2}(R) &\equiv& \bar{K}_{2}(r)  = \bar{K}_{2}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,,
 \\
M(R) &\equiv& \bar{M}(r)  = \bar{M}(R + \epsilon^{2} \xi_{2} \alpha Y_{2m})\,.
\end{aligned}$$

## Matter Representation

We here consider NSs that are uniformly rotating, and thus, we model them as a perfect fluid. Uniform rotation should be a reasonable approximation unless one considers newly-born NSs. The stress energy-momentum tensor of the matter field $T_{\mu\nu}^\text{mat}$ is then given by

$$
T_{\mu\nu}^\text{mat} = (\rho + p ) u_\mu u_\nu + p \; g_{\mu\nu}\,,
$$

where the four-velocity $u^\mu$ is given by

$$
u^\mu = (u^0, 0,0, \epsilon \Omega_* u^0)\,,
$$

and $\Omega_*$ is the constant angular velocity of the NS. By using the normalization condition $u_\mu u^\mu = -1$, we obtain the time component of the four-velocity $u^0$ as

$$\begin{aligned}
u^0 &=& e^{-\nu/2} +  \epsilon^2\frac{e^{-3\nu/2}}{2} [ \omega_1^2 P_1’{}^2 R^2 \sin^2\theta  \\
& & -  e^{\nu} (2 h_2 +\nu’ \xi_2) \alpha Y_{2m} ] + \mathcal{O}(\epsilon^4)\,.
\end{aligned}$$

We here consider 4 realistic EoSs: APR [APR], SLy [SLy,shibata-fitting], Lattimer-Swesty with nuclear incompressibility of 220MeV (LS220) [LS, ott-EOS] and Shen [Shen1,Shen2,ott-EOS], the latter two with temperature of 0.1MeV and an electron fraction determined by the neutrino-less, beta-equilibrium condition. All of the EoS described above are “realistic” in that they allow NSs with masses larger than 1.93$M_\odot$, the lower bound of the recently found massive pulsar J1614-2230 [1.97NS]. For comparison purposes, we also consider polytropic EoSs, i.e. EoSs of the form

$$
p=K \rho^{1+1/n}\,,

$$

where $K$ is an amplitude constant and $n$ is the constant polytropic index. One can approximate the NS EoS with polytropes in the range $n \approx 0.5-1$ [flanagan-hinderer-love,lattimer-prakash-2001]. No single polytrope, however, is believed to be an accurate representation of a realistic EoS.

The APR EoS is constructed by using the variational chain summation methods, which is expected to include all leading many-body correlation effects. The APR EoS uses Hamiltonians that include a three-nucleon interaction, which predicts that a transition exists from NS matter to a phase with neutral pion condensation at a baryon number density of $\sim 0.2 \ \mathrm{fm^{-3}}$. The SLy EoS is calculated from a non-relativistic mean field theory approach, with a new set of Skyrme-type effective nucleon-nucleon interactions, suitable for describing very neutron rich matter. Unlike the APR EoS that describes only the NS’s liquid core, the SLy EoS is a “unified EoS” in the sense that it is supposed to describe also the NS crust. The LS220 EoS is constructed from a finite-temperature compressible liquid-droplet model with a Skyrme nuclear force. Such an EoS is derived within the single heavy nucleus approximation and the assumption of nuclear statistical equilibrium. The Shen EoS uses a relativistic mean-field theory model and assumes nuclear statistical equilibrium. Nuclear incompressibility of the Shen EoS occurs at 281MeV.