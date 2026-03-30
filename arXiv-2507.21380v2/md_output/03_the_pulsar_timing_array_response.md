# The pulsar timing array response

Here we give a brief overview of the observed response of a PTA to a gravitational wave source. For this paper, we will restrict out attention to purely monochromatic sources. The metric perturbation of such a source propagating from the $\hat{n}$-direction, with some frequency $f_0$, is given by
$$\begin{aligned}
    h_{ij}(t, x) = h \epsilon_{ij} e^{-2 \pi i f_0 (\hat{n} \cdot \mathbf{x} + t - t_0)}.
\end{aligned}$$
Note that here $\hat{n}$ points in the direction of the source, not in the direction of propagation. The polarization tensor, $\epsilon_{ij}$, is a normalized linear combination of the $+$ and $\times$ modes:
$$\begin{aligned}
    \epsilon^+ = \hat{\phi}_i \hat{\phi}_j - \hat{\theta}_i \hat{\theta}_j, \\
    \epsilon^\times = \hat{\theta}_i \hat{\phi}_j + \hat{\phi}_i \hat{\theta}_j,
\end{aligned}$$
where $\hat{\theta}$ and $\hat{\phi}$ are the spherical-coordinate basis vectors at $\hat{n}$. Alternatively, we can describe the polarization states in terms of the left- and right-circular polarized modes: $\epsilon^L_{ij} = (\epsilon^+_{ij} + i \epsilon^\times_{ij}) / \sqrt{2}$, $ \epsilon^R_{ij} = (\epsilon^+_{ij} - \epsilon^\times_{ij}) / \sqrt{2}$. This is the more natural basis for gravitational wave polarization, as astrophysical sources of gravitational radiation will primarily be circularly polarized. For a ring of test particles perpendicular to the direction of propagation, the left-handed polarization mode rotates the particles in a clockwise sense, as viewed from the origin.

The effect of the metric perturbation is to produce a relative delay in the arrival time of the pulse from a pulsar located at a position, $\mathbf{r}$, known as the timing residual. The timing residual is given by [Maiorano2021]
$$\begin{aligned}
    \tau(t; \mathbf{r}) &= \int df \tilde{\tau}(\mathbf{r}) e^{-2 \pi i f t}, \\
    \tilde{\tau}(f; \mathbf{r}) &= \frac{i \tilde{h}_{ij}(f) \hat{r}^i \hat{r}^j}{2 \pi f (1 + \hat{n} \cdot \hat{r})} \left( 1 - \mathcal{P}(f; \mathbf{r})\right),
\end{aligned}$$
where $\tilde{h}_{ij}(f)$ is the Fourier transform of the wave, and $\mathcal{P}(f; \mathbf{r}) = e^{2 \pi i f r (1 + \hat{n} \cdot \hat{r})}$ arises from the so-called “pulsar term" of the timing residual. Essentially, the total timing residual is the sum of the effect of the wave at Earth and at the pulsar. The pulsar sees the wave at a relative phase of $2 \pi i f r (1 + \hat{n} \cdot \hat{r})$ compared to the earth.

We will be interested in monochromatic waves for which $\tilde{h}(f)$ is a delta function. We will also assume that the distances in the PTA are poorly determined, so that we neglect the oscillatory pulsar term, $\mathcal{P}$. For PTAs with unknown pulsar distances, the pulsar term may be treated as a source of noise in the limit where $2 \pi f r$ is large. In this limit, $\mathcal{P}$ is highly oscillatory and any two-point correlations of the timing residuals due to the pulsar term effectively vanish.

Now, the Earth term only depends on the angular position $\hat{r}$ of the pulsar on the sky. Thus, we can compute the PTA response of a pulsar at position $\hat{r}$ on the sky for a circularly polarized wave propagating in the $\hat{n}$ direction:
$$\begin{aligned}
    &\tilde{\tau}^{R,L}(\hat{r}) = -\frac{h}{4 \pi f_0} \zeta^{R,L}(\hat{r};\hat{n},t_0), \\
    \begin{split}
    &\zeta^{R,L}(\hat{r};\hat{n},t_0) = e^{2 \pi i f_0 t_0} \,
    \times \\
    &\frac{\left(\cos\theta \sin \theta_{\rm gw} - \cos \theta_{\rm gw} \cos \left( \phi - \phi_{\rm gw} \right) \sin \theta \pm i \sin \theta \sin \left( \phi - \phi_{\rm gw}\right) \right)^2}{\cos\theta \cos \theta_{\rm gw} + \cos\left(\phi - \phi_{\rm gw} \right) \sin\theta \sin \theta_{\rm gw} - 1},
    \end{split}
\end{aligned}$$
where $(\theta, \phi)$ are the angular coordinates of the pulsar on the sky, and $(\theta_{\rm gw}, \phi_{\rm gw})$ are the angular coordinates of the gravitational wave source. We have also defined the dimensionless angular response, $\zeta^{R,L}(\hat{r};\hat{n},t_0)$, which is independent of the amplitude of the wave.

Now, let us consider a gravitational wave source with an arbitrary polarization. The polarization state is described by the angles, $(2\chi, 2\psi)$, which define the polarization ellipse shown in the diagram, Figure [ref:fig:Stokes_diagram]. The angle $\psi$ defines the direction of the linear polarization, and $\chi$ describes the degree of circular polarization. An alternative description of polarization is the Poincaré sphere, also shown in Figure [ref:fig:Stokes_diagram]. A point on the Poincaré sphere, defined by the unit vector $\hat{n}_{\rm p} = (\sin\theta_p \cos\phi_p, \sin\theta_p \sin\phi_p, \cos\theta_p )$, is related to the polarization angles by $\theta_p = \frac{\pi}{2} - 2 \chi$ and $\phi_p = 2\psi$. The axes of the Poincaré sphere correspond to the polarization states described by the Stokes parameters: $Q$, $U$, and $V$. Using the polarization angles, we can write the general PTA response for a polarized, monochromatic gravitational wave as:
$$\begin{aligned}
\begin{split}
    \zeta(\hat{r};\hat{n},\hat{n}_p, t_0) = &\sin(\theta_p/2) e^{-i \phi_p / 2} \cdot \zeta^L(\hat{r}; \hat{n}, t_0) \\
    &+ \cos(\theta_p / 2) e^{i \phi_p / 2} \cdot \zeta^R(\hat{r}; \hat{n}, t_0).
\end{split}
\end{aligned}$$
The timing residuals are obtained by $\tilde{\tau} = \zeta h / 2 \pi f_0$. From this, we can see that the gravitational wave source is fully described by its position $\hat{n}$, its polarization state $\hat{n}_p$, and its initial phase $2 \pi f_0 t_0$. For the purposes of map making, the phase will turn out to be irrelevant, and so, the full state-space is $S_2 \times S_2$. Our goal is to construct a quadratic estimate of the power at each point across this state space.

![](Figures/Stokes_angles.pdf)