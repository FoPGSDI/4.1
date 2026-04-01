# APPENDIX A
## THE SPLIT INTO TENSORIAL SPHERICAL HARMONICS AND THE REGGE-WHEELER CHOICE OF GAUGE

*(cf. Regge and Wheeler 1957)*

At the beginning of the analysis we expand the arbitrary metric of an equilibrium configuration in spherical harmonics. All quantities which transform as scalar fields under rotations are expanded in scalar spherical harmonics, $Y^l_m(\theta, \phi)$, which have *even parity*, $\pi = (-1)^l$. All fields which transform as vectors with respect to rotations are expanded in vector spherical harmonics, which are of two types: *even-parity* harmonics ($\pi = (-1)^l$)

$$\psi^l_{Mj} = \partial_j Y^l_M(\theta, \phi) \tag{A1a}$$

and *odd-parity* harmonics ($\pi = (-1)^{l+1}$)

$$\phi^l_{Mj} = \epsilon_j{}^{kl}\, \partial_k Y^l_M(\theta, \phi) \,. \tag{A1b}$$

Here $j$ and $k$ run over $x^3 = \theta$ and $x^4 = \phi$, and $\epsilon_j^{\ kl}$ is the antisymmetric tensor with respect to rotations

$$\epsilon_{34} = -1/\sin\theta \,, \qquad \epsilon_{43} = \sin\theta \,, \qquad \epsilon_{33} = \epsilon_{44} = 0 \,. \tag{A2}$$

All fields which transform as tensors with respect to rotations are expanded in tensor spherical harmonics, which have *even parity*

$$\Psi^l_{Mjk} = Y^l_M \gamma_{jk} \,, \qquad \hat{\Psi}^l_{Mjk} = \partial_j \partial_k Y^l_M \,, \tag{A3a}$$

or *odd parity*

$$\chi^l_{Mjk} = \tfrac{1}{2}\bigl(\epsilon_j{}^{kl}\, \partial_k \psi^l_{Ml} + \epsilon_k{}^{lm}\, \partial_l \psi^l_{Mj}\bigr) \,. \tag{A3b}$$

Here $\gamma_{jk}$ is the metric for rotations

$$\gamma_{33} = 1 \,, \qquad \gamma_{23} = \gamma_{32} = 0 \,, \qquad \gamma_{44} = \sin^2\theta \,; \tag{A4}$$

and the slash in equation (A3a) denotes covariant differentiation with respect to the metric $\boldsymbol{\gamma}$.

Consider first the *odd-parity* parts of the fluid displacement and metric perturbation. $\xi_r$, $h_{00}$, $h_{0r}$, and $h_{rr}$ are all scalars under the rotation group and therefore have vanishing odd-parity parts. $(\xi_\theta, \xi_\phi)$ and $(h_{0\theta}, h_{0\phi})$ are vectors under the rotation group; and $h_{\theta\theta}$, $(h_{\theta\phi},\, h_{\phi\phi})$ is a tensor under the rotation group. Consequently, the *odd-parity* $(l,\, M,\, \varpi = [-1]^{l+1})$ harmonics of the perturbation are

$$\begin{aligned}
\xi_r &= 0 \,,
&\xi_\theta &= U(r_f)\Psi^l{}_{M} \,,
&\xi_\phi &= V(r_f)\Psi^l{}_{M} \,;
&h_{0r} &= h_0(r_f)\Psi^l{}_{M} \,; \\
h_{0j} &= h_0(r_f)e^l{}_{jM} \,,
\quad
h_{1j} &= h_1(r_f)e^l{}_{jM} \,;
&&
h_{j k} &= h_2(r_f)\chi^l{}_{jkM} \,.
\end{aligned} \tag{A5}$$

Similarly, the general *even-parity* $(l,\, M,\, \varpi = [-1]^l)$ harmonics of the perturbation are

$$\begin{aligned}
\xi_r &= X(r_f)Y^l{}_{M} \,,
&\xi_\theta &= V(r_f)\Psi^l{}_{M} \,,
&\xi_\phi &= V(r_f)\Psi^l{}_{M} \,; \\
h_{00} &= e^{2\Phi}H_0 Y^l{}_{M} \,,
&h_{0r} &= H_1 Y^l{}_{M} \,,
&h_{rr} &= e^{2\Lambda}H_2 Y^l{}_{M} \,; \\
h_{0j} &= h_0(r_f)\Psi^l{}_{jM} \,,
\quad
h_{1j} &= h_1(r_f)\Psi^l{}_{jM} \,;
&&
h_{jk} &= r^2 K(r_f)g^l{}_{jkM} + r^2 G(r_f)g^l{}_{jkM} \,.
\end{aligned} \tag{A6}$$

In order to simplify the above expressions we introduce, for each spherical harmonic, a small coordinate transformation similar to that which Regge and Wheeler (1957) used in studying perturbations of the Schwarzschild geometry:

$$x^{\mu'} = x^\mu + \varphi^\mu(x) \,. \tag{A7}$$

In the new coordinate system (Regge--Wheeler gauge) the metric perturbation has the new form

$$h_{\mu'\nu'} = h_{\mu\nu} + \varphi_{\mu;\nu} + \varphi_{\nu;\mu} \,. \tag{A8}$$

For the *odd-parity* harmonic $(l,\, M,\, \varpi = [-1]^{l+1})$, we take

$$\eta_0 = 0 \,, \qquad \eta_1 = \Lambda(r_f)\Psi^l{}_{M} \,, \tag{A9}$$

where $\Lambda$ is chosen so as to annul the function $h_2(r_f)$ in equations (A5):

$$h_2(r_f) = 0 \quad \text{for odd-parity harmonic.} \tag{A10}$$

For the *even-parity* $(l,\, M,\, \varpi = [-1]^l)$, we take

$$\eta_0 = M_0(r_f)Y^l{}_{M} \,, \qquad \eta_1 = M_1(r_f)\Psi^l{}_{M} \,, \tag{A11}$$

where $M_0$, $M_1$, and $M_2$ are chosen so as to annul the functions $h_0(r_f)$, $h_1(r_f)$, and $G(r_f)$ in equations (A6):

$$h_0(r_f) = h_1(r_f) = G(r_f) = 0 \,. \tag{A12}$$

When we specialize to components $M = 0$ and change notation slightly, the general odd-parity gauge (A5), (A10) in the Regge--Wheeler gauge takes on the forms (7) used in this paper; and the general even-parity gauge (A6), (A12) takes on the form (7).
