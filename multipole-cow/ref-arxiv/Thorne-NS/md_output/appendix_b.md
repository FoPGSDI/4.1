# APPENDIX B
## EQUATIONS OF MOTION FOR ODD PARITY

The fluid 4-velocity corresponding to the odd-parity displacements of equation (6) is, to first order in the displacement,

$$
u_0 = e^{\Phi} \,,
\quad
u_r = 0 \,,
\quad
u_\theta = 0 \,,
\quad
u_\phi = e^{-\Phi} U_{,t} \sin\theta\, \partial_\theta P_l(\cos\theta) \,.
\tag{B1}
$$

For odd parity the pressure and density are unchanged since they are scalars under the rotation group. Consequently, the perturbation in the stress-energy tensor, $T_{\mu\nu} = (\rho + p)\,u_\mu u_\nu - p\,g_{\mu\nu}$, has as its only non-vanishing components

$$\begin{aligned}
\delta T_{t\phi} = \delta T_{\phi t} &= [(\rho+p)\,U_\phi - p\,h_t]
    \sin\theta\,\partial_\theta P_l(\cos\theta) \,, \\
\delta T_{r\phi} = \delta T_{\phi r} &= -p\,h_r
    \sin\theta\,\partial_\theta P_l(\cos\theta) \,.
\end{aligned}\tag{B2}$$

We have calculated by hand the perturbations in the Ricci curvature tensor, $\delta R_{\mu\nu}$, associated with the odd-parity metric perturbation of equation (6b). Our result, which has been checked against similar hand computations of Edelstein and Vishveshwara (1966) and against IBM 7094 computations using a program of Thorne and Zimmerman (1967), is

$$\begin{aligned}
\delta R_{t\phi} &=
    e^{-\lambda}\bigl[{-h_{0,rr}'} + \tfrac{1}{2}(\nu'+\lambda')\,h_{0,r}'
    + [-l(l+1)\,e^{-\nu} + 2r^{-1}e^{-\lambda}]\,r^{-1}h_0'\bigr] \\
&\quad
    - h_{0,t}'\bigl[-2r^{-1}e^{-\nu} + \tfrac{1}{2}(\nu'+\lambda')\,e^{-\nu}\bigr] \,, \\[4pt]
\delta R_{r\phi} &=
    \tfrac{1}{2}\bigl[e^{-\nu-\lambda}h_{0,tr}'' + 2r^{-1}e^{-\lambda}h_{0,r}'\bigr] \,, \\[4pt]
\delta R_{\theta\phi} &=
    \tfrac{1}{2}\bigl[e^{-\nu-\lambda}h_{0,t\theta}'' + 2e^{-\lambda}r^{-1}h_{0,\theta}'\bigr] \,, \\[4pt]
\delta R_{\phi\phi} &=
    \tfrac{1}{2}\bigl[e^{-\nu}h_{0,tt}'' + 2r^{-2}e^{-\lambda}h_0''
    - 2\cot\theta\,\partial_\theta P_l(\cos\theta)\bigr] \,.
\end{aligned}\tag{B3}$$

All other components vanish.

The Einstein field equations,

$$
\delta R_{\mu\nu} = 8\pi\,\delta\!\left(T_{\mu\nu} - \tfrac{1}{2}g_{\mu\nu}T\right) \,,
\tag{B4}
$$

which govern the odd-parity perturbation can be put into the following form by combining equations (B3), (B2), (6b), and (3):

$$\begin{aligned}
h_{0,rr} &- r^{-1}[e^{(\nu-\lambda)/2}]^{-1}[l(l+1)+1]^{1/2}[l(l+1)-2]\,r^{-1}h_0
    - 2r^{-2}e^{-\lambda}h_0 = 0 \,,
\end{aligned}\tag{B5a}$$

$$\begin{aligned}
h_{0,t} &= e^{(\nu+\lambda)/2}\bigl[l(l+1)/r^2\bigr]\,h_1 \,,
\end{aligned}\tag{B5b}$$

$$\begin{aligned}
16\pi\,(\rho+p)\,U &=
    -\bigl[{-h_{1,t}'} + (-P-l+2)\,r^{-1}e^{-\lambda} - 2/r^2\bigr] \\
&\quad
    - \tfrac{1}{2}(\lambda'+\nu')\,h_1'
    + (r^{-1}e^{-\lambda})\,h_0' = F(t) \,.
\end{aligned}\tag{B5c}$$

Note that equations (B5a) and (B5b) guarantee that the right-hand side of equation (B5c) is time-independent and, therefore, the fluid does not pulsate:

$$
U_{,tt} = 0 \,.
\tag{B6}
$$

Rather, the fluid undergoes a continuous, non-varying differential rotation.

In physical terms equation (B5a) is the propagation equation for the odd-parity gravitational waves, which do not couple to the star in any way. Equations (B5b) and (B5c) together fix $h_0$, $h_1$, and, thence, the dragging of inertial reference frames, once the differential rotation and the decoupled gravitational waves have been specified.

When gravitational waves are absent, $U$, $h_0$, $h_1$ are all constant in time; the star rotates rigidly (eq. [B5b]), and the local inertial frames are dragged along by the star with angular velocity as measured by a distant observer

$$
\Omega_\mathrm{star} = \bigl[(h_0 - U_\phi)/(r^2\sin\theta)\bigr]
    \partial_\theta P_l(\cos\theta) \,;
\tag{B7}
$$

and the local inertial frames are dragged along with the star with angular velocity as measured by a distant observer

$$
\Omega_\mathrm{drag} = \Omega_\mathrm{local} =
    (h_0/r^2\,e^\nu)\sin\theta\,\partial_\theta P_l(\cos\theta) \,.
\tag{B8}
$$
