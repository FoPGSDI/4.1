# APPENDIX C

## EQUATIONS OF MOTION FOR EVEN PARITY

The fluid 4-velocity corresponding to the even-parity equation of (7) in, to first order in the displacement,

$$\begin{aligned}
u^{0} &= e^{-\nu/2}\bigl[1 - \tfrac{1}{2}H_{0}\,P_{l}(\cos\theta)\bigr] , &
u^{r} &= e^{-\lambda/2} r^{-1} e^{-(\lambda+\nu)/2} \dot{W}_{l}\,P_{l}(\cos\theta) ,
\\
u^{\theta} &= -r^{-1} e^{-\nu/2} \dot{V}_{l}\,\partial_{\theta}P_{l}(\cos\theta) , &
u^{\phi} &= 0 \;.
\end{aligned} \tag{C1}$$

In the pulsating configuration the Lagrangian change in the number density of baryons is

$$\Delta n = -n\,\Theta^{k}{}_{;k} - \tfrac{1}{2}n\,\dot{h}^{l}{}_{l}\,({}^{(3)}g)^{1/2} \;. \tag{C2}$$

Here $\Theta^{k}{}_{;k}$ is the divergence of the fluid displacement with respect to the 3-geometry at constant time, $t$; and ${}^{(3)}g$ is the determinant of the metric of that 3-geometry. In the Regge-Wheeler choice of gauge, equation (C2) reduces to[^3]

$$\Delta n/n = \bigl[-r^{-2}e^{(\lambda-\nu)/2}\dot{W} - l(l+1)r^{-2}\dot{V} + \tfrac{1}{2}\dot{H}_{1} + \dot{K}\bigr]P_{l}(\cos\theta) \;. \tag{C3}$$

The corresponding Eulerian changes in pressure and in density of mass-energy are

$$\begin{aligned}
\delta p &= (p+\rho)\,(\Delta n/n) - p'\,r^{-1}e^{-(\lambda+\nu)/2}W\,P_{l}(\cos\theta) ,
\\
\delta\rho &= \gamma p\,(\Delta n/n) - \rho'\,r^{-1}e^{-(\lambda+\nu)/2}W\,P_{l}(\cos\theta) \;;
\end{aligned} \tag{C4}$$

and the Eulerian changes in the stress-energy tensor are

$$\begin{aligned}
\delta T^{0}{}_{0} &= -\delta\rho , &
\delta T^{1}{}_{1} &= \delta T^{2}{}_{2} = \delta T^{3}{}_{3} = \delta p , &
\delta T^{0}{}_{1} &= (p+\rho)\,u^{1} ,
\\
\delta T^{1}{}_{0} &= -(p+\rho)\,u^{1} , &
\delta T^{2}{}_{0} &= (p+\rho)\,u^{0}\mu^{2} , &
\delta T^{3}{}_{0} &= (p+\rho)\,u^{0}\mu^{3} \;.
\end{aligned} \tag{C5}$$

All other components of $\delta T^{\mu}{}_{\nu}$ vanish.

The way we had the perturbations in the Einstein field equations, $\delta G^{\mu}{}_{\nu} = 8\pi\,\delta T^{\mu}{}_{\nu}$, and in the equations of motion of the fluid, $\delta(T^{\mu}{}_{\nu;\mu}) = 0$, associated with the perturbed stress-energy tensor (C5) and the perturbed metric (7b); and we have checked our result on the IBM 7094 using the computer programs of Thorne and Zimmerman (1967).[^4] The origins of equations (8) and (9) are as follows. Equation (8a) is $\delta G^{t}{}_{t} - \delta G^{r}{}_{r} = 8\pi\,\delta T^{t}{}_{t} - 8\pi\,\delta T^{r}{}_{r}$, and it has been used to eliminate $H_{2}$ from all equations. Equation (8b) is $\delta G^{r}{}_{t} = 8\pi\,\delta T^{r}{}_{t}$, and it has been used to eliminate $H_{1}$ from all equations. Equation (8c) is $\delta G^{r}{}_{r} = 8\pi\,\delta T^{r}{}_{r}$, and it has been used to eliminate $\dot{K}$ from all other equations. Equation (8d) is $\delta G^{\theta}{}_{\theta} = 8\pi\,\delta T^{\theta}{}_{\theta}$, and it has been used to eliminate $\dot{H}_{1}$ from all other equations. Equation (9a) is $\delta G^{r}{}_{\theta} = 0$. Equation (9b) is $\delta(T^{r}{}_{\theta;\mu}) = 0$. Combined with equation (8a), equation (9b) is $\delta G^{\theta}{}_{r} = 0$. Equation (9c) is $\delta(T^{r}{}_{\mu;\nu}) = 0$. Equation (9d) is $\delta(T^{r}{}_{\mu;\nu}) = 0$.

Several other useful dynamical equations for the perturbed configuration, which can be derived from equations (2), (3), (8), and (9) with a large amount of effort, are these: $\delta G^{r}{}_{t} = 8\pi\,\delta T^{r}{}_{t}$ is

$$H_{1}'' + \bigl[\tfrac{1}{2}r^{-1}(\lambda' - \nu') - 2r^{-1}\bigr]H_{1}' + \bigl[2r\nu' e^{-\nu} p - \nu'' - \tfrac{1}{4}(\nu')^{2}\bigr]H_{1} = e^{\lambda}(K_{1} + H_{00}) - 16\pi(p+\rho)V \;. \tag{C6}$$

and $8(G^1 + 2G\gamma) = 8\pi(T_r^0 - T_\phi^0 + 2T_r^r)$, when combined with equation (8c), is

$$\begin{aligned}
H_{0,tt}
  &= -e^{-2\Lambda} H_0''
     + e^{-2\Lambda} \bigl[1 - 5m/r + 2\pi r^2 - 3p\bigr] H_0'
\\
  &\quad
     + 2\pi^{-1} r^{-1} \bigl\{4m/r - 10m^2/r^2 - 3p\bigr\} p
\\
  &\quad
     + 2\pi r^{-1} \bigl[1 - \gamma + (\gamma - 3)\,2m/r\bigr]
     - 32\pi r^2 e^{-2\Lambda} p^3
\\
  &\quad
     + l[l+1]\,e^{-2\Lambda} H_0
     + 2r^{-2} \bigl[-2 + 6m/r + 8\pi r^2 \rho\bigr] K'
\\
  &\quad
     + 8\pi e^2 (\rho + \gamma p) K
     - 8\pi e^2 (\rho + \gamma p)\,r^{e^{-\Lambda}} W''
\\
  &\quad
     - 8\pi (r^2 - p)^2\,r^{-1} e^{-\Lambda} W
     - 8\pi l(l+1)\,e^2 (\rho + \gamma p)\,r^{-1} V
\\
  &= 0 \,.
\end{aligned} \tag{C7}$$

---

[^3]: The peculiar form (7a) in which we chose to define $\delta n$ was motivated by a desire to take on as simple a form as possible and to thereby simplify the eigensys. (14).

[^4]: $G^{r}{}_{r}$, $K_{r}$, $T^{r}{}_{r}$, and $T^{\phi}{}_{\phi}$ for even-parity pulsations are presented in Thorne and Zimmerman (1967) precisely as they were output by the computer.
