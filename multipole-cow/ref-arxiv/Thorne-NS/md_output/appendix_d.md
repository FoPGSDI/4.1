# APPENDIX D

## BOUNDARY CONDITIONS FOR EVEN-PARITY EIGENFUNCTIONS

Near the center of the star ($r \ll M$, $r \ll r_s^{-1/3}$), the eigenfunctions (14) reduce to the simpler form

$$\begin{aligned}
H &= \tfrac{1}{2} r^l K'' + \bigl[1 - l(l+1)/2\bigr] K \,, \tag{D1a} \\[4pt]
H' + \bigl[1 + l(l+1)/2\bigr] H/r
     &= K'' + 3K' + \bigl[1 - l(l+1)/2\bigr] K/r \,, \tag{D1b} \\[4pt]
W' &= -l(l+1)\,V \,, \tag{D1c} \\[4pt]
K' - H' &= \omega^2 e^{-\nu} \bigl(V + W/r\bigr) \,, \tag{D1d}
\end{aligned}$$

where $r_s$ is the value of $r$ at the center of the star. Equation (D1a) arises from the leading terms in (14c) and (14b). Equation (D1d), the remaining content of (14c) and (14d), is the sum of (14c) and (14d) used to eliminate $W$ and $V$.

Equations (D1) are a 4th-order system of linear differential equations, so there are four linearly independent solutions. Only two of the independent solutions---those of equations (15a)---are physically acceptable. Two unacceptable solutions correspond to the two arbitrary constants $C$ and $D$ in

$$\begin{aligned}
K &= C\,r^{-(l+1)} + \cdots \,, & H &= C\,r^{-(l+1)} + \cdots \,, \\
W &= D\,r^{-(l+1)} + \cdots \,, & V &= \tfrac{1}{2} D/(l+1)\cdot r^{-(l+1)} + \cdots \,;
\tag{D2}
\end{aligned}$$

and the third unacceptable solution corresponds to the arbitrary constant, $E$, in

$$\begin{aligned}
K &= E\,r^{2-l(l+1)/4} + \cdots \,, \\
H &= \tfrac{1}{2}\bigl\{1 - l(l+1)/2 - l(l+1)/2\bigr\} E\,r^{2-l(l+1)/4} + \cdots \,, \\
W &= \omega^2 e^{-\nu}\{1+l\}\bigl[1 + l(l+1)/2\bigr] E\,r^{l(l+1)(l+1)/2} + \cdots \,, \\
V &= -\omega^2 e^{-\nu}\bigl\{-1 - l(l+1)/2\bigr\} E\,r^{l(l+1)/2} + \cdots \,.
\tag{D3}
\end{aligned}$$

The solutions (D2) and (D3) are unacceptable for $l \geq 1$ because they lead to perturbations in the density, the pressure, and the geometry of spacetime which diverge at $r = 0$. (Cf. eqs. [C3], [C4], and [7b].)

At the surface of the star we demand that the Lagrangian change in the pressure (eq. [160]) vanish and that the fluid displacement and the perturbation in the spatial geometry not diverge.

As for radial perturbations (see Bardeen, Thorne, and Meltzer 1966), so also here, the precise analytic conditions which these demands place on the eigenfunctions depend upon the distribution of pressure, density, and adiabatic index near the surface of the equilibrium configuration. In all cases the eigenfunctions will have to assume the regular form (15b) at the surface, but the relation of the constants $k_j$, $h_j$, $W_j$, and $v_j$ to each other will vary from case to case.

All cases of current interest fall into two classes: (1) *Absolute zero temperature at the surface*; example---the Harrison-Wakano-Wheeler configurations (see Harrison et al. 1965). In this case the pressure, density, and adiabatic index near the surface have the form

$$p \sim p_0(R - r) \,, \qquad \rho \sim p_0 + p_1(R - r) \,, \qquad \gamma p \sim \gamma_0 + \gamma_1(R - r) \,; \tag{D4}$$

and, consequently, all five solutions to the eigenequations (14) have the regular form (15b) near the surface. However, only four of the five solutions are acceptable because only four have vanishing Lagrangian change in pressure.

(2) *Polytropic pressure-density relation near the surface*; example---any hot stellar model. In this case the pressure, density, and adiabatic index near the surface are related by

$$p = Q\rho^{1+1/N} \,, \qquad \gamma = \gamma_0 \,; \tag{D5}$$

and the equation of hydrostatic equilibrium then sets up the pressure and density distributions

$$p = a(R - r)^{N+1} \,, \qquad \rho = b(R - r)^N \,. \tag{D6}$$

Here $a$, $b$, $\gamma_0$, and $Q$ are constants. An examination of the eigenequations (14) near the surface for such configurations reveals that only four of the five solutions have the regular form (15b). These four solutions are physically acceptable because equations (D5) and (D6) guarantee that $\gamma p$ vanishes at the surface, and therefore---for regular solutions---that $\Delta p$ vanishes (see eq. [16]). The fifth, unacceptable, solution has the surface behavior

$$\begin{aligned}
K &= k_0 + k_1(R - r) + \cdots \,, \quad
H = h_0 + h_1(R - r) + \cdots \,, \quad
W = G(R - r)^{-N} \,, \\[4pt]
V &= -\frac{G(1 - 2M/R)^{1/2}\,(N\gamma_0 - N - 1)}{\omega^2 R}\,(R - r)^{-N} \,.
\tag{D7}
\end{aligned}$$

Here $k_0$, $k_1$, $h_0$, $h_1$ are constants freely adjustable by the addition of varying amounts of the regular solutions (15b); while $G$ is the arbitrary constant which characterizes the solution.
