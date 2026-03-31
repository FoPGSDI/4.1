# Calculation of the Love number

## Equilibrium configuration

The
geometry of spacetime of a spherical, static star can be described
by the line element [mtw]
$$
ds_0^2=g^{(0)}_{\alpha \beta}dx^\alpha
dx^\beta=-e^{\nu(r)}dt^2+e^{\lambda(r)}dr^2+r^2
\left(d\theta^2+\sin^2\theta d \phi^2 \right).
$$ The star’s stress-energy
tensor is given by $$ T_{\alpha
\beta}=\left(\rho+p\right)u_\alpha u_\beta+pg^{(0)}_{\alpha
\beta},$$ where
$\vec{u}=e^{-\nu/2}\partial_t$ is the fluid’s four-velocity and
$\rho$ and $p$ are the density and pressure. Numerical integration
of the Tolman-Oppenheimer-Volkhov equations (see e.g. [mtw])
for neutron star models with a polytropic pressure-density relation
$$ P=K\rho^{1+1/n},$$ where $K$ is a constant and $n$ is the
polytropic index, gives the equilibrium stellar model with radius
$R$ and total mass $M=m(R)$.

## Static linearized perturbations due to an external
tidal field

We examine the behavior of the
equilibrium configuration under linearized perturbations due to an
external quadrupolar tidal field following the method of [tc].
The full metric of the spacetime is given by $$
g_{\alpha\beta}=g^{(0)}_{\alpha \beta} +h_{\alpha
\beta},$$ where $h_{\alpha\beta}$ is a
linearized metric perturbation. We analyze the angular dependence of
the components of $h_{\alpha\beta}$ into spherical harmonics as in
[rw]. We restrict our analysis to the $l=2$, static,
even-parity perturbations in the Regge-Wheeler gauge [rw].
With these specializations, $h_{\alpha\beta}$ can be written as
[rw, tc]:
$$ h_{\alpha \beta} ={\rm diag}
\left[e^{-\nu(r)}H_0(r),   e^{\lambda(r)} H_2(r),   r^2 K(r),
r^2 \sin^2\theta K(r)\right] Y_{2m}(\theta, \varphi).
$$
The nonvanishing components of the perturbations of the
stress-energy tensor ([ref:tab]) are $\delta T^0_0=-\delta
\rho=-(dp/d\rho)^{-1}\delta p$ and $\delta T^i_i=\delta p$. We
insert this and the metric metric perturbation ([ref:hrw]) into the
the linearized Einstein equation $\delta G_{\alpha}^{\beta} = 8 \pi
\delta T_{\alpha}^{\beta}$ and combine various components. From
$\delta G^\theta_\theta-\delta G^\phi_\phi=0$ it follows that that
$H_2=H_0\equiv H$, then $\delta G^r_\theta=0$ relates $K’$ to $H$,
and after using $\delta G_\theta ^\theta+\delta G_\phi^\phi=16\pi
\delta p$ to eliminate $\delta p$, we finally subtract the $r-r$
component of the Einstein equation from the $t-t$ component to
obtain the following differential equation for $H_0\equiv H$ (for
$l= 2$):
$$\begin{aligned}
&&H{”}+H{’} \left[{2 \over r} + e^{\lambda} \left( {2m(r)\over
r^2}
+ 4 \pi r \left(p-\rho\right)\right) \right]\\
&&+H\left[ -{6 e^{\lambda} \over r^2 } + 4 \pi e^{\lambda}\left( 5
\rho + 9 p +
 {\rho + p \over \left(dp/d\rho\right)} \right)
 - \nu{’}^2 \right]=0,
\end{aligned}$$ where the prime denotes $d/dr$. The
boundary conditions for Eq. ([ref:Hdiffeq]) can be obtained as
follows. Requiring regularity of $H$ at $r=0$ and solving for $H$
near $r=0$ yields $$ H(r)=a_0
r^2\left[1-\frac{2\pi}{7}\left(5\rho(0)+9p(0)+\frac{\rho(0)+p(0)}{(dp/d\rho)(0)}\right)r^2+
O(r^3)\right], $$ where $a_0$ is a
constant. To single out a unique solution from this one-parameter
family of solutions parameterized by $a_0$, we use the continuity of
$H(r)$ and its derivative across $r=R$. Outside the star, Eq.
([ref:Hdiffeq]) reduces to
$$
H{”} + \left( {2\over r} - \lambda{’}\right) H{’} -\left( {6
e^{\lambda} \over r^2 } + \lambda{’}^2 \right) H= 0,

$$
and changing variables to $x =( r/M -1)$ as in [tc] transforms
Eq. ([ref:Hdiffeqouts1]) to a form of the associated Legendre
equation with $l=m=2$:
$$
\left(x^2 -1 \right) H{”} + 2 x H{’} - \left( 6 + {4\over x^2 -1 }
\right)H=0.
$$
The general solution to Eq. ([ref:houtsideeq]) in terms of the
associated Legendre functions $Q_2{\,}^2(x)$ and $P_2{\,}^2(x)$ is
given by
$$H = c_1 Q_2{\,}^2 \left( {r \over M}-1 \right) + c_2
P_2{\,}^2\left( {r \over M}-1 \right),
$$
where $c_1$ and $c_2$ are coefficients to be determined.
Substituting the expressions for $Q_2{\,}^2(x)$ and $P_2{\,}^2(x)$
from [as] yields for the exterior solution $$\begin{aligned}
H&=&c_1 \left({r \over M}\right)^2\left(1-{2M \over r} \right)
\left[ -{M(M-r)(2M^2+6Mr-3r^2)\over r^2(2M-r)^2} + {3 \over 2}
\log{\left({r \over r-2M}\right)} \right]\\
&&+  3 c_2  \left({r \over M}\right)^2\left(1-{2M \over r} \right).
\end{aligned}$$ The asymptotic behavior of the solution
([ref:Hout]) at large $r$ is
$$ H = {8 \over 5} \left({M \over
r}\right)^3 c_1+O\left(\left(\frac{M}{r}\right)^4\right) + 3\left({r
\over M}\right)^2
c_2+O\left(\left(\frac{r}{M}\right)\right),$$
where the coefficients $c_1$ and $c_2$ are determined by matching
the asymptotic solution ([ref:asympH]) to the expansion ([ref:Qdef])
and using Eq. ([ref:Lovedef1]):
$$ c_1 = {15 \over 8}{1 \over M^3} \lambda {\cal{E}},
{\;}{\;}{\;}{\;}{\;}c_2 = {1 \over 3} M^2
{\cal{E}}.$$ We now solve for
$\lambda$ in terms of $H$ and its derivative at the star’s surface
$r=R$ using Eqs. ([ref:coefficients]) and ([ref:Hout]), and use the
relation ([ref:k2def]) to obtain the expression:
 $$\begin{aligned}
&& k_2 = \frac{8C^5}{5}\left(1-2C\right)^2
\left[2+2C\left(y-1\right)-y\right]\times\\
&&\bigg\{2C\left(6-3 y+3 C(5y-8)\right)+4C^3\left[13-11y+C(3 y-2)+2
C^2(1+y)\right] \\
&&
+3(1-2C)^2\left[2-y+2C(y-1)\right]\log\left(1-2C\right)\bigg\}^{-1},
\end{aligned}$$ where we have defined the star’s compactness parameter $C\equiv
M/R$ and the quantity $y\equiv RH{’}(R)/ H(R)$, which is obtained by
integrating Eq. ([ref:Hdiffeq]) outwards in the region $0<r<R$.

## Newtonian limit

The first term in the expansion of the expression ([ref:k2expr]) in
$M/R$ reproduces the Newtonian result:
$$ k_2^N=\frac{1}{2 }\left(\frac{2- y}{ y + 3
}\right),$$ where the superscript $N$
denotes "Newtonian". In the Newtonian limit, the differential
equation ([ref:Hdiffeq]) inside the star becomes $$
H”+\frac{2}{r}H’+\left(\frac{4\pi\rho}{dp/d\rho}-\frac{6}{r^2}\right)H=0.
$$ For a polytropic index of $n=1$, Eq. ([ref:Newtoniandiff]) can
be transformed to a Bessel equation with the solution that is
regular at $r=0$ given by $ H=A \sqrt{{r/R}}   J_{5/2}(\pi r/ R)$,
where $A$ is a constant. At $r=R,$ we thus have $y=RH{’}/ H =(\pi^2
- 9 )/3,$ and from Eq. ([ref:k2expr]) it follows that
$$ k_2^N(n=1)=\left( -{1 \over 2} + {15 \over 2 \pi^2}\right) \approx
0.25991,$$ which agrees with the known result of
[bo].