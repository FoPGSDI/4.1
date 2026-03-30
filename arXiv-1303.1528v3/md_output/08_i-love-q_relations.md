# I-Love-Q Relations

Now that the moment of inertia, quadrupole moment and Love numbers have been calculated, let us present the universal I-Love-Q relations. We first show numerical results and a fitting curve through these. Then, we obtain analytic I-Love-Q relations for the $n=0$ and 1 polytropic EoSs in the Newtonian limit.

## Numerical Results

![](I-Love-Newton-Dup-maxC-PRD.eps)
Figure [ref:fig:I-Love-Q] shows universal relations between dimensionless quantities, $\bar{I}$, $\bar{Q}$, $\bar{\lambda}^\mathrm{(tid)}$ and $\bar{\lambda}^\mathrm{(rot)}$ for various EoSs. Notice that these (barred) 4 quantities are essentially {*independent of the NS spin*}, to second order in the slow-rotation approximation[^1].  The parameter varied along each curve is the NS central density, or equivalently the NS compactness. Therefore, for the polytropic EoS the I-Love-Q relations are independent of the polytropic amplitude coefficient $K$ in Eq. (Eq. polytropic). The bottom part of each panel shows the relative fractional difference between each of the curves and curve corresponding to the $n=1$ polytropic EoS. For reference, the top axes show the NS mass with the APR EoS. The vertical dashed lines correspond to $M_*=1M_\odot$ for the APR EoS and points to the left of these lines correspond to more massive NSs with higher compactness. Observe that, for realistic EoSs with $M_*>1M_\odot$, the fractional relative differences are $\mathcal{O}(1)%$. Observe also that the polytropic I-Love-Q relations deviate from those with realistic EoSs as one increases $n$, ie. as the NS becomes more centrally-concentrated. We see this as evidence that the I-Love-Q trio is most sensitive to the NS outer layers, where realistic EoSs mostly agree with each other. Curiously, the fractional relative difference between the $n=1$ and $n=0$ polytrope (constant density NS star) is also of ${\cal{O}}(1)%$, but this case will be studied analytically in Sec. [ref:analytic-exp].

Observe that, in general, the dependence of the I-Love-Q relations on any EoS becomes weaker as the NS compactness $C$ increases, ie. from right to left in any of the panels of Fig. [ref:fig:I-Love-Q]. This may be evidence that at least part of the universality observed is due to the NS sequence approaching a BH as $C \to 0.5$, where the latter has no internal-structure dependence by the no-hair theorems. As $C \to 0.5$, the asymptotic values of $\bar{I}$, $\bar{\lambda}^\mathrm{(tid)}$ and $\bar{Q}$ for a BH are $\bar{I} \to 4$ [membrane], $\bar{\lambda}^\mathrm{(tid)} \to 0$ [binnington-poisson] [see also Eq. (Eq. k2)] and $\bar{Q} \to 1$ [geroch,hansen,poisson-quadrupole] respectively. In each panel of Fig. [ref:fig:I-Love-Q], we show the BH limit of either $\bar{I}=4$, $\bar{Q} =1$ or $\bar{\lambda}^\mathrm{(rot)} = 16$ as a horizontal dashed line. Observe that the I-Love-Q relations asymptote to such BH values as $C$ increases. However, one cannot quite reach this limit, as one can never construct a BH solution by increasing the central density of a NS solution by a finite amount. We think that this is why the relative fractional differences shown in the bottom panels of Fig. [ref:fig:I-Love-Q] do not decrease to zero as one increases the central density.

Unlike the other relations, the $\bar{\lambda}^\mathrm{(rot)}$--$\bar{\lambda}^\mathrm{(tid)}$ relation (bottom, right panel of Fig. [ref:fig:I-Love-Q]) depends very weakly on the EoS even when the NS compactness is relatively small (as one approaches the Newtonian limit). In fact, one can show that the relation $\bar{\lambda}^\mathrm{(rot)}=\bar{\lambda}^\mathrm{(tid)}$ holds exactly in the Newtonian limit for any EoS, as we will discuss in Sec. [ref:sec:Love-Newton].

{\renewcommand{1.2}{1.2}
| \noalign{\smallskip}
$y_i$ | $x_i$ |  | \multicolumn{1}{c}{$a_i$} | \multicolumn{1}{c}{$b_i$} | \multicolumn{1}{c}{$c_i$} | \multicolumn{1}{c}{$d_i$} | \multicolumn{1}{c}{$e_i$} |
|---|---|---|---|---|---|---|---|
| \noalign{\smallskip}
$\bar{I}$ | $\bar{\lambda}^\mathrm{(tid)}$ |  | 1.47 | 0.0817 | 0.0149 | $2.87\times 10^{-4}$ | $-3.64\times 10^{-5}$ |
| $\bar{I}$ | $\bar{Q}$ |  | 1.35 | 0.697 | -0.143 | $9.94\times 10^{-2}$ | $-1.24\times 10^{-2}$ |
| $\bar{Q}$ | $\bar{\lambda}^\mathrm{(tid)}$ |  | 0.194 | 0.0936 | 0.0474 | $-4.21\times 10^{-3}$ | $1.23\times 10^{-4}$ |
| \noalign{\smallskip} |
}

![](I-Q-fit-Dup-maxC-PRD.eps)

Given the universality of the I-Love-Q relations, one can fit them all with a single curve:

$$
\ln y_i = a_i + b_i \ln x_i + c_i (\ln x_i)^2 + d_i (\ln x_i)^3+ e_i (\ln x_i)^4\,,

$$

where the coefficients are summarized in Table [ref:table:coeff].
Figures [ref:fig:I-Love-fit]  and [ref:fig:I-Q-fit] show the fitting curves for the I-Love, Q-Love and I-Q relations, together with the relative fractional difference between the fitting curves and all other EoS curves. For the polytropic EoSs, we do not show the results when $n=0$, 2, 2.5 and $3$ since such EoSs do not model NSs well; instead, we add results when $n=0.6$ and 0.8. As one can see, the fitting curves are accurate to within $\mathcal{O}(1)%$ accuracy.

[^1]: Formally, the barred quantities depend on the mass $M_{*}$, which is not the observed mass. The two are related by $M_{\rm ons} = M_{*} \left(1 + \chi^{2} \delta M\right)$. Rapidly rotating NS calculations indicate that $\delta M = {\cal{O}}(0.3)$, and thus, for $\chi < 0.1$, such spin dependence introduces corrections of ${\cal{O}}(0.003)$.

## Analytical Explanations

The universal I-Love-Q relations presented in the previous subsection are quite intriguing, and thus, they beg for an analytic explanation. We will attempt one here, by investigating these relations for certain EoSs that allow for an analytical treatment. In particular, we will study the $n=0$ and $n=1$ polytropic EoSs in the Newtonian limit, for which the moment of inertia, the quadrupole moment and the Love numbers can be computed fully analytically. We can then derive the I-Love-Q relations analytically as well to try to obtain an analytical explanation for these relations.

### Newtonian $\lambda$ for Generic EoSs

In Newtonian theory, the curl of the equation of hydrostatic equilibrium, $\nabla p = \rho \nabla \Psi$ where $\Psi$ is the total gravitational potential, with the Newtonian force vanishes, ie. $\nabla \rho \times \nabla \Psi =0$. Thus, surfaces of constant $\rho$ and $\Psi$ coincide. One can express such surfaces in terms of a radial parameter $a$ as [mora-will]

\allowdisplaybreaks
$$\begin{aligned}
r(a, \theta, \phi) &=& a \left[ 1 + \sum_{\ell, m} f_{\ell} Y_{\ell m} (\hat{\Omega}) Y_{\ell m} (\hat{n}) \right]\,,  \\
&=& a \left[ 1 + \frac{5}{4\pi} \sum_{\ell} f_{\ell} P_\ell (\hat{\Omega} \cdot \hat{n}) \right]\,,
\end{aligned}$$

where $f_\ell$ is the dimensionless distortion function of the constant $a$. This function is related to the $\ell =2$ tidal apsidal constant $k_2^\mathrm{(tid)}{}^{,\text{N}}$ by

$$
k_2^\mathrm{(tid)}{}^{,\text{N}} = \frac{3-\eta_2(a_*)}{2[2 + \eta_2(a_*)]}\,,

$$

where $a=a_*$ denotes the surface of the star and

$$
\eta_2 (a) \equiv \frac{d \ln f_2}{d \ln a}\,.
$$

This function can be obtained by solving the Clairaut-Radau equation [brooker-olle,mora-will]

$$
a\frac{d\eta_2}{d a} + 6 \mathcal{D} (\eta_2 +1) + \eta_2 (\eta_2 -1) -6 =0\,,

$$

with the boundary condition $\eta_2 (0) =0$, where $\mathcal{D}(a) \equiv \rho(a)/\bar{\rho}(a)$ with $\bar{\rho}$ representing the mean density of the star.

In the Newtonian limit, both rotational and tidal apsidal constants can be calculated from Eq. (Eq. Clairaut-Radau) with the same boundary condition [mora-will], and hence,

$$
\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} = \bar{\lambda}^\mathrm{(rot)}{}^{,\text{N}}\,.

$$

This shows that the rotationally-induced and tidally-induced NS deformabilities are exactly the same in the Newtonian limit. In GR, non-linear effects modify this relation and break the equality. Equation (Eq. Love-Love-Newton) is shown as a black, thin, solid line in the bottom, right panel of Fig. [ref:fig:I-Love-Q]. Notice that all the curves approach this Newtonian result as one decreases the compactness, as expected.

One can also calculate the $\ell=2$ tidal apsidal constant $k_2^\mathrm{(tid)}{}^{,\text{N}}$ by taking the Newtonian limit ($R \gg M(R)$, $\rho \gg p$) of Eq. (Eq. k2) [hinderer-love], where Eq. (Eq. h2RR) in the Newtonian limit is given by

$$
\frac{d^2 h_2}{dR^2} + \frac{2}{R} \frac{d h_2}{dR} + \left( 4 \pi \rho \frac{d\rho}{dp} - \frac{6}{R^2} \right)h_2 =0\,.

$$

### Polytropic I-Love-Q relations: $n=0$

Now, let us investigate the I-Love-Q relations in the Newtonian limit for specific EoSs. First, we focus on the polytropic EoS with $n=0$, which corresponds to the incompressible EoS with

$$
\rho = \rho_c \; \Theta(R_{*} - R)= \frac{3}{4\pi} \frac{M_*}{\mathcal{R}_*^3} \Theta(R_{*} - R)\,,

$$

where $\Theta(R _{*}- R)$ is the Heaviside function, which is unity inside the star
and zero outside. By substituting Eq. (Eq. rho-n0) into Eq. (Eq. I-Newton), one obtains

$$
I^\text{N} = \frac{2}{5} M_* \mathcal{R}_*^2\,, \quad \bar{I}^\text{N} = \frac{2}{5} \frac{1}{C^2}\,.

$$

Not surprisingly, this is the Newtonian moment of inertia for a sphere of constant density.

Next, we solve Eq. (Eq. h2RR-Newton) to obtain $\lambda^\mathrm{(tid)}{}^{,\text{N}}$. As pointed out in [damour-nagar], one must be careful when solving Eq. (Eq. h2RR-Newton) for the incompressible EoS. This is because $\rho$ can be expressed as a step-function with a discontinuity at the NS surface. Thus, $d\rho/dp$ in Eq. (Eq. h2RR-Newton) gives a delta-function centered at the NS surface. To be more precise, by using $d\rho/dR = - \rho_c \delta (R_*-R)$, $p(\mathcal{R}_*)=0$, $M_*=(4\pi/3) \mathcal{R}_*^3 \rho_c$ and the hydrostatic equilibrium equation [or the Newtonian limit of Eq. (Eq. TOV-zeroth)] at the NS surface, $dp(\mathcal{R}_*)/dR = -M_* \rho(\mathcal{R}_*)/\mathcal{R}_*^2$, the coefficient of $h_2$ in Eq. (Eq. h2RR-Newton) that is proportional to $d\rho/dp$ becomes $4\pi \rho d\rho/dp = 4\pi \rho (d\rho/dR) (dp/dR)^{-1} = (3/\mathcal{R}_*) \delta (R_{*}-R)$ near the surface. By taking such term into account, one can obtain the correct $k_2^\mathrm{(tid)}{}^{,\text{N}}$ by first solving Eq. (Eq. h2RR-Newton) with $d\rho/dp=0$ and then shifting $y^\text{N}$ by -3 [damour-nagar]. By taking all of this into account, one obtains $y^\text{N} = -1$, and thus, $k_2^\mathrm{(tid)}{}^{,\text{N}} = 3/4$ [damour-nagar]. This agrees with the classic result in [brooker-olle]. By using Eq. (Eq. tid-Love-def), $\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}}$ becomes

$$
\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} = \frac{1}{2}\frac{1}{C^5}\,.

$$

Let us now move on to the rotationally-induced quadruple moment. Rotating configurations of constant density stars can be described by Maclaurin spheroids [shapiro-teukolsky]. The quadrupole moment in Newtonian theory is given by [laarakkers]

$$
Q^\mathrm{(rot)}{}^{,\text{N}} = 2\pi \int_0^{\pi} \int_0^{r_* (\theta)} \rho(r,\theta) r^4 P_2 (\cos \theta) \sin \theta \ dr \ d\theta\,.

$$

The surface of the star $r=r_* (\theta)$ is in turn given by

$$
r_* (\theta) = \left( \frac{\sin^2\theta}{b^2} + \frac{\cos^2\theta}{c^2} \right)^{-1/2}\,,

$$

where $b$ and $c$ are the semi-major and semi-minor axes, respectively. By substituting Eq. (Eq. RNS-n0) and $\rho=\rho_c$ into Eq. (Eq. Q-integ), one obtains

$$
Q^\mathrm{(rot)}{}^{,\text{N}} = -\frac{4\pi}{15} \rho_c b^2 c (b^2-c^2)\,.

$$

From the equation of hydrostatic equilibrium,

$$
\frac{d\boldsymbol{v}}{dt} = - \frac{1}{\rho} \boldsymbol{\nabla}p - \boldsymbol{\nabla} \Psi\,,
$$

where boldfaced quantities refer to three-dimensional Euclidean vectors,
with $\boldsymbol{v} = \boldsymbol{\Omega} \times \boldsymbol{r}$ and $\boldsymbol{\Omega}$ the NS angular velocity vector, one obtains [shapiro-teukolsky]

$$\begin{aligned}
\Omega_* &=& \left\{ 2\pi \rho_c \left[ \frac{\sqrt{1-e^2}(3-2e^2)}{e^3} \sin^{-1}e - \frac{3(1-e^2)}{e^2} \right] \right\}^{1/2}  \\
&=& \sqrt{\frac{8\pi}{15}\rho_c} \ e + \mathcal{O}(e^3)\,,

\end{aligned}$$

where $\Omega_* = |\boldsymbol{\Omega}|$  and $e$ is the eccentricity defined by

$$
e \equiv \sqrt{1- \frac{c^2}{b^2}}\,.

$$

Using Eqs. (Eq. Omega-n0) and (Eq. e-n0), we can eliminate $c$ from Eq. (Eq. Q-n0) and substitute $b= \mathcal{R}_* + \mathcal{O}(\Omega_*^2)$ to obtain

$$
Q^\mathrm{(rot)}{}^{,\text{N}} = - \frac{1}{2} \mathcal{R}_*^5 \Omega_*^2 + \mathcal{O}(\Omega_*^4)\,.
$$

Keeping only the leading term, one obtains [laarakkers]

$$
\bar{Q}^{\text{N}} = \frac{25}{8} \frac{1}{C}\,.

$$

The dimensionless rotational Love number $\bar{\lambda}^\mathrm{(rot)}{}^{,\text{N}}$ can be calculated as

$$
\bar{\lambda}^\mathrm{(rot)}{}^{,\text{N}} = (\bar{I}^\text{N})^2 \bar{Q}^{\text{N}} = \frac{1}{2} \frac{1}{C^5}\,,
$$

which agrees with $\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}}$ given in Eq. (Eq. lambda-bar-n0). This then verifies Eq. (Eq. Love-Love-Newton).

From Eqs. (Eq. I-bar-n0), (Eq. lambda-bar-n0) and (Eq. Q-bar-n0), one obtains the *I-Love-Q* relations in the Newtonian limit for the incompressible EoS:

$$\begin{aligned}
\bar{I}^\text{N} &=& C_{\bar{I} \bar{\lambda}}^{(n=0)} \left[ \bar{\lambda}^\mathrm{(rot)}{}^{,\text{N}} \right]^{2/5}\,, \quad \bar{I}^\text{N} = C_{\bar{I} \bar{Q}}^{(n=0)} \left[ \bar{Q}^{\text{N}} \right]^{2}\,,  \\
\bar{Q}^{\text{N}} &=& C_{\bar{Q} \bar{\lambda}}^{(n=0)} \left[ \bar{\lambda}^\mathrm{(rot)}{}^{,\text{N}} \right]^{1/5}\,,

\end{aligned}$$

with

$$\begin{aligned}

C_{\bar{I} \bar{\lambda}}^{(n=0)} &=& \frac{2^{7/5}}{5} \approx 0.528\,, \\
C_{\bar{I} \bar{Q}}^{(n=0)} &=& \frac{128}{3125} \approx 0.0410\,, \\

C_{\bar{Q} \bar{\lambda}}^{(n=0)} &=& \frac{25}{2^{14/5}} \approx 3.59\,.
\end{aligned}$$

Of course, universality would be established if the constants $C_{A}$ are independent of the EoS, with $A$ any pair in the I-Love-Q trio. We will compute the same relations for the $n=1$ polytrope next, and thus, we will verify the degree of universality quantitatively.

### Polytropic I-Love-Q relations: $n=1$

Let us now concentrate on the $n=1$ polytrope and look first at the moment of inertia.
Equation (Eq. TOV-zeroth) in the Newtonian limit gives the equation of hydrostatic equilibrium:

$$
\frac{dp}{dR} = - \frac{\rho M}{R^2}\,.

$$

From Eqs. (Eq. tt-zeroth), (Eq. TOV-zeroth-Newton) and $p=K \rho^2$, one can solve this equation to obtain

$$
\rho = \frac{1}{4} \frac{M_*}{\mathcal{R}_*^2} \frac{1}{R} \sin \left( \frac{\pi R}{\mathcal{R}_*} \right)\,.
$$

One can then calculate $I^\text{N}$ (and $\bar{I}^\text{N}$) by substituting the above equation in Eq. (Eq. I-Newton) to find

$$
I^\text{N} = \frac{2 (\pi^2-6)}{3\pi^2} M_* \mathcal{R}_*^2\,, \quad \bar{I}^\text{N} = \frac{2 (\pi^2-6)}{3\pi^2} \frac{1}{C^2}\,.

$$

Let us now consider the tidal apsidal constant. The solution to Eq. (Eq. h2RR-Newton) can be written in terms of Bessel functions,
as $h_2^\text{N} \propto (R/\mathcal{R}_*)^{-1/2} J_{5/2} (\pi R/\mathcal{R}_*)$ [hinderer-love,damour-nagar]. With this, we find that the apsidal constant is

$$
k_2^\mathrm{(tid)}{}^{,\text{N}} = -\frac{1}{2} + \frac{15}{2\pi^2}\,.
$$

This constant agrees with the numerical results of [brooker-olle]. The dimensionless tidal Love number, $\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}}$,
is then

$$
\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} = \frac{15-\pi^2}{3\pi^2} \frac{1}{C^5}\,.

$$

Finally, let us look at the rotationally-induced quadrupole moment. From Eq. (Eq. Love-Love-Newton), one easily finds that

$$
\bar{Q}^{\text{N}} = \frac{\bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}}}{(\bar{I}^\text{N})^2} = \frac{3\pi^2 (15-\pi^2)}{4(\pi^2-6)^2} \frac{1}{C}\,.

$$

We now have all the necessary ingredients to compute the I-Love-Q relations in the Newtonian limit for an
$n=1$ polytrope. From Eqs. (Eq. I-bar-n1), (Eq. lambda-bar-n1) and (Eq. Q-bar-n1), one finds

$$\begin{aligned}
\bar{I}^\text{N} &=& C_{\bar{I} \bar{\lambda}}^{(n=1)} \left[ \bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} \right]^{2/5}\,, \quad \bar{I}^\text{N} = C_{\bar{I} \bar{Q}}^{(n=1)} \left[\bar{Q}^{\text{N}} \right]^{2}\,,  \\
\bar{Q}^{\text{N}} &=& C_{\bar{Q} \bar{\lambda}}^{(n=1)} \left[ \bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} \right]^{1/5}\,,

\end{aligned}$$

with

$$\begin{aligned}

C_{\bar{I} \bar{\lambda}}^{(n=1)} &=& \frac{ 2(\pi^2 -6)}{3^{3/5} \pi^{6/5} (15-\pi^2)^{2/5}} \approx 0.527\,, \\
C_{\bar{I} \bar{Q}}^{(n=1)} &=& \frac{32 (\pi^2-6)^5}{27 \pi^6 (\pi^2 -15)^2} \approx 0.0406\,, \\

C_{\bar{Q} \bar{\lambda}}^{(n=1)} &=& \frac{3^{6/5} \pi^{12/5} (15-\pi^2)^{4/5}}{4(\pi^2-6)^2} \approx 3.60\,.
\end{aligned}$$

Observe that the numbers shown in Eqs. (Eq. coeff-I-Love-n1)--(Eq. coeff-Q-Love-n1) are almost identical to those in Eqs. (Eq. coeff-I-Love-n0)--(Eq. coeff-Q-Love-n0).

The I-Love-Q relations for the $n=0$ and $n=1$ polytropic EoS in the Newtonian limit are shown as black thin solid lines in Fig. [ref:fig:I-Love-Q]. Notice that the relations for the $n=0$ and $n=1$ polytropic EoSs in GR approach the Newtonian ones as the compactness decreases. Notice, however, that we have here analytically shown that the I-Love-Q relations are very similar for the $n=0$ and $n=1$ polytropic EoSs *only*. This does not mean that the dependence of the I-Love-Q relations on the EoSs is weak in the Newtonian limit for all EoSs. Indeed, Fig. [ref:fig:I-Love-Q] shows that these relations for some EoSs, such as APR and SLy,  do not approach the Newtonian limit of the $n=0$ and 1 polytrope.

### Analytical Reasoning

The universality of the I-Love-Q relations rests on two ingredients. The first ingredient is that the functional form of the relations must
be the same for different EoSs. For the $n=0$ and $n=1$ polytropes in the Newtonian limit, this is verified by comparing Eqs. (Eq. I-Love-Q-n0) to (Eq. I-Love-Q-n1), and noting that regardless of the EoS, $\bar{I}^\text{N} \propto  \left[ \bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} \right]^{2/5}$, $\bar{I}^\text{N} \propto \left[ \bar{Q}^{\text{N}} \right]^{2}$ and $\bar{Q}^{\text{N}} \propto \left[ \bar{\lambda}^\mathrm{(tid)}{}^{,\text{N}} \right]^{1/5}$. This fact is perhaps expected, since all multipole moments must be proportional to the product of a dimensionless constant and the compactness to some power. The power is determined by the Newtonian dimensional structure of the particular multipole moment, e.g. $I \propto \mathcal{R}_*^{2}$ and thus $\bar{I} \propto C^{-2}$. This power will be the same regardless of the EoS, and thus the power exponent in the I-Love-Q relations will also be EoS independent.

The second ingredient, and perhaps the most difficult to understand, is the requirement that the constants of proportional  (ie. the $C_{A}$’s) be the same regardless of the EoS. For the $n=0$ and $n=1$ polytropes in the Newtonian limit, this is again verified by noting that the coefficients in Eqs. (Eq. coeff-I-Love-n0)--(Eq. coeff-Q-Love-n0) are almost identical to those in Eqs. (Eq. coeff-I-Love-n1)--(Eq. coeff-Q-Love-n1); their ratios are

$$\begin{aligned}
\frac{C_{\bar{I} \bar{\lambda}}^{(n=0)}}{C_{\bar{I} \bar{\lambda}}^{(n=1)}} &=& \frac{2^{2/5} 3^{3/5} \pi^{6/5} (15-\pi^2)^{2/5}}{5 \pi^2-30} \approx 1.002\,, \\
\frac{C_{\bar{I} \bar{Q}}^{(n=0)}}{C_{\bar{I} \bar{Q}}^{(n=1)}} &=& \frac{108 \pi^6 (\pi^2-15)^2}{3125 (\pi^2-6)^5} \approx 1.008\,, \\
\frac{C_{\bar{Q} \bar{\lambda}}^{(n=0)}}{C_{\bar{Q} \bar{\lambda}}^{(n=1)}} &=& \frac{25 (\pi^2-6)^2}{ 2^{4/5} 3^{6/5} \pi^{12/5} (15-\pi^2)^{4/5}} \approx 0.997\,.
\end{aligned}$$

In principle, there is no reason to expect that these coefficients should be equal to each other regardless of the EoS. Rather, one expects them to depend on the NS internal structure. One possible explanation is to argue that these coefficients depend on integrals of the energy density that are more heavily weighted toward the NS’s outer layers, i.e. the structure of the NS in its outer layers is what is mostly determining these coefficients. But it is precisely in the outer layers that nuclear physics uncertainties are lowest. Therefore, the EoSs in this regime are more similar to each other than in the core, thus leading to some degree of universality.

We have found some evidence to support this interpretation, shown in Fig. [ref:fig:I-Love-Q]. Focus on the I-Love-Q relations for the polytropic EoSs with $n=2$, $2.5$ and $3$, which greatly modify the internal structure in the NS’s outer layers, far from the core. In fact, these polytropes lead to essentially no energy density near the NS surface, with most of it concentrated near the core. We see that, indeed, when we choose an EoS that affects the NS structure far from the core, we significantly lose universality in the I-Love-Q relations, as one can see from the $n=2$, $2.5$ and $3$ curves in the bottom part of each panel in Fig. [ref:fig:I-Love-Q].

![](drhodp-PRD.eps)

Further evidence can be found by investigating a few of the terms that control the behavior of the moment of inertia, the quadrupole moment and the tidal Love number as a function of $C$. First, in the Newtonian limit, both $I^\text{N}$ [Eq. (Eq. I)] and $Q^\mathrm{(rot)}{}^{,\text{N}}$ [Eq. (Eq. Q-integ)] can be written in integral form, where the radial dependence of the integrand is proportional to $\rho(R) R^4$. The top panel of Fig. [ref:fig:drhodp] shows  $(\rho/\rho_c) (R/\mathcal{R}_*)^4$ as a function of  $R/\mathcal{R}_*$ for $C=0.17$ with various EoSs. $I$ and $Q^\mathrm{(rot)}$ are proportional to the area under the curves in this panel. Observe that the curves are similar, even for different realistic EoSs, and the dominant contribution comes from the NS outer layer, somewhere between $R/\mathcal{R}_* \approx 0.7-0.9$. This may explain the similar behavior of the I-C and Q-C curves with different realistic EoSs in the high-compactness regime (see Figs. [ref:fig:IbarMC] and [ref:fig:QbarMC]), as well as the universal I-Q behavior. Observe also that as one increase the polytropic index $n$, the NS becomes more centrally-condensed, and $I$ and $Q$ are not dominated just by the NS outer layers.

Similarly, we can study the behavior of the structure-dependent term that determines the tidal Love number in the Newtonian limit. Equation (Eq. h2RR-Newton) shows that there is only one such term and it is proportional to $\rho d\rho/dp = \rho (d\rho/dR) (dp/dR)^{-1}$. The bottom panel of Fig. [ref:fig:drhodp] plots $\rho (d\rho/dp)/(\rho_c^2 /p_c)$ as a function of $R/\mathcal{R}_*$ for $C=0.17$ with various EoSs. Sudden changes in the slope corresponds to nuclear phase transitions. Similar to the top panel, the behavior of  $\rho d\rho/dp$ is similar among realistic EoSs and the dominant contribution comes from the NS outer layers. This partially explains the similar behavior observed in the Love-C curves with different realistic EoSs in the high-compactness regime (see Fig. [ref:fig:lambdaMC]), as well as the universal Q-Love and I-Love relations.

Another possible explanation for the universality of the I-Love-Q relations involves the behavior of this trio as one approaches the BH limit. The no-hair theorems [robinson,israel,israel2,hawking-uniqueness0,hawking-uniqueness,carter-uniqueness] of GR state that the exterior multipolar structure of an isolated, stationary, axisymmetric BH solution in GR is completely determined by its mass and its spin angular momentum. Therefore, the quadrupole moment, for example, is completely determined by the spin angular momentum through Eq. (Eq. quadrupole) with $A = 0$ [geroch,hansen]. For NSs, such a result does not exist, but one might still expect the I-Q relation to become less structure dependent as one approaches the BH limit ($C \to 0.5$). Indeed, Fig. [ref:fig:I-Love-Q] shows that the loss of universality (the relative fractional difference shown in the bottom part of each panel) decreases logarithmically as the mass increases from $0.3 M_{\odot}$ to $1.4 M_{\odot}$.

Of course, one can never increase the compactness enough by a finite amount to turn a NS into a BH, ie. the NS sequence of varying compactness does not terminate in a BH for finite central density.  Still, it is interesting to see that as $C$ increases, the I-Love-Q relations are indeed approaching the BH limit, as explicitly shown in Fig. [ref:fig:I-Love-Q]. Moreover, universality in the I-Q curve suggests a universal relation between the NS spin and the NS quadrupole moment that is almost independent of the internal structure. Such a relation is similar to the no-hair relations for BHs [geroch,hansen].

Before proceeding, let us point out that this {*effacing of internal structure*} is not the same as what is discussed in the effacement principle [damour-effacement] in GR. The latter states that the equations of motion of compact objects of any size and structure depend only on integral parameters, like the mass and spin, and it is independent of its actual shape and internal structure. Of course, this principle holds in GR but only for BHs because of the no-hair theorems. The effacement principle is violated for NSs, but the violation is small, with corrections to the acceleration entering at 5PN order for a binary of non-spinning compact objects. Since the effacement principle deals with the motion of a body only, and not on the multipolar structure of its exterior gravitational field, the effacement we find here is not a consequence of the standard effacement principle.