# Quasi-Normal Modes of Black Holes

One of the most interesting aspects of gravitational wave detection
will be the connection with the existence of black
holes [Thorne98]. Although there are presently several indirect
ways of identifying black holes in the universe, gravitational waves
emitted by an oscillating black hole will carry a unique fingerprint
which would lead to the direct identification of their existence.

As we mentioned earlier, gravitational radiation from black hole
oscillations exhibits certain characteristic frequencies which are independent of the
processes giving rise to these oscillations. These “quasi-normal”
frequencies are directly connected to the parameters of the black hole
(mass, charge and angular momentum) and for stellar mass black holes
are expected to be inside the bandwidth of the constructed
gravitational wave detectors.

The perturbations of a Schwarzschild black hole reduce to a simple
wave equation which has been studied extensively. The wave equation
for the case of a Reissner-Nordström black hole is more or less
similar to the Schwarzschild case, but for Kerr one has to solve a
system of coupled wave equations (one for the radial part and one for
the angular part). For this reason the Kerr case has been studied less
thoroughly. Finally, in the case of Kerr-Newman black holes we face
the problem that the perturbations cannot be separated in their
angular and radial parts and thus apart from special
cases [KDK93] the problem has not been studied at all.

## Schwarzschild Black Holes

The study of perturbations of Schwarzschild black holes assumes a
small perturbation $h_{\mu\nu}$ on a static spherically symmetric
background metric

$$
  ds^2 =g_{\mu \nu}^0 dx^{\mu}dx^{\nu} = -e^{v(r)}dt^2 +
  e^{\lambda(r)}dr^2 + r^2 \left( d\theta^2 + \sin^2 \theta d\phi^2
  \right),

$$

with the perturbed metric having the form

$$
  g_{\mu\nu} = g_{\mu\nu}^0 + h_{\mu\nu},

$$

which leads to a variation of the Einstein equations i.e.

$$
  \delta G_{\mu\nu} = 4 \pi \delta T_{\mu\nu}.
$$

By assuming a decomposition into tensor spherical harmonics for each
$h_{\mu\nu}$ of the form

$$
  \chi (t,r,\theta, \phi) = \sum_{\ell m}{\chi_{\ell m}(r, t)\over r}
  Y_{\ell m}(\theta, \phi),

$$

the perturbation problem is reduced to a single wave equation, for
the function $\chi_{\ell m}(r,t)$ (which is a combination of the
various components of $h_{\mu\nu}$). It should be pointed out that
equation ([ref:harmonics]) is an expansion for scalar quantities only. From the 10 independent components of the $h_{\mu\nu}$ only
$h_{tt}$, $h_{tr}$, and $h_{rr}$ transform as scalars under rotations.
The $h_{t \theta}$, $h_{t\phi}$, $h_{r\theta}$, and $h_{r\phi}$
transform as components of two-vectors under rotations and can be
expanded in a series of vector spherical harmonics while the
components $h_{\theta\theta}$, $h_{\theta\phi}$, and $h_{\phi\phi}$
transform as components of a $2\times2$ tensor and can be expanded in
a series of tensor spherical harmonics (see [TC67, Zer70, Monc74c]
for details). There are two classes of vector spherical harmonics
(\mathbf{polar} and \mathbf{axial}) which are build out of combinations of
the Levi-Civita volume form and the gradient operator acting on the
scalar spherical harmonics. The difference between the two families is
their parity. Under the parity operator $\pi$ a spherical harmonic
with index $\ell$ transforms as $(-1)^\ell$, the polar class of
perturbations transform under parity in the same way, as $(-1)^\ell$,
and  the axial perturbations as $(-1)^{\ell+1}$[^1]. Finally, since we are
dealing with spherically symmetric spacetimes the solution will be
independent of $m$, thus this subscript can be omitted.

The radial component of a perturbation outside the event horizon
satisfies the following wave equation,

$$
  {{\partial^2}\over {\partial t^2}} \chi_\ell +
  \left( -{{\partial^2}\over {\partial r_*^2}} + V_\ell(r)
  \right)\chi_\ell = 0,

$$

where $r_*$ is the “tortoise” radial coordinate defined by

$$
  r_*=r+2M \log (r/2M-1),
$$

and $M$ is the mass of the black hole.

For “axial” perturbations
$$
  V_\ell(r)=\left(1 - {2M\over r}\right) \left[ {{\ell(\ell+1)}\over
  r^2} + {{2\sigma M} \over r^3} \right]

$$

is the effective potential or (as it is known in the literature)
Regge-Wheeler potential [RW57], which is a single potential
barrier with a peak around $r=3M$, which is the location of the
unstable photon orbit. The form ([ref:rwpot]) is true even if we
consider scalar or electromagnetic test fields as perturbations.
The parameter $\sigma$ takes the values 1 for scalar perturbations, 0
for electromagnetic perturbations, and $-3$ for gravitational
perturbations and can be expressed as $\sigma=1-s^2$, where $s=0, 1,
2$ is the spin of the perturbing field.

For “polar” perturbations the effective potential was derived by
Zerilli [Zer70] and has the form

$$
  V_\ell(r)=\left( 1 - \frac{2M}r \right) \frac{2n^2(n+1)r^3 + 6n^2Mr^2 +
  18nM^2r +18M^3}{r^3(nr+3M)^2},

$$

where

$$
  2 n= (\ell-1)(\ell+2).
$$

Chandrasekhar [Chandra75] has shown that one can transform the
equation ([ref:qnmwave]) for “axial” modes to the corresponding
one for “polar” modes via a transformation involving differential
operations. It can also be shown that both forms are connected to the
Bardeen-Press [BP73] perturbation equation derived via the
Newman-Penrose formalism. The potential $V_\ell(r_*)$ decays
exponentially near the horizon, $r_*\to-\infty$, and as $r_*^{-2}$ for
$r_*\to +\infty$.

From the form of equation ([ref:qnmwave]) it is evident that the study
of black hole perturbations will follow the footsteps of the theory
outlined in section [ref:section_2].

Kay and Wald [KW87] have shown that solutions with data of
compact support are bounded. Hence we know that the time independent
Green function $G(s,r_*,r_*’)$ is analytic for $Re(s)>0$. The
essential difficulty is now to obtain the solutions $f_\pm$
(cf. equation ([ref:bs1.9])) of the equation

$$
  s^2\hat{\chi} -\hat{\chi}”+ V \hat{\chi} = 0,

$$

(prime denotes differentiation with respect to $r_*$)
which satisfy for real, positive $s$:

$$
  f_+\sim e^{-s r_*} \quad \text{for }  r_*\to\infty,
  \quad \quad f_-\sim e^{+r_* x} \quad \text{for } r_*\to -\infty.

$$

To determine the quasi-normal modes we need the analytic continuations
of these functions.

As the horizon ($r_*\to\-\infty$) is a regular singular point of
([ref:wave1]), a representation of $f_-(r_*,s)$ as a converging series
exists. For $M={1\over 2}$ it reads:

$$
  f_-(r,s)=(r-1)^s\sum_{n=0}^\infty a_n(s)(r-1)^n.

$$

The series converges for all complex $s$ and $|r-1|<1$ [Persides].
(The analytic extension of $f_-$ is investigated in [JC85].)
The result is that $f_-$ has an extension to the complex $s$ plane
with  poles only at negative real integers. The representation of
$f_+$ is more complicated: Because infinity is a singular point no
power series expansion like ([ref:exp_hor]) exists. A representation
coming from the iteration of the defining integral equation is given
by Jensen and Candelas [JC85], see also [NS92]. It turns out
that the continuation of $f_+$ has a branch cut $Re(s)\leq 0$ due to
the decay $r^{-2}$ for large $r$ [JC85].

The most extensive mathematical investigation of quasi-normal modes of
the Schwarzschild solution is contained in the paper by Bachelot and
Motet-Bachelot [BB91]. Here the existence of an infinite number
of quasi-normal modes is demonstrated. Truncating the
potential ([ref:rwpot]) to make it of compact support leads to the
estimate ([ref:bs1.15]).

The decay of solutions in time is not exponential because of the weak
decay of the potential for large $r$. At late times, the quasi-normal
oscillations are swamped by the radiative tail [Price72a, Price72b].
This tail radiation is of interest in its own right since it originates
on the background spacetime. The first authoritative study of nearly
spherical collapse, exhibiting radiative tails, was performed by
Price [Price72a, Price72b].

Studying the behavior of a massless scalar field propagating on a
fixed Schwarzschild background, he showed that the field dies off with
the power-law tail,

$$
  \chi(r,t) \sim t^{-(2\ell+P+1)},
$$

at late times, where $P=1$ if the field is initially static, and $P=2$
otherwise. This behavior has been seen in various calculations, for
example the gravitational collapse simulations by Cunningham, Price
and Moncrief [CPMa, CPMb, CPMc]. Today it is apparent in any
simulation involving evolutions of various fields on a black hole
background including Schwarzschild,
Reissner-Nordström [GPP94a], and Kerr [KLP96, KLPA97].
It has also been observed in simulations of axial oscillations of
neutron stars [ak96], and should also be present for polar
oscillations. Leaver [Leaver86] has studied in detail these tails
and associated this power low tail with  the branch-cut integral
along the negative imaginary $\omega$ axis in the complex
$\omega$ plane. His suggestion that there will be radiative tails
observable at ${\cal J}^+$ and ${\cal H}^+$ has been verified by
Gundlach, Price, and Pullin [GPP94a]. Similar results were
arrived at recently by Ching et al. [suen2] in a more extensive
study of the late time behavior. In a nonlinear study Gundlach, Price,
and Pullin [GPP94b] have shown that tails develop even when the
collapsing field fails to produce a black hole.
Finally, for a study of tails in the presence of a cosmological constant
refer to [BCKL97], while for a recent study, using analytic methods,
of the late-time tails of linear scalar fields outside Schwarzschild and
Kerr black holes refer to [Barak99, BO99].

Using the properties of the waves at the horizon and infinity
given in equation ([ref:bconditions]) one can search for the
quasi-normal mode frequencies since practically the whole problem has
been reduced to a boundary value problem with $s=i\omega$ being the
complex eigenvalue. The procedure and techniques used to solve the
problem will be discussed later in section [ref:section_6], but it is
worth mentioning here a simple approach to calculate the QNM
frequencies proposed by Schutz and Will [SW85]. The approach is
based on the standard WKB treatment of wave scattering on the peak of
the potential barrier, and it can be easily shown that the complex
frequency can be estimated from the relation

$$
  (M\omega_n)^2 = V_\ell(r_0) - i \left(n+{1\over 2}\right) \left[-2
    {{d^2V_\ell(r_0)} \over {dr_*^2}}\right]^{1/2},

$$

where $r_0$ is the peak of the potential barrier. For $\ell=2$ and
$n=0$ (the fundamental mode) the complex frequency is $M\omega \approx
(0.37, -0.09)$, which for a $10 M_\odot$ black hole corresponds to a
frequency of $1.2$ kHz and damping time of $0.55$ ms. A few more QNM
frequencies for $\ell = 2, 3$ and 4 are listed in table [ref:table1].

| n | $\ell=2$ |  | $\ell=3$ |  | $\ell=4$ | \cr

      0 | 0.37367 | -0.08896 i | 0.59944 | -0.09270 i | 0.80918 | -0.09416 i |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.34671 | -0.27391 i | 0.58264 | -0.28130 i | 0.79663 | -0.28443 i |
| 2 | 0.30105 | -0.47828 i | 0.55168 | -0.47909 i | 0.77271 | -0.47991 i |
| 3 | 0.25150 | -0.70514 i | 0.51196 | -0.69034 i | 0.73984 | -0.68392 i |

Figure [ref:bhmodes] shows some of the modes of the Schwarzschild
black hole. The number of modes for each harmonic index $\ell$ is
infinite, as was mathematically proven by Bachelot and
Motet-Bachelot [BB91]. This was also implied in an earlier work
by Ferrari and Mashhoon [FM84b], and it has been seen in the
numerical calculations in [AL92, Nollert93]. It can be also seen
that the imaginary part of the frequency grows very quickly. This
means that the higher modes do not contribute significantly in the
emitted gravitational wave signal, and this is also true for the
higher $\ell$ modes (octapole etc.).
<!-- TODO: figure content -->
As is apparent in figure [ref:bhmodes] that there is a special
purely imaginary QNM frequency. The existence of “algebraically
special” solutions for perturbations of Schwarzschild,
Reissner-Nordström and Kerr black holes were first pointed out by
Chandrasekhar [Chandra84]. It is still questionable whether these
frequencies should be considered as QNMs [Leaver90] and there is
a suggestion that the potential might become transparent for these
frequencies [Nils94]. For a more detailed discussion refer
to [LM96].

As a final comment we should mention that as the order of the modes
increases the real part of the frequency remains constant, while the
imaginary part increases proportionally to the order of the mode.
Nollert [Nollert93] derived the following approximate formula for
the asymptotic behavior of QNMs of a Schwarzschild black hole,

$$
  M \omega_n \approx 0.0437 + {{\gamma_1}\over {(2n+1)^{1/2}}} + \dots -
  i \left[ -{1\over 8} (2n+1) + {{\gamma_1}\over {(2n+1)^{1/2}}} + \dots
  \right],

$$

where $\gamma_1 = 0.343$, $0.7545$ and $2.81$ for $\ell=$ 2, 3 and 6,
correspondingly, and $n \rightarrow \infty$. The above relation was
later verified in [nils93] and [Liu95].

For large values of $\ell$ the distribution of QNMs is given
by [press71, FM84a, FM84b, Iyer87]

$$
  3 \sqrt{3} M \omega_n \approx \ell+ {1\over2} - i \left( n+ {1\over
  2}\right).
$$
For a mathematical proof refer to [BZ97].

The perturbations of Reissner-Nordström black holes, due to the
spherical symmetry of the solution, follow the footsteps of the
analysis that we have presented in this section. Most of the work was
done during the seventies by
Zerilli [Zer74], Moncrief [Monc74a, Monc74b] and later by
Chandrasekhar and Xanthopoulos [Chandra80, Xanth81]. For an
extensive discussion refer to [Chandra83]. We have again wave
equations of the form ([ref:qnmwave]), one for each parity with
potentials which are like ([ref:rwpot]) and ([ref:zerpot]) plus extra
terms which relate to the charge of the black hole. An interesting
feature of the charged black holes is that any perturbation of the
gravitational (electromagnetic) field will also induce electromagnetic
(gravitational) perturbations. In other words, any perturbation of the
Reissner-Nordström spacetime will produce both electromagnetic and
gravitational radiation. Again it has been shown that the solutions
for the odd parity oscillations can be deduced from the solutions for
even parity oscillations and vice versa [Chandra80]. The QNM
frequencies of the  Reissner-Nordström black hole have been
calculated by Gunter [Gunter80], Kokkotas, and
Schutz [KS88], Leaver [Leaver90], Andersson [nils93b],
and lately for the nearly extreme case by Andersson and
Onozawa [AO96].

[^1]: In the
literature the {\em polar} perturbations are also called {\em
even-parity} because they are characterized by their behavior under
parity operations as discussed earlier, and in the same way the {\em
axial} perturbations are called {\em odd-parity}. We will stick to the
polar/axial terminology since there is a confusion with the
definition of the parity operation, the reason is that to most people,
the words “even” and “odd” imply that a mode transforms under
$\pi$ as $(-1)^{2n}$ or $(-1)^{2n+1}$ respectively (for $n$ some
integer). However only the polar modes with even $\ell$ have even
parity and only axial modes with even $\ell$ have odd parity. If
$\ell$ is odd, then polar modes have odd parity and axial modes have
even parity. Another terminology is to call the polar perturbations
{\em spheroidal} and the axial ones {\em toroidal}. This definition is
coming from the study of stellar pulsations in Newtonian theory and
represents the type of fluid motions that each type of perturbation
induces. Since we are dealing both with stars and black holes we will
stick to
 the polar/axial terminology.

## Kerr Black Holes

The Kerr metric represents an axisymmetric, black hole solution to the
source free Einstein equations. The metric in $(t,r,\theta,\varphi)$
coordinates is

$$\begin{aligned}
  ds^2&=&-\left(1-{2Mr\over \Sigma}\right)dt^2 -{4Mar\sin^2\theta\over
    \Sigma} dtd\varphi +{\Sigma\over \Delta} dr^2  + \Sigma d\theta^2
   \\
  &&+ \left(r^2+a^2 +{2Ma^2r\sin^2\theta\over \Sigma}\right)\sin^2\theta
  d\varphi^2,

\end{aligned}$$

with

$$
  \Delta=r^2-2Mr+a^2, \quad \quad \Sigma=r^2+a^2\cos^2\theta.
$$

$M$ is the mass and $0\leq a\leq M$ the rotational parameter of the
Kerr metric. The zeros of $\Delta$ are

$$
  r_\pm=M\pm (M^2-a^2)^{1/2},
$$

and determine the horizons. For $r_+\leq r<\infty$ the spacetime
admits locally a timelike Killing vector. In the ergosphere region

$$
  r_+\leq r<M+(m^2-a^2\cos^2\theta)^{1/2},
$$

the Killing vector $\partial / \partial t$ which is timelike at
infinity, becomes spacelike. The scalar wave equation for the Kerr
metric is

$$\begin{aligned}
  & &\left[ {(r^2+a^2)^2\over\Delta} - a^2\sin^2\theta \right]
  {\partial^2\chi\over\partial t^2}+{4Mar\over\Delta}
  {\partial^2\chi\over\partial t\partial\varphi}+\left[{a^2\over\Delta}
  - {1\over\sin^2\theta}\right] {\partial^2\chi\over\partial \varphi^2}
   \\
  &&- \Delta^{-\sigma} {\partial \over {\partial r}}
  \left(\Delta^{\sigma+1} {\partial \chi \over \partial r} \right)
  -{1\over\sin\theta}{\partial\over\partial\theta}
  \left(\sin\theta{\partial \chi \over\partial\theta}\right)
  - 2 \sigma \left[ {{a(r-M)}\over \Delta} + {i\cos\theta \over
  \sin^2\theta} \right]{\partial \chi \over \partial \varphi}  \\
  &&- 2 \sigma \left[ {M(r^2-a^2)\over \Delta} -r - i a \cos \theta\right]
  {\partial \chi \over \partial t} + \left( \sigma^2 \cot^2\theta -
  \sigma \right)\chi=0,

\end{aligned}$$

where
$\sigma = 0, \pm 1, \pm 2$ for scalar, electromagnetic or gravitational
perturbations, respectively. As the Kerr metric outside the  horizon
$(r>r_+)$ is globally hyperbolic, the Cauchy problem for the scalar
wave equation ([ref:teu1]) is well posed for data on any Cauchy surface.
However, the coefficient of $\partial^2\chi / \partial\varphi^2$
becomes negative in the ergosphere. This implies that the time
independent equation we obtain after the Fourier or Laplace
transformation is not elliptic!

For linear hyperbolic equations with time independent coefficients, we
know that solutions determined by  data with compact support are
bounded by $ce^{\gamma t}$, where $\gamma$ is independent of the
data. It is not known whether all such solutions are bounded in time,
i.e.\ whether they are stable.

Assuming harmonic time behavior $\chi=e^{i\omega
t}\hat{\chi}(r,\theta,\varphi)$, a separation in angular and radial
variables was found by Teukolsky [Teu72]:

$$
  \hat{\chi}(r,\theta,\varphi)=R(r,\omega) S(\theta,\omega)e^{im\varphi}.
$$

Note that in contrast to the case of spherical harmonics, the
separation is $\omega$-dependent. To be a solution of the wave
equation ([ref:teu1]), the functions $R$ and $S$ must satisfy

$$
  {1\over \sin\theta}{d \over d \theta}\left[\sin\theta {d S\over d
  \theta} \right] \!+\! \left[a^2\omega^2\cos^2\theta + 2a\omega \sigma
  \cos\theta - {m^2 \!\! + \! \sigma^2 \!\! + \! 2m \sigma
  \cos\theta \over \sin^2\theta} \! + \! E \right] S = 0,

$$

$$
  \Delta^{-\sigma}{d\over {d r}}\left[\Delta^{\sigma+1}{d R\over {d
  r}}\right]+ {1\over\Delta}\left[K^2 + 2 i \sigma(r-1)K -\Delta(4 i
  \sigma r\omega+\lambda)\right] R = 0,

$$

where $K = (r^2+a^2)\omega + am$, $\lambda = E + a^2\omega^2 +
2am\omega - \sigma(\sigma+1)$, and $E$ is the separation constant.

For each complex $a^2\omega^2$ and positive integer $m$,
equation ([ref:teu1a]) together with the boundary conditions of
regularity at the axis, determines a singular Sturm-Liouville
eigenvalue problem. It has solutions for eigenvalues
$E(\ell,m^2,a^2\omega^2)$, $|m| \leq \ell$. The eigenfunctions are the
spheroidal (oblate) harmonics $S_{\ell|m|}(\theta)$. They exist for
all complex $\omega^2$. For real $\omega^2$ the spheroidal harmonics
are complete in the sense that any function  of $z=\cos\theta$,
absolutely integrable over the interval $[-1, 1]$, can be expanded
into spheroidal harmonics of fixed $m$ [Seidel89]. Furthermore,
functions $A(\theta,\varphi)$ absolutely integrable over the sphere
can be expanded into

$$
  A(\theta,\varphi)=\sum_{\ell=0}^\infty\sum_{m=-\ell}^{+\ell}
  A(\ell,m,\omega)S_{\ell|m|}(\theta)e^{im\varphi}.
$$

For general complex $\omega^2$ such an expansion is not possible.
There is a countable number of “exceptional values” $\omega^2$
where no such expansion exists [MS54].

Let us pick one such solution $S_{\ell|m|}(\theta,\omega)$ and
consider some solutions $R(r,\omega)$ of ([ref:teu1b]) with the
corresponding $E(\ell,|m|,\omega)$. Then $R_{\ell|m|}(r,\omega)
S_{\ell|m|}e^{im\varphi}e^{-i\omega t}$ is a solution
of ([ref:teu1]). Is it possible to obtain “all” solutions by summing
over $\ell,m$ and integrating over $\omega$? For a solution in
spacetime for which a Fourier transform in time exists at any space
point (square integrable in time), we can expand the Fourier transform
in spheroidal harmonics because $\omega$ is real. The coefficient
$R(r,\omega,m)$ will solve equation ([ref:teu1b]). Unfortunately, we
only know that a solution determined by data of compact support is
exponentially bounded. Hence we can only perform a Laplace
transformation.

We proceed therefore as in section [ref:section_2]. Let
$\hat{\chi}(s,r,\theta,\varphi)$ be the Laplace transform of a
solution determined by data $\chi(t,r,\theta,\varphi)=0$ and
$\partial_t\chi(t,r,\theta,\varphi)=\rho$, while $\chi$ is analytic in
$s$ for real $s>\gamma\geq 0$, and has an analytic continuation onto
the half-plane $Re(s)\geq\gamma$. For real $s$ we can expand
$\hat{\chi}$ into a converging sum of spheroidal harmonics [MS54]

$$
  \hat{\chi}(s,r,\theta,\varphi) = \sum_{\ell=0}^\infty
  \sum_{m=-\ell}^{+\ell} R(\ell,m,s)
  S_{\ell|m|}(-s^2,\theta)e^{im\varphi}.

$$

$R(\ell,m,s)$ satisfies the radial equation ([ref:bs1.8])
with $i\omega=s$. This representation of $\hat{\chi}$ does, however,
not hold for all complex values in the half-plane on which it is
defined. Nevertheless is it true that for all values of $s$ which are
not exceptional an expansion of the form ([ref:expan1])
exists[^1].

To define quasi-normal modes we first have to define the correct Green
function of ([ref:teu1b]) which determines $R(\ell,m,s)$ from the data
$\rho$. As usual, this is done by prescribing decay for real $s$ for
two linearly independent solutions $^\pm R(\ell,m,s)$ and analytic
continuation. Out of the Green function for $R(\ell,m,s)$ and
$S_{\ell|m|}(-s^2,\theta)$ for real $s$ we can build the Green
function $G(s,r,r’,\theta,\theta’,\phi,\phi’)$ by a series
representation like ([ref:expan1]). Analytic continuation defines $G$
on the half plane $Re(s)>\gamma$. For non exceptional $s$ we have a
series representation. (We must define $G$ by this complicated
procedure because the partial differential operator corresponding
to ([ref:teu1a]), ([ref:teu1b]) is not elliptic in the ergosphere.)

Normal and quasi-normal modes appear as poles of the analytic
continuation of $G$. Normal modes are determined by poles with
$Re(s)>0$ and quasi-normal modes by $Re(s)<0$. Suppose all such values
are different from the exceptional values. Then we have always the
series expansion of the Green function near the poles and we see that
they appear as poles of the radial Green function of $R(\ell,m,s)$.

To relate the modes to the asymptotic behavior in time we study the
inverse Laplace transform and deform the integration path to include
the contributions of the poles. The decay in time is dominated either
by the normal mode with the largest (real) eigenfrequency or the
quasi-normal mode with the largest negative real part.

It is apparent that the calculation of the QNM frequencies of the
Kerr black hole is more involved than the Schwarzschild and
Reissner-Nordström cases. This is the reason that there have only
been a few attempts [Leaver85, SI90, Kokkotas91,Onoz97] in this
direction.

The quasi-normal mode frequencies of the Kerr-Newman black hole have
not yet been calculated, although they are more general than all other
types of perturbation. The reason is the complexity of the
perturbation equations and, in particular, their non-separability. This
can be understood through the following analysis of the perturbation
procedure. The equations governing a perturbing massless field of spin
$\sigma$ can be written as a set of $2\sigma+1$ wavelike equations in
which the various different helicity components of the perturbing
field are coupled not only with each other but also with the curvature
of the background space, all with four independent variables as
coordinates over the manifold. The standard problem is to decouple the
$2\sigma+1$ equations or at least some physically important subset of
them and then to separate the decoupled equations so as to obtain
ordinary differential equations which can be handled by one of the
previously stated methods. For a discussion and estimation of the QNM
frequencies in a restrictive case refer to [KDK93].

[^1]: The definition of $S_{\ell|m|}(-s^2,\theta)$ on the
  complex plane is made unique by fixing certain conventions about the
  branch cuts. The exceptional points are the beginnings of branch
  cuts. $R(l,m,s)$ as a solution of ([ref:teu1b]) is defined also at
  the exceptional points; just the series does not exist there.

## Stability and Completeness of Quasi-Normal Modes

From the normal modes one can learn a lot about stability. Take as an
example linear stellar oscillation within the framework of Newton’s
theory of gravity. As outlined  at the beginning of
section [ref:section_2], we have a sequence of normal modes $\omega_n$
with $\omega^2$ real. The general solution is a convergent linear
combination of the corresponding eigenfunctions. Hence all
eigenfunctions are bounded if  $\omega$ is purely imaginary,
i.e.\ $\omega^2<0$.

Therefore the spectrum contains all the information about
stability. To discuss stability for systems with quasi-normal modes,
let us consider a case like equation ([ref:bs1.5]) with the assumption
that $V$ is of compact support but not necessarily positive.

Data of compact support define solutions which grow at most
exponentially in time

$$
  |\phi(t,x)|< c e^{at},

$$

where $a$ is independent of the data. As outlined in
appendix [ref:appendix_1], eigenvalues necessarily have $s_n>0$ and
the eigenfunctions determine solutions growing exponentially in time. If
no eigenvalues exist, the solution can not grow
exponentially. Polynomial growth is still possible and related to the
properties of the Laplace transform of the Green function at $s=0$. As
the potential has compact support, the functions $f_\pm(s,x)$ are
analytic for all $s$. Hence, the Green function can at most have a
pole at $s=0$. A pole of order two and higher implies polynomial
growth in time. If the potential is positive, energy conservation
shows that the field can grow at most linearly in time and therefore
we can have at most a pole of order 2 at $s=0$.

If we define stability as boundedness in time for all solutions
with data of compact support, properties of quasi-normal modes can not
decide the stability issue. However, the appearance of a normal mode
proves instability. If the support of the potential is not compact
everything becomes more complicated. In particular, it is a non
trivial problem to obtain the behavior of the Green function at $s=0$.

In the case of the Schwarzschild black hole, stability is demonstrated
by Kay and Wald [KW87] who showed the boundedness of all solutions
with data of compact support.

The issue is more subtle for Kerr. There is a conserved energy, but
because of the ergoregion its integrand is not positive definite,
hence the conserved energy could be finite while the field still might
grow exponentially in parts of the spacetime. Papers by Press and
Teukolsky [PT73], Hartle and Wilkins [HW74], and
Stewart [Stewart] try to exclude the existence of an
exponentially growing normal mode. Their work makes the stability very
plausible but is not as conclusive as the Wald-Kay result. However
this is a delicate issue as we see if, for example, we multiply the
Regge-Wheeler potential by a factor $\epsilon$: For any $\epsilon > 0$
we obtain an infinite number of QNMs, for $\epsilon=0$, however there
is no QNM! Whiting [Whiting] has proven that there are no
exponentially growing modes, and in his proof he showed that the
growth of the  modes is at most linear. Recent numerical evolution
calculations [KLP96, KLPA97] for slowly and fast rotating Kerr
black holes pick up all the expected features (QNM ringing, tails) and
show no sign of exponential growth. It should be noted that the
massive scalar perturbations of Kerr are known to be
unstable [DDR76, ZE79, Det80]. These unstable modes are known to be
very slowly growing (with growth times similar to the age of universe).

Let us finally turn to the “completeness of QNMs”. A general
mathematical theorem (spectral theorem) implies that for systems like
strings or membranes the general solution can be expanded into a
converging sum of normal modes. A similar result can not be expected
for QNMs, the reasons are given  in section [ref:section_2]. There is,
however, the possibility that an infinite sum of the
form ([ref:bs1.14]) will be a representation of a solution for late
times. This property has been shown by Beyer [Horst] for the
Pöschl-Teller potential which has a similar form as the potential
on Schwarzschild ([ref:rwpot]). The main difference is its exponential
decay at both ends. In [NP99] Nollert and Price propose a
definition of completeness and show its adequateness for a particular
model problem. There are also systematic studies [suen1] about
the relation between the structure of the QNM’s of the Klein-Gordon
equation and the form of the potential. In these studies there is a
discussion on both the requirements for QNMs to form a
complete set and the definition of completeness.