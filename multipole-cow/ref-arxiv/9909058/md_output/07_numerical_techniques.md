# Numerical Techniques

For each candidate source of gravitational waves, gravitational wave
astronomy needs answers to two questions. Firstly, how much energy
will be carried away by the emitted gravitational waves? Secondly,
which are the “preferred” frequencies at which a 10 $M_\odot$ black
hole or a neutron star will oscillate? The answer to the first
question is that the energy will depend on the degree of asymmetry
that the process generates, and it will depend critically on the
initial data. In the previous section we have tried to provide some
guesstimates for the energy emitted during the oscillation phase of
black holes and neutron stars. The answer to the second question is
related to the numerical solution of the perturbation equations. The
numerical schemes developed for this purpose will be described in this
section.

Let us describe why the numerical calculation of quasi-normal mode
frequencies is delicate. Consider again the case treated in
section [ref:section_2] of the wave equation with a potential with
compact support. We try to find a complex number $s$ with negative
real part such that the solution which is $e^{-sx}$ for large positive
$x$, is $e^{sx}$ for large negative $x$. Note that these solutions
grow exponentially with $|x|$ and therefore one has to be very careful
to make sure that there is no exponentially decaying part in the
solution. The situation becomes even more complicated if we do not
know $f^{\pm}$ explicitly because one can not characterize the correct
solution by some growth property. This is for example the case for the
Schwarzschild solution.

## Black Holes

We would like to point out here that the attempts to calculate the QNM
frequencies date back to the beginning of 70s. More specifically in
their study of the black hole oscillations excited by an infalling
particle Davis et al. [DRPP71] found that the peak of the
spectrum is (for a Schwarzschild black hole) at around $M\omega =
0.32$ (geometrical units). This number is very close to the one
calculated by much more accurate methods later. In what will follow we
will describe the various methods used to calculate the QNM
frequencies.

### Evolving the Time Dependent Wave Equation

This approach was actually the first one used to study the QNM
excitation by Vishveshwara [vishu] but it has only been recently
revived thanks to increased power of computers. In general, one does
not need to Fourier decompose the perturbed Einstein equations but
instead evolve them for given sets of initial data and at the end to
Fourier analyze the resulting waveform. This procedure has certain
advantages since one does not need to be so careful in considering the
appropriate boundary conditions on the horizon, at infinity, on the surface or
at the center of the star. This does not mean that these boundaries are not
important for the evolution schemes. The difference is that in the
time independent case one formulates a boundary value problem, and the
eigenvalues (quasi-normal frequencies) depend critically on the correct
conditions on the various boundaries. A major disadvantage of the
evolution schemes is that one cannot get the complete spectrum of the
QNMs neither for a star nor for a black hole. The reason is that
although any perturbation is the sum of the harmonics involved, in
practice only a few of them will be observed and in the best case one
can succeed in getting a few extra modes by “playing” around with
the initial data.

To be more specific, by evolving a perturbation on a black hole
background one can get a QNM signal as the one shown in figure [ref:bhwave], but
in this signal the Fourier transform will show that there are present
at most two frequencies (the slowest damped ones); then by doing
“matched filtering” of this signal we will get the right frequencies
and damping times, but then all the extra modes of the spectrum are
missing! See for example the work by Bachelot and
Motet-Bachelot [BB91, BB92], and Krivan et al [KLP96,
  KLPA97]. This is true also for stars; in recent evolutions of axial
perturbations of stellar backgrounds Andersson and
Kokkotas [ak96] saw only a few of the $w$-modes (the ones that
damped slowest), while in similar calculations for even parity stellar
perturbations Allen et al. [aaks] have seen only the $f$-mode,
$2-3$ $p$-modes and two of the $w$-modes.

Of course with more detailed studies for various sets of initial data
one might be successful to get a few more modes, but more important
than deriving extra modes is understanding the physical situation
which generates the appropriate initial data.

Finally, the evolutions of the time dependent perturbation equations
can be extremely useful (and probably will be the only way) for the
calculation of the QNM frequencies and waveforms for the perturbations
of the Kerr-Newman black hole and for slowly and fast rotating
relativistic stars.

### Integration of the Time Independent Wave Equation

This technique was used by Chandrasekhar and Detweiler [CD75] and
is based on the definition of QNMs given in
section [ref:section_2]. They assumed that a QNM is a solution
corresponding to incoming waves on the horizon and outgoing at
infinity. Then by taking a series expansion of the Zerilli wave
equation ([ref:qnmwave], [ref:rwpot], [ref:zerpot]) at both limits (horizon and
infinity) of the form given by ([ref:bconditions]) they found initial
values for the numerical integration of the equation. Their
integration goes from both limits towards a common point which was set
close to the peak of the potential i.e.\ around $r = 3M$. The values
of $\omega$ for which the Wronskian of the two numerically taken
solutions vanishes are the quasi�eigensolutions of the problem. In
this way they managed to calculate the first $2-3$ QNM frequencies of
the Schwarzschild black hole for various harmonic indices. The
accuracy of the method improves for increasing $\ell$. Later,
Gunter [Gunter80] and Kokkotas with Schutz [KS88] used the
same approach to calculate the QNMs of the Reissner-Nordström black
hole.

The approach used by Nollert and Schmidt [Nollert90, NS92] is
more elaborate and based on a better estimate of the values of the
quasi-eigenfunctions on both boundaries ($\pm \infty$); this leads to
a more accurate estimate of frequencies and one also finds frequencies
which damp extremely fast. Andersson [nils92] suggested an
alternative integration scheme. The key idea is to separate ingoing
and outgoing wave solutions by numerically calculating their analytic
continuations to a place in the complex $r$-coordinate plane where
they have comparable amplitudes. This method is extremely accurate.

### WKB Methods

This technique, originally in the form suggested by Schutz and
Will [SW85], based on elementary quantum mechanical arguments,
was later developed into a powerful technique with which accurate
results have been derived. The idea is that one can reduce the QNM
problem into the standard WKB treatment of scattering of waves on the
peak of the potential barrier. The simplest way to find the QNM
frequencies is to use the well known Bohr-Sommerfeld (BS) rule. Using
this rule it is possible to reproduce not only the Schutz-Will formula
([ref:qnmsw]) but also to give a way to extend the accuracy of that
formula by taking higher order terms [IW87,GWKS90]. The classical form of the BS rule
for equations like ([ref:qnmwave]) is

$$
  \int_{r_A}^{r_B}\left[\omega^2 - V(r)\right]^{1/2} d r= (n+{1\over 2})\pi.

$$

where $r_A, r_B$ are the two roots (turning points) of $\omega^2
-V(r)=0$. A more general form can be found [Dunh32, BO78] which
is valid for complex potentials. This form can be extended to the
complex $r$ plane where the contour encircles the two turning points
which are connected by a branch cut. In this way one can calculate the
eigenfrequencies even in the case of complex potentials, as it is the
case for Kerr black holes.

This method has been used for the calculation of the eigenfrequencies
of the Schwarzschild [Iyer87], Reissner-Nordström [KS88],
Kerr [SI90, Kokkotas91] and  Kerr-Newman [KDK93] black holes
(restricted case). In general with this approach one can calculate
quite accurately the low-lying (relatively small imaginary part) QNM
modes, but it fails to give accurate results for higher-order modes.

This WKB approach was improved considerably when the phase integral
formalism of Fröman and Fröman [FF65] was introduced. In a
series of papers [FFAH92, AL92, AAS93a, ANS93] the method was
developed and a great number of even extremely fast damped QNMs of the
Schwarzschild black hole have been calculated with a remarkable
accuracy. The application of the method for the calculation of QNM
frequencies of the Reissner-Nordström black hole [AAS94] has
considerably improved earlier results [KS88] for the QNMs which
damp very fast. Close to the logic of this WKB approach were the
attempts of Blome, Mashhoon and Ferrari [BM84, FM84a, FM84b] to
calculate the QNM modes using an inversion of the black hole
potential. Their method was not very accurate but stimulated future
work using semi�analytic methods for estimating QNMs.

### The Method of Continued Fractions

In 1985, Leaver [Leaver85] presented a very accurate method for
calculating the QNM frequencies. His  method can be applied to the
calculation of the QNMs of Schwarzschild, Kerr and, with some
modifications, Reissner-Nordström black holes [Leaver90]. This
method is very accurate also for the high-order modes.

His approach was analogous to the determination of the eigenvalues of
the H$^+$ ion developed in [BH35]. A series representation of the
solution $f^-$ is assumed to represent also $f^+$ for the value of the
quasi-normal mode frequency. For normal modes the method may work
because $f^+$ is certainly bounded at infinity. In the case of
quasi-normal modes this is not so clear because $f^+$ grows
exponentially. Nevertheless, the method works very well numerically
and was improved by Nollert [Nollert93] such that he was able to
calculate very high mode numbers (up to 100,000!). In this way he
obtained the asymptotic distributions of modes described in
([ref:asymqnm]). An alternative way of using the recurrence relations
was suggested in [MP89].

Nollert [Nollert90] explains in his PhD thesis why the method
works. As initiated by Heisenberg et al. [Heisen] he considers
potentials depending analytically on a parameter $\lambda$ such that
for $\lambda = 1$ the potential has bound states -- normal modes --
and for for $\lambda = -1$ just quasi-normal modes. This is, for
example,  the case if we multiply the Regge-Wheeler potential
([ref:rwpot]) by $-\lambda$. Assuming that the modes depend
continuously on $\lambda$, one can try to relate normal modes to
quasi-normal modes and their methods of calculation.

In the case of QNMs of the Kerr black hole, one has to deal in
practice with two coupled equations, one which governs the radial part
([ref:teu1b]) and another which governs the angular dependence of the
perturbation ([ref:teu1a]). For both of them one can construct
recurrence relations for the coefficients of the series expansion of
their solutions, and through them calculate the QNM frequencies.

For the case of the QNMs of the Reissner-Nordström black hole, the
asymptotic form of the solutions is similar to that shown in equation
([ref:exp_hor]) but the coefficients $a_n$ are determined via a four
term recurrence relation. This means that the nice properties of
convergence of the three term recurrence relations have been lost and
one should treat the problem with great caution. Nevertheless,
Leaver [Leaver90] has overcome this problem and showed how to
calculate the QNMs for this case.

As a final comment on this excellent method we should point out that
it has a disadvantage compared to the WKB based methods in that it is
a purely numerical method and it cannot provide much intuition about
the properties of the QNM spectrum.

## Relativistic Stars

Although the perturbation equations in the exterior of a star are
similar to those of the black-hole and the techniques described
earlier can be applied here as well, special attention must be given
to the interior of the star where the perturbation equations are more
complicated.

For the time independent case the system of
equations ([ref:starper1], [ref:starper2], [ref:starper3]) inside the star reduces to a
4th order system of ODEs [DL85]. One can then even treat it as
two coupled time independent wave equations[^1]. The first equation will correspond to the fluid and
the second equation will correspond to spacetime perturbations. In
this way one can easily work in the Cowling approximation (ignore the
spacetime perturbations) if the aim is the calculation of the QNM
frequencies of the fluid modes ($f$, $p$, $g$, \dots) or the Inverse
Cowling Approximation [AKS96] (ignore the fluid perturbations) if
the interest is in $w$-modes. The integration procedure inside the
star is similar to those used for Newtonian stars and involves
numerical integration of the equations from the center towards the
surface in such a way that the perturbation functions are regular at
the center of the star and the Lagrangian variation of the pressure is
zero on the surface (for more details refer to [lindblom83,
  kokkotas92]). The integrations inside the star should provide the
values of the perturbation functions on the surface of the star where
one has to match them with the perturbations of the spacetime
described by Zerilli’s equation ([ref:qnmwave], [ref:rwpot], [ref:zerpot]).

In principle the integrations of the wave equation outside the star
can be treated as in the case of the black holes. Leaver’s method of
continued fraction has been used in  [kojima88, leins93],
Andersson’s technique of integration on the complex $r$ plane was used
in [AKS95] while a simple but effective WKB approach was used by
Kokkotas and Schutz [kokkotas92, YEF94].

Finally, there are a number of additional approaches used in the past
which improved our understanding of stellar oscillations in GR. In the
following paragraphs they will be discussed briefly.

- \mathbf{Resonance Approach.} This method was developed by
  Thorne [Thorne69a], the basic assumption being that there are
  no incoming or outgoing waves at infinity, but instead standing
  waves. Then by searching for resonances one can identify the QNM
  frequencies. The damping times can be estimated from the half-width
  of each resonance. This is a simple method and can be used for the
  calculations of the fluid QNMs. In a similar fashion Chandrasekhar
  and Ferrari [chandra91a] have calculated the QNM frequencies
  from the poles of the ratio of the amplitudes of the ingoing and
  outgoing waves.
- \mathbf{Direct Numerical Integration.} This method was used by
  Lindblom and Detweiler [lindblom83] for the calculation of the
  frequencies and damping times of the $f$-modes for various stellar
  models for thirteen different equations of state. In this case,
  after integration of the perturbation equations inside the star, one
  gets initial data for the integration of the Zerilli equation
  outside. The numerical integration is extended up to “infinity”
  (i.e.\ at a distance where the solutions of Zerilli equations
  become approximately simple sinusoidal ingoing and outgoing waves),
  where this solution is matched with the asymptotic solutions of the
  Zerilli equations which describe ingoing and outgoing waves. The QNM
  frequencies are the ones for which the amplitude of the incoming
  waves is zero. This method is more accurate than the previous one at
  least in the calculation of the damping times as has been verified
  in [AK98], but still is not appropriate for the calculation of
  the $w$-modes.
- \mathbf{Variational Principle Approach.} Detweiler and
  Ipser [DI73] derived a variational principle for non-radial
  pulsational modes. Associated with that variational principle is a
  conservation law for the pulsational energy in the star. The time
  rate of change of that pulsational energy, as given by the
  variational principle, is equal to minus the power carried off by
  gravitational waves. This method was used widely for calculating
  the $f$ [Det75a] and $g$-modes [Finn86] and in studies of
  stability [DI73, Det75b].
- \mathbf{WKB.} This is a very simple method but quite accurate, and
  contrary to the previous three methods it can be used for the
  calculation of the QNM frequencies of the $w$-modes (this was the
  first method used for the derivation of these modes). In practice
  one substitutes the numerical solutions of the Zerilli equation with
  their approximate WKB solutions and identifies the QNM frequencies
  as the values of the frequency for which the amplitude of the
  incoming waves is zero [kokkotas92, YEF94].
- \mathbf{WKB-Numerical.} This is a combination of direct integration
  of the Zerilli equation and the WKB method. The trick is that
  instead of integrating outwards, one integrates inwards (using
  initial data at infinity for the incoming wave solution); this
  procedure is numerically more stable than the outward
  integration. On the stellar surface one needs of course not the
  solution for incoming waves but the one for outgoing
  ones. But through WKB one can derive approximately the value of the
  outgoing wave solution from the value of the incoming wave
  solution. In this way errors are introduced, nevertheless the
  results are quite accurate [kokkotas92].

[^1]: Chandrasekhar and
  Ferrari [chandra91a] have also reduced the time independent
  perturbation equations (using a different gauge) into a 5th order
  system which involves only the spacetime perturbations with the
  fluid perturbations being calculated via algebraic relations from
  the spacetime perturbations. It was later proven by Ipser and
  Price [IP91, PI91] that this system of ODEs can be reduced to
  the standard equations in the Regge-Wheeler gauge. That this is
  possible is apparent from equations
  ([ref:starper1], [ref:starper2], [ref:starper3], [ref:starconstr]) and it is discussed
  in [aaks]