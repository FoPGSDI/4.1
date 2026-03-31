# Introduction and Motivation

A key challenge of current astrophysical research is to obtain
information about the equation of state (EoS) of the ultra-dense
nuclear matter making up neutron stars (NSs). The observational
constraints on the internal structure of NSs are weak: the observed
range of NS masses is $M\sim 1.1 - 2.2 M_\odot$ [Lattimer],
and there is no current method to directly measure the radius. Some
estimates using data from X-ray spectroscopy exist, but those are
highly model-dependent  (e. g. [astroph]). Different
theoretical models for the NS internal structure predict, for a
neutron star of mass $M\sim 1.4 M_\odot$, a central density in the
range of $ \rho_c\sim 2-8 \times 10^{14}{\rm g}{\rm cm}^{-3}$ and a
radius in the range of $R\sim 7-16$km [Lattimer]. Potential
observations of pulsars rotating at frequencies above $1400$Hz could
be used to constrain the EoS if the pulsar’s mass could also be
measured (e. g. [Bejger]).

Direct and model-independent constraints on the EoS of NSs could be
obtained from gravitational wave observations. Coalescing binary
neutron stars are one of the most important sources for ground-based
gravitational wave detectors [Cutler:2002me]. LIGO
observations have established upper limits on the coalescence rate
per comoving volume [Abbott:2007xi], and at design sensitivity
LIGO II is expected to detect inspirals at a rate of $\sim 2/$day
[2004ApJ...601...L179].

In the early, low frequency part of the inspiral ($f\leq 100$Hz,
where $f$ is the gravitational wave frequency), the waveform’s phase
evolution is dominated by the point-mass dynamics and finite-size
effects are only a small correction. Toward the end of the inspiral
the internal degrees of freedom of the bodies start to appreciably
influence the signal, and there have been many investigations of how
well the EoS can be constrained using the last several orbits and
merger, including constraints from the gravitational wave energy
spectrum [faber] and from the NS tidal disruption signal for
NS-black hole binaries [vallisneri]. Several numerical
simulations of the hydrodynamics of NS-NS mergers have studied the
dependence of the gravitational wave spectrum on the radius and EoS
(see, e.g. [baumgarte] and references therein). However,
trying to extract EoS information from this late time regime
presents several difficulties: (i) the highly complex behavior
requires solving the full nonlinear equations of general relativity
together with relativistic hydrodynamics; (ii) the signal depends on
unknown quantities such as the spins and angular momentum
distribution inside the stars, and (iii) the signals from the
hydrodynamic merger are outside of LIGO’s most sensitive band.

During the early regime of the inspiral the signal is very clean and
the influence of tidal effects is only a small correction to the
waveform’s phase. However, signal detection is based on matched
filtering, i. e. integrating the measured waveform against
theoretical templates, where the requirement on the templates is
that the phasing remain accurate to $\sim 1$ cycle over the
inspiral. If the accumulated phase shift due to the tidal
corrections becomes of order unity or larger, it could corrupt the
detection of NS-NS signals or alternatively, detecting a phase
perturbation could give information about the NS structure. This has
motivated several analytical and numerical investigations of tidal
effects in NS binaries
[bc,ks,kochanek,ts,mw,shibata,gualteri,pons,berti]. The
influence of the internal structure on the gravitational wave phase
in this early regime of the inspiral is characterized by a single
parameter, namely the ratio $\lambda$ of the induced
 quadrupole to the perturbing tidal field. This ratio $\lambda$ is related to
the star’s tidal Love number $k_2$ by $k_2 = 3G \lambda R^{-5}/2$,
where $R$ is the star’s radius. [tidal] have shown that for an
inspiral of two non-spinning $1.4M_{\odot}$ NSs at a distance of 50
Mpc, LIGO II detectors will be able to constrain $\lambda$ to
$\lambda \leqslant 2.01 \times 10^{37} {\rm g}{\,}{\rm cm}^2{\rm
s}^2 $ with $90%$ confidence. This number is an upper limit on
$\lambda$ in the case that no tidal phase shift is observed. The
corresponding constraint on radius would be $R \leqslant 13.6
{\,}{\rm km} {\;}\left(15.3{\,} {\rm km} \right)$ for a $n=0.5$
$\left(n=1.0\right)$ fully relativistic polytrope, for
$1.4M_{\odot}$ NSs [tidal].

 Because neutron stars are compact
objects with strong internal gravity, their Love numbers could be
very different from those for Newtonian stars that have been
computed previously, e. g. by [bo].

 Knowledge of Love number values could also be useful for
comparing different numerical simulations of NS binary inspiral by
focusing on models with the same masses and values of $\lambda$.

In [tidal], the $l=2$ tidal Love numbers for fully
relativistic neutron star models with polytropic pressure-density
relation $P=K\rho^{1+1/n}$, where $K$ and $n$ are constants, were
computed for the first time. The present paper will give details of
this computation. Using polytropes allows us to explore a wide range
of stellar models, since most realistic models can be reasonably
approximated as a polytrope with an effective index in the range
$n\sim 0.5-1.0$ [Lattimer]. Our prescription for computing
$\lambda$ is valid for an arbitrary pressure-density relation and
not restricted to polytropes. In Sec. [ref:sec:lovedef], we start by
defining $\lambda$ in the fully relativistic context in terms of
coefficients in an asymptotic expansion of the metric in the star’s
local asymptotic rest frame and discuss the extent to which it is
uniquely defined. In Sec. [ref:sec:lovecalc], we discuss our method
of calculating $\lambda$, which is based on static linearized
perturbations of the equilibrium configuration in the Regge-Wheeler
gauge as in [tc]. Section [ref:sec:results] contains the
results of the numerical computations together with a discussion.
Unless otherwise specified, we use units in which $c=G=1$.