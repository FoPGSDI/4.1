# What can we learn from observations?

Our present understanding of neutron stars comes mainly from
X-ray and radio-timing observations.
These observations provide some insight into the
structure of these objects and the properties of supranuclear
matter. The most commonly and accurately observed parameter
is the  rotation period, and we know that radio pulsars can spin
very fast (the shortest observed period being the 1.56 ms of PSR
1937+21).
Another
basic observable, that can be obtained (in a few cases)
with some accuracy from todays observations,
is the mass of the neutron star.
As Finn \shortcite{finn} has shown, the observations of radio pulsars
indicate that $1.01 < {M/M_\odot} < 1.64$. Similarly, van der Kerkwijk
et al \shortcite{Kerk95} find that data for X-ray pulsars  indicate
$1.04 < {M/M_\odot} < 1.88$. The data used in these two studies
is actually consistent with (if one includes error bars)
$M<1.44 M_\odot$.
We now recall that the various
EOS that have been proposed by theoretical physicists
can be divided into two major categories: i) the “soft” EOS which
typically lead to neutron star models with maximum masses around
$1.4 M_\odot$ and radii usually smaller than 10 km, and ii)
  the “stiff” EOS
with the maximum values  $M\sim 1.8 M_\odot$ and $R\sim 15$ km
[ab77].
From this one can deduce that, even though the constraint put on the
neutron star mass by present observations seems strong,
it does not rule out many
of the proposed EOS.
In order to arrive at a more useful result we are
likely to need detailed observations also of the stellar radius.
Unfortunately, available data provide
little information about the radius. The recent observations of
quasiperiodic oscillations in low mass X-ray binaries indicate
that $R<6M$, but again this is not a severe constraint.
Although a number of attempts have been made, using either X-ray observations
[lewin] or the limiting spin period of neutron stars [fip],
to
put constraints on the mass-radius relation,
we do not yet have a method
which can provide the desired answer.

## A detection scenario

In view of this situation, any method that can be used to infer
neutron star parameters is a welcome addition.
Of specific interest may be the new possibilities
offered once gravitational-wave observations become reality.
An obvious question is to what extent one can solve the inverse
problem in gravitational-wave astronomy. In this paper we address
this issue for the case of waves from pulsating neutron stars.
In Appendix A we provide extended tables with frequencies and
damping times for the most relevant pulsation modes (as far as
gravitational waves are concerned) of various stellar models
created from a range of realistic EOS.
We will suppose that these modes can be detected by some
future generation of gravitational-wave detectors, and
investigate to what precision we can hope to calculate the
parameters of the source from observed data.

Let us suppose that a nearby supernova explodes, say
in the Local Group of galaxies, and is
followed by a core collapse that leads to the formation of
a compact object. As the dust from the collapse settles the
compact object pulsates wildly in its various oscillation modes,
generating a gravitational-wave signal which
is composed of an overlapping of different frequencies..
We will assume that the results of Allen et al \shortcite{aaks}
can be brought to bear on this situation, i.e. that  most of the energy
is radiated through the $f$-mode, a few $p$-modes and the first
$w$-mode. Our detector picks up this signal, and a subsequent
Fourier analysis of the data stream yields the frequencies
and the energy content in each mode.

The first question to be answered by the gravitational-wave
astronomer concerns what kind of compact object could produce the
detected signal. Is it a
black hole or a neutron star? The pulsation of  these  objects
lead to qualitatively similar gravitational waves, eg. exponentially
damped oscillations, but the
question should nevertheless be relatively easy to answer.
If more than one of the stellar pulsation modes is observed the answer
is clear, but even if we only observe only one single mode the two
cases should be easy to distinguish.
The fundamental (quadrupole)
quasinormal mode frequency of a Schwarzschild black
hole follows from
$$
f \approx 12 {\rm kHz } \left( {M_\odot \over M} \right) ,
$$
while the associated e-folding time is
$$
\tau \approx 0.05 {\rm ms } \left( {M \over M_\odot} \right) \ .
$$
That is, the oscillations of a 10 $M_\odot$ black hole lie in
the frequency range of the $f$-mode for a typical neutron star (see
Appendix A). But the two signals will differ greatly in
 the damping time, the e-folding time of
the black hole being nearly three orders of magnitude shorter than
that of the neutron star $f$-mode.

Having excluded the possibility that our signal came from a
black hole, we want to know  the mass and the radius of
the newly born neutron star. We also want to decide which
of the proposed EOS that best represents
this star. To address these questions we can use a set of empirical
relations deduced from the
data of  Appendix A: Relations that can be used to
estimate the mass,
the radius and the EOS of the neutron star with good precision.

## Empirical relations

Let us first consider the frequency of the $f$-mode.
It is well known that the characteristic time-scale of any dynamical
process
is related to the mean density of the mass involved
(see Misner, Thorne and Wheeler \shortcite{MTW} ch. 36.2).
This notion should be relevant for the
fluid oscillation modes of a star, and  we consequently expect that
$\omega_f \sim \bar{\rho}^{1/2}$. That is,
we should normalize the
$f$-mode frequency with the average density of the star.
The result of doing this is shown
in Figure [ref:ffw].
From this Figure it is apparent that the relation between
the $f$-mode frequencies and the mean density is almost linear,
and a  linear fitting leads to the following simple relation:
$$
\omega_f (kHz)
\approx 0.78 + 1.635 \left(  {{\bar M} \over {{\bar R}^3}} \right)^{1/2}
\ ,

$$
where we have introduced the dimensionless variables
$$
{\bar M} = {M \over {1.4 M_\odot}} \quad \text{and} \quad
{\bar R} = {R \over 10 {\rm km} } \ .
$$
>From equation ([ref:rfw]) follows that the typical  $f$-mode frequency
is around 2.4 kHz.

<!-- TODO: figure content -->

To deduce a corresponding relation for the damping rate
of the $f$-mode, we can use the rough estimate given by the
 quadrupole formula. That is, the damping time should follow from
$$
\tau_f \sim {\text{oscillation energy} \over{ \text{power emitted in GWs}}}
\sim R\left({R\over M}\right)^3 \ .
$$
Using this normalization we plot the functional
$(\tau_f M^3/R^4)^{-1}$ as a function of the stellar compactness, cf.
Figure [ref:fftau].
The data shown in this Figure leads to a
relation between the damping time of the $f$-mode and the stellar parameters
$M$ and $R$:
$$
{1 \over \tau_f ({\rm s})} \approx {{\bar M}^3 \over {\bar R}^4}
\left[ 22.85 - 14.65 \left( {{\bar M}\over {\bar R}} \right)\right] \ .

$$
The small deviation of the numerical data
from the above formula is apparent in
Figure [ref:fftau], and one can easily see that a typical value for
the damping time of the $f$-mode is a tenth of a second.

<!-- TODO: figure content -->

<!-- TODO: figure content -->

For the damping rate of the $p$-modes the situation is not so
favorable. This is because the damping  of the $p$-modes
is more sensitive to changes in the modal distribution inside the
star. Thus, different EOS lead to rather different $p$-mode
damping rates, cf. Figure [ref:fptau].
Previous evidence for polytropes [akpoly] actually indicate that this
would be the case. Clearly, an empirical relation based on the data
in Figure [ref:fptau] would not be very robust.

The situation is slightly better if we consider the oscillation
frequency
of the $p$-mode. From the data graphed in
Figure [ref:fpw]
we can deduce a relation
between the $p$-mode frequency and the parameters of the star:
$$
\omega_p ({\rm kHz})
\approx {1\over {{\bar M}}}\left(1.75 + 5.59 {{\bar M} \over {{\bar R}}} \right)
\ ,
$$
and we see that a typical
 $p$-mode frequency is around 7 kHz. Although
the data for several EOS deviate significantly from ([ref:rpw])
it is still a useful result. Stellar masses and radii deduced from it
will not be as accurate as ones based on $f$-mode
data, but on the other hand, if
$M$ and $R$ are obtained in some other way (say, from a combination
of observed $f$- and $w$-modes) the $p$-mode can be used to
deduce the relevant EOS.

<!-- TODO: figure content -->

That empirical relations based on $p$-mode data would be less
robust and useful than those for the $f$-mode was expected, since
the $p$-modes are sensitive to changes in the
matter distribution inside the star. In contrast, the
gravitational-wave $w$-modes  should lead to very robust results.
It is well known
[ks92,akk96] that the $w$-modes do not excite a significant
fluid motion. Thus, they are
more or less independent of the characteristics of the fluid:
The frequencies do not depend on the sound speed and the damping
times cannot be modeled by the quadrupole formula.
But we can nevertheless deduce the appropriate
normalization
for the $w$-mode data listed in Appendix A.
Analytic results for model problems
for the $w$-modes [ks86,nils96],
show that the frequency
of the $w$-mode is inversely proportional to the size of the star.
This is clear from the data in Figure [ref:fww].
Meanwhile, the
damping time is related to
the compactness of the star, i.e. the more relativistic the star is the longer
the $w$-mode oscillation lasts. This is shown in Figure [ref:fwtau].
These properties have already been discussed in some detail by
Andersson, Kojima and Kokkotas \shortcite{akk96}
for uniform density stars.

<!-- TODO: figure content -->

<!-- TODO: figure content -->

For the present numerical data (shown in Figures [ref:fww] and
[ref:fwtau])
we find
the following relations for the frequency and damping of the first $w$-mode:
$$
\omega_w ({\rm kHz}) \approx {1\over {\bar R}}
         \left[ 20.92 - 9.14 \left( {\bar M} \over {\bar R} \right)
\right] \ ,

$$
and
$$
{1 \over \tau_w ({\rm ms})} \approx {1 \over {\bar M}}
\left[ 5.74 + 103 \left( {{\bar M}\over {\bar R}} \right)
- 67.45 \left( {{\bar M}\over {\bar R}} \right)^2 \right] \ .

$$
We see that a typical value for the $w$-mode frequency is 12 kHz, but
since the frequency depends strongly on the radius of the star
it varies greatly for different EOS.
For example, for a very stiff EOS (L) the
$w$-mode frequency is around 6 kHz while for the
softest  EOS in our set (G) the typical frequency is around 14 kHz.
The $w$-mode damping time is comparable to that of an oscillating
black hole with the same mass,  i.e. it is typically less than a tenth
of a millisecond.

## A simple experiment

In principle,
the relations we have deduced between the various pulsation modes and the
stellar parameters can be used to infer $M$ and $R$ (or some
combination thereof) from detected mode-data.
The five relations ([ref:rfw]) --- ([ref:rwtau])
form an over-determined system of five equations
for the two unknown quantities $R$ and $M$. One would expect
this system to
 provide an accurate characterization of the star
in the ideal case when the gravitational-wave
signal carries energy in all modes
($f$, $p$ and $w$).

This idea is promising and simple enough, but we need to examine
how well it can work in
practice. To do this we have constructed a set of independent
polytropic stellar models ($p = K \rho^{1+1/N}$) with varying polytropic
indices ($N=0.8;1;1.2$). We have determined the $f$-mode, the first
$p$-mode and
the slowest
damped  $w$-mode for each of these models.
We let this data
represent “observed”
gravitational-wave signals, and use
 various combinations of the
relations ([ref:rfw]) -- ([ref:rwtau]) to extract the values of the
masses and radii of the stellar models.

In Figure [ref:fig5] we show the
result of a combination of the relations for $f$- and $w-$modes
[([ref:rfw]), ([ref:rftau]),
([ref:rww]) and ([ref:rwtau])] for one of the polytropic models.
In the Figure, a filled
circle represents the true parameters of the star, and it is clear that
estimates based on the above relations can be very accurate.
More detailed results are listed in
Table 2. The typical errors of a parameter estimation based on the
oscillation frequencies of the $f$- and the $w$-mode
[the combination of ([ref:rfw]) and ([ref:rww])]
are  (5% ,2%) where
the first number is the error in the radius and the second
is the error in the mass.
Combining the frequency and damping of the $f$-mode [([ref:rfw]) and
 ([ref:rftau])] we find  (6.5% ,17.6%). The
$f$-mode frequency and the $w$-mode damping rate
[([ref:rfw]) and  ([ref:rwtau])] lead to (5.6% ,1.4%),
while the $w$-mode frequency and the $f$-mode damping
[([ref:rww]) and ([ref:rftau])] yield  (3.2% ,1.9%).
A combination of the $w$-mode frequency and damping rate
[([ref:rww]) - ([ref:rwtau])] leads to  (3.9% ,1.6%).
Finally, by combining the damping rates of the $f$- and the $w$-mode
[([ref:rwtau]) and ([ref:rftau])] we get  (2.1% ,6.3%).
These results are rather impressive. The robustness of our empirical
relations for $f$- and $w$-modes,
and the precision with which they can be used to deduce stellar
masses and radii, is surprising.
The errors are notably larger when we replace either
the $f$- or the $w$-mode with the $p$-mode. For example,
for
the combination ([ref:rfw]) and  ([ref:rpw]), the oscillation frequencies
of the $f$- and $p$-modes, the method estimates the stellar
parameters to within (23%,123%).
That is, from this combination we can at best get upper and lower bounds of
the parameters of the observed object.

<!-- TODO: figure content -->

| lrlcccccc@{}}
  $N$ | $R (km)$ | $M(M_\odot)$ | ([ref:rfw]) and  ([ref:rww]) | ([ref:rfw]) and ([ref:rftau]) | ([ref:rfw]) and ([ref:rwtau]) | ([ref:rww]) and ([ref:rftau]) | ([ref:rww]) and ([ref:rwtau]) | ([ref:rwtau]) and ([ref:rftau]) |
|---|---|---|---|---|---|---|---|---|
| 0.8 | 10.03 | 1.08 | (-3.5 , 0.1) | ( -5.3 ,- 5.4) | (-3.9 ,-1.3) | (-2.6 ,-1.0) | (-2.2 ,-1.3) | ( -2.8 , -1.3) |
| 0.8 | 9.49 | 1.35 | (-0.02,-1.8) | (  4.7 , 11.8) | ( 0.8 , 0.6) | (-4.6 ,-0.6) | (-2.9 ,-1.1) | ( -1.5 , -0.5) |
| 0.8 | 8.99 | 1.50 | ( 1.3 ,-1.1) | (  0.1 , -5.1) | ( 2.1 , 1.1) | (-6.8 ,-1.6) | (-1.6 ,-1.4) | (  0.6 ,  0.02) |
| 1.0 | 9.65 | 1.13 | ( 3.4 ,-4.0) | (  9.7 , 15.2) | ( 4.4 ,-0.6) | (-0.5 ,-1.5) | (-0.3 ,-1.6) | ( -0.5 , -1.7) |
| 1.0 | 8.86 | 1.27 | ( 6.6 ,-3.0) | ( -3.5 ,-40.0) | ( 7.9 , 1.5) | (-1.7 ,-2.5) | ( 1.9 ,-2.0) | ( -0.8 ,-28.9) |
| 1.0 | 7.42 | 1.35 | ( 9.7 , 2.7) | ( 10.9 ,  6.6) | (10.2 , 4.2) | (-4.4 ,-1.4) | ( 7.7 , 1.9) | (  4.3 , -8.7) |
| 1.2 | 12.77 | 1.24 | (-1.6 , 1.6) | (-12.0 ,-32.0) | (-2.5 ,-0.9) | ( 3.0 ,-4.2) | ( 0.7 ,-1.5) | (  4.4 , -1.8) |
| 1.2 | 10.48 | 1.44 | ( 7.4 ,-3.5) | ( -2.2 ,-38.7) | ( 7.8 ,-1.9) | ( 3.8 ,-3.0) | ( 5.6 ,-3.2) | ( -3.2 , -4.7) |
| 1.2 | 8.97 | 1.46 | (11.4 , 0.3) | ( 10.3 , -3.4) | (11.2 ,-0.2) | ( 1.8 ,-1.5) | (12.0 , 0.5) | (  1.3 , -8.9) |

That the $p$-mode relation is less useful for an inversion to yield
the stellar mass and radius is, however, not completely bad news.
Once we  have estimated the mass and the radius we want to
identify which of the proposed EOS that best fits the observed data.
When combined with data deduced from the other modes the $p$-modes
can provide the answer to this question, eg. via the results in
Figure [ref:fptau]. If we  observe a $p$-mode we should at least
be able to exclude the unsuitable EOS.

The most suitable EOS can, of course, also be deduced from the
mass and radius of the star.
As we show in Figure [ref:fig6] the mass-radius relation
is characteristic for each EOS in our sample.
From this data it should not be difficult to infer which EOS
that can lead to a mass and radius obtained via the empirical
relations.
Alternatively, one can use the approach suggested by Lindblom
\shortcite{lind92}. He has shown how one can reconstruct the
 density-pressure relation in the interior of
a neutron star from a sample of observed masses and radii.

<!-- TODO: figure content -->