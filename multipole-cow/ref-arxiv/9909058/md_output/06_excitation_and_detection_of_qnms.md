# Excitation and Detection of QNMs

A critical issue related to the discussions of the previous sections
is the excitation of the QNMs. The truth is that, although the QNMs
are predicted by our perturbation equations, it is not always clear
which ones will be excited and under what initial conditions. As we
have already mentioned in the introduction there is an excellent
agreement between results obtained from perturbation theory and full
nonlinear evolutions of Einstein equations for head-on collisions of
two black holes. Still there is a degree of arbitrariness in the
definition of initial data for other types of stellar or black hole
perturbations. This due to the arbitrariness in specifying the
gravitational wave content in the initial data.

The construction of acceptable initial data for the evolution of
perturbation equations is not a trivial task. In order to specify
astrophysically relevant initial data one should first solve the fully
nonlinear 3-dimensional initial value problem for (say) a newly formed
neutron star that settles down after core collapse or two colliding
black holes or neutron stars. Afterwards, starting from the Cauchy
data on the initial hypersurface, one can evolve forward in time with
the linear equations of perturbation theory instead of the full
nonlinear equations. Then most of the long-time evolution problems of
numerical relativity (throat stretching when black holes form,
numerical instabilities or effects due to the approximate outer
boundary conditions) are avoided. Additionally, the interpretations of
the computed fields in terms of radiation is
immediate [AP96]. This scheme has been used by  Price and
Pullin [PP94] and Abrahams and Cook [AC94] with great
success for head-on colliding black holes. The success was based on
the fact that the bulk of the radiation is generated only in the very
strong-field interactions around the time of horizon formation and the
radiation generation in the early dynamics can be practically ignored
(see also the discussion in [AST95]). The extension of this
scheme to other cases, like neutron star collisions or supernovae
collapse, is not trivial. But if one can define even numerical data on
the initial hypersurface then the perturbation method will be probably
enough or at least a very good test for the reliability of fully
numerical evolutions. For a recent attempt towards applying the above
techniques in colliding neutron stars see [ARAKLP, AKLP98].

Before going into details we would like to point out an important
issue, namely the effect of the potential barrier on the QNMs of
black holes. That is, for any set of initial data that one can impose,
the QNMs will critically depend on the shape of the potential barrier,
and this is the reason that the close limit approximation of the two
black-hole collision used by Price and Pullin [PP94] was so
successful, because whatever initial data you provide inside the $r<
3M$ region (the peak of the potential barrier is around $3M$) the
barrier will “filter” them and an outside observer will observe only
the QNM ringing (see for example recent studies  by Allen, Camarda
and Seidel [ACS98]). This point of view is complementary to the
discussion earlier in this section, since roughly speaking even before
the creation of the final black hole the common potential barrier has
been created and anything that was to be radiated had to be
“filtered” by this common barrier.

## Studies of Black Hole QNM Excitation

In the study of QNMs of black holes the attention, in most cases, was
focused on estimating the spectrum and its properties. But
there is limited work in the direction of understanding what details
of the perturbation determine the strength of the QNM ringing.

In 1977 Detweiler [Det77] discussed the resonant oscillations of
a rotating black hole, and after identifying the QNMs as “resonance
peaks” in the emitted spectrum he showed that the modes formally
correspond to poles of a Green function to the inhomogeneous Teukolsky
equation [Teu73]. This idea has been extended in a more
mathematically rigorous way by Leaver [Leaver86]. Leaver extracts
the QNM contribution to the emitted radiation as a sum over
residues. This sum arises when the inversion contour of the Laplace
transform, which was used to separate the dependence on the spatial
variables from the time dependence, is deformed analytically in the
complex frequency plane. In this way the contribution from the QNM can
be accounted for. Sun and Price [SP88,SP90] discussed in detail
the way that QNM are excited by given Cauchy data based to some extent
on numerical results obtained by Leaver [Leaver86]. Lately,
Andersson [Nils95] used the phase-integral method to determine
some characteristics of the QNM excitation.

<!-- TODO: figure content -->

## Studies of Stellar QNM Excitation

The discussion in the previous sections has shown that we have an
acceptable knowledge of stellar pulsations and that we can extract
information from the detection of such oscillations. But as was clear
from the discussion earlier in this section, the energy released as
gravitational radiation during the stellar collapse is basically
unknown. Additionally, it is not known how the energy is distributed
in the various fluid and spacetime modes. Both uncertainties depend
strongly on the details of gravitational collapse, or in general the
mechanism that excites the modes. Unless full 3D general relativistic
codes are generated for the gravitational collapse or the final stages
of binary coalescence we will never be able to give definite answers
to the above questions. In the meantime, perturbation theory is a
reliable way to get some first hints and indications. A survey of the
literature reveals several indications that the pulsation modes are
present in the gravitational wave signals from coalescing
stars. Waveforms obtained by Nakamura and Oohara [nakamura91]
show clear mode presence. Ruffert et al. [ruffert] have also
obtained gravitational wave signals from coalescing stars that show
late-time oscillations. Their waveforms and spectra show oscillations
at frequencies between 1.5 and 2 kHz.

The situation is similar for rotating core collapse. Most available
studies use Newtonian hydrodynamics and account for gravitational wave
emission through the quadrupole formula. The collapse of a
non-rotating star is expected to bounce at nuclear densities, but if
the star is rotating the collapse can also bounce at subnuclear
densities because of the centrifugal force. In each case, the emerging
gravitational waves are dominated by a burst associated with the
bounce. But the waves that follow a centrifugal bounce can also show
large amplitude oscillations that may be associated with pulsations in
the collapsed core. Such results have been obtained by Mönchmeyer et
al. [monch91]. Some of their models show the presence of modes
with different angular dependence superimposed. Typically, these
oscillations have a period of a few ms and damp out in 20 ms. The
calculations also show that the energy in the higher multipoles is
roughly three orders of magnitude smaller than that of the
quadrupole. More recent simulations by Yamada and Sato [yamada]
and Zwerger and Müller [zwerger] also show post-bounce
oscillations. The cited examples are encouraging, and it seems
reasonable that besides the fluid modes the spacetime modes should
also be excited in a generic case. To show that this is the case one
must incorporate general relativity in the simulations of collapse and
coalescence. As yet there have been few attempts to do this, but an
interesting example is provided by the core collapse studies of Seidel
et al. [ed2, ed1, seidel90, seidel91]. They considered axial (odd
parity) or polar (even parity) perturbations of a time-dependent
background (that evolves according to a specified collapse
scenario). The extracted gravitational waves are dominated by a sharp
burst associated with the bounce at nuclear densities. But there are
also features that may be related to the fluid and the
$w$-modes. Especially for the axial case, since there are no axial
fluid modes for a non-rotating star, it is plausible that this mode
corresponds to one of the axial $w$-modes of the core. Furthermore,
the power spectrum for one of the simulations discussed in [ed2]
(cf.\ their Fig.\ 2) shows some enhancement around 7 kHz.

Recently, Andersson and Kokkotas [ak96] have studied the
excitation of axial modes by sending gravitational wave pulses to hit
the star (see figure [ref:nswave1]). The results of this study are
encouraging because they provide the first indication that the
$w$-modes can be excited in a dynamical scenario. Similar results have
recently been obtained by Borelli [BF96] for particles falling
onto a neutron star. The excitation of the axial QNMs  by test
particles scattered by a neutron star have been recently studied in
detail [TSM99, FGB99, AP99] and some interesting conclusions can
be drawn. For example, that the degree of excitation depends on the
details of the particle’s orbit, or that in the cases that we have
strong excitation of the quadrapole oscillations there is a comparable
excitation of higher multipole modes.

These results have prompted the study of an astrophysically more
relevant problem, the excitation of even parity (or polar) stellar
oscillations. In recent work Allen et al. [aaks] have studied the
excitation of the polar modes using two sets of initial data. First,
as previously for the axial case, they excited the modes by an
incoming gravitational wave pulse. In the second case, the initial
data described an initial deformation of both the fluid and the
spacetime. In the first case the picture was similar to that of the
axial modes i.e.\ the star was excited and emitted gravitational waves
in both spacetime and the fluid modes. Nearly all of the energy was
radiated away in the $w$-modes.

<!-- TODO: figure content -->

This should be expected since the incoming pulse does not have enough
time to couple to the fluid (which has a much lower speed of
propagation of information). The infalling gravitational waves are
simply affected by the spacetime curvature associated with the star,
and the outgoing radiation contains an unmistakable $w$-mode
signature. In the second case one has considerable freedom in choosing
the degree of initial excitation of the fluid and the spacetime. In
the study the choice was some “plausible” initial data inspired by
the treatment of the problem in the Newtonian limit. Then the energy
emitted as gravitational waves was more evenly shared between the
fluid and the spacetime modes, and one could see the $f$, $p$, and $w$
modes in the signal. The characteristic signal was (as expected) a
short burst ($w$-modes) followed by a slowly damped sinusoidal wave
(fluid modes).

The previous picture emerged as well in a recent study in the close
limit of head-on collision of two neutron stars [ARAKLP]. There
the excitation of both fluid and spacetime modes is apparent but most
of the energy is radiated via the fluid modes. These new results show
the importance of general relativity for the calculation of the energy
emission in violent processes. Up to now the general way of
calculating the energy emission has been the quadrupole formula. But
this formula only accounts for the fluid deformations and motions, and
as we have seen this accounts only for a part of the total energy.

\vskip 4cm

<!-- TODO: figure content -->

## Detection of the QNM Ringing

It is well known in astrophysics that many stars will end their lives
with a violent supernova explosion. This will leave behind a compact
object which will oscillate violently in the first few seconds. Huge
amounts of gravitational radiation will be emitted and the initial
oscillations will consequently damp out. The gravitational waves will
carry away information about the compact object. If the supernova
remnant is a black hole we will observe a short monochromatic burst
lasting a few tenths of a millisecond (see figures [ref:qnm1]
and [ref:bhwave]). If it is a neutron star we will observe a more
complicated signal which will be the overlapping of several
frequencies (see figure [ref:nswave2]). The stellar signal will
consist of a short burst (similar to that from a black hole) followed
by a long lasting sinusoidal wave. Supernovae are quite frequent, the
expected event rate is 5.8 $\pm$ 2.4 events per century per
galaxy [BT91]. The amount of energy emitted in such an event
depends on the details of the collapse. A spherical collapse will not
produce any gravitational waves at all, while as much as  $10^{-3}
M_\odot c^2$ [bfs97] can be radiated away in a highly
non�spherical one. The collapse releases an enormous amount of energy,
at least equal to the binding energy of a neutron star, about $0.15
M_\odot c^2$. Most of this energy should be carried away by neutrinos,
and this is supported by the neutrino observations at the time of
supernova SN1987A. But even if only $1%$ of the energy released in
neutrinos is radiated in gravitational waves then the above number
makes sense. The present numerical codes used to simulate collapse
predict that the energy emitted as gravitational waves will be of the
order of $10^{-4}-10^{-7} M_\odot c^2$ [BM94]. However, most of
these codes are based on  Newtonian dynamics and the few fully general
relativistic ones are not 3-dimensional. Modern computers are still
not able to perform realistic simulations of gravitational collapse in
3D, including all the important nuclear reactions and neutrino and
photon transport. For example, most of the codes fail to explain the
high average pulsar velocity which is believed to be a result of a
boost that the neutron star gets during the collapse due to anisotropy
in the neutrino distribution [Burrows].

Although collapse may be the most frequent source for excitation of
black hole and stellar oscillations there are other situations in
which significant pulsations take place. For example, after the merger
of two coalescing black holes or neutron stars it is natural to expect
that the final object will oscillate. Thus the well known waveform for
inspiralling binaries [BIWW96] will be followed by a short, but
not yet properly known, period (the merger phase) and will end with
the characteristic quasi-normal signal (ringing) of the newly created
neutron star or black hole. During the inspiralling phase the stellar
oscillations can be excited by the tidal fields of the two
stars [KS95]. A detailed description of the gravitational wave
emission and detection from binary black hole coalescences can be
found in two recent articles by Flanagan and Hughes [FH98a,
  FH98b]. In the same way smaller bodies falling on a neutron star or
black hole will excite oscillations. Stellar or black hole
oscillations can also be excited by a close encounter with another compact
object [TSM99, FGB99, AP99].

Another potential excitation mechanism for stellar pulsation is a
starquake, e.g., associated with a pulsar glitch. The typical energy
released in this process may be of the order of $10^{-10} M_\odot
c^2$. This is an interesting possibility considering the recent
discovery of so-called magnetars: Neutron stars with extreme magnetic
fields [DT92]. These objects are sometimes seen as soft gamma-ray
repeaters, and it has been suggested that the observed gamma rays are
associated with starquakes. If this is the case, a fraction of the
total energy could be released through nonradial oscillations in the
star. As a consequence, a burst from a soft gamma-ray repeater may be
associated with a gravitational wave signal.

Finally, a phase-transition could lead to a sudden contraction during
which a considerable part of the stars gravitational binding energy
would be released, and it seems inevitable that part of this energy
would be channeled into pulsations of the remnant. Transformation of a
neutron star into a strange star is likely to induce pulsations in a
similar fashion.

One way of calibrating the sensitivity of detectors is to calculate
the amplitude of the gravitational wave that would be produced if a
certain fraction of the released energy were converted into
gravitational waves. To obtain rough estimates for the typical
gravitational wave amplitudes from a pulsating star we use the
standard relation for the gravitational wave flux which is valid far
away from the star [bfs96]

$$
  F = {c^3 \over 16\pi G} \vert \dot{h} \vert =
  {1\over 4\pi r^2 } {dE \over dt},
$$

where $h$ is the gravitational wave amplitude and $r$ the distance of
the detector from the source. Combining this with i) ${dE/dt} =
E/2\tau$ where $\tau$ is the damping time of the pulsation and $E$
is the available energy, ii) the assumption that the signal is
monochromatic (with frequency $f$), and iii) the knowledge that the
effective amplitude achievable after matched filtering scales as the
square root of the number of observed cycles, $t_{\rm eff} = h\sqrt{n}
= h\sqrt{f\tau}$, we get the estimates [ak96, AK98]

$$
  h_{\rm eff} \sim 2.2 \times 10^{-21} \left( {E\over 10^{-6} M_\odot c^2}
  \right)^{1/2} \left(  { 2 \mbox { kHz} \over f_{gw} } \right)^{1/2}
  \left( {50 \text{ kpc} \over r} \right)

$$

for the $f$-mode, and

$$
  h_{\rm eff}\sim 9.7\times 10^{-22} \left( { E \over 10^{-6} M_\odot c^2}
  \right)^{1/2} \left( { 10 \text{ kHz} \over f_{gw} } \right)^{1/2}
  \left( { 50 \text{ kpc} \over r } \right)

$$

for the fundamental $w$-mode. Here we have used typical parameters for
the pulsation modes, and the distance scale used is that to
SN1987A. In this volume of space one would not expect to see more than
one event per ten years or so. However, the assumption that the energy
release in gravitational waves in a supernova is of the order of
$10^{-6} M_\odot c^2$ is very conservative [bfs97].

Similar relations can be found for black holes [bfs96]:

$$
  h_{\rm eff} \sim 5\times 10^{-22}
  \left( {E\over 10^{-3} M_\odot c^2} \right)^{1/2}
  \left(  { 1 \text{ kHz} \over f_{gw} } \right)^{1/2}
  \left( {15 \text{ Mpc} \over r} \right),

$$

for stellar black holes, and

$$
  h_{\rm eff}\sim 3\times 10^{-18}
  \left( { E \over 10^3 M_\odot c^2} \right)^{1/2}
  \left( { 1 \text{ mHz} \over f_{gw} } \right)^{1/2}
  \left( { 3 \text{ Gpc} \over r } \right),

$$

for galactic black holes.

An important factor for the detection of gravitational waves are the
pulsation mode frequencies. Existing resonant gravitational wave
detectors, as well as laser interferometric ones which are under
construction, are only sensitive in a certain bandwidth. The spherical
and bar detectors are typically tuned to $0.6-3$ kHz, while the
interferometers are sensitive within $10-2000$ Hz. The initial part of
the QNM waveform, which carries away whatever deformation a collapse
left in the spacetime, is expected to be for a neutron star in the
frequency range of $5-12$ kHz ($w$-mode). The subsequent part of the
waveform is constructed from combination of the $f$- and
$p$-modes. Still the present gravitational wave detectors are
sensitive only in the frequencies of the $f$-mode. For a black hole
the frequency will depend on the mass and rotation rate[^1], thus for a 10 solar mass black hole the frequency of
the signal will be around 1 kHz, around 100 Hz for a 100 $M_\odot$
black hole and around 1 mHz for galactic black holes.

[^1]: For neutron stars the frequencies depend not only on the mass and rotation rate, but also on the radius and the equation of state.

## Parameter Estimation

For astronomy it is important not only to observe various astronomical
phenomena but also to try to mine information from these
observations. From the observations of solar and stellar oscillations
(of normal stars) astronomers have managed to get details of the
internal structure of stars. In our days the GONG program [gong]
for detailed observation of the solar seismology is well
underway. This has suggested that, in a similar way, information
about neutron star parameters (mass, radius), and internal structure or
the mass and the rotation rate of black holes can be found, using the theory
of QNMs. It will be instructive to briefly examine the
case of oscillating black holes since they are much “cleaner”
objects than stars. From the normal mode analysis of black hole
oscillations we can get a spectrum which is related to the parameters
of the black hole (mass $M$ and angular momentum $a$). In particular,
for the frequency of the first quasi-normal mode (which as we have
stated previously is the most important one for the gravitational wave
detection) the following approximate relations have been
suggested [Ech89, Finn92]:

$$
  M \omega \approx \left[ 1 - {{63} \over {100}}(1-a)^{3/10} \right]
  \approx (0.37 +0.19 a),

$$

$$
  \tau \approx {{4 M}\over { (1-a)^{9/10}}} \left[ 1 - {{63} \over
  {100}}(1-a)^{3/10} \right]^{-1} \approx M(1.48+2.09 a).

$$

These two relations can be inverted and thus from the “observed”
frequency and the damping time we can derive the parameters of the
oscillating black hole. In practice, the noise of the detector will
contaminate the signal but still (depending on the signal to noise
ratio) we will get a very accurate estimate of the black hole
parameters. A similar set of empirical relations cannot be derived in
the case of neutron star oscillations since the stars are not as
“clean” as black holes, since more than one frequency
contributes. Although we expect that most of the dynamical energy
stored in the fluid oscillations will be radiated away in the
$f$-mode, some of the $p$-modes may be excited as well and a
significant amount of energy could be radiated away through these
modes [aaks]. As far as the spacetime modes are concerned we
expect that only the curvature modes (the standard $w$-modes) will be
excited, but it is possible that the radiated energy can be shared
between the first two $w$-modes [aaks]. Nevertheless, Andersson
and Kokkotas [AK98], using the properties of the various families
of modes ($f$, $p$, and $w$), managed to create a series of empirical
relations which can provide quite accurate estimates of the mass,
radius and equation of state of the oscillating star, if the $f$ and
the first $w$-mode can be observed. In figure [ref:ffreq] one can see
an example of the relation between the stellar parameters and the
frequencies of the $f$ and the first $w$-mode for various equations of
state and various stellar models. There it is apparent that the
relation between the $f$-mode frequencies and the mean density is
almost linear, and a linear fitting leads to the following simple
relation:

$$
  \omega_f (\text{kHz}) \approx 0.78 + 1.635 \left[ \left( {M
  \over{1.4M_\odot}} \right) \left( {{10 \text{ km}}\over R}
  \right)^3 \right]^{1/2}.

$$

We can also find the following relation for the frequency of the first
$w$-mode:

$$
  \omega_w (\text{kHz}) \approx \left({{10 \text{ km}}\over
  {R}}\right) \left[ 20.92 - 9.14 \left( {M \over {1.4
  M_\odot}}\right) \left({{10 \text{ km}} \over R} \right)
  \right].

$$

<!-- TODO: figure content -->

From tests performed using polytropic stars to provide data for the
above relations it was seen that these equations  predict the masses
and the radii of the polytropes usually with an error less than
10%. There is work underway towards extracting the parameters of the
star from a noisy signal [KAA98].