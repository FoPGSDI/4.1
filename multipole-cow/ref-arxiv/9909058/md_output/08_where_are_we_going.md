# Where Are We Going?

From the previous discussions we can draw some general conclusions and
suggestions for future work. We believe that independent of the
advancement of numerical relativity and computer power there remains
much future work in perturbation theory, in particular  parallel to
the numerical relativity efforts.

## Synergism Between Perturbation Theory and Numerical Relativity

Fifteen years ago with the advancements in computer power, we could
not think that perturbation theory would continue to play such an
important role in many problems. Today it is widely accepted that
there is a need for the development of new approximation schemes to
accompany large scale simulations. Fully relativistic computer
simulations  give only numerical answers to problems; often these
answers do not provide physical understanding of which principles are
important, or even what principles govern a given process. Even more, in
some cases, simulation results can be simply incorrect or
misleading. By closely coupling various perturbation schemes it is
possible to interpret and confirm simulation results.

For example, during the late stages of black hole or stellar
coalescence or supernovae collapse, the system settles down to a
slightly perturbed black hole or neutron star. Numerical codes,
evolving the full nonlinear Einstein equations, should be able to
accurately compute the waveforms  required for gravitational wave
detection. Parallel to this, it should be possible to evolve the
perturbations of both black holes and stars (governed by their own
linear evolution equations) for the same set of initial data. This is
an important check of the fully nonlinear codes. An excellent example
is the head-on collision of two black holes. This sounds like an
impossible task for perturbation theory but it can be achieved if the
two black holes are close together and they can be considered as
having already merged into a single perturbed black hole (close
limit); see more details in a recent review by
Pullin [Pullin98]. Much work has already been performed for
head-on collisions of two non-rotating black holes but the more
realistic problem for the inspiral collision of rotating black holes
is still open.

Following in the same spirit are more recent
calculations [ARAKLP] of the close limit for two identical
neutron stars. The results show the excitation of stellar QNMs and it
remains for numerical relativity to verify the results.

Finally, perturbation theory can be also used as a tool to construct a
gauge invariant measure of the gravitational radiation in a
numerically generated perturbed black hole spacetime [ACS98,
  Seidel98].

## Second Order Perturbations

A basic feature of linearized perturbation theory is that there is no
“built-in” indication of how good the approximation is. But if one
has results for a physical quantity to second order in a perturbation
parameter, then the difference between the results of the first and
second order theory is a quantitative indication of the error in the
perturbation calculation. From this point of view second order
perturbation theory is a practical tool for
calculations. Nevertheless, there remain many technical problems in
this approach, for example the gauge that should be
used [GNPP98c, BMMS97] and the  amount of calculations. Thus
second order perturbation theory has turned out to be much more
difficult than linearized theory, but if one overcomes these
difficulties second order calculations will be a great deal easier
than numerical relativity. From this point of view we encourage work
in this direction. In particular, the perturbations of stellar
oscillations should be extended to second order, and the study of the
second order perturbations of black holes extended to the Kerr case.

## Mode Calculations

Although the modes have been well studied there remain a few open
issues. In the black hole case, the Kerr-Newman spectrum is known only
in a restricted case [KDK93], so the general case should be
studied, probably by evolving the perturbation equations. For
non-rotating stars, there remain certain technical questions, i.e.\ to
find if there exist a class of $w$-modes with even larger imaginary
part or the existence of extra interface modes. But for rotating
stars, even slowly rotating ones, there remains much
work [Nick98]. There should be estimates of the fluid and
spacetime modes for slowly rotating stars for various realistic
equations of state, and the results should be incorporated into the
method suggested in [AK98] for the estimates of the stellar
parameters. There should be work towards understanding the possible
interaction of $r$-modes with $g$-modes [schutz, friedman]. And
finally the time dependent perturbation equations for
rotating stars should be evolved, because in this way we expect to see most of the new features that
rotation induces in the spectra (splitting, instabilities etc).

## The Detectors

As we have mentioned earlier in section [ref:section_5], there are
techniques available for the extraction of the QNM signal from the
noise of the detectors [Ech89, FH98a, FH98b]. But there are still
issues related to the sensitivity of the planned detectors in parts of
the spectrum. For example, the QNM frequencies of stellar black holes
($10 - 100 M_\odot$) will be of around $100-1000$ Hz, i.e.\ in the
frequencies where the laser interferometers are sensitive. The QNM
frequencies of galactic size black holes will be detectable only from
space (LISA) since their frequencies will be in the mHz regime. For
the QNM frequencies of stars, there is a lack of appropriate
detectors. Laser interferometers will be sensitive enough in the
frequency regime of the $f$-mode but it will be very hard to detect
signals in the frequencies of the $p$- and $w$-modes. Nevertheless,
from the discussion in section [ref:section_5] it is apparent that
there is a wealth of information in the signal of oscillating neutron
stars, and in order to extract this information we need the $p$-
or/and $w$-modes. This suggests that the ideas considering detectors,
or arrays of detectors, in this high frequency regime [papa]
should be considered more seriously.