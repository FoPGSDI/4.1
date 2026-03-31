# 9. Discussion and Conclusions

We have developed a comprehensive framework for the oscillation modes, tidal response, and gravitational wave signatures of non-spherical compact objects, using the rubber duck as the maximally asymmetric worked example.  The results organize into three main findings, which we summarize in turn before discussing future directions.


## 9.1 Summary of Results

**Gravitational Zeeman splitting.** --- The duck's complete breaking of spherical symmetry lifts the $(2\ell+1)$-fold degeneracy of every QNM multiplet.  The mechanism is entirely analogous to the Zeeman effect in atomic physics: the deformation spectrum $\varepsilon_{\ell' m'}$ plays the role of the external field, coupling to the mode structure through Gaunt integrals and Wigner $3j$ symbols.  The resulting splitting is first-order in $\varepsilon$ and follows strict selection rules from angular momentum coupling --- only even-$\ell'$ deformations contribute, with the quadrupole ($\ell' = 2$) driving the primary splitting and the hexadecapole ($\ell' = 4$) producing finer structure.

Physically, the bill and body contribute to different channels of the deformation spectrum: the axial elongation of the body dominates $\varepsilon_{2,0}$, splitting $m = 0$ from $m = \pm 1, \pm 2$, while the bill-to-tail asymmetry and equatorial ellipticity populate $\varepsilon_{2,\pm 1}$ and $\varepsilon_{2,\pm 2}$, respectively.  The neck and bill tip contribute primarily at $\ell' = 4$.  For the $\ell = 2$ $f$-mode, the predicted five-fold splitting has magnitude $\Delta\omega / \omega \sim \mathcal{O}(\varepsilon) \sim 0.3$ --- a 30% effect that dwarfs the $\sim 10^{-3}$ rotational splitting in the fastest-spinning pulsars.  The perturbative predictions were validated against FEM eigenvalue computations (Section 6), with agreement at the expected $\mathcal{O}(\varepsilon^2)$ level.

**Tidal Love tensor and I-Love-Q breakdown.** --- For a spherical body, the tidal response is fully characterized by the scalar Love number $k_2$.  The duck demands more: its response to an external tidal field $\mathcal{E}_{ij}$ is governed by a fourth-rank Love tensor $\lambda_{ijkl}$ with up to 15 independent components, reflecting the 15-dimensional space of traceless symmetric $3 \times 3$ matrices.  The direction-dependent Love number $k_2(\theta, \phi)$ varies by a factor of $\sim 2$ across the duck surface --- the bill region, with its high local compactness, is tidally stiff, while the rounded body is comparatively soft.

The I-Love-Q universality relations, which hold to $\mathcal{O}(1\%)$ across all spherical equations of state, break down at $\mathcal{O}(100\%)$ for the duck.  The duck maximally violates universality because its moment of inertia $\bar{I}$ and tidal deformability $\bar{\Lambda}$ receive shape corrections that scale differently with the deformation amplitudes --- $\bar{I}$ is corrected at $\mathcal{O}(\varepsilon^2)$ via the principal moments, while $\bar{\Lambda}$ receives a geometric $\mathcal{O}(\varepsilon^2)$ enhancement from the $R^{-5}$ scaling.  The magnitude of the discrepancy constitutes a "sphericity test" for compact objects: any observed violation of I-Love-Q at the $\gtrsim 10\%$ level would be evidence for significant non-spherical deformation.

**Observational indistinguishability.** --- Despite the rich phenomenology predicted by the duck framework, current gravitational wave observations cannot distinguish duck-shaped compact objects from spherical ones.  In the PTA band, the duck quadrupole correction to the stochastic GW background is $\mathcal{O}(10^{-5})$ or smaller, some 4--9 orders of magnitude below the $\sim 30\%$ measurement uncertainty on $A_{\mathrm{GWB}}$.  The population-level degeneracy theorem (Section 8.1.4) formalizes this: any SMBHB model consistent with NANOGrav 15-year data has a binary duck counterpart that reproduces the observed spectrum to better than $10^{-5}$.  The Hellings-Downs angular correlation is identical for ducks and spheres, as it depends only on the spin-2 tensor structure of the radiation.

In the LIGO band, the duck tidal deformability correction is $\mathcal{O}(\varepsilon^2) \sim 10\%$ --- comparable to the equation-of-state uncertainties that already limit the interpretation of GW170817.  The GW170817 constraint $\tilde{\Lambda} = 300^{+420}_{-230}$ accommodates the full range of duck orientations and nuclear equations of state.  QNM spectroscopy offers the most promising avenue, as the $\Delta\omega/\omega \sim 0.3$ splitting exceeds the spectral resolution $1/Q \sim 0.05$--$0.2$, but current LIGO ringdown SNR is marginal ($d_{\max} \sim 15\,\mathrm{Mpc}$ at O4 sensitivity).


## 9.2 Future Directions

Several avenues for extending the present work suggest themselves, spanning a range of ambition and ornithological scope.

**Duck spectroscopy with next-generation detectors.** --- Cosmic Explorer and the Einstein Telescope will achieve order-of-magnitude improvements in sensitivity at $f \sim 2\,\mathrm{kHz}$, pushing the horizon for ringdown detection to $d \sim 200$--$300\,\mathrm{Mpc}$.  At these distances, $\sim 5$--$10$ events per year would yield sufficient ringdown SNR to resolve the $m$-splitting of the $\ell = 2$ $f$-mode.  A systematic Bayesian analysis comparing the Kerr ringdown template against the "Zeeman duck" template would constitute a definitive test of compact object sphericity.  We defer the construction of such a template bank to future work.

**Multi-messenger constraints on duck morphology.** --- While neither gravitational waves nor electromagnetic imaging can independently constrain the duck shape at current sensitivity, a combined multi-messenger analysis could in principle break the degeneracy.  For instance, a simultaneous measurement of the tidal deformability (from the inspiral) and the $m$-splitting pattern (from the ringdown) of a single event would overdetermine the deformation spectrum $\varepsilon_{\ell m}$, enabling a consistency test.  Whether nature provides duck-shaped compact objects at distances where both channels are accessible remains an open question.

**Anatidae survey.** --- The present work has focused exclusively on the DASSL rubber duck.  A comprehensive survey of QNMs across the family Anatidae --- including the mallard (*Anas platyrhynchos*), the common teal (*Anas crecca*), and the domestic goose (*Anser anser domesticus*) --- would map out the landscape of avian deformation spectra and their gravitational wave consequences.  We note that the goose, with its elongated neck, is expected to exhibit anomalously large $\varepsilon_{3,m}$ coefficients, potentially enabling odd-parity QNM splitting that is suppressed in the more compact duck.  This survey is left for future work.

**The duck matter problem.** --- Throughout this paper, we have assumed that the duck is composed of nuclear matter at extreme density ($\rho \sim 10^{14}$--$10^{15}\,\mathrm{g\,cm^{-3}}$), supported against gravitational collapse by a standard equation of state.  However, the duck shape requires anisotropic stresses $\sigma \sim G \rho^2 R^2 \varepsilon \sim 10^{33}\,\mathrm{dyn\,cm^{-2}}$ to maintain static equilibrium --- roughly $10^5$ times the isotropic nuclear pressure.  The microphysical origin of such stresses is unclear.  Exotic candidates include a crystalline QCD phase with preferential alignment along the bill axis, magnetically confined pasta phases with toroidal topology, or a hitherto undiscovered "anatine" phase of dense matter.  The equation of state of duck matter remains an open problem in nuclear astrophysics.


## 9.3 Closing Remarks

We have shown that the duck hypothesis --- the proposition that compact astrophysical objects may possess the geometry of a rubber duck --- is internally consistent, theoretically rich, and observationally unconstrained.  The formalism developed here generalizes straightforwardly to any non-spherical geometry and provides a systematic framework for quantifying the observational consequences of compact object morphology.

Until next-generation gravitational wave detectors achieve the sensitivity required for duck spectroscopy, the duck hypothesis remains consistent with all current gravitational wave data.  We encourage the community to keep an open mind --- and an open bill.


---

**References cited in this section:**

References are shared with the preceding sections; see the bibliography for the full list.
