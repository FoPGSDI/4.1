# 1. Introduction

The spherical cow approximation (SCA) is among the most celebrated simplifications in theoretical physics.
From Newtonian gravity to quantum chromodynamics, the replacement of a complicated object with a sphere of equivalent mass has enabled generations of physicists to obtain closed-form solutions where none would otherwise exist.
The SCA is so deeply embedded in the culture of physics that it has become both a pedagogical device and a punchline --- a shorthand for the productive audacity of theorists who privilege tractability over realism.

Recently, Lehmann [1] undertook a rigorous investigation of the SCA's validity in the context of gravitational wave (GW) emission, computing the full multipole structure of a Holstein dairy cow and evaluating the resulting gravitational radiation.
The central result is devastating: the SCA fails catastrophically for GW physics.
A perfect sphere, possessing only a monopole moment, radiates *zero* gravitational waves.
The cow's gravitational wave luminosity is entirely determined by its higher multipole moments --- the very features that the SCA discards.
This observation elevates the multipole structure of non-spherical bodies from a curiosity to a matter of fundamental importance for gravitational wave astronomy.

In the present work, we extend this program from the bovine to the anatine.
We develop the general theory of oscillation modes, tidal deformability, and gravitational wave signatures for non-spherical self-gravitating bodies, using a rubber duck as the worked example.
The resulting framework --- which we term the theory of **quacky normal modes** (QNMs) --- reveals a rich phenomenology that is entirely absent in the spherical limit.


## 1.1 Why the Duck

The choice of test geometry is not arbitrary.
While the cow is roughly ellipsoidal, retaining an approximate $\mathbb{Z}_2$ symmetry about the sagittal plane and a near-axial symmetry along the spine, the rubber duck breaks *all* continuous symmetries of $SO(3)$.
The duck possesses three geometrically distinct substructures --- bill, neck, and tail --- each of which contributes independently to the spherical harmonic content of the surface deformation.
In the language of the deformation spectrum,

$$R_{\mathrm{duck}}(\theta, \phi) = R_0 \left[ 1 + \sum_{\ell \geq 1} \sum_{m=-\ell}^{\ell} \varepsilon_{\ell m}\, Y_\ell^m(\theta, \phi) \right], \tag{1}$$

the duck exhibits significant power at all multipole orders $\ell \leq 8$, with deformation amplitudes $|\varepsilon_{\ell m}| \sim 0.1$--$0.3$ extending to high $\ell$.
By contrast, the cow's deformation spectrum is dominated by $\ell = 2$ (the body ellipsoid) with rapid falloff at higher $\ell$.

The duck is therefore a **maximally non-spherical** test case: it probes the regime where the SCA fails most dramatically, and where the full machinery of degenerate perturbation theory is required to capture the physics.
As we shall demonstrate, this rich multipole content has observable consequences for every aspect of compact object physics --- from the mode spectrum to the tidal response to the gravitational wave signature.


## 1.2 Three Threads

This paper develops three interconnected threads, each of which reveals new physics arising from the duck's non-spherical geometry.

**Quacky normal modes and gravitational Zeeman splitting.** --- The normal modes of a spherical star are characterized by quantum numbers $(n, \ell, m)$, where each $(n, \ell)$ multiplet is $(2\ell + 1)$-fold degenerate.
The duck shape lifts this degeneracy completely: the deformation $\varepsilon_{\ell' m'}$ couples to the mode spectrum through Gaunt integrals involving Wigner $3j$ symbols, producing a perturbation matrix

$$V_{m_1 m_2}^{(n,\ell)} = -|R_{n\ell}'(R_0)|^2 R_0^3 \sum_{\ell', m'} \varepsilon_{\ell' m'} \int_{S^2} Y_\ell^{m_1 *}\, Y_\ell^{m_2}\, Y_{\ell'}^{m'}\, d\Omega \tag{2}$$

whose eigenvalues give the split frequencies.
We call this **gravitational Zeeman splitting**, in analogy with the lifting of magnetic quantum number degeneracy by an external field --- except here the "field" is the duck itself.
The bill, tail, and body contribute differently to the splitting: the quadrupole deformation $\varepsilon_{2,m}$ drives the primary splitting (analogous to the linear Zeeman effect), while the hexadecapole $\varepsilon_{4,m}$ produces finer structure from the neck and bill tip.
For the $\ell = 2$ $f$-mode, we predict a five-fold splitting with $\Delta\omega / \omega \sim \mathcal{O}(\varepsilon) \sim 0.3$ --- a 30% effect, far larger than the $\sim 10^{-3}$ rotational splitting in millisecond pulsars.

**Tidal Love numbers and the Love tensor.** --- For a spherical body, the tidal response is characterized by a single dimensionless number $k_2$.
For the duck, the response to an external tidal field $\mathcal{E}_{ij}$ is described by a fourth-rank **Love tensor** $\lambda_{ijkl}$:

$$Q_{ij}^{\mathrm{induced}} = -\lambda_{ijkl}\, \mathcal{E}_{kl}, \tag{3}$$

a symmetric $5 \times 5$ matrix (in the space of traceless symmetric tensors) with up to 15 independent components.
The direction-dependent tidal deformability $k_2(\theta, \phi)$ reveals the anisotropic stiffness of the duck: the bill, with its high local compactness, is tidally stiffer than the body.
More dramatically, the celebrated I-Love-Q universality relations of Yagi and Yunes [2, 3] --- which hold to $\mathcal{O}(1\%)$ for all spherical neutron star equations of state --- break down catastrophically for the duck.
The duck-shaped compact object maximally violates I-Love-Q universality, with deviations of $\mathcal{O}(100\%)$ for $\varepsilon \sim 0.3$.
The magnitude of the violation is itself a probe of non-sphericity: a "sphericity test" for compact objects.

**Observational indistinguishability.** --- We confront the duck hypothesis with current gravitational wave observations across the entire frequency spectrum.
In the nanohertz band, a cosmological population of supermassive binary ducks ($M \sim 10^6$--$10^9\, M_\odot$) produces a stochastic gravitational wave background with characteristic strain

$$h_c(f) = A_{\mathrm{duck}} \left( \frac{f}{f_{\mathrm{yr}}} \right)^{-2/3}, \tag{4}$$

where $A_{\mathrm{duck}} = 2.4 \times 10^{-15}$ --- a value indistinguishable from the signal reported by the NANOGrav 15-year dataset [4], the European Pulsar Timing Array [5], the Parkes Pulsar Timing Array [6], and the Chinese Pulsar Timing Array.
The duck correction to the GW power scales as $\delta\dot{E}/\dot{E} \sim \kappa_{\mathrm{duck}} (\pi M_c f / c^3)^{4/3} \sim 10^{-5}$--$10^{-9}$ at PTA frequencies, which is 4--9 orders of magnitude below the $\sim 30\%$ measurement uncertainty on the GWB amplitude.
The spectral index $\gamma = 13/3$ and the Hellings-Downs angular correlation [7] are identical to the SMBHB prediction, as they depend on the inspiral dynamics and the quadrupolar nature of gravitational radiation, not on the shape of the inspiraling bodies.
We prove a **population-level degeneracy theorem**: for any SMBHB population model consistent with current PTA data, there exists a binary duck population with identical masses and merger rates that reproduces the observed signal to arbitrary precision.

In the audio-frequency band, stellar-mass duck mergers ($M \sim 1$--$3\, M_\odot$) produce ringdown waveforms with anomalous QNM splitting.
The gravitational Zeeman splitting of the $\ell = 2$ $f$-mode, with $\Delta\omega / \omega \sim 0.3$, is in principle resolvable by next-generation detectors such as Cosmic Explorer [8] and the Einstein Telescope [9].
The duck's anisotropic tidal deformability $\tilde{\Lambda}_{\mathrm{duck}}$, which differs from the spherical value by $\mathcal{O}(10$--$100\%)$, lies within the broad constraints from GW170817 ($\tilde{\Lambda} = 300^{+420}_{-230}$) [10, 11].
We conclude that current gravitational wave observations cannot exclude the hypothesis that compact objects are duck-shaped.


## 1.3 Paper Organization

The paper is organized in three parts.

**Part I: General Formalism** (Sections 2--5) develops the theoretical framework.
Section 2 characterizes the duck geometry and its spherical harmonic decomposition.
Section 3 derives the quacky normal mode spectrum using degenerate perturbation theory and the Hadamard domain perturbation formula.
Section 4 develops the Love tensor formalism and establishes the breakdown of I-Love-Q universality.
Section 5 constructs the post-Newtonian dynamics of binary duck systems, including the quadrupole-monopole interaction and modified gravitational wave phase evolution.

**Part II: Numerical Simulations** (Sections 6--7) validates and extends the analytic results.
Section 6 presents finite-element solutions of the Helmholtz eigenvalue problem on the duck domain, providing numerical QNM frequencies for cross-validation with the perturbative predictions.
Section 7 describes full numerical relativity simulations using AthenaK [12, 13, 14], evolving a duck-shaped TOV star and extracting gravitational waveforms from the Newman-Penrose scalar $\Psi_4$.
The visual centerpiece of the paper is the evolution of the duck-shaped star as it relaxes to spherical equilibrium: the bill melts first (high-$\ell$ features damp fastest), followed by the tail and neck, until only the fundamental $\ell = 2$ breathing mode remains.

**Part III: Observational Confrontation** (Section 8) places constraints on the duck hypothesis using gravitational wave data.
Section 8.1 demonstrates the indistinguishability of the duck GWB from the NANOGrav 15-year signal.
Section 8.2 assesses the detectability of duck QNM splitting and anomalous tidal deformability with current and future GW detectors.

Section 9 presents our conclusions and discusses prospects for definitively testing the duck hypothesis with future multi-messenger observations.


## 1.4 Notation and Conventions

Throughout this paper we employ geometric units $G = c = 1$ unless otherwise stated; factors of $G$ and $c$ are restored in expressions intended for numerical evaluation.
The metric signature is $(-,+,+,+)$.
Spherical harmonics $Y_\ell^m(\theta, \phi)$ follow the Condon-Shortley phase convention, with the normalization

$$\int_{S^2} Y_\ell^{m*}\, Y_{\ell'}^{m'}\, d\Omega = \delta_{\ell\ell'}\, \delta_{mm'}. \tag{5}$$

Greek indices $\mu, \nu, \ldots$ run over spacetime coordinates $\{0, 1, 2, 3\}$; Latin indices $i, j, \ldots$ run over spatial coordinates $\{1, 2, 3\}$.
The compactness parameter is $C = M/R$ (in geometric units), or equivalently $C = GM/(Rc^2)$.
The dimensionless tidal deformability is $\Lambda = (2/3) k_2 C^{-5}$.
We use $M_\odot$ for the solar mass and define $f_{\mathrm{yr}} \equiv 1/\mathrm{yr} \approx 31.7\,\mathrm{nHz}$ as the reference frequency for PTA analyses.


---

**References cited in this section:**

[1] L. Lehmann, "Higher multipoles of the cow," arXiv:2504.00506 (2025).

[2] K. Yagi and N. Yunes, "I-Love-Q," Science **341**, 365 (2013).

[3] K. Yagi and N. Yunes, "I-Love-Q relations in neutron stars and their applications to astrophysics, gravitational waves, and fundamental physics," Phys. Rev. D **88**, 023009 (2013).

[4] G. Agazie *et al.* (NANOGrav Collaboration), "The NANOGrav 15 yr data set: Evidence for a gravitational-wave background," Astrophys. J. Lett. **951**, L8 (2023).

[5] J. Antoniadis *et al.* (EPTA Collaboration), "The second data release from the European Pulsar Timing Array," Astron. Astrophys. **678**, A50 (2023).

[6] D. J. Reardon *et al.* (PPTA Collaboration), "Search for an isotropic gravitational-wave background with the Parkes Pulsar Timing Array," Astrophys. J. Lett. **951**, L6 (2023).

[7] R. W. Hellings and G. S. Downs, "Upper limits on the isotropic gravitational radiation background from pulsar timing analysis," Astrophys. J. **265**, L39 (1983).

[8] Cosmic Explorer Collaboration, https://cosmicexplorer.org.

[9] Einstein Telescope Collaboration, https://www.et-gw.eu.

[10] B. P. Abbott *et al.* (LIGO Scientific Collaboration and Virgo Collaboration), "GW170817: Observation of gravitational waves from a binary neutron star inspiral," Phys. Rev. Lett. **119**, 161101 (2017).

[11] B. P. Abbott *et al.* (LIGO Scientific Collaboration and Virgo Collaboration), "GW170817: Measurements of neutron star radii and equation of state," Phys. Rev. Lett. **121**, 161101 (2018).

[12] J. M. Stone *et al.*, AthenaK (2024).

[13] H. Zhu *et al.* (2024).

[14] C. E. Fields *et al.* (2024).
