# 8. Confrontation with Gravitational Wave Data

Having developed the theoretical framework (Part I) and validated it numerically (Part II), we now confront the duck hypothesis with observational data across the gravitational wave spectrum.
The central question is: **can current or near-future gravitational wave observations distinguish duck-shaped compact objects from spherical ones?**

The answer, as we shall demonstrate, is no.


## 8.1 PTA Band --- Supermassive Binary Ducks (nHz)

### 8.1.1 Population of Supermassive Binary Ducks

Consider a cosmological population of supermassive binary ducks with component masses $M \sim 10^6$--$10^9\, M_\odot$, evolving via GW-driven inspiral in the nanohertz band ($f \sim 1$--$100\,\mathrm{nHz}$).
The orbital velocities at these frequencies are highly sub-relativistic: at $f = 10\,\mathrm{nHz}$ for a chirp mass $\mathcal{M}_c = 10^9\,M_\odot$,

$$\frac{v}{c} = \left( \frac{\pi G \mathcal{M}_c f}{c^3} \right)^{1/3} \approx 4.7 \times 10^{-2}, \tag{51}$$

implying that the binary dynamics are deep in the Newtonian regime.  The finite size of the compact objects --- whether spherical or duck-shaped --- enters the dynamics only through tidal and spin-multipole couplings suppressed by powers of $(R/r)^2 \sim (v/c)^4$.

### 8.1.2 The Phinney Formula and Duck Correction

The characteristic strain spectrum from an astrophysical population of circular GW-driven binaries is given by the Phinney formula [15]:

$$h_c^2(f) = \frac{4G}{3\pi c^2} \frac{1}{f} \int dz\, d\mathcal{M}_c\, \frac{d^2 n}{dz\, d\mathcal{M}_c} \frac{(\mathcal{M}_c^{(z)})^{5/3}}{(1+z)^{1/3}} \frac{(\pi f)^{2/3}}{d_L^2(z)} \frac{dE/df_r}{dE_0/df_r}, \tag{52}$$

where the last factor is the ratio of the actual GW energy spectrum to the point-particle prediction.  For circular, GW-driven inspirals of spherical objects, $dE/df = dE_0/df$ and the integral yields the familiar power-law spectrum:

$$h_c(f) = A_{\mathrm{GWB}} \left( \frac{f}{f_{\mathrm{yr}}} \right)^{-2/3}, \tag{53}$$

with spectral index $\gamma = 13/3$.

For duck-shaped binaries, the quadrupole-monopole interaction (Section 5.1) modifies the orbital energy and GW luminosity.  The fractional correction to the energy spectrum is:

$$\frac{\delta(dE/df)}{dE/df} \sim \kappa_{\mathrm{duck}} \left( \frac{\pi \mathcal{M}_c f}{c^3} \right)^{4/3}, \tag{54}$$

where $\kappa_{\mathrm{duck}} = Q^C_{\mathrm{eigenvalue}} / (M R_{\mathrm{eq}}^2) \sim 0.3$ is the shape-induced quadrupole parameter defined in Section 2.  This correction enters at 2PN order and is proportional to the fourth power of the orbital velocity.

### 8.1.3 Numerical Evaluation

We evaluate the duck correction at representative PTA frequencies for a fiducial binary with $\mathcal{M}_c = 10^9\,M_\odot$.  Table 11 presents the orbital velocity and duck correction at three frequencies spanning the PTA band.

**Table 11.** Duck quadrupole correction to the GW energy spectrum for $\mathcal{M}_c = 10^9\,M_\odot$, $\kappa_{\mathrm{duck}} = 0.3$.

| $f$ (nHz) | $v/c$ | $\delta_{\mathrm{duck}}$ | NANOGrav uncertainty |
|:---:|:---:|:---:|:---:|
| 1   | $0.025$ | $1.2 \times 10^{-7}$ | $\sim 30\%$ |
| 10  | $0.054$ | $2.5 \times 10^{-6}$ | $\sim 30\%$ |
| 100 | $0.12$  | $5.4 \times 10^{-5}$ | $\sim 30\%$ |

The duck correction is 4--7 orders of magnitude below the current measurement uncertainty.  Even at the most optimistic upper end of the PTA band ($f \sim 100\,\mathrm{nHz}$), the correction is $\delta_{\mathrm{duck}} \sim 10^{-5}$, which would require a fractional strain measurement precision of $\sim 10^{-7}$ --- far beyond any foreseeable PTA capability.

### 8.1.4 Population-Level Degeneracy Theorem

We now state the central result of this subsection.

> **Theorem (Population-level degeneracy).** *For any SMBHB population model $\{M_i, q_i, z_i\}$ consistent with the NANOGrav 15-year dataset, there exists a binary duck population with identical masses, mass ratios, and redshift distribution that reproduces the observed $h_c(f)$ to better than $10^{-5}$ at all PTA frequencies.*

*Proof.* --- The duck correction (Eq. 54) modifies $h_c(f)$ by a factor $[1 + \mathcal{O}(\delta_{\mathrm{duck}})]^{1/2}$.  Since $\delta_{\mathrm{duck}} < 10^{-6}$ for all binaries with $f < 300\,\mathrm{nHz}$ and $\mathcal{M}_c < 10^{10}\,M_\odot$, the integrated strain spectrum from any population of binary ducks agrees with the corresponding spherical binary population to better than $10^{-5}$ in $h_c$.  The spectral index $\gamma = 13/3$ is unchanged, as it is determined by the GW-driven frequency evolution $\dot{f} \propto f^{11/3}$, which receives the same negligible correction.  $\square$

### 8.1.5 NANOGrav 15-Year Comparison

The NANOGrav 15-year dataset [4] reports a gravitational wave background with amplitude and spectral index:

$$A_{\mathrm{GWB}} = (2.4 \pm 0.7) \times 10^{-15}, \qquad \gamma = 13/3, \tag{55}$$

at the reference frequency $f_{\mathrm{yr}} = 1/\mathrm{yr} \approx 31.7\,\mathrm{nHz}$.  The signal is consistent with the Hellings-Downs angular correlation at $\sim 3$--$4\sigma$ significance.

The duck prediction is:

$$A_{\mathrm{duck}} = A_{\mathrm{SMBHB}} \left[1 + \mathcal{O}(10^{-7})\right] = 2.4 \times 10^{-15}, \tag{56}$$

which is, by the degeneracy theorem, indistinguishable from the SMBHB prediction.

### 8.1.6 Hellings-Downs Correlation

The Hellings-Downs curve [7] describes the expected angular correlation of timing residuals between pulsar pairs as a function of their angular separation $\theta$:

$$\Gamma(\mu) = \frac{1}{3} \left\{ 1 + \frac{3}{2}(1 - \mu) \left[ \ln\!\left(\frac{1-\mu}{2}\right) - \frac{1}{6} \right] \right\}, \tag{57}$$

where $\mu = \cos\theta$.  This correlation arises from the quadrupolar nature of gravitational radiation and the geometry of the Earth-term response [16].  Crucially, $\Gamma(\mu)$ depends only on the spin-2 tensor structure of the metric perturbation --- it is independent of the internal structure or shape of the GW sources.  A population of inspiraling ducks produces exactly the same Hellings-Downs correlation as a population of inspiraling black holes or neutron stars, since all produce quadrupolar ($\ell = 2$) gravitational radiation in the transverse-traceless gauge.

As demonstrated in Ref. [16], the Hellings-Downs curve is equivalent to summing the total power in a PTA sky map --- the "optimal statistic" for an isotropic stochastic background.  The sky-map formalism encodes additional information about anisotropy and polarization in the higher multipoles of the GW power distribution, but neither the isotropic monopole nor the angular correlation function can distinguish ducks from black holes.

### 8.1.7 Electromagnetic Indistinguishability

One might hope to resolve the duck shape directly through electromagnetic imaging.  At mass $M = 10^9\,M_\odot$, the Schwarzschild radius is $R_s \approx 3 \times 10^{12}\,\mathrm{m}$.  At redshift $z = 0.1$ (luminosity distance $d_L \approx 470\,\mathrm{Mpc}$), the angular size of the duck is:

$$\theta_{\mathrm{duck}} \approx \frac{R_s}{d_L} \approx \frac{3 \times 10^{12}\,\mathrm{m}}{1.4 \times 10^{25}\,\mathrm{m}} \approx 2 \times 10^{-13}\,\mathrm{rad} \approx 4 \times 10^{-8}\,\mathrm{arcsec}. \tag{58}$$

For comparison, the Event Horizon Telescope achieves an angular resolution of $\sim 20\,\mu\mathrm{as} = 2 \times 10^{-5}\,\mathrm{arcsec}$ [17].  The duck is $\sim 5$ orders of magnitude too small to resolve.  Even a hypothetical space-based VLBI at millimeter wavelengths with Earth-Moon baseline ($\sim 4 \times 10^{8}\,\mathrm{m}$) would achieve only $\sim 0.1\,\mu\mathrm{as}$, still $\sim 3$ orders of magnitude too coarse.

The duck is electromagnetically indistinguishable from a sphere.

### 8.1.8 Mock PTA Skymap

To visualize the observational landscape, we construct a mock PTA skymap using the methodology of Ref. [16].  We generate a Mollweide projection of the nanohertz GW sky containing:

1. An isotropic stochastic background with $A_{\mathrm{GWB}} = 2.4 \times 10^{-15}$, represented as a smooth monopole contribution.
2. Five "loud" individual supermassive binary duck sources at random sky positions, with strains $h \sim (1$--$5) \times 10^{-15}$ at $f = 10\,\mathrm{nHz}$.

The resulting skymap (Fig. 16) is, by construction, identical to the skymap that would be produced by five supermassive binary black holes with the same masses and sky positions.  The PTA angular resolution, limited to $N_{\mathrm{res}} \sim 2N_{\mathrm{pulsar}}$ independent resolution elements [16], cannot resolve the shape of the sources.

### 8.1.9 The Punchline

We arrive at the central observational result of this paper:

> **The nanohertz gravitational wave background observed by NANOGrav, EPTA, PPTA, and CPTA is equally consistent with a cosmological population of inspiraling ducks.**

This is not a failure of the observations --- it is a fundamental consequence of the separation of scales between the body size ($R \sim R_s \sim GM/c^2$) and the orbital separation ($r \sim c/\pi f \gg R$) in the PTA band.  The gravitational wave emission depends on the *orbital* quadrupole moment, not the *body* quadrupole moment, and the two are related by factors of $(R/r)^2 \sim (v/c)^4 \ll 1$.


## 8.2 LIGO Band --- Stellar-Mass Ducks (10--1000 Hz)

In the audio-frequency band accessible to ground-based interferometers, stellar-mass ducks ($M \sim 1$--$3\,M_\odot$) present two observational channels: tidal effects during inspiral and QNM spectroscopy during ringdown.

### 8.2.1 Inspiral --- Tidal Deformability

**(a) The tidal phase correction.** --- Tidal effects enter the gravitational wave phasing at 5PN order (relative to the leading quadrupole formula) [18, 19]:

$$\delta\Psi_{\mathrm{tidal}}(f) = -\frac{39}{2}\,\tilde{\Lambda}\,\left(\frac{\pi \mathcal{M}_c f}{c^3}\right)^{5/3} + \cdots, \tag{59}$$

where $\tilde{\Lambda}$ is the mass-weighted tidal deformability of the binary:

$$\tilde{\Lambda} = \frac{16}{13} \frac{(M_1 + 12 M_2) M_1^4\,\Lambda_1 + (M_2 + 12 M_1) M_2^4\,\Lambda_2}{(M_1 + M_2)^5}. \tag{60}$$

**(b) Duck tidal deformability.** --- For the duck, the tidal deformability $\Lambda_{\mathrm{duck}}$ differs from the spherical value $\Lambda_{\mathrm{sphere}}$ by the perturbative correction derived in Section 4.4:

$$\Lambda_{\mathrm{duck}} = \Lambda_{\mathrm{sphere}} \left[ 1 + \sum_{\ell \geq 1} \alpha_\ell \sum_m |\varepsilon_{\ell m}|^2 + \mathcal{O}(\varepsilon^3) \right]. \tag{61}$$

The geometric contribution alone gives $\delta\Lambda/\Lambda \approx 15 \sum_\ell \sum_m |\varepsilon_{\ell m}|^2 \sim 150\%$ for the duck (Section 4.4), while the physical correction from $\delta k_2/k_2$ adds a comparable contribution.  In total, $\Lambda_{\mathrm{duck}}$ can differ from $\Lambda_{\mathrm{sphere}}$ by $\mathcal{O}(10$--$100\%)$, depending on the equation of state and the duck orientation.

**(c) Comparison with GW170817.** --- The binary neutron star merger GW170817 [10] yielded the constraint [11]:

$$\tilde{\Lambda} = 300^{+420}_{-230} \qquad (90\%\text{ credible interval}). \tag{62}$$

A nuclear duck with $\Gamma = 2$, $C \approx 0.15$ has $\Lambda_{\mathrm{sphere}} \sim 290$--$880$ (Section 4.1), and the duck correction shifts this by $\mathcal{O}(100\%)$.  Given the broad uncertainty on $\tilde{\Lambda}$, the duck-shaped compact object is consistent with GW170817 for a wide range of equations of state and orientations.

**(d) Anisotropic Love tensor.** --- The duck's Love tensor $\lambda_{ijkl}$ (Section 4.2) introduces a qualitatively new effect: the effective tidal deformability $\tilde{\Lambda}$ depends on the orientation of the duck relative to the orbital plane.  For a duck with its bill aligned along the orbital angular momentum, $\tilde{\Lambda}$ is maximized (the body is softest along the tidal axis); for a bill-in-plane orientation, $\tilde{\Lambda}$ is minimized.  This orientation dependence is absent for spherical neutron stars and constitutes a genuine observable signature of non-sphericity.  However, the GW170817 uncertainty band is broad enough to accommodate the full range of duck orientations.

### 8.2.2 Ringdown --- QNM Spectroscopy

**(a) Duck QNM spectrum in the ringdown.** --- Following a duck-duck merger, the post-merger remnant rings down through its quasinormal mode spectrum (Section 3).  For a spherical remnant, the $\ell = 2$ fundamental mode has a single complex frequency $\omega_0 + i/\tau_0$.  For the duck, the five-fold degeneracy is lifted by gravitational Zeeman splitting, producing five distinct frequencies with splittings (Section 3.5):

$$\frac{\Delta\omega}{\omega_0} \sim \mathcal{O}(\varepsilon) \sim 0.3. \tag{63}$$

The ringdown waveform of the duck is:

$$h_+(t) + i h_\times(t) = \frac{1}{D_L} \sum_{n,\ell,m} A_{n\ell m}\, e^{-t/\tau_{n\ell m}}\, e^{i\omega_{n\ell m} t}\, {}_{-2}Y_{\ell m}(\iota, \varphi), \tag{64}$$

where each $(n,\ell,m)$ substate has a distinct frequency $\omega_{n\ell m}$ and damping time $\tau_{n\ell m}$ determined by the duck's deformation spectrum.

**(b) Resolvability criterion.** --- Two QNM frequencies can be resolved in the ringdown signal if the splitting exceeds the spectral width of each mode:

$$\frac{\Delta\omega}{\omega_0} > \frac{1}{Q_{\mathrm{mode}}}, \tag{65}$$

where $Q_{\mathrm{mode}} = \omega_R / (2\omega_I)$ is the quality factor.  For the $\ell = 2$ $f$-mode of a typical neutron star, $Q \sim 5$--$20$ [20, 21], giving $1/Q \sim 0.05$--$0.2$.  Since $\Delta\omega/\omega \sim 0.3$ for the duck, the splitting satisfies $\Delta\omega/\omega > 1/Q$ --- the duck QNMs are in principle spectrally resolvable.

**(c) SNR requirements.** --- Resolving $m$-split QNMs requires sufficient signal-to-noise ratio in the ringdown portion of the waveform.  The ringdown SNR for a single damped sinusoid scales as [22]:

$$\rho_{\mathrm{ring}} \sim \frac{h_0}{S_n(f_{\mathrm{QNM}})^{1/2}} \sqrt{\frac{\tau}{2}}, \tag{66}$$

where $h_0$ is the initial ringdown amplitude, $S_n$ is the detector noise spectral density, and $\tau$ is the damping time.  To resolve the splitting, we require $\rho_{\mathrm{ring}} \gtrsim 8$ per resolved mode [23].  For the $\ell = 2$ $f$-mode at $f \sim 2\,\mathrm{kHz}$:

**Table 12.** Detectability of duck QNM $m$-splitting at various distances for a $1.4\,M_\odot$ duck remnant with $\Delta\omega/\omega = 0.3$.

| Detector | $d_{\max}$ (single mode, $\rho > 8$) | $m$-splitting resolvable? |
|:---:|:---:|:---:|
| LIGO O4 (current) | $\sim 15\,\mathrm{Mpc}$ | Marginal --- requires exceptionally loud event |
| Cosmic Explorer | $\sim 300\,\mathrm{Mpc}$ | Yes --- $\sim 10$ events/year |
| Einstein Telescope | $\sim 200\,\mathrm{Mpc}$ | Yes --- $\sim 5$ events/year |

Current LIGO sensitivity is insufficient to resolve the duck $m$-splitting except for the most nearby events.  However, next-generation detectors (Cosmic Explorer [8], Einstein Telescope [9]) will achieve order-of-magnitude improvements in sensitivity at $\sim 2\,\mathrm{kHz}$, bringing duck QNM spectroscopy into the realm of feasibility.

**(d) The diagnostic power of QNM spectroscopy.** --- If a post-merger ringdown signal exhibits anomalous frequency splitting inconsistent with the Kerr spectrum or standard neutron star QNMs, the duck hypothesis provides a natural explanation.  The gravitational Zeeman pattern --- five split lines from the $\ell = 2$ $f$-mode, with splittings proportional to the deformation amplitudes $\varepsilon_{2,m}$ --- constitutes a unique fingerprint of non-spherical geometry.  Extracting this pattern from ringdown data would simultaneously measure the duck's deformation spectrum, providing a direct probe of compact object morphology.


## 8.3 Breaking the Degeneracy (Speculative)

We briefly discuss four approaches that could, in principle, distinguish ducks from spherical compact objects.  All are currently beyond observational reach.

**1. Gravitational wave memory.** --- The nonlinear GW memory from a binary merger encodes the full multipole structure of the system, including contributions from body multipole moments.  The duck-to-sphere difference in the memory signal scales as $\delta h_{\mathrm{mem}} / h_{\mathrm{mem}} \sim \kappa_{\mathrm{duck}} (v/c)^4$, requiring SNR $\gtrsim 100$ in the memory channel --- beyond current capabilities but potentially accessible to third-generation detector networks [24].

**2. Multi-band observations.** --- Combining PTA observations ($f \sim \mathrm{nHz}$) with LISA ($f \sim \mathrm{mHz}$) and ground-based detectors ($f \sim 10$--$10^3\,\mathrm{Hz}$) for a single source would track the inspiral from the Newtonian regime through merger, where duck effects grow as $(v/c)^4$.  The transition from the PTA to LISA band occurs over $\sim 10^4$--$10^6$ years, during which the velocity increases from $v/c \sim 0.05$ to $v/c \sim 0.3$, amplifying the duck correction by a factor of $\sim 10^3$.  Multi-band parameter estimation could in principle constrain $\kappa_{\mathrm{duck}}$ to $\mathcal{O}(10^{-3})$.

**3. Tidal disruption events.** --- A duck-shaped tidal field differs from a spherical one, potentially producing asymmetric tidal disruption event (TDE) lightcurves.  However, at the tidal disruption radius $r_T \sim R_\star (M/M_\star)^{1/3}$, the monopole contribution to the tidal field dominates the quadrupole by $(R/r_T)^2 \ll 1$, and the duck correction to the TDE rate is negligible.

**4. Direct imaging.** --- As shown in Section 8.1.7, EHT angular resolution is insufficient by $\sim 5$ orders of magnitude.  Even hypothetical space-based VLBI at sub-millimeter wavelengths with Solar System baselines ($\sim 1\,\mathrm{AU}$) would achieve $\sim 10^{-4}\,\mu\mathrm{as}$ --- still marginal for resolving supermassive ducks at cosmological distances.

In summary, no currently available or near-future observation can definitively break the degeneracy between duck-shaped and spherical compact objects.  The duck hypothesis remains unfalsified.


---

**Figures referenced in this section:**

- **Fig. 15:** Duck GWB spectrum ($h_c$ vs $f$) overlaid on NANOGrav 15-year data.  The duck and SMBHB predictions are visually identical.
- **Fig. 16:** Mock PTA skymap (Mollweide projection) showing five loud supermassive binary duck sources against an isotropic stochastic background.
- **Fig. 17:** Injected duck ringdown waveform in simulated LIGO noise, showing $m$-split QNM frequencies.

---

**References cited in this section:**

[4] G. Agazie *et al.* (NANOGrav Collaboration), "The NANOGrav 15 yr data set: Evidence for a gravitational-wave background," Astrophys. J. Lett. **951**, L8 (2023).

[5] J. Antoniadis *et al.* (EPTA Collaboration), "The second data release from the European Pulsar Timing Array," Astron. Astrophys. **678**, A50 (2023).

[6] D. J. Reardon *et al.* (PPTA Collaboration), "Search for an isotropic gravitational-wave background with the Parkes Pulsar Timing Array," Astrophys. J. Lett. **951**, L6 (2023).

[7] R. W. Hellings and G. S. Downs, "Upper limits on the isotropic gravitational radiation background from pulsar timing analysis," Astrophys. J. **265**, L39 (1983).

[8] Cosmic Explorer Collaboration, https://cosmicexplorer.org.

[9] Einstein Telescope Collaboration, https://www.et-gw.eu.

[10] B. P. Abbott *et al.* (LIGO Scientific and Virgo Collaborations), "GW170817: Observation of gravitational waves from a binary neutron star inspiral," Phys. Rev. Lett. **119**, 161101 (2017).

[11] B. P. Abbott *et al.* (LIGO Scientific and Virgo Collaborations), "GW170817: Measurements of neutron star radii and equation of state," Phys. Rev. Lett. **121**, 161101 (2018).

[15] E. S. Phinney, "A practical theorem on gravitational wave backgrounds," arXiv:astro-ph/0108028 (2001).

[16] L. C. Stein and N. Yunes, "Mapping the gravitational-wave sky with pulsar timing arrays," arXiv:2507.21380 (2025).

[17] Event Horizon Telescope Collaboration, "First M87 Event Horizon Telescope results. I.," Astrophys. J. Lett. **875**, L1 (2019).

[18] E. E. Flanagan and T. Hinderer, "Constraining neutron-star tidal Love numbers with gravitational-wave detectors," Phys. Rev. D **77**, 021502 (2008).

[19] T. Hinderer *et al.*, "Tidal deformability of neutron stars with realistic equations of state and their gravitational wave signatures in binary inspiral," Phys. Rev. D **81**, 123016 (2010).

[20] K. D. Kokkotas and B. G. Schmidt, "Quasi-normal modes of stars and black holes," Living Rev. Rel. **2**, 2 (1999).

[21] N. Andersson and K. D. Kokkotas, "Towards gravitational wave asteroseismology," MNRAS **299**, 1059 (1998).

[22] E. Berti, V. Cardoso, and C. M. Will, "Gravitational-wave spectroscopy of massive black holes with the space interferometer LISA," Phys. Rev. D **73**, 064030 (2006).

[23] E. Berti, J. Cardoso, V. Cardoso, and M. Cavaglia, "Matched-filter searches for gravitational waves from ringdown of astrophysical black holes," Phys. Rev. D **76**, 104044 (2007).

[24] M. Favata, "Nonlinear gravitational-wave memory from binary black hole mergers," Astrophys. J. **696**, L159 (2009).
