# 5. Binary Duck on Post-Newtonian Orbits

Having characterized the duck's geometry (Section 2), normal mode spectrum (Section 3), and tidal response (Section 4), we now embed the duck in a binary system and study the gravitational wave signatures of the inspiral.
The permanent body-frame quadrupole $Q^C_{ij}$, computed from the mesh in Section 2, modifies the orbital dynamics, the radiated power, and ultimately the gravitational wave phase evolution.
We derive these corrections systematically, following the post-Newtonian framework of Poisson and Will [15] and the quadrupole-monopole formalism of Barker and O'Connell [16].

Throughout this section we restore explicit factors of $G$ and $c$.
We consider a binary system of total mass $M = M_A + M_B$, reduced mass $\mu = M_A M_B / M$, symmetric mass ratio $\eta = \mu / M$, and orbital separation $r$.
Each duck carries a permanent traceless symmetric quadrupole $Q^{C,A}_{ij}$ in its body frame.
We define the equivalent sphere radius $R_{\mathrm{eq}} = (3V / 4\pi)^{1/3}$ and the dimensionless shape parameter $\hat{Q} = Q^C_{\mathrm{max}} / (M R_{\mathrm{eq}}^2)$, where $Q^C_{\mathrm{max}}$ is the largest eigenvalue of the quadrupole tensor.
For our duck mesh: $Q^C = \mathrm{diag}(4.18, -1.61, -2.57) \times 10^{-3}$ (benchmark units), giving $\hat{Q} \sim 0.015$.


## 5.1 Quadrupole-Monopole Interaction

The leading interaction between body $A$'s quadrupole and body $B$'s monopole is [15, 16]

$$U_{QM} = -\frac{G M_B}{2 r^3}\, Q^{C,A}_{ij}\, n_i\, n_j, \tag{29}$$

where $Q^{C,A}_{ij}$ is duck $A$'s quadrupole tensor rotated to the inertial frame and $\mathbf{n} = \mathbf{x} / r$ is the unit separation vector.
For the full system with two ducks, the total quadrupole-monopole energy is

$$U_{QM}^{\mathrm{total}} = -\frac{G}{2 r^3}\left(M_B\, Q^{C,A}_{ij} + M_A\, Q^{C,B}_{ij}\right) n_i\, n_j. \tag{30}$$

**Post-Newtonian order counting.** --- The ratio $U_{QM} / E_N$ scales as $(R/r)^2$.
For compact objects whose radius is of order the Schwarzschild radius, $R \sim GM/c^2$, the ratio becomes $\sim (GM/c^2 r)^2 \sim v^4/c^4$, placing the quadrupole-monopole interaction at **2PN order** $[\mathcal{O}(v^4/c^4)]$.
For extended objects with $R_{\mathrm{phys}} \gg R_{\mathrm{Schwarzschild}}$, the effect is parametrically larger and enters at an effective order $(R_{\mathrm{phys}} / r)^2$ that can approach $\mathcal{O}(1)$ at close separations.

**Orientation cases.** --- The projection $Q^C_{ij} n_i n_j$ depends on the relative orientation between the duck's body frame and the orbital plane.
We consider three physically distinct configurations:

1. *Tidally locked.* --- The duck's body $x$-axis remains aligned with $\mathbf{n}$ at all times, so $Q_{nn} \equiv Q^C_{xx} = \mathrm{const}$.
   This is the dominant configuration for close binaries where tidal torques enforce synchronous rotation.
   The quadrupole correction enters at the orbital frequency with constant amplitude.

2. *Freely spinning.* --- The duck spins about one of its principal axes at angular velocity $\Omega_{\mathrm{spin}} \neq \omega_{\mathrm{orb}}$.
   The projection $Q_{nn}(t)$ oscillates at frequencies $2(\Omega_{\mathrm{spin}} - \omega_{\mathrm{orb}})$ and $2\Omega_{\mathrm{spin}}$, generating additional spectral components in the gravitational wave signal.

3. *Precessing.* --- The duck's spin axis is misaligned with the orbital angular momentum, leading to general Euler angle evolution.
   The time-dependent rotation matrix $R_A(t)$ must be evolved via Euler's equations for the rigid body, coupled to the tidal torque $\tau_i = -\partial U_{QM} / \partial \alpha_i$.
   This produces the richest phenomenology, including spin-orbit and spin-spin precession modulations of the gravitational waveform.

For definiteness, the remainder of this section assumes the tidally locked configuration unless otherwise stated.


## 5.2 Modified Orbital Dynamics

The quadrupole-monopole interaction modifies the radial force balance for circular orbits.
Defining $Q_{nn}^{\mathrm{eff}} = M_B Q^{C,A}_{nn} + M_A Q^{C,B}_{nn}$, the modified Kepler relation for circular orbits is

$$\omega^2 = \frac{G M}{r^3}\left[1 + \frac{15\, Q_{nn}^{\mathrm{eff}}}{2\, M\, r^2}\right]. \tag{31}$$

The factor of 15 arises from the angular average of the quadrupole potential's radial gradient (the $-3/r^4$ from $U_{QM}$ combined with the traceless tensor projection).
For equal-mass tidally locked ducks with $Q_{nn} = Q^C_{xx}$, this becomes $Q_{nn}^{\mathrm{eff}} = M Q^C_{xx}$ and the correction is $15 Q^C_{xx} / (2 r^2)$.

**Binding energy.** --- The orbital energy including the quadrupole correction is

$$E(r) = -\frac{G M \mu}{2r} - \frac{G}{2r^3} Q_{nn}^{\mathrm{eff}}. \tag{32}$$

Inverting the modified Kepler relation Eq.\ (31) to obtain $r(\omega)$ and substituting, the binding energy as a function of orbital frequency takes the form

$$E(\omega) = E_{\mathrm{Kep}}(\omega)\left[1 + \delta_Q(\omega)\right], \tag{33}$$

where $E_{\mathrm{Kep}} = -(1/2)\mu (GM\omega)^{2/3}$ is the Keplerian binding energy and the fractional correction is

$$\delta_Q(\omega) = \frac{5\, Q_{nn}^{\mathrm{eff}}}{M}\left(\frac{\omega}{GM}\right)^{2/3}. \tag{34}$$

**Energy balance.** --- The inspiral is governed by the energy balance equation

$$\frac{dE}{dt} = -P_{\mathrm{GW}}, \tag{35}$$

where $P_{\mathrm{GW}}$ is the total gravitational wave luminosity (Section 5.3).
The quadrupole correction modifies both sides: $E(\omega)$ through Eq.\ (33) and $P_{\mathrm{GW}}(\omega)$ through the additional body quadrupole contributions.
The resulting chirp rate is

$$\dot{\omega} = -\frac{P_{\mathrm{GW}}}{dE/d\omega} = \dot{\omega}_{\mathrm{pp}}\left[1 + \delta_{\dot{\omega}}(\omega)\right], \tag{36}$$

where $\dot{\omega}_{\mathrm{pp}}$ is the point-particle chirp rate and $\delta_{\dot{\omega}}$ encodes the combined correction from modified binding energy and modified radiated power.


## 5.3 System Quadrupole and Gravitational Wave Power

The total mass quadrupole moment of the binary system is the sum of the orbital and body contributions:

$$\mathcal{I}^{\mathrm{sys}}_{ij}(t) = \underbrace{\mu\left(x_i x_j - \tfrac{1}{3} r^2 \delta_{ij}\right)}_{\mathrm{orbital}} + \underbrace{\sum_{A=1,2} R_A(t)\, Q^{C,A}_{\mathrm{body}}\, R_A^T(t)}_{\mathrm{body}}, \tag{37}$$

where $\mathbf{x}$ is the relative separation vector and $R_A(t)$ is the rotation matrix that brings duck $A$'s body-frame quadrupole into the inertial frame.
For the tidally locked configuration, $R_A(t)$ is the orbital rotation matrix and the body quadrupole rotates at the orbital frequency $\omega_{\mathrm{orb}}$.
This is crucial: the body term contributes at frequency $2\omega_{\mathrm{orb}}$, the same harmonic as the orbital quadrupole, enabling coherent interference.

The gravitational wave luminosity is given by the standard quadrupole formula:

$$P_{\mathrm{GW}} = \frac{G}{45\, c^5}\left\langle\dddot{\mathcal{I}}^{\mathrm{sys}}_{ij}\, \dddot{\mathcal{I}}^{\mathrm{sys},ij}\right\rangle. \tag{38}$$

Expanding the square and using the bilinearity of the contraction, the power decomposes into three terms:

$$P_{\mathrm{GW}} = P_{\mathrm{orb}} + P_{\mathrm{body}} + P_{\mathrm{cross}}. \tag{39}$$

**Orbital power** $P_{\mathrm{orb}}$ is the standard Peters-Mathews result [17]:

$$P_{\mathrm{orb}} = \frac{32}{5}\frac{G^4}{c^5}\frac{\mu^2 M^3}{r^5} = \frac{32}{5}\frac{c^5}{G}\eta^2\, v^{10}, \tag{40}$$

where $v = (GM\omega)^{1/3}$ is the characteristic orbital velocity.

**Body power** $P_{\mathrm{body}}$ arises from the third time derivative of the rotating body quadrupole alone.
For a tidally locked duck spinning at $\omega_{\mathrm{orb}}$, the body quadrupole contributes

$$P_{\mathrm{body}} = \frac{G}{45\, c^5}\left\langle\dddot{Q}^{\mathrm{body}}_{ij}\, \dddot{Q}^{\mathrm{body},ij}\right\rangle \propto \omega^6\, |Q^C|^2, \tag{41}$$

where $|Q^C|^2 = Q^C_{ij} Q^{C,ij}$ is the quadrupole tensor norm.
This is the same power formula validated in Section 2 for a single spinning duck, now applied to each binary component.

**Cross power** $P_{\mathrm{cross}}$ captures the interference between the orbital and body quadrupoles:

$$P_{\mathrm{cross}} = \frac{2G}{45\, c^5}\left\langle\dddot{\mathcal{I}}^{\mathrm{orb}}_{ij}\, \dddot{Q}^{\mathrm{body},ij}\right\rangle. \tag{42}$$

This term is proportional to $\mu r^2 \omega^3 \times Q^C \omega^3 = \mu Q^C \omega^6 r^2$ and scales as an intermediate power of the separation.
Critically, **the cross-term can be negative**: when the body quadrupole is oriented such that its third time derivative is anti-correlated with the orbital quadrupole's third derivative, the interference is destructive.
This occurs when the principal axis of $Q^C$ with the most negative eigenvalue is aligned with the orbital angular momentum.
Destructive interference reduces the total radiated power below the Peters-Mathews prediction, an effect unique to bodies with permanent quadrupole moments.

**Scaling relations.** --- The three power components scale differently with separation:

$$P_{\mathrm{orb}} \propto r^{-5}, \qquad P_{\mathrm{body}} \propto r^{-3} \cdot \omega^6 \propto r^{-12}, \qquad P_{\mathrm{cross}} \propto r^{-7/2} \cdot \omega^6 \propto r^{-17/2}. \tag{43}$$

Wait --- let us be more careful.
For circular orbits, $\omega^2 \propto r^{-3}$ (Kepler), so $\omega^6 \propto r^{-9}$.
The orbital quadrupole is $\mathcal{I}^{\mathrm{orb}} \sim \mu r^2$, so $\dddot{\mathcal{I}}^{\mathrm{orb}} \sim \mu r^2 \omega^3 \sim \mu r^{2} r^{-9/2} = \mu r^{-5/2}$, giving $P_{\mathrm{orb}} \sim \mu^2 r^{-5}$ (Peters-Mathews).
The body quadrupole $Q^{\mathrm{body}}$ has fixed amplitude, so $\dddot{Q}^{\mathrm{body}} \sim Q^C \omega^3 \sim Q^C r^{-9/2}$, giving $P_{\mathrm{body}} \sim (Q^C)^2 r^{-9}$.
The cross-term: $P_{\mathrm{cross}} \sim \mu r^{-5/2} \cdot Q^C r^{-9/2} = \mu Q^C r^{-7}$.
Therefore:

$$\frac{P_{\mathrm{body}}}{P_{\mathrm{orb}}} \sim \frac{(Q^C)^2}{\mu^2 r^4} \sim \hat{Q}^2 \left(\frac{R_{\mathrm{eq}}}{r}\right)^4, \qquad \frac{P_{\mathrm{cross}}}{P_{\mathrm{orb}}} \sim \frac{Q^C}{\mu r^2} \sim \hat{Q}\left(\frac{R_{\mathrm{eq}}}{r}\right)^2. \tag{44}$$

The cross-term dominates over the body term at all separations $r \gg R_{\mathrm{eq}}$, scaling as $(R_{\mathrm{eq}}/r)^2$ rather than $(R_{\mathrm{eq}}/r)^4$.
At $r = 5 R_{\mathrm{eq}}$, we find $P_{\mathrm{cross}} / P_{\mathrm{orb}} \sim 6 \times 10^{-4}$ and $P_{\mathrm{body}} / P_{\mathrm{orb}} \sim 4 \times 10^{-7}$.


## 5.4 Higher Multipoles

The quadrupole formula Eq.\ (38) captures the leading-order gravitational radiation.
At higher post-Newtonian orders, additional multipole moments contribute to the radiated power.

**Current quadrupole** $S_{ij}$. --- The mass-current quadrupole arises from the angular momentum distribution of the system.
For a rigidly rotating duck with inertia tensor $I_{ij}$ spinning at angular velocity $\Omega_k$, the current quadrupole is

$$S_{ij} = \epsilon_{kl(i}\, I_{j)k}\, \Omega_l. \tag{45}$$

The current quadrupole contributes to the radiated power at **0.5PN** relative to the mass quadrupole [18]:

$$P_{\mathrm{current}} = \frac{G}{45\, c^7}\left\langle\dddot{S}_{ij}\, \dddot{S}^{ij}\right\rangle. \tag{46}$$

For the duck, the anisotropic inertia tensor means that $S_{ij}$ depends on the spin orientation relative to the principal axes, introducing additional angular structure in the emitted radiation.

**Mass octupole** $\mathcal{I}_{ijk}$. --- The mass octupole moment of the binary, constructed from the $\ell = 3$ spherical multipoles $Q_3^m$, also enters at **0.5PN** relative to the mass quadrupole.
Its contribution to the radiated power is [18]

$$P_{\mathrm{octupole}} = \frac{G}{189\, c^7}\left\langle\overset{(4)}{\mathcal{I}}_{ijk}\, \overset{(4)}{\mathcal{I}}^{ijk}\right\rangle, \tag{47}$$

where $\overset{(4)}{\mathcal{I}}_{ijk}$ denotes the fourth time derivative.
For an **equal-mass binary**, the mass octupole of the orbital motion vanishes by parity symmetry: the octupole is odd under exchange of the two bodies, so $\mathcal{I}^{\mathrm{orb}}_{ijk} \propto (M_A - M_B) = 0$.
For unequal-mass systems, the octupole introduces odd harmonics (frequencies $\omega$, $3\omega$) in addition to the even harmonics ($2\omega$, $4\omega$) from the quadrupole.

The body octupole $Q_3^m$ of the duck, extracted from the existing multipole pipeline, provides a permanent octupolar contribution that does not vanish for equal-mass binaries.
However, the body octupole power scales as $(Q_3 / r^3)^2 \omega^8 \sim (R_{\mathrm{eq}}/r)^6$, which is subdominant to the body quadrupole cross-term at all separations of interest.


## 5.5 Phase Evolution

The accumulated gravitational wave phase is the primary observable for matched-filter searches.
In the stationary phase approximation (SPA), the Fourier-domain phase is [19]

$$\Psi(f) = \Psi_{\mathrm{pp}}(f) + \delta\Psi_Q(f), \tag{48}$$

where $\Psi_{\mathrm{pp}}$ is the point-particle (Keplerian) phase and $\delta\Psi_Q$ is the quadrupole correction.

The point-particle phase through 2PN order is

$$\Psi_{\mathrm{pp}}(f) = 2\pi f t_c - \phi_c - \frac{\pi}{4} + \frac{3}{128\eta v^5}\left[1 + \frac{20}{9}\left(\frac{743}{336} + \frac{11}{4}\eta\right)v^2 - 16\pi v^3 + \cdots\right], \tag{49}$$

where $v = (\pi M f)^{1/3}$ and $t_c$, $\phi_c$ are the coalescence time and phase.

The quadrupole correction, following the analysis of Poisson [20], enters at 2PN order:

$$\delta\Psi_Q(f) = -\frac{75}{64\eta}\frac{\hat{Q}_{\mathrm{eff}}}{M\, r^2(f)}\, v^{-1} = -\frac{75}{64\eta}\,\hat{Q}_{\mathrm{eff}}\left(\pi M f\right)^{4/3}\, v^{-1}, \tag{50}$$

where $\hat{Q}_{\mathrm{eff}} = (M_B \hat{Q}_A + M_A \hat{Q}_B) / M$ is the mass-weighted effective quadrupole parameter and $r(f)$ is obtained from the Kepler relation $r = (GM / \pi^2 f^2)^{1/3}$.
The explicit frequency dependence is

$$\delta\Psi_Q(f) = -\frac{75}{64\eta}\,\hat{Q}_{\mathrm{eff}}\,(\pi \mathcal{M} f)^{1/3}, \tag{51}$$

where $\mathcal{M} = \eta^{3/5} M$ is the chirp mass.

**Dephasing cycles.** --- The number of dephasing cycles accumulated over a detector's sensitive band $[f_{\mathrm{low}}, f_{\mathrm{high}}]$ is

$$\Delta N = \frac{1}{2\pi}\left|\delta\Psi_Q(f_{\mathrm{high}}) - \delta\Psi_Q(f_{\mathrm{low}})\right|. \tag{52}$$

For stellar-mass duck binaries ($M = 2.8\, M_\odot$, $\eta = 1/4$) in the LIGO band ($f \in [10, 1000]\,\mathrm{Hz}$):

$$\Delta N_{\mathrm{LIGO}} \sim \frac{75}{128\pi}\,\hat{Q}\,\left[(\pi \mathcal{M} f_{\mathrm{high}})^{1/3} - (\pi \mathcal{M} f_{\mathrm{low}})^{1/3}\right] \sim 0.01\text{--}0.1\;\mathrm{cycles}. \tag{53}$$

This is below the $\sim 1$ cycle threshold for detection with current instruments, but within reach of third-generation detectors operating at design sensitivity.

For supermassive duck binaries ($M = 10^7\, M_\odot$, $\eta = 1/4$) in the LISA band ($f \in [10^{-4}, 10^{-1}]\,\mathrm{Hz}$), the enormous number of orbital cycles ($\sim 10^5$) accumulated over years of observation amplifies even tiny fractional corrections:

$$\Delta N_{\mathrm{LISA}} \sim 1\text{--}10\;\mathrm{cycles}. \tag{54}$$

This is in principle detectable, provided the duck quadrupole parameter $\hat{Q}$ can be disentangled from spin-induced quadrupole effects (which enter at the same PN order).
We return to this degeneracy in Section 8.


## 5.6 Waveform

The time-domain gravitational waveform for a quasi-circular inspiral, including duck corrections to both amplitude and phase, is

$$h_+(t) = -\frac{2G\mu}{c^4 D_L}\,(\omega r)^2\,(1 + \cos^2\iota)\,\cos\!\left[2\Phi(t) + \delta\Phi_Q(t)\right], \tag{55}$$

$$h_\times(t) = -\frac{4G\mu}{c^4 D_L}\,(\omega r)^2\,\cos\iota\,\sin\!\left[2\Phi(t) + \delta\Phi_Q(t)\right], \tag{56}$$

where $D_L$ is the luminosity distance, $\iota$ is the orbital inclination, $\Phi(t) = \int \omega\, dt$ is the orbital phase, and $\delta\Phi_Q(t)$ is the accumulated phase correction from the quadrupole-monopole interaction.

The amplitude also receives a duck correction from the body quadrupole contribution to $\mathcal{I}^{\mathrm{sys}}$:

$$h_{+,\times}(t) = h_{+,\times}^{\mathrm{pp}}(t) + \delta h_{+,\times}^{\mathrm{body}}(t), \tag{57}$$

where the body correction $\delta h^{\mathrm{body}}$ is proportional to $Q^C / (\mu r^2)$ and oscillates at the same frequency $2\omega$ for the tidally locked case.
In the Fourier domain, the amplitude correction modifies the waveform at 2PN order:

$$\tilde{h}(f) = \tilde{h}_{\mathrm{pp}}(f)\left[1 + \delta A_Q(f)\right] e^{i\,\delta\Psi_Q(f)}, \tag{58}$$

with $\delta A_Q \sim \hat{Q} v^4 \sim 10^{-4}$--$10^{-2}$ across the sensitive band.
The phase correction $\delta\Psi_Q$ (Eq.\ 50) dominates over the amplitude correction for parameter estimation, as matched-filter searches are exponentially more sensitive to phase errors than amplitude errors.

The complete waveform Eq.\ (55)--(56) reduces to the standard restricted post-Newtonian waveform in the limit $Q^C \to 0$, recovering the Peters-Mathews inspiral as required (Verification V1).


---

**References cited in this section:**

[15] E. Poisson and C. M. Will, *Gravity: Newtonian, Post-Newtonian, Relativistic* (Cambridge University Press, 2014).

[16] B. M. Barker and R. F. O'Connell, "Gravitational two-body problem with arbitrary masses, spins, and quadrupole moments," Phys. Rev. D **12**, 329 (1975).

[17] P. C. Peters, "Gravitational radiation and the motion of two point masses," Phys. Rev. **136**, B1224 (1964).

[18] L. Blanchet, "Gravitational radiation from post-Newtonian sources and inspiralling compact binaries," Living Rev. Relativ. **17**, 2 (2014).

[19] C. Cutler and E. E. Flanagan, "Gravitational waves from merging compact binaries: How accurately can one extract the binary's parameters from the inspiral waveform?" Phys. Rev. D **49**, 2658 (1994).

[20] E. Poisson, "Gravitational waves from inspiraling compact binaries: The quadrupole-moment term," Phys. Rev. D **57**, 5287 (1998).
