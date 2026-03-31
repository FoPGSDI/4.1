# Section 3: Quacky Normal Modes --- Perturbative Splitting Theory

The oscillation modes of a spherical compact object are classified by
radial order $n$ and angular harmonic $(\ell,m)$, with each $(n,\ell)$
level $(2\ell+1)$-fold degenerate.  A nonspherical deformation lifts
this degeneracy, producing a *gravitational Zeeman splitting* of the
quasinormal mode (QNM) spectrum.  In this section we develop
first-order perturbation theory for the mode frequencies of a
duck-shaped object, starting from the Hadamard variation formula for
eigenvalues under domain deformation and reducing the splitting
computation to Gaunt integrals expressible via Wigner 3j symbols.  The
resulting selection rules and explicit matrix elements are the
foundation for both the analytical predictions tested against FEM
eigenvalues in Sec. 6 and the ringdown waveform templates constructed
in Sec. 7.


## 3.1 Problem Setup

We model the oscillation modes of a compact duck by the scalar
Helmholtz eigenvalue problem on the interior domain $D$ bounded by the
duck surface $\partial D$, with Dirichlet boundary conditions:

$$-\nabla^2 \psi_{n\ell m} = \omega_{n\ell}^2 \, \psi_{n\ell m}, \qquad \psi \big|_{\partial D} = 0.$$
<!-- Eq. (3.1) -->

This is the Cowling-approximation analog of the fluid f-mode and
p-mode problem, capturing the essential spectral structure while
admitting exact unperturbed solutions [1,2].

**Unperturbed spectrum.**
On the sphere of radius $R_0$, the eigenfunctions separate as
$\psi_{n\ell m}(r,\theta,\phi) = R_{n\ell}(r) \, Y_\ell^m(\theta,\phi)$,
where $R_{n\ell}(r) = j_\ell(\omega_{n\ell} r)$ is a spherical Bessel
function.  The Dirichlet condition $j_\ell(\omega_{n\ell} R_0) = 0$
fixes the eigenfrequencies:

$$\omega_{n,\ell} = \frac{z_{n,\ell}}{R_0},$$
<!-- Eq. (3.2) -->

where $z_{n,\ell}$ is the $n$-th positive zero of $j_\ell$.  Each
$(n,\ell)$ level is $(2\ell+1)$-fold degenerate: the frequency is
independent of the azimuthal order $m$.  The first few zeros are
$z_{1,0} = \pi$, $z_{1,1} \approx 4.493$, $z_{1,2} \approx 5.763$,
$z_{1,3} \approx 6.988$.

**Duck boundary.**
We parametrize the duck surface as a deformation of the
equivalent-radius sphere (Sec. 2):

$$R_{\mathrm{duck}}(\theta,\phi) = R_0 \left[ 1 + \sum_{\ell' \geq 1} \sum_{m'=-\ell'}^{\ell'} \varepsilon_{\ell' m'} \, Y_{\ell'}^{m'}(\theta,\phi) \right].$$
<!-- Eq. (3.3) -->

The deformation coefficients $\varepsilon_{\ell' m'}$ are extracted
from the duck mesh via regularized least-squares SH fitting (Sec. 2);
the expansion is real-valued with $\varepsilon_{\ell',-m'} =
(-1)^{m'} \varepsilon_{\ell',m'}^*$.  The normal displacement from the
sphere is $\delta n(\theta,\phi) = R_0 \sum_{\ell',m'}
\varepsilon_{\ell' m'} Y_{\ell'}^{m'}(\theta,\phi)$ to first order.


## 3.2 Hadamard Formula for Domain Perturbation

The Hadamard variational formula gives the first-order sensitivity of
a Dirichlet eigenvalue to a normal displacement $\delta n$ of the
boundary [3,4].  For a simple eigenvalue $\omega^2$ with eigenfunction
$\psi$ normalized to unit $L^2$ norm on $D$, the variation is

$$\delta \omega^2 = -\int_{\partial D} \left| \frac{\partial \psi}{\partial \mathbf{n}} \right|^2 \delta n \, dS,$$
<!-- Eq. (3.4) -->

where $\partial \psi / \partial \mathbf{n}$ is the outward normal
derivative on the boundary.  This classical result dates to
Hadamard [3]; a modern rigorous treatment including the degenerate
case and higher-order corrections is given by Suzuki and
Tsuchiya [5], who establish the formula via shape calculus on
Sobolev spaces and prove analyticity of the eigenvalue branches
with respect to the deformation parameter.

**Application to the duck deformation.**
Substituting the duck boundary displacement
$\delta n = R_0 \sum_{\ell',m'} \varepsilon_{\ell' m'} Y_{\ell'}^{m'}$
into Eq. (3.4), and evaluating on the sphere $\partial D_0 = S^2(R_0)$
where $dS = R_0^2 \, d\Omega$, we obtain

$$\delta \omega_{n\ell m}^2 = -R_0^3 \sum_{\ell',m'} \varepsilon_{\ell' m'} \int_{S^2} \left| \frac{\partial \psi_{n\ell m}}{\partial r} \bigg|_{r=R_0} \right|^2 Y_{\ell'}^{m'}(\theta,\phi) \, d\Omega.$$
<!-- Eq. (3.5) -->

The radial derivative on the sphere is
$\partial \psi_{n\ell m} / \partial r |_{R_0} = R'_{n\ell}(R_0) \, Y_\ell^m(\theta,\phi)$,
where $R'_{n\ell}(R_0) = j'_\ell(z_{n,\ell}) \cdot (1/R_0)$.  Equation (3.5)
therefore involves the product of three spherical harmonics integrated
over the sphere --- the Gaunt integral.


## 3.3 Degenerate Perturbation Theory

Since the unperturbed $(n,\ell)$ eigenvalue is $(2\ell+1)$-fold
degenerate, we must apply degenerate perturbation theory.  The
first-order frequency shifts within the multiplet are determined by
diagonalizing the $(2\ell+1) \times (2\ell+1)$ perturbation matrix

$$V_{m_1, m_2}^{(n,\ell)} = -\left| R'_{n\ell}(R_0) \right|^2 R_0^3 \sum_{\ell',m'} \varepsilon_{\ell' m'} \int_{S^2} Y_\ell^{m_1 *} \, Y_\ell^{m_2} \, Y_{\ell'}^{m'} \, d\Omega.$$
<!-- Eq. (3.6) -->

The angular integral is the *Gaunt integral*, which admits a
closed-form expression in terms of Wigner 3j symbols [6]:

$$\int_{S^2} Y_\ell^{m_1 *} \, Y_\ell^{m_2} \, Y_{\ell'}^{m'} \, d\Omega = (-1)^{m_1} \sqrt{\frac{(2\ell+1)^2 (2\ell'+1)}{4\pi}} \begin{pmatrix} \ell & \ell & \ell' \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} \ell & \ell & \ell' \\ -m_1 & m_2 & m' \end{pmatrix}.$$
<!-- Eq. (3.7) -->

Substituting into Eq. (3.6), the perturbation matrix becomes

$$V_{m_1,m_2}^{(n,\ell)} = -\left| R'_{n\ell}(R_0) \right|^2 R_0^3 \sum_{\ell',m'} \varepsilon_{\ell' m'} \, (-1)^{m_1} \sqrt{\frac{(2\ell+1)^2(2\ell'+1)}{4\pi}} \begin{pmatrix} \ell & \ell & \ell' \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} \ell & \ell & \ell' \\ -m_1 & m_2 & m' \end{pmatrix}.$$
<!-- Eq. (3.8) -->

The eigenvalues $\{ \delta\omega^2_\alpha \}_{\alpha=1}^{2\ell+1}$ of
$V$ give the first-order frequency shifts, and the eigenvectors give
the linear combinations of $Y_\ell^m$ that diagonalize the
perturbation --- the "good quantum numbers" of the duck.

The perturbed frequencies are then

$$\omega_{n\ell\alpha}^2 = \omega_{n\ell}^{(0)\,2} + \delta\omega^2_\alpha + \mathcal{O}(\varepsilon^2),$$
<!-- Eq. (3.9) -->

or equivalently, to first order in $\varepsilon$,

$$\omega_{n\ell\alpha} = \omega_{n\ell}^{(0)} + \frac{\delta\omega^2_\alpha}{2\,\omega_{n\ell}^{(0)}} + \mathcal{O}(\varepsilon^2).$$
<!-- Eq. (3.10) -->


## 3.4 Selection Rules

The Wigner 3j symbols in Eq. (3.8) impose stringent selection rules on
which deformation harmonics $\varepsilon_{\ell' m'}$ contribute to the
splitting of a given $(n,\ell)$ multiplet.

**(i) Triangle inequality.**
The 3j symbol
$\bigl(\begin{smallmatrix} \ell & \ell & \ell' \\ 0 & 0 & 0 \end{smallmatrix}\bigr)$
vanishes unless $|\ell - \ell| \leq \ell' \leq \ell + \ell$, i.e.,

$$0 \leq \ell' \leq 2\ell.$$
<!-- Eq. (3.11) -->

**(ii) Parity.**
The same 3j symbol with all-zero bottom row vanishes unless
$\ell + \ell + \ell' = 2\ell + \ell'$ is even, which requires

$$\ell' \in \{0, 2, 4, \ldots, 2\ell\} \quad \text{(even only)}.$$
<!-- Eq. (3.12) -->

Odd-$\ell'$ deformation harmonics do not contribute to first-order
splitting.  This is physically expected: odd-$\ell'$ deformations
(such as the dipole $\ell'=1$, which displaces the center of mass)
cannot split eigenvalues at first order because they do not break the
left-right symmetry of the mode pattern.

**(iii) Azimuthal selection.**
The second 3j symbol requires

$$m' = m_1 - m_2.$$
<!-- Eq. (3.13) -->

The diagonal elements ($m_1 = m_2$) receive contributions only from
$m' = 0$ deformations, while off-diagonal coupling ($m_1 \neq m_2$)
requires $m' \neq 0$.

**Consequences for the duck.**
The hierarchy of contributions is:

- $\varepsilon_{0,0}$: produces an overall shift (rescaling of $R_0$), no splitting.
- $\varepsilon_{2,m}$: *primary splitting* --- the quadrupolar deformation dominates the lifting of degeneracy.
- $\varepsilon_{4,m}$: *secondary splitting* --- hexadecapole corrections refine the pattern.
- $\varepsilon_{6,m}$ and higher: contribute only to $\ell \geq 3$ multiplets and are progressively smaller.

For the duck, the measured deformation spectrum (Tab. 2) shows
$|\varepsilon_{2,m}| \sim 0.1$--$0.3$ and $|\varepsilon_{4,m}| \sim
0.03$--$0.05$, confirming that the quadrupole dominates.


## 3.5 Explicit $\ell = 2$ f-Mode Splitting

The astrophysically most important case is the $\ell = 2$ f-mode
(fundamental quadrupole oscillation, $n = 1$), which is the dominant
channel for gravitational wave emission during ringdown.  The 5-fold
degenerate multiplet splits via the $5 \times 5$ perturbation matrix.

From the selection rules (Sec. 3.4), only $\ell' = 0, 2, 4$
contribute, and since $\ell' = 0$ gives a uniform shift, the splitting
is controlled by $\varepsilon_{2,m'}$ (primary) and
$\varepsilon_{4,m'}$ (secondary).

**Restricting to $\ell' = 2$ (dominant contribution).**
The perturbation matrix is

$$V_{m_1,m_2} = -\left| R'_{1,2}(R_0) \right|^2 R_0^3 \sum_{m'=-2}^{2} \varepsilon_{2,m'} \, (-1)^{m_1} \sqrt{\frac{25 \cdot 5}{4\pi}} \begin{pmatrix} 2 & 2 & 2 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} 2 & 2 & 2 \\ -m_1 & m_2 & m' \end{pmatrix},$$
<!-- Eq. (3.14) -->

with the constraint $m' = m_1 - m_2$ from the second 3j symbol.

The key numerical value of the all-zero 3j symbol is

$$\begin{pmatrix} 2 & 2 & 2 \\ 0 & 0 & 0 \end{pmatrix} = \sqrt{\frac{2}{35}}.$$
<!-- Eq. (3.15) -->

Defining the overall prefactor

$$\mathcal{A} \equiv -\left| R'_{1,2}(R_0) \right|^2 R_0^3 \sqrt{\frac{25 \cdot 5}{4\pi}} \cdot \sqrt{\frac{2}{35}} = -\left| R'_{1,2}(R_0) \right|^2 R_0^3 \sqrt{\frac{50}{28\pi}},$$
<!-- Eq. (3.16) -->

the matrix elements become

$$V_{m_1,m_2} = \mathcal{A} \, (-1)^{m_1} \, \varepsilon_{2,\, m_1-m_2} \begin{pmatrix} 2 & 2 & 2 \\ -m_1 & m_2 & m_1 - m_2 \end{pmatrix},$$
<!-- Eq. (3.17) -->

where $\varepsilon_{2,m'}$ is understood to vanish for $|m'| > 2$.
The 3j symbols
$\bigl(\begin{smallmatrix} 2 & 2 & 2 \\ -m_1 & m_2 & m_1-m_2 \end{smallmatrix}\bigr)$
are tabulated in standard references [6,7]; their explicit values
give a $5 \times 5$ Hermitian matrix that is straightforwardly
diagonalized numerically.

For an axisymmetric deformation ($\varepsilon_{2,0} \neq 0$,
$\varepsilon_{2,m\neq 0} = 0$), the matrix is diagonal with entries
proportional to

$$\begin{pmatrix} 2 & 2 & 2 \\ -m & m & 0 \end{pmatrix},$$

giving splitting proportional to $m^2$: modes with $|m| = 2$ are
shifted most, $|m| = 1$ intermediate, and $m = 0$ least, producing the
pattern $\omega_{|m|=2} > \omega_{|m|=1} > \omega_{m=0}$ (or reversed,
depending on the sign of $\varepsilon_{2,0}$).  Non-axisymmetric
deformations ($\varepsilon_{2,\pm 1}$, $\varepsilon_{2,\pm 2}$)
introduce off-diagonal couplings and mix the $m$-substates.

**Numerical example.**
For the duck with $\varepsilon_{2,0} = -0.3$, $\varepsilon_{2,2} =
0.15$, and $\varepsilon_{4,0} = 0.05$ (representative values from the
SH fit), diagonalization of the $5 \times 5$ matrix yields five
distinct frequencies.  The fractional splitting $\Delta\omega/\omega
\sim |\varepsilon_{2,0}| \sim 0.3$ is large --- comparable to or
exceeding the linewidth of the modes --- ensuring that the split lines
are spectrally resolved.


## 3.6 Physical Interpretation

The mapping between duck morphological features and deformation
harmonics provides a direct physical interpretation of the QNM
splitting pattern.

**Axial elongation ($\varepsilon_{2,0}$).**
The $m = 0$ quadrupole deformation measures the prolate/oblate
distortion of the duck body along the symmetry axis.  This is the
dominant deformation coefficient and produces the largest splitting,
separating modes by $|m|$: the $m = 0$ mode (oscillation along the
body axis) has a different frequency from $|m| = 1$ (tilting) and
$|m| = 2$ (equatorial squeezing).  For the duck, the elongated body
gives $\varepsilon_{2,0} < 0$ (prolate), pushing $|m| = 2$ to higher
frequencies.

**Bill-to-tail tilt ($\varepsilon_{2,\pm 1}$).**
These components measure the asymmetry between the bill and tail
ends.  They couple $m$ and $m \pm 1$ substates, mixing the
axisymmetric pattern and tilting the nodal lines of the
eigenfunctions away from the body axis.

**Equatorial ellipticity ($\varepsilon_{2,\pm 2}$).**
The $m = \pm 2$ quadrupole deformation measures the difference between
the duck's width and depth.  These components couple $m$ and $m \pm 2$
substates, breaking the $m \leftrightarrow -m$ degeneracy that persists
for axisymmetric deformations.

**Higher-order features ($\varepsilon_{4,m}$).**
The hexadecapole corrections encode finer morphological structure:
the neck constriction, bill tip geometry, and tail curl.  These
provide secondary splitting corrections and are most important for
$\ell \geq 3$ modes.

**Gravitational Zeeman analogy.**
The splitting pattern is precisely analogous to the Zeeman effect in
atomic physics.  There, a magnetic field $\mathbf{B}$ breaks rotational
symmetry and lifts the $(2\ell+1)$-fold degeneracy of hydrogen
energy levels, with the splitting proportional to $m \cdot B$.  Here,
the duck shape deformation $\varepsilon_{\ell' m'}$ plays the role of
the symmetry-breaking field, the Wigner 3j symbols replace the
Clebsch-Gordan coefficients of the magnetic coupling, and the
Hadamard formula replaces the first-order energy perturbation.  The
analogy extends to selection rules: just as only $\Delta m = 0, \pm 1$
transitions are allowed for electric dipole radiation in the Zeeman
effect, the 3j selection rule $m' = m_1 - m_2$ constrains which
deformation harmonics couple which substates.  We therefore term
this phenomenon the *gravitational Zeeman effect*.


## 3.7 Gravitational Wave Emission from Duck QNMs

Each duck QNM radiates gravitational waves at its own (shifted)
frequency, producing a richer ringdown signal than a spherical
remnant.

**Power per mode.**
For a mode of angular order $(\ell, m)$ oscillating at frequency
$\omega_{\ell m}$ with mass multipole amplitude $\delta Q_\ell^m$,
the radiated power is [8]

$$\dot{E}_{\ell,m} = \frac{G}{c^{2\ell+1}} \, \omega_{\ell m}^{2\ell+2} \, |\delta Q_\ell^m|^2 \times \frac{4\pi(\ell+1)(\ell+2)}{\ell(\ell-1)[(2\ell+1)!!]^2}.$$
<!-- Eq. (3.18) -->

For the duck, the key difference from the sphere is that different
$m$-substates now have *different* frequencies $\omega_{\ell m}$,
different multipole amplitudes $|\delta Q_\ell^m|$ (because the
eigenfunctions are rotated by the non-axisymmetric perturbation), and
different damping rates $\tau_{\ell m}^{-1}$ (since damping depends on
the mode frequency and coupling to radiation).

**Ringdown waveform.**
The gravitational wave signal from a ringing duck, observed at
luminosity distance $D_L$ and inclination $(\iota, \varphi)$, is

$$h_+(t) + i h_\times(t) = \frac{1}{D_L} \sum_{n,\ell,m} \mathcal{A}_{n\ell m} \, e^{-t/\tau_{n\ell m}} \, e^{i \omega_{n\ell m} t} \, {}_{-2}Y_{\ell m}(\iota, \varphi),$$
<!-- Eq. (3.19) -->

where $\mathcal{A}_{n\ell m}$ are the complex excitation amplitudes
(determined by the initial perturbation and duck geometry),
$\omega_{n\ell m}$ are the split frequencies from Sec. 3.5, and
${}_{-2}Y_{\ell m}$ are the spin-weighted spherical harmonics encoding
the angular radiation pattern.

For a spherical remnant, all $m$-substates within a given $\ell$ have
the same $\omega$ and $\tau$, and the sum over $m$ can be absorbed
into the angular dependence.  For the duck, the $(2\ell+1)$
substates produce $(2\ell+1)$ distinct damped sinusoids that beat
against one another, creating a characteristic *amplitude modulation*
in the ringdown waveform.  This modulation pattern encodes the
deformation coefficients $\varepsilon_{\ell' m'}$ and constitutes a
direct observable signature of the duck's non-spherical geometry.

**Observational strategy.**
QNM spectroscopy of the ringdown [9,10] can in principle resolve
the $m$-splitting if the fractional frequency difference
$\Delta\omega/\omega \sim \varepsilon$ exceeds the inverse quality
factor $1/Q = 2\omega_I/\omega_R$ of the modes.  For the duck,
$\varepsilon \sim 0.3$ while typical neutron star f-mode quality
factors are $Q \sim 5$--$20$ [1], so $\Delta\omega/\omega \gg 1/Q$:
the splitting is easily resolved.  We construct explicit template
waveforms from Eq. (3.19) for injection studies in Sec. 8.2.


---

**References for this section:**

[1] K. D. Kokkotas and B. G. Schmidt, Living Rev. Rel. **2**, 2 (1999).

[2] T. G. Cowling, MNRAS **101**, 367 (1941).

[3] J. Hadamard, *Lecons sur le Calcul des Variations* (Hermann, Paris, 1910).

[4] D. Henry and A. Douanla, J. Math. Anal. Appl. **430**, 995 (2015).

[5] T. Suzuki and S. Tsuchiya, arXiv:2309.00273 (2023).

[6] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum* (World Scientific, 1988).

[7] DLMF, NIST Digital Library of Mathematical Functions, Chap. 34.

[8] K. S. Thorne, Rev. Mod. Phys. **52**, 299 (1980).

[9] E. Berti, V. Cardoso, and A. O. Starinets, Class. Quant. Grav. **26**, 163001 (2009).

[10] M. Isi, M. Giesler, W. M. Farr, M. A. Scheel, and S. A. Teukolsky, Phys. Rev. Lett. **123**, 111102 (2019).
