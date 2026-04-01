# THE REAL PULSATIONS OF A RELATIVISTIC STELLAR MODEL

The physical pulsations of a relativistic stellar model are *real* linear combinations of
the complex normal modes discussed above and of their complex conjugates. We shall
describe these physical pulsations, first in words, and then in terms of the mathematics
of complex normal modes.

When an equilibrium configuration is perturbed, all of the perturbation energy is
concentrated initially in the motion of the fluid; there is no gravitational radiation.
However, as time passes the pulsation energy is gradually converted to gravitational
waves and radiated away. At any particular moment of time during the radiating phase,
the external gravitational waves have an amplitude which increases with radius until the
wave front is reached and which vanishes (no gravitational radiation) beyond the wave
front. The wave front moves forward with the speed of light carrying information about
the existence of the perturbation with it.

The mathematical description which accompanies these words is conveniently divided
into two regions: The region far behind the wave front both in distance and time, where
a very simple mathematical description suffices; and the region near the wave front,
where a more complicated description is necessary.

## a) *The Region Far behind the Wave Front*

Far behind the wave front the pulsation and the gravitational waves consist of a
linear combination of discrete, complex normal modes with purely outgoing waves
("*outgoing normal modes*"):

$$
  K(r,t) = \sum_n A_n \bigl[ K_n^{(0)}(r)\, e^{i\omega_n t} + K_n^{(0)*}(r)\, e^{-i\omega_n^* t} \bigr]
  \quad \text{and similarly for } H_0, W, V \;.
\tag{21}
$$

Here $K_n^{(0)}$ is the eigenfunction for the $n$th outgoing normal mode, $\omega_n$ is the
corresponding eigenfrequency, $A_n$ is a real amplitude, and the sum is over complex
conjugate pairs. Note that the pulsations described by equation (21) are far from the
most general ones possible: They are restricted to a particular spherical harmonic
$l$, $M \neq 0$, $\tau =$ $[-1]^l$). In order to construct the most general pulsation, one must: (1) combine these
individual harmonic pulsations with those for $M \neq 0$, which are obtained from these by
rotations about the center of the star; (2) transform the resultant harmonic pulsations
for each value of $l$ to some common gauge (recall that the Regge-Wheeler gauge depends
upon $l$); (3) take linear combinations of the transformed harmonic pulsations.

### i) *Quasi-normal Modes*

Of considerable interest are the "quasi-normal modes" of pulsation which, far behind
the wave front, are constructed from one single outgoing complex normal mode

$$
  K_n(r,t) = K_n^{(0)}(r)\, e^{i\omega_n t} + K_n^{(0)*}(r)\, e^{-i\omega_n^* t}
  = \bigl[ K_n^{(0)}(r)\, e^{i\sigma_n t} + K_n^{(0)*}(r)\, e^{-i\sigma_n t} \bigr] e^{-t/\tau_n}
\tag{22}
$$

and similarly for $H_0, W, V$.

In the wave zone, but still far behind the wave front, the gravitational waves from such
a quasi-normal pulsation have the form (cf. eq. [15]):

$$\begin{aligned}
  K(r,t) &= 2 |C_n^{(0)}| \exp\bigl[-(t - r - 2M\ln r)/\tau_n\bigr]
             \cos\bigl[\sigma_n(t - r - 2M\ln r) + \delta_n\bigr] \\
  H(r,t) &= -2\sigma_n |C_n^{(0)}| r \exp\bigl[-(t - r - 2M\ln r)/\tau_n\bigr]
             \sin\bigl[\sigma_n(t - r - 2M\ln r) + \delta_n\bigr]
             \quad \text{for } r \gg M \text{ and } r \gg 1/\sigma_n \;.
\end{aligned}\tag{23}$$

Here $\delta_n$ is the phase of $C_n^{(0)}$.

Let us examine some of the properties of the quasi-normal modes (22) and (23).

### ii) *Stability of the Quasi-normal Modes*

It is clear from equations (22) and (23) that the imaginary part, $1/\tau_n$, of the complex
eigenfrequency of the outgoing complex normal mode becomes, for the physical
quasi-normal mode, the damping rate or growth rate of pulsations. If $\tau_n$ is positive, the
quasi-normal pulsations are damped by the emission of gravitational waves. But if $\tau_n$
is negative, the quasi-normal "pulsations" are unstable, and the perturbation grows
exponentially in time. Stability versus instability is thus determined by where in the
complex plane $\omega_n$ lies: The region below the real axis is unstable; the region above the
real axis is stable.

It is actually more useful to restate this criterion in terms of $\omega_n{}^2$---the quantity appearing in the eigenequations (14). Since we have adopted the convention that $\sigma_n \geq 0$
(eq. [13]), there is a one-to-one relation between $\omega_n$ and $\omega_n{}^2$, and there is no ambiguity
about stability in the complex $\omega^2$ plane: *For any complex outgoing normal mode* ($C_n^{(0)} \neq 0$),
*the corresponding quasi-normal pulsation is stable if $\omega_n{}^2$ lies on the positive real axis
or in the upper half of the complex plane; it is unstable if $\omega_n{}^2$ lies on the negative real axis
or in the lower half of the complex plane.*

### iii) *Gravitational Radiation from the Quasi-normal Modes*

Turn now from stability to the form of the gravitational waves of the quasi-normal
pulsations. From equation (23) it is clear that the gravitational radiation emitted by
the $n$th quasi-normal mode has, as seen by a man at radius $r$, the frequency

$$
  f(r) = \frac{\sigma_n}{2\pi\,(1 - 2M/r)^{1/2}}
         \quad \text{as} \quad r \to \infty \;.
\tag{24}
$$

The same gravitational redshift is present here as for $r \to \infty$. Note that the logarithmic
term in the sinusoidal part of the gravitational waves (23) is unrelated to the gravitational redshift; it is related instead to the coordinate-time delay for the propagation of radiation through the Schwarzschild geometry.

The power radiated in gravitational waves by the quasi-normal pulsations is equal
to the time rate of change of the fluid's pulsation energy. (Here we refer to energies
measured at infinity with the energy redshift effect taken into account; cf. Harrison,
Thorne, Wakano, and Wheeler [1965], pp. 19–21; also RSSD § 3.1.2.) It is easy to
write down an expression for the kinetic energy of fluid motion but very difficult to
construct an expression for the compressional and gravitational potential energy.
Therefore, we shall calculate only the kinetic energy, and we shall then argue that the
potential energy must be equal to the kinetic energy but must pulsate out of phase with
it.

### (iv) Pulsation Energy and Power Radiated

The kinetic energy of pulsation is given by (cf. RSSD, eq. [4.26'])

$$
E_{\rm kin} = \int_{\rm star}
\underbrace{\left(\rho + p\right)}_{\text{inertial mass}\atop\text{per unit volume}}
\underbrace{\left(\frac{1}{2}v^2\right)}_{\text{physical velocity}\atop\text{of fluid}}
\underbrace{e^{(\nu - \lambda)/2}}_{\text{gravitational}\atop\text{redshift}}
\underbrace{d\mathcal{V}}_{\text{proper}\atop\text{volume}}
= \int_0^R \frac{1}{4}(\rho+p)\,e^{(\lambda-\nu)/2}\bigl[r^{-2}\xi_r^2(r) + (r^{-1}\xi_\perp(r))^2\bigr]\,e^{\nu/2}\,4\pi r^2\,dr\,d\theta\,d\phi \,.
\tag{25}
$$

Upon substituting equations (7a), (11), (12), and (22) into this expression and
integrating over $\theta$, we obtain for the fluid kinetic energy of the $n$th quasi-normal mode

$$
E_{\rm kin} = 2\pi\,(2l+1)^{-1}\,\sigma_n^2\,e^{-2i\sigma_n t}
\int_0^R (\rho+p)\,e^{(\lambda-\nu)/2}
\Bigl\{
r^{-2}\bigl|W_n^{(l)}(r)\bigr|^2
+ l(l+1)\,r^{-2}\bigl|V_n^{(l)}(r)\bigr|^2
\Bigr\}
\Bigl[
1 + l(l+1)\,\bigl|V_n^{(l)}(r)\bigr|^2
\sin^2\bigl(\sigma_n t + \delta_n(r)\bigr)
+ \bigl\{l + \delta_n(r)\bigr\}
\Bigr]\,dr \,,
\tag{26}
$$

Here $\delta_n(r)$ and $\delta_n(r)$ are the phases of $W_n^{(l)}(r)$ and $V_n^{(l)}(r)$, and we have assumed that
$1/\tau_n \ll \sigma_n$.

For situations of physical interest---e.g., neutron stars formed in a supernova explosion or neutron stars undergoing relaxation oscillations---rough estimates suggest
that $1/\tau_n \ll \sigma_n$ (see Wheeler 1966; Zee and Wheeler 1967; Misner and Zapolsky 1967;
Chau 1967). When this is true, the eigenfunctions (14) couple the real and imaginary
parts of the eigenfunctions only very slightly; and, consequently, one has

$$
\delta_n(r) \ll 1, \quad \delta_n(r) \ll 1 \,.
\tag{27}
$$

and

$$
E_{\rm kin} = E_{\rm pulse}^n \sin^2(\sigma_n t) \,.
\tag{28a}
$$

Here

$$
E_{\rm pulse}^n = 2\pi\,(2l+1)^{-1}\,\sigma_n^2
\int_0^R (\rho+p)\,e^{(\lambda-\nu)/2}
\Bigl[
r^{-2}\bigl|W_n^{(l)}(r)\bigr|^2
+ l(l+1)\,r^{-2}\bigl|V_n^{(l)}(r)\bigr|^2
\Bigr]\,dr \,.
\tag{28b}
$$

When $1/\tau_n \ll \sigma_n$, the total fluid pulsation energy---kinetic plus potential---is very
nearly constant over a total pulsation period, $\Delta t = 2\pi/\sigma_n$. This is possible only if the
potential energy of pulsation is given by

$$
E_{\rm pot} \approx E_{\rm pulse}^n \cos^2(\sigma_n t) \,,
\tag{29a}
$$

so that

$$
E_{\rm kin} + E_{\rm pot} = E_{\rm pulse}^n \,.
\tag{29b}
$$

is nearly time independent. Equations (29a) and (29b) for the total pulsation energy
should be roughly correct even when $1/\tau_n$ is not small compared to $\sigma_n$.

Since the power radiated as gravitational waves is the negative of the time rate of
change of the fluid pulsion energy, the *radiated power as measured far from the star is*

$$
  P = -2\,\tau_{n}^{-1} E_{n,\text{puls}}^{-k/\tau_{n}} \;.
\tag{30}
$$

## b) The Region near the Wave Front

Thus far we have described mathematically the quasi-normal pulsations only in the
region far behind the wave front. Near the wavefront the description is more complicated because one must superimpose a large number of normal modes with various
amounts of outgoing and incoming radiation in order to obtain a sharp "turning-on" of
the perturbation at the wave front.

We shall confine our analysis near the wavefront to those quasi-normal modes with
$|f/\tau_{n}| \ll \sigma_{n}$ ($\sigma_{n}$ near the positive real axis); and we shall pattern our analysis after the
classic analysis of nuclear particle decay by Gamow (1931) and by Breit and Wert
(1935).

We build the wave front and associated quasi-normal pulsations out of a linear
combination of real standing-wave normal modes with eigenfrequencies near $\sigma_{n}$---and
hence also near $\omega_{n} = \sigma_{n} + i/\tau_{n}$. Each of these standing-wave modes is normalized to
have the same real amplitude, $B_{n}$, at the center of the star (cf. eq. [18]). In the wave
zone reality of the eigenfunctions guarantees that the ingoing and outgoing amplitudes
are complex conjugates of each other (cf. eq. [18]):

$$
  C^{(\mathrm{I})} = C^{(\mathrm{O})*} \;.
\tag{31}
$$

These wave-zone amplitudes are conveniently expressed as functions of the real frequency, $\omega$, by means of a power series expansion in the complex frequency plane about
the purely outgoing frequency $\omega_{n} = \sigma_{n} + i/\tau_{n}$:

$$
  C_{\omega}^{(\mathrm{O})} = C_{\Omega}^{(\mathrm{O})} = \int dC^{(\mathrm{O})}/d\omega\big|_{\Omega}\,(\omega - \omega_{n} - i/\tau_{n}) \;.
\tag{32}
$$

The particular combination of real, standing-wave normal modes which yields the desired wave front and subsequent quasi-normal pulsations is this:

$$
  K(r,t) = \int_{-\infty}^{+\infty}
    \frac{K_{n}(r)}{(\omega - \sigma_{n})^{2} + \tau_{n}^{-2}}
    \frac{1}{2\pi}\,e^{i\omega t}\,d\omega
  + \begin{pmatrix} \text{complex} \\ \text{conjugate} \end{pmatrix} ,
\tag{33}
$$

and similarly for $H_{0}, W, V$.

In the wave zone this becomes, for times[^1] $t > 0$ [combine eqs. (18), (32), and (33); integrate $\omega \to \infty$ since the resonant denominator in eq. (33) is huge for negative
$\omega$; close the integration contours at complex infinity; and use the method of residues to
evaluate the integrals]:

$$
  K(r,t)
  = \begin{cases}
      0 & \text{if} \quad t < r + 2M\ln r \;, \\[4pt]
      \dfrac{1}{2}\,dC^{(\mathrm{O})}/d\omega\big|_{\Omega}\,
        \exp\!\bigl[-(t - r - 2M\ln r)/\tau_{n} + i\sigma_{n}(t - r - 2M\ln r) + \delta\bigr]
      & \text{if} \quad t > r + 2M\ln r \;.
    \end{cases}
\tag{34a}
$$

$$
H_0(r,t) =
\begin{cases}
0 & \text{if} \quad t < r + 2M \ln r \,, \\[6pt]
-2\pi_r \bigl|dC^{(1)}/du\bigr|_{\omega} \exp\bigl[-\bigl(t - r - 2M\ln r\bigr)/\tau_r\bigr]
\sin\bigl[\omega_r\bigl(t - r - 2M\ln r\bigr) + \delta\bigr]
& \text{if} \quad t > r + 2M \ln r \,.
\end{cases}
\tag{34b}
$$

Here $\delta + \pi/2$ is the phase of $[dC^{(1)}/du]_{\omega}$.

Note that in the wave zone ($r \gg 2M$, beyond the wave front) there is no gravitational radiation; information that the equilibrium configuration was perturbed has not yet had time to reach radius $r$. However, immediately behind the wave front ($r > r + 2M$ in the original), an observer will measure gravitational waves corresponding to the $\omega$ quasi-normal mode of pulsation of harmonic $(l,\, M = 0,\, \tau \equiv [-\frac{1}{2}])$ (cf. eqs. [23] with eqs. [34]).

[^1]: For times $t < 0$ (eq. (30))---incoming gravitational waves with a trailing wave front. We choose to assume $t = 0$. The perturbation is first turned on inside the star at time $t = 0$, thereby discarding the incoming solution for $t < 0$.
