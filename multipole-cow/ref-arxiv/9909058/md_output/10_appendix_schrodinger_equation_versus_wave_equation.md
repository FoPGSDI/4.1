# Appendix: Schrödinger Equation Versus Wave Equation

The purpose of this appendix is a brief comparison of properties of
the Schrö\-dinger equation

$$
  i\dot\Psi = \left(-{d^2\over {dx^2}}+ V(x)\right)\Psi.
$$

and the wave equation

$$
  -\ddot\Phi = \left(-{d^2\over {dx^2}}+ V(x)\right)\Phi
$$

for the same potential $V(x)$.

The ansatz of “stationary states”

$$
  \Psi=e^{-iEt}\psi (x)
$$

leads to the time independent Schrödinger equation

$$
  -\psi”+ V\ \psi=E\psi.

$$

The ansatz (corresponding to Laplace transformation)

$$
  \Phi=e^{st}\phi (x)
$$

gives for the wave equation

$$
  -\phi”+ V\ \phi =-s^2\phi .

$$

The equations ([ref:a_4]) and ([ref:a_6]) are the same if we set $
E=-s^2$. Hence time independent scattering theory -- the theory of the
operator $-d^2 / dx^2 + V(x) $ -- does not know whether we deal with
the wave equation or the Schrödinger equation.

Consider first a positive  potential of compact support. From quantum
mechanics we know that there are no bound states. This is easy to
understand: outside the support of the potential the solutions are
$\exp (\pm x \sqrt{-E})$, hence only negative values of $E$ are
possible eigenvalues. However, if we multiply ([ref:a_4]) by $\psi$
and integrate over all $x$ we obtain a contradiction because of the
positivity of $V$. Hence the operator $-d^2 / dx^2 + V(x)$ has no
eigenfunctions and there are only scattering states, the continuous
spectrum. The operator is selfadjoint on the Hilbert space of square
integrable functions on the real line. Its resolvent, $\mathbf{R}_E$ is
defined for all complex $E$ outside the continuous spectrum, which
consists of the non negative real numbers. $\mathbf{R}_E$ is an integral
operator whose kernel is the Green function constructed in
section [ref:section_1]. There $G(s,x,x’)$ was considered as a
function of $s$. As a function of $E$, $G$ is analytic on the whole
complex plane with the exception of  real $E\leq 0$, the continuous
spectrum. In terms of $s$, the Green function -- and the resolvent --
is defined on the “physical half space” $Re(s)>0$.

Resonances of Schrödinger operators can be defined as the
poles of the analytic extension of the Green function or the
resolvent. There is a huge amount of literature on the subject, in
particular mathematical papers. A convenient starting point may be the
proceedings of a conference on resonances in
1984 [albererio]. Existence of resonances, asymptotic
distribution and also their interpretation is treated. There is in
particular the recent development of “Geometric scattering Theory",
where the use of “pseudo differential operators” gives strong and
interesting results [Melrose]. We describe one such result on the
asymptotic distribution of quasi-normal mode frequencies of the
Schwarzschild spacetime in section [ref:section_2]. It is amusing to
note that in the field of quantum mechanics the same difficulties in
defining the notion of a resonance occurred as in relativity in the
context of “normal modes of black holes”! [Simon] is a good
reference explaining this point.

Historically, quasi-normal modes appeared the first time in Gamow’s
treatment of the $\alpha$-decay. The model he studies is a potential
with two positive square potentials. Radioactive decay is exponential
in experiments. However, even the decay in time of solutions of the
free Schrödinger equation is a power law decay ($1 / t$ for
1-dimensional systems), similarly for potentials with compact
support. So we face the difficulty to characterize that part of the
time evolution, in which exponential decay is a good approximation,
knowing that the final decay is polynomial! In [Skibsted] such
estimates are derived.

Let us finally consider general potentials of
compact support. If the potential well is deep enough a finite number
of negative eigenvalues $E_n$ may exist describing  bound states of
the quantum system. What are the properties of the corresponding
solutions of the wave equation? Let $\psi_n(x)$ be an eigenfunction of
$- d^2 / dx^2 + V(x)$ with eigenvalue $E_n<0$. Then ($ E=-s^2$)

$$
  \Phi = e^{\sqrt{-E_n}\ t}\psi_n(x)
$$

is a solution of the wave equation. For large positive $x$ we have

$$
  \psi_n= e^{-\sqrt{-E_n}\ x}
$$

as the only square integrable solution of the Schrödinger equation
with vanishing potential. Thus the solution of the wave equation for
large $x$ is

$$
  \Phi = e^{\sqrt{-E_n}\ t}e^{-\sqrt{-E_n}\ x} = e^{\sqrt{-E_n}\
  (t-x)}.
$$

This solution grows exponentially in time and falls off exponentially
in space for $x \to\infty$. This apparently strange behaviour is
possible because the energy density of the conserved energy of the
wave equation

$$
  \epsilon={1\over 2}\left((\dot\Phi )^2+(\Phi’)^2 +V\Phi^2\right)
$$

is not positive definite, if the potential is somewhere negative. In
this situation the solution can grow in time and nevertheless have
conserved finite energy.

In the complex $s$ plane the eigenvalues appear as poles of the Green
function with values $s_n=\sqrt{-E_n}$. To the right of the largest
eigenvalue we have analyticity of the Green function.

Let us close this section by remarking that functional analysis
techniques can also be used to develop existence theory: In the case
discussed above, $- d^2 / dx^2 + V(x) $ is selfadjoint on the space of
square-integrable function. “Functional calculus" can be used to show
the existence and uniqueness of solutions of the time dependent
Schrödinger and wave equation, given appropriate initial data. From
this point of view the time independent Green function is the primary
object, which is unique because there is a unique selfadjoint operator
$- d^2 / dx^2 + V(x) $ on the real line.