# Normal Modes -- Quasi-Normal Modes -- Resonances

Before discussing quasi-normal modes it is useful to remember what
normal modes are!

Compact classical linear oscillating systems such as finite strings,
membranes, or cavities filled with electromagnetic radiation have
preferred time harmonic states of motion ($\omega$ is real):

$$
  \chi_n(t,x)=e^{i\omega_n t}\chi_n(x), \quad n = 1, 2, 3 \dots,

$$

if dissipation is neglected. (We assume $\chi$ to be some complex
valued field.) There is generally an infinite collection of such
periodic solutions, and the “general solution” can be expressed as a
superposition,

$$
  \chi(t,x)=\sum_{n=1}^\infty a_n e^{i\omega_n t}\chi_n(x),

$$
of such normal modes. The simplest example is a string of length $L$
which is fixed at its ends. All such systems can be described by
systems of partial differential equations of the type ($\chi$ may be
a vector)

$$
  {\partial\chi\over\partial t} = \mathbf{A}\chi,

$$

where $\mathbf{A}$ is a linear operator acting only on the spatial
variables. Because of the finiteness of the system the time evolution
is only determined if some boundary conditions are prescribed.
The search for solutions periodic in time leads to a boundary value
problem in the spatial variables. In simple cases it is of the
Sturm-Liouville type. The treatment of such boundary value problems
for differential equations played an important role in the development
of Hilbert space techniques.

A Hilbert space is chosen such that the differential operator
becomes symmetric. Due to the boundary conditions dictated by the
physical problem, $\mathbf{A}$ becomes a self-adjoint operator on the
appropriate Hilbert space and has a pure point spectrum. The
eigenfunctions and eigenvalues determine the periodic
solutions ([ref:bs1.1]).

The definition of self-adjointness is rather subtle from a physicist’s
point of view since fairly complicated “domain issues” play an
essential role. (See [BS95] where a mathematical exposition for
physicists is given.) The wave equation modeling the finite string has
solutions of various degrees of differentiability. To describe all
“realistic situations”, clearly $C^\infty$ functions should be
sufficient. Sometimes it may, however, also be convenient to consider
more general solutions.

From the mathematical point of view the collection of all smooth
functions is not a natural setting to study the wave equation because
sequences of solutions exist which converge to non-smooth solutions.
To establish such powerful statements like ([ref:bs1.2]) one has to
study the equation on certain subsets of the Hilbert space of square
integrable functions. For “nice” equations it usually  happens that
the eigenfunctions are in fact analytic. They can then be used to
generate, for example, all smooth solutions by a pointwise converging
series ([ref:bs1.2]). The key point is that we need some mathematical
sophistication to obtain the “completeness property” of the
eigenfunctions.

This picture of “normal modes”  changes when we consider “open
systems” which can lose energy to infinity. The simplest case are
waves on an infinite string. The general solution of this problem is

$$
  \chi(t,x)=A(t-x) + B(t+x)

$$

with “arbitrary” functions $A$ and $B$. Which solutions should we
study? Since we have all solutions, this is not a serious question.
In more general cases, however, in which the general solution is not
known, we have to select a certain class of solutions which we
consider as relevant for the physical problem.

Let us consider for the following discussion, as an example, a wave
equation with a potential on the real line,

$$
  {\partial^2\over{\partial t^2}}\chi
  +\left(-{\partial^2\over{\partial x^2}} + V(x)\right)\chi = 0.

$$
Cauchy data $\chi(0,x), \partial_t\chi(0,x)$ which have two
derivatives determine a unique twice differentiable solution.
No boundary condition is needed at infinity to determine the
time evolution of the data! This can be established by fairly simple
PDE theory [John].

There exist solutions for which the support of the fields are
spatially compact, or -- the other extreme -- solutions with infinite
total energy for which the fields grow at spatial infinity in a quite
arbitrary way!

From the point of view of physics smooth solutions with spatially
compact support should be the relevant class -- who cares what happens
near infinity! Again it turns out that mathematically it is more
convenient to study all solutions of finite total energy. Then the
relevant operator is again self-adjoint, but now its spectrum is
purely “continuous”. There are no eigenfunctions which are square
integrable. Only “improper eigenfunctions” like plane waves exist.
This expresses the fact that we find a solution of the form ([ref:bs1.1])
for any real $\omega$ and by forming appropriate superpositions one
can construct solutions which are “almost eigenfunctions”. (In the
case $V(x)\equiv 0$ these are wave packets formed from plane waves.)
These solutions are the analogs of normal modes for infinite systems.

Let us now turn to the discussion of “quasi-normal modes” which are
conceptually different to normal modes. To define quasi-normal modes
let us consider the wave equation ([ref:bs1.5]) for  potentials with
$V\geq 0$ which vanish for $|x|>x_0$. Then in this case all solutions
determined by data of compact support are bounded: $|\chi(t,x)|< C$.
We can use Laplace transformation techniques to represent such
solutions. The Laplace transform $\hat{\chi}(s,x)$ ($s>0$ real) of a
solution $\chi(t,x)$ is

$$
  \hat{\chi}(s,x)=\int_0^\infty e^{-st} \chi(t,x) dt,

$$

and satisfies the ordinary differential equation

$$
  s^2\hat{\chi} -\hat{\chi}”+ V \hat{\chi} = + s \chi(0,x)+\partial_t
  \chi (0,x),

$$

where

$$
  s^2\hat{\chi} -\hat{\chi}”+ V \hat{\chi} = 0

$$

is the homogeneous equation. The boundedness of $\chi$ implies that
$\hat{\chi}$ is analytic for positive, real $s$, and has an analytic
continuation onto the complex half plane $Re(s)>0$.

Which solution $\hat{\chi}$ of this inhomogeneous equation gives the
unique solution in spacetime determined by the data? There is no
arbitrariness; only one of the Green functions for the inhomogeneous
equation is correct!

All Green functions can be constructed by the following well known
method. Choose any two linearly independent solutions of the
homogeneous equation $f_-(s,x)$ and $f_+(s,x)$, and define
$$
  G(s,x,x’) = {1 \over W(s)}
  \begin{array}{ll}
    f_-(s,x’) f_+(s,x) \quad & (x’ < x), \\
    f_-(s,x) f_+(s,x’) \quad & (x’ > x),
  \end{array}

$$

where $W(s)$ is the Wronskian of $f_-$ and $f_+$. If we denote the
inhomogeneity of ([ref:bs1.7]) by $j$, a solution of ([ref:bs1.7]) is

$$
  \hat{\chi}(s,x) = \int_{-\infty}^\infty G(s,x,x’) j(s,x’) dx’.

$$

We still have to select a unique pair of solutions $f_-, f_+$. Here
the information that the solution in spacetime is bounded can be
used. The definition of the Laplace transform implies that
$\hat{\chi}$ is bounded as a function of $x$. Because the potential
$V$ vanishes for $|x|>x_0$, the solutions of the homogeneous
equation ([ref:bs1.7b]) for $|x|>x_0$ are

$$
  f = e^{\pm sx}.

$$

The following pair of solutions

$$
  f_+=e^{-sx} \quad \text{for } x>x_0 , \quad \quad f_-=e^{+sx} \quad
  \text{for } x<-x_0,

$$

which is linearly independent for $Re(s)> 0$, gives the unique Green
function which defines a bounded solution for $j$ of compact support.
Note that for $Re(s)>0$ the solution $f_+$ is exponentially decaying
for large $x$ and $f_-$ is exponentially decaying for small $x$. For
small $x$ however, $f_+$ will be a linear combination $a(s)e^{-s x} +
b(s)e^{s x}$ which will in general grow exponentially. Similar
behavior is found for $f_-$.

Quasi-Normal mode frequencies $s_n$ can be defined as those complex
numbers for which

$$
  f_+(s_n,x)=c(s_n)f_-(s_n,x),

$$

that is the two functions become linearly dependent, the Wronskian
vanishes and the Green function is singular! The corresponding
solutions $f_+(s_n,x)$ are called quasi eigenfunctions.

Are there such numbers $s_n$? From the boundedness of the solution in
spacetime we know that the unique Green function must exist for
$Re(s)>0$. Hence $f_+,f_-$ are linearly independent for those values
of $s$. However, as solutions $f_+,f_-$ of the homogeneous
equation ([ref:bs1.7b]) they have a unique continuation to the complex
$s$ plane. In [BB91] it is shown that for positive potentials
with compact support there is always a countable number of zeros of
the Wronskian with $Re(s)<0$.

What is the mathematical and physical significance of the quasi-normal
frequencies $s_n$ and the corresponding quasi-normal functions $f_+$?
First of all we should note that because of $Re(s)<0$ the function
$f_+$ grows exponentially for small and large $x$! The corresponding
spacetime solution  $e^{s_n t} f_+(s_n, x) $ is therefore not a
physically relevant solution, unlike the normal modes.

If one studies the inverse Laplace transformation and expresses $\chi$
as a complex line integral ($a>0$),

$$
  \chi (t,x) = {1\over {2\pi i}}\int_{-\infty}^{+\infty}
  e^{(a+is)t}\hat{\chi}(a+is,x)ds,

$$

one can deform the path of the complex integration and show that the
late time behavior of solutions can be approximated in finite parts of
the space by a finite sum of the form

$$
  \chi(t,x)\sim\sum_{n=1}^N\ a_n e^{(\alpha_n+i\beta_n)t}
  f_+(s_n,x).

$$

Here we assume that $Re(s_{n+1})<Re(s_n)<0$, $s_n=\alpha_n+i\beta_n$.
The approximation $\sim$ means that if we choose $x_0$, $x_1$,
$\epsilon$ and $t_0$ then there exists a constant
$C(t_0,x_0,x_1,\epsilon)$ such that

$$
  \left \vert \chi(t,x)-\sum_{n=1}^N a_n e^{(\alpha_n+i\beta_n)t}
    f_+(s_n,x) \right \vert \le Ce^{(-|\alpha_{N+1}|+\epsilon )t}

$$

holds for $t>t_0$, $x_0<x<x_1$, $\epsilon>0$ with
$C(t_0,x_0,x_1,\epsilon)$ independent of $t$. The constants $a_n$
depend only on the data [BB91]! This implies in particular that
all solutions defined by data of compact support decay exponentially
in time on spatially bounded regions. The generic leading order decay
is determined by the quasi-normal mode frequency with the largest real
part $s_1$, i.e.\ slowest damping. On finite intervals and for late
times the solution is approximated by a finite sum of quasi
eigenfunctions ([ref:bs1.14]).

It is presently unclear whether one can strengthen ([ref:bs1.15]) to a
statement like ([ref:bs1.2]), a pointwise expansion of the late time
solution in terms of quasi-normal modes. For one particular potential
(Pöschl-Teller) this has been shown by Beyer [Horst].

Let us now consider the case where the potential is positive for all
$x$, but decays near infinity as happens for example for the wave
equation on the static Schwarzschild spacetime. Data of compact
support determine again solutions which are bounded [KW87].
Hence we can proceed as before. The first new point  concerns the
definitions of $f_\pm$. It can be shown that the homogeneous
equation ([ref:bs1.7b]) has for each real positive $s$ a unique
solution $f_+(s,x)$ such that $\lim_{x\to\infty} (e^{sx}f_+(s,x))=1$
holds and correspondingly for $f_-$. These functions are uniquely
determined, define the correct Green function and have analytic
continuations onto the complex half plane $Re(s)>0$.

It is however quite complicated to get a good representation of these
functions. If the point at infinity is not a regular singular point, we
do not even get converging series expansions for $f_\pm$. (This is
particularly serious for values of $s$ with negative real part because
we expect exponential growth in $x$).

The next new feature is that the analyticity properties of $f_\pm$ in
the complex $s$ plane depend on the decay of the potential. To obtain
information about analytic continuation, even use of analyticity
properties of the potential in $x$ is made! Branch cuts may occur.
Nevertheless in a lot of cases an infinite number of quasi-normal mode
frequencies exists.

The fact that the potential never vanishes may, however, destroy the
exponential decay in time of the solutions and therefore the essential
properties of the quasi-normal modes. This probably happens if the
potential decays slower than exponentially. There is, however, the
following way out: Suppose you want to study a solution determined by
data of compact support from $t=0$ to some large finite time $t=T$.
Up to this time the solution is -- because of domain of dependence
properties -- completely independent of the potential for sufficiently
large $x$. Hence we may see an exponential decay of the
form ([ref:bs1.14]) in a time range $t_1<t<T$. This is the behavior
seen in numerical calculations. The situation is similar in the case
of $\alpha$-decay in quantum mechanics. A comparison of quasi-normal
modes of wave equations and resonances in quantum theory can be found
in the appendix, see section [ref:appendix_1].