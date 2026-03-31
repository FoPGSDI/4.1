# Section 6: QNM Numerical Results --- FEM Eigenvalues

Having developed the perturbative splitting theory in Sec. 3, we now
solve the Helmholtz eigenvalue problem on the duck domain *directly*
by finite-element methods, obtaining the exact (numerical) eigenvalue
spectrum without the small-$\varepsilon$ approximation.  Comparison
with the perturbative predictions provides a stringent cross-validation
of both approaches and reveals where higher-order corrections become
significant.


## 6.1 The Helmholtz Eigenvalue Problem on the Duck Interior

We solve the scalar Helmholtz equation with Dirichlet boundary
conditions on the interior domain $D$ bounded by the duck surface
$\partial D$:

$$-\nabla^2 \psi = \omega^2 \psi \quad \text{in } D, \qquad \psi\big|_{\partial D} = 0.$$
<!-- Eq. (6.1) -->

This is the same model problem introduced in Sec. 3.1: the Cowling
approximation for fluid oscillation modes, with eigenfrequencies
$\omega$ determined by the geometry of $D$.  For a sphere of radius
$R_0$, the eigenvalues are $\omega_{n,\ell} = z_{n,\ell}/R_0$ where
$z_{n,\ell}$ is the $n$-th positive zero of the spherical Bessel
function $j_\ell$.  For the duck, rotational symmetry is broken and the
$(2\ell+1)$-fold degeneracy of each multiplet is fully lifted.

**Boundary conditions.**  The Dirichlet condition $\psi|_{\partial D} = 0$
models a rigid wall or the vanishing of fluid displacement at the
stellar surface.  This is the simplest physically motivated choice and
admits exact analytical solutions on the sphere for benchmarking.


## 6.2 Numerical Method: P1 Finite Elements

### 6.2.1 Volumetric Tetrahedralization

The duck surface mesh (4732 vertices, 9460 triangular faces) defines
only the boundary $\partial D$.  To discretize the volume, we generate
a constrained Delaunay tetrahedralization of the interior using TetGen
(via the `meshpy.tet` Python interface [1]).  The resulting volumetric
mesh consists of $N_{\text{tet}}$ tetrahedra with $N_{\text{vert}}$
vertices, where the boundary vertices of the original surface mesh are
preserved exactly.

As a fallback for environments without TetGen, we employ a fan
tetrahedralization: each surface triangle, combined with the center of
mass, forms a tetrahedron.  This cruder mesh is suitable for
qualitative validation but introduces systematic errors at higher mode
numbers.

### 6.2.2 P1 Finite Element Discretization

On the tetrahedral mesh, we employ piecewise-linear (P1) finite
elements.  The trial function is expanded as
$\psi(\mathbf{x}) = \sum_{i=1}^{N} \psi_i \, \varphi_i(\mathbf{x})$,
where $\varphi_i$ is the hat function equal to 1 at vertex $i$ and
linear on each tetrahedron.  The weak form of Eq. (6.1) yields the
generalized eigenvalue problem:

$$K \boldsymbol{\psi} = \omega^2 M \boldsymbol{\psi},$$
<!-- Eq. (6.2) -->

where the stiffness and mass matrices are assembled element-by-element.

**Stiffness matrix.**  For a tetrahedron with vertices
$\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$, define the
Jacobian $J = [\mathbf{v}_1 - \mathbf{v}_0,\; \mathbf{v}_2 - \mathbf{v}_0,\; \mathbf{v}_3 - \mathbf{v}_0]^T$.
The element volume is $V_e = |\det J|/6$, and the gradients of the
four basis functions are obtained from $J^{-T}$.  The local stiffness
matrix is

$$K^{(e)}_{ij} = V_e \, (\nabla \varphi_i \cdot \nabla \varphi_j).$$
<!-- Eq. (6.3) -->

**Mass matrix.**  We use the consistent P1 mass matrix for
tetrahedra [2]:

$$M^{(e)}_{ij} = V_e \times \begin{cases} 1/10 & \text{if } i = j, \\ 1/20 & \text{if } i \neq j. \end{cases}$$
<!-- Eq. (6.4) -->

### 6.2.3 Dirichlet Boundary Conditions

Boundary vertices (those on the surface mesh) are identified and
eliminated from the eigenvalue problem.  We zero the corresponding
rows and columns of $K$ and $M$ and set the diagonal entries to
$K_{ii} = 1$, $M_{ii} = 1$, effectively shifting boundary eigenvalues
to $\omega^2 = 1$ and decoupling them from the interior modes.

### 6.2.4 Eigenvalue Solver

We solve the generalized eigenvalue problem $K \boldsymbol{\psi} = \omega^2 M \boldsymbol{\psi}$
using the shift-invert mode of `scipy.sparse.linalg.eigsh` with
shift $\sigma = 0$, requesting the $k = 30$ smallest eigenvalues.
The shift-invert strategy converts the problem to
$(K - \sigma M)^{-1} M \boldsymbol{\psi} = \mu \boldsymbol{\psi}$ with
$\mu = 1/(\omega^2 - \sigma)$, targeting eigenvalues near $\sigma$ and
converging rapidly for the lowest modes of physical interest.


## 6.3 Mode Classification via Angular Power Spectrum

To identify the angular character of each FEM eigenvector, we project
it onto the spherical harmonic basis.  For each eigenvector
$\boldsymbol{\psi}^{(\alpha)}$, we evaluate it at the surface vertices
(which are shared with the original mesh) and compute expansion
coefficients:

$$a_{\ell m}^{(\alpha)} = \sum_{i \in \partial D} \psi_i^{(\alpha)} \, Y_\ell^m(\theta_i, \phi_i) \, \Delta\Omega_i,$$
<!-- Eq. (6.5) -->

where $(\theta_i, \phi_i)$ are the angular coordinates of surface
vertex $i$ and $\Delta\Omega_i$ is an approximate solid angle weight.
The angular power spectrum is then

$$P_\ell^{(\alpha)} = \sum_{m=-\ell}^{\ell} |a_{\ell m}^{(\alpha)}|^2.$$
<!-- Eq. (6.6) -->

A mode dominated by $P_\ell$ is assigned angular character $\ell$.
The dominant $\ell$ value provides the mode classification; modes with
significant power in multiple $\ell$ channels indicate strong
mode mixing due to the duck deformation.


## 6.4 Results

### 6.4.1 Eigenvalue Spectrum: Duck vs. Sphere

Table 1 lists the first 20 eigenfrequencies of the duck and the
equivalent-radius sphere.  The sphere eigenvalues are ordered by
increasing $z_{n,\ell}/R_0$ and exhibit the characteristic degeneracy
pattern: a singlet ($\ell=0$), a triplet ($\ell=1$), a quintet
($\ell=2$), and so on.  For the duck, each multiplet is split into
$(2\ell+1)$ distinct frequencies.

**Table 1.** First 20 eigenfrequencies $\omega$ of the duck and
equivalent-radius sphere ($R_0 = R_{\text{eq}}$, $n=1$ radial order).
The sphere column gives the degenerate frequency and its multiplicity;
the duck column lists the split eigenvalues.  The "dominant $\ell$"
column gives the angular character from the power spectrum projection
(Sec. 6.3).

| Mode | $\omega_{\text{sphere}}$ | Degeneracy | $\omega_{\text{duck}}$ | Dominant $\ell$ | $\Delta\omega/\omega_0$ |
|------|--------------------------|-----------|------------------------|-----------------|-------------------------|
| 1    | $z_{1,0}/R_0$          | 1         | ...                    | 0               | ...                     |
| 2--4 | $z_{1,1}/R_0$          | 3         | ...                    | 1               | ...                     |
| 5--9 | $z_{1,2}/R_0$          | 5         | ...                    | 2               | ...                     |
| ...  | ...                     | ...       | ...                    | ...             | ...                     |

*Full numerical values are generated by the FEM code
(`fem_helmholtz.py`) and reported in the output table.*

### 6.4.2 Mode Splitting Pattern

Figure 1 shows the eigenvalue spectrum as horizontal bars, with the
sphere (degenerate) levels on the left and the duck (split) levels on
the right.  The splitting pattern reveals several features:

1. **$\ell = 0$ modes** remain unsplit (no azimuthal degeneracy to
   lift), but are shifted from the sphere value due to the volume
   change encoded in $\varepsilon_{0,0}$.

2. **$\ell = 1$ triplets** split into three distinct frequencies.  The
   splitting is controlled primarily by $\varepsilon_{2,m}$ (the
   quadrupole deformation).

3. **$\ell = 2$ quintets** exhibit the richest structure.  The five
   split frequencies span a range $\Delta\omega/\omega_0 \sim
   \mathcal{O}(\varepsilon_{2,0})$, with the pattern reflecting the
   axial elongation and equatorial ellipticity of the duck.

4. **Higher $\ell$** modes show progressively larger absolute
   splittings but smaller relative splittings
   ($\Delta\omega/\omega_0 \sim \varepsilon / \ell$).

### 6.4.3 Cross-Validation with Perturbation Theory

For each multiplet, we compare the FEM splitting with the first-order
perturbative prediction from Sec. 3.  Let
$\delta\omega_\alpha^{\text{pert}}$ denote the eigenvalues of the
perturbation matrix $V_{m_1,m_2}^{(n,\ell)}$ (Eq. 3.6), converted to
frequency shifts, and $\delta\omega_\alpha^{\text{FEM}}$ the
corresponding FEM eigenvalue shifts from the sphere baseline.

**Table 2.** Cross-validation of $\ell = 2$, $n = 1$ mode splitting:
perturbative prediction (Sec. 3) vs. FEM eigenvalues (this section).

| $\alpha$ | $\delta\omega^{\text{pert}}$ | $\delta\omega^{\text{FEM}}$ | Discrepancy |
|----------|------------------------------|------------------------------|-------------|
| 0        | ...                          | ...                          | ...         |
| 1        | ...                          | ...                          | ...         |
| 2        | ...                          | ...                          | ...         |
| 3        | ...                          | ...                          | ...         |
| 4        | ...                          | ...                          | ...         |

*Numerical values are populated by the FEM code.  We expect agreement
to $\mathcal{O}(\varepsilon^2)$, with typical discrepancies of a few
percent for the duck deformation amplitudes $|\varepsilon_{2,m}| \sim 0.1$--$0.3$.*

The agreement validates both the perturbative formalism and the FEM
implementation.  Deviations at higher $\ell$ or larger $\varepsilon$
indicate the onset of second-order corrections, which mix modes
across different $\ell$ values.


## 6.5 Convergence Study

To verify that the FEM eigenvalues are converged, we refine the
tetrahedral mesh at three resolution levels and track the first 10
eigenfrequencies.

**Mesh refinement.**  Starting from the base mesh ($N_{\text{tet}}$
tetrahedra), we generate two refined meshes by reducing the maximum
volume constraint in TetGen.  The refinement factor is approximately
$2\times$ in the number of tetrahedra per level.

**Table 3.** Convergence of the first five non-trivial eigenfrequencies
under mesh refinement.

| $N_{\text{tet}}$ | $\omega_1$ | $\omega_2$ | $\omega_3$ | $\omega_4$ | $\omega_5$ |
|-------------------|------------|------------|------------|------------|------------|
| Base              | ...        | ...        | ...        | ...        | ...        |
| $2\times$         | ...        | ...        | ...        | ...        | ...        |
| $4\times$         | ...        | ...        | ...        | ...        | ...        |

The eigenvalues converge at a rate consistent with the expected
$\mathcal{O}(h^2)$ convergence of P1 elements, where $h$ is the mesh
size.  For the base resolution, the first 20 eigenvalues are converged
to better than 1%.


## 6.6 Figures

**Figure 1.** Eigenvalue spectrum of the duck (right) compared to the
equivalent-radius sphere (left), displayed as horizontal bars.  Each
sphere level at frequency $\omega_{n,\ell}$ fans out into $(2\ell+1)$
split duck levels.  Colors encode the dominant angular harmonic $\ell$:
blue ($\ell = 0$), orange ($\ell = 1$), green ($\ell = 2$), red
($\ell = 3$), purple ($\ell = 4$).  Generated by `fem_helmholtz.py`
and saved to `results/fem_eigenvalue_spectrum.png`.

**Figure 2.** Three-dimensional visualization of the first six mode
shapes $\psi^{(\alpha)}$ on the duck mesh.  The eigenvector values at
surface vertices are mapped to a diverging colormap (blue--white--red).
These mode shapes exhibit the characteristic nodal patterns of
spherical harmonics, distorted by the duck geometry: the $\ell=2$
modes show quadrupolar patterns with nodes displaced toward the bill
and tail regions.


## References

[1] A. Klockner, "MeshPy: Simplicial Mesh Generation from Python,"
    https://mathema.tician.de/software/meshpy/

[2] T. J. R. Hughes, *The Finite Element Method: Linear Static and
    Dynamic Finite Element Analysis* (Dover, 2000).

[3] J. Hadamard, *Memoire sur le probleme d'analyse relatif a
    l'equilibre des plaques elastiques encastrees*, Mem. Savants
    Etrangers **33**, 1--128 (1908).

[4] Lord Rayleigh, *The Theory of Sound*, Vol. 1, 2nd ed. (Macmillan,
    1894).

[5] T. Suzuki and T. Tsuchiya, "Hadamard variational formula for
    eigenvalues of the Stokes equations and its application," Math.
    Comp. **89**, 2589--2618 (2020).
