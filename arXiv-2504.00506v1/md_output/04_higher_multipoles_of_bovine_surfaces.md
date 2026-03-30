# Higher multipoles of bovine surfaces

In the previous section, we used the multipole expansion of the mass distribution and gravitational potential to study bovine phenomenology induced by the nonspherical geometry. This is a paradigmatic example of how a particular problem might be solved beyond the spherical cow approximation, with contributions ordered by their asphericity. However, it gives us little intuition for the shape of the cow itself at successive terms in the series. The monopole term clearly corresponds to the spherical cow due to its spherical symmetry, but what is the shape of the cow that corresponds to inclusion of the dipole or quadrupole? In this section, we answer this question by defining a multipole expansion of the cow’s *geometry* based on the multipole expansion of functions on the 2-sphere.

We begin by defining the boundary surface of the cow, $\partial\mathcal C$, in a form that is amenable to a direct multipole expansion. To that end, our goal is to represent $\partial\mathcal C$ in terms of a function on $S^2$ that can be expanded in spherical harmonics. Clearly, a typical cow is topologically equivalent to a sphere,[^1] so there must exist many homeomorphisms $\mathcal F\colon S^2\to\partial\mathcal C$. Such a map $\mathcal F$ can be written in terms of a set of maps $f_i\colon S^2\to\mathbb R$ (e.g., the vector components of points in the image), and the $f_i$ admit a representation in spherical harmonics as
$$
    \boldsymbol{\mathrm{f}}(\Omega) = \sum_{\ell, m} \boldsymbol{\mathrm{f}}_\ell^m(\Omega),
    \qquad
    \boldsymbol{\mathrm{f}}_\ell^m(\Omega)\equiv
        \langle \boldsymbol{\mathrm{f}}|Y_\ell^m\rangle_{S^2} Y_\ell^m(\Omega)
    ,
$$
where $\Omega$ denotes the angular coordinates on the sphere, and $\langle\cdot|\cdot\rangle_{S^2}$ denotes the $L^2$ inner product on $S^2$,
$$
    \langle u | v\rangle = \int_{S^2}\mathrm{d}\Omega\,uv^*
    .
$$
Here, we are simply using the fact that the spherical harmonics provide a complete orthonormal system for real-valued functions on the sphere. The monopole term, $\boldsymbol{\mathrm{f}}_0^0(\Omega) = \langle \boldsymbol{\mathrm{f}} |Y_0^0\rangle Y_0^0$, is spherically symmetric in the $\boldsymbol{\mathrm{f}}$ coordinate system, and can thus be used to define the spherical cow. Higher-order terms give corrections to the cow boundary $\partial\mathcal C$.

While the decomposition of $\mathcal F$ into real-valued component maps $\boldsymbol{\mathrm{f}}$ is in principle arbitrary, there is a simple and natural choice. In order to respect the spherical symmetry of the monopole term, let us define $\boldsymbol{\mathrm{f}}(\theta, \phi) = (f_r, f_{\Delta\theta}, f_{\Delta\phi})$, where $f_r(\Omega)$ gives the radial component of $\mathcal F(\Omega)$, and $f_{\Delta\theta}$ and $f_{\Delta\phi}$ give the difference in the polar and azimuthal angles produced by application of $\mathcal F$. Specifically, let us define
$$
    \bar{f}_{\Delta\theta}(\theta, \phi) \equiv
        [\mathcal F(\theta, \phi)]_\theta - \theta,
$$
where $[\boldsymbol{\mathrm{x}}]_\theta$ denotes the $\theta$ coordinate of the point $\boldsymbol{\mathrm{x}}$. Then let $f_{\Delta\theta}(\theta,\phi) \equiv \bar{f}_{\Delta\theta}(\theta,\phi) - (4\pi)^{-1/2}\bar{f}_{\Delta\theta}^{(0)}$, where $\bar{f}_{\Delta\theta}^{(0)}$ denotes the monopole coefficient of $\bar{f}_{\Delta\theta}$. This definition of $f_{\Delta\theta}$ effects a rotation of the coordinates such that the monopole of $f_{\Delta\theta}$ vanishes. An analogous definition can be made for $f_{\Delta\phi}$.

This choice of $\boldsymbol{\mathrm{f}}$ guarantees that the image of the monopole term, with constant $(r, \Delta\theta, \Delta\phi)$, is truly a sphere, and not merely a point, or a surface homeomorphic to $S^2$, and involves no spurious rotation of the coordinate system. In particular, this means that a perfectly spherical cow is always absorbed entirely by the monopole term of the radial coordinate, and has vanishing higher multipole contributions, as desired.

All that remains is to implement this prescription is to identify the map $\mathcal F$. This is an example of a “surface matching” problem, which has been well studied in computational geometry. In the remainder of this section, I will describe two surface matching algorithms that give rise to multipole expansions of the cow surface.

[^1]: Here I ignore the alimentary canal, since it is well known that including this only results in pure bovine waste.

## Distance gradient flow method

First, let us solve the surface matching problem in a geometrically-intuitive way. The method follows from the realization that if the cow were simply rounder, every point on the boundary would have a unique angular coordinate, meaning that it would clearly be trivial to project onto the sphere. For such a round cow, $\partial\mathcal C$ could by defined by its radial coordinate via a single function $r_{\partial\mathcal C}\colon S^2\to\mathbb R$, which would then admit a multipole expansion. While this is not possible for a general cow, any cow can be first fattened up into a rounder shape and then projected to the sphere, which defines the map $\mathcal F^{-1}$.

More precisely, suppose there exists an automorphism $\varphi\colon\mathbb R^3\to\mathbb R^3$ such that $\varphi(\mathcal C)$, the image of the cow, is a star-shaped domain. Then it is possible to choose the origin such that for every point $\boldsymbol{\mathrm{x}}\in\varphi(\partial\mathcal C)$, the ray from the origin through $\boldsymbol{\mathrm{x}}$ has no other intersections with $\partial\mathcal C$, meaning that $\boldsymbol{\mathrm{x}}$ has unique angular coordinates. Then the map $p\colon\varphi(\partial\mathcal C)\to S^2$ given by $p(\boldsymbol{\mathrm{x}}) = \boldsymbol{\mathrm{x}}/\|\boldsymbol{\mathrm{x}}\|$ is smooth and bijective, so we can define $\mathcal F \equiv (p\circ\varphi)^{-1}$.

Such an inflated cow can in fact be produced computationally by flowing points in the direction of increasing distance from the cow. To make this precise, we define the signed distance function $d_{\mathcal C}\colon\mathbb R^3 \to \mathbb R$ by
$$
    d_{\mathcal C}(\boldsymbol{\mathrm{x}}) =
        \left(\min_{\boldsymbol{\mathrm{y}}\in\partial\mathcal C}\|\boldsymbol{\mathrm{x}} - \boldsymbol{\mathrm{y}} \|\right)
        \times
        \begin{cases}
            -1 & \boldsymbol{\mathrm{x}}\in\mathcal C \\
            1 & \boldsymbol{\mathrm{x}}\notin\mathcal C.
        \end{cases}
$$
Now, observe that $\lim_{\|\boldsymbol{\mathrm{x}}\|\to\infty}d_{\mathcal C}(\boldsymbol{\mathrm{x}})/\|\boldsymbol{\mathrm{x}}\| = 1$, so at large distances, the level sets of $d_{\mathcal C}$ approach a dilated $S^2$, and thus they become trivially star-shaped. This behavior is demonstrated in [ref:fig:puffy-cow]. Now consider the gradient flow of $d_{\mathcal C}$, i.e., the function $\varphi\colon{\mathbb R^3\times\mathbb R}\to\mathbb R^3$ satisfying $\partial_t\varphi(\boldsymbol{\mathrm{x}}, t) = \nabla d_{\mathcal C}(\varphi(\boldsymbol{\mathrm{x}}, t))$ with $\varphi(\boldsymbol{\mathrm{x}}, 0) = \boldsymbol{\mathrm{x}}$. This function translates the point $\boldsymbol{\mathrm{x}}$ along integral curves associated with the gradient of $d_{\mathcal C}$, orthogonal to the level sets. At large $t$, the function $\varphi(\cdot, t)$ maps $\mathcal C$ to a star-shaped domain, and the inverse is given by $\varphi(\cdot, -t)$.

Using this map to define $\mathcal F$, we can then expand the components $(f_r, f_{\Delta\theta}, f_{\Delta\phi})$ in spherical harmonics. The first few multipoles are enumerated in [ref:tab:distance-gradient-flow-components]. The magnitude of the monopole term relative to higher-order terms gives us the first quantitative assessment of the reliability of the spherical cow approximation. The resulting geometry $\partial\mathcal C$ is shown at progressively higher multipole order in [ref:fig:gradient-flow-multipoles].

While the distance gradient flow method has an appealingly simple geometric interpretation, it is far from the simplest approach to implement. The flow tends to squeeze points along the “sutures” of the inflated cow (see [ref:fig:puffy-cow]), meaning that the integral curves must be computed with very high precision for many points. All that is really needed is for the cow to become *rounder.* While the distance flow is an extremely simple way to accomplish that, there are more sophisticated techniques to more smoothly round out a cow, although they suffer from the clear deficiency that they do not produce a cow that looks as though it has been inflated with a bicycle pump. For example, one could use a smoothing operator to asymptotically map $\partial\mathcal C$ towards a surface with constant curvature, i.e., the sphere. Ultimately, however, a more established approach is provided by the method of harmonic maps.

## Harmonic map method

Another way to define a map from the cow to the sphere is to identify a set of harmonic coordinates on $\partial\mathcal C$. A major advantage of this strategy is that there exist powerful theorems that guarantee the features of the resulting map. In the previous section, we relied on ample handwaving to justify the smoothness and bijectivity of the map generated by the distance gradient flow, but in the harmonic map method, these properties can be easily proven.

A harmonic coordinate system on a manifold $M$ is a set of coordinate maps $x_i$ on $M$ such that each $x_i$ is harmonic, meaning that $\Delta x_i = 0$. Here $\Delta$ is the Laplace-Beltrami operator, a generalization of the Laplacian to Riemannian manifolds: just as the ordinary Laplacian $\nabla^2$ can be written as the divergence of the gradient of a scalar function on $\mathbb R^n$, the Laplace-Beltrami operator $\Delta$ is the Riemannian divergence of the Riemannian gradient of a scalar function on $M$. Explicitly, given a coordinate chart with a metric $g_{ij}$, the the Laplace-Beltrami operator takes the form
$$
    \Delta f = \frac{1}{\sqrt{|g|}}\partial_i\left[
        \sqrt{|g|}g^{ij}f
    \right]
    .
$$

It will turn out that if the map $\mathcal F^{-1}\colon \partial\mathcal C\to S^2$ is harmonic, then for sufficiently fine discretizations of $\partial\mathcal C$ and $S^2$, $\mathcal F$ is bijective without self-intersections of the mapped triangles. I refer the reader to Ref. [floater:2005] for the details. In particular, let us assume that $\partial\mathcal C$ is a genus-0 surface (again, omitting the alimentary canal in order to avoid BS), and that the discretizations of the surfaces take the form of triangulations. For genus-0 surfaces, harmonic maps are equivalent to conformal maps [Gu:2002], and there exists a substantial literature on computing discrete conformal maps between meshes, notably Refs. [Kharevych:2006,Springborn:2008].

In practice, as discussed in the foregoing references, we can construct a harmonic map to the sphere by first deleting one vertex $v_0$ and the connecting faces, leaving a mesh with the topology of the disk. Efficient routines for computing harmonic maps between meshes with disk topology are widely available, and indeed, one such algorithm is implemented in ‘libigl.harmonic‘. Once mapped harmonically to the circular disk, the mesh can be further mapped to the sphere via stereographic projection, which is also conformal. Finally, the deleted vertex $v_0$ is mapped to the pole of the stereographic projection, and the corresponding faces are restored. This fully specifies the homeomorphism $\mathcal F^{-1}\colon\partial\mathcal C\to\mathcal S^2$, and with it the inverse $\mathcal F$.

[ref:tab:harmonic-map-components] shows multipole coefficients for the map $\mathcal F$ constructed by this route. I selected the vertex for deletion in the middle of the cow’s back, where the geometry is relatively smooth.[^1] The basic structure of the multipole coefficients is very similar to those obtained by the distance gradient flow method ([ref:tab:distance-gradient-flow-components]).

| \ell | m | f_r | f_{\Delta\theta} | f_{\Delta\phi} |
|---|---|---|---|---|
| 0 | 0 | 1.0750 | 0 | 0 |
| 1 | 0 | 0.0010 | 0.7340 | 0.0016 |
| 1 | 1 | 0.0567 + 0.0145i | -0.0008 + 0.0045i | 0.0236 - 0.0683i |
| 2 | 0 | -0.2681 | 0.0051 | 0.2863 |
| 2 | 1 | -0.0007 - 0.0006i | -0.0085 + 0.0474i | 0.0006 - 0.0045i |
| 2 | 2 | 0.1178 - 0.0438i | 0.0010 - 0.0006i | -0.1602 - 0.0209i |
| 3 | 0 | 0.0008 | -0.1178 | -0.0017 |
| 3 | 1 | 0.0182 - 0.0307i | -0.0020 + 0.0028i | 0.0542 - 0.0645i |
| 3 | 2 | -0.0003 + 0.0025i | 0.0489 + 0.0091i | 0.0001 + 0.0002i |
| 3 | 3 | 0.0392 + 0.0763i | -0.0010 - 0.0011i | -0.1240 + 0.0255i |
| \ell | m | f_r | f_{\Delta\theta} | f_{\Delta\phi} |
|---|---|---|---|---|
| 0 | 0 | 0.5825 | 0 | 0 |
| 1 | 0 | -0.0996 | 2.5793 | -0.4817 |
| 1 | 1 | 0.0237 + 0.0000i | -0.3161 + 1.6600i | 0.4750 + 0.4784i |
| 2 | 0 | 0.0106 | -0.0125 | -0.2309 |
| 2 | 1 | 0.0193 + 0.0061i | 0.0425 - 0.2846i | 0.1581 + 1.5432i |
| 2 | 2 | 0.0591 + 0.0193i | 0.1197 - 0.1800i | 0.8938 - 0.0098i |
| 3 | 0 | 0.0005 | 0.1995 | 0.2132 |
| 3 | 1 | -0.0020 + 0.0003i | 0.0226 - 0.0683i | 0.1065 - 0.3695i |
| 3 | 2 | -0.0358 - 0.0127i | -0.0657 + 0.0821i | -0.2082 - 0.0496i |
| 3 | 3 | -0.0083 - 0.0094i | -0.0272 + 0.0286i | -0.1859 + 0.1900i |

<!-- TODO: figure content -->

[^1]: Experimental implementation of this procedure is strongly discouraged, as the deletion of vertices from live cows is highly unethical.