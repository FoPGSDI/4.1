# Higher multipoles of bovine potentials

Let us first consider the simplest definition of the spherical cow, as it appears in the context of gravitational problems. Of course, the same methods apply to the study of the electrostatic potential in the case of a charged cow. However, charging cows are rarely observed in nature relative to charging bulls [lott:1981,douglass:1999]. (Still, caution is advised in experimental settings, per Ref. [murphy:2010].)

When computing the gravitational potential $\phi(\boldsymbol{\mathrm{x}})$ at a point $\boldsymbol{\mathrm{x}}$ outside the cow, the Laplace equation can be expanded in spherical harmonics, which yields the multipole expansion
$$

    \phi(\boldsymbol{\mathrm{x}}) = -G\sum_{\ell = 0}^{\infty}
        \frac{
            \left[\frac{4\pi}{2\ell+1}\right]^{1/2}
        }{\|\boldsymbol{\mathrm{x}}\|^{\ell+1}}
        \sum_{m=-\ell}^\ell (-1)^m Y_\ell^{-m}(\boldsymbol{\mathrm{\hat x}})
        Q_\ell^m
        ,
$$
where the spherical multipole moment $Q_\ell^m$ is defined by
$$
    Q_\ell^m\equiv \int_{\mathcal C}\,\mathrm{d}^3\boldsymbol{\mathrm{x}}\,
    \rho(\boldsymbol{\mathrm{x}})\,\|\boldsymbol{\mathrm{x}}\|^\ell Y_\ell^m(\boldsymbol{\mathrm{\hat x}})
    ,
$$
given a mass density $\rho(\boldsymbol{\mathrm{x}})$. I will take $\rho$ to be constant over the cow, and I will work in “cow coordinates,” where the $x$ axis is aligned with the forward direction of the cow, the $y$ axis is orthogonal to the ground, and the positive $z$ axis points to the cow’s right. I also retain the scaling of the benchmark cow, so the length units correspond to a bounding box for the cow with dimensions $(1.044, 0.6397, 0.3403)$ in $x$, $y$, and $z$, respectively. I call these units “benchmark units.” The first few multipole moments for the benchmark cow are given in these units in [ref:tab:mass-moments]. The monopole corresponds to the spherical cow, so the higher multipole coefficients indicate the size of corrections to the SCA.

| \ell | m | Q_\ell^m | \ell | m | Q_\ell^m |
|---|---|---|---|---|---|
| 0 | 0 | 0.0539 | 3 | 3 | -0.0003 + 0.0004i |
| 2 | 0 | -0.0029 | 4 | 0 | 0.0003 |
| 2 | 1 | -6.43×10^{-6} + 3.57×10^{-7}i | 4 | 1 | 5.62×10^{-7} + 1.40×10^{-6}i |
| 2 | 2 | 0.0027 - 0.0008i | 4 | 2 | -0.0003 + 0.0001i |
| 3 | 0 | -5.79×10^{-6} | 4 | 3 | -1.42×10^{-6} + 5.35×10^{-7}i |
| 3 | 1 | 0.0003 - 0.0001i | 4 | 4 | 0.0001 - 0.0003i |
| 3 | 2 | -1.00×10^{-6} - 2.41×10^{-6}i | 5 | 0 | 1.05×10^{-6} |

We can now use these coefficients to study gravitational phenomena in the bovine potential. There are several potential applications. For instance, the quadrupole encodes the equivalent of Earth’s equatorial bulge for a heavy cow. Additionally, a cow with a nonvanishing quadrupole moment will experience a torque in a nonuniform gravitational field, so a cow falling from a high altitude in Earth’s gravitational field will favor a particular orientation. Thus, the quadrupole moment is crucial to answering the question of whether a dropped cow tends to land on its feet. For the same reason, the quadrupole plays a key role in the tidal locking of cows in orbit.[^1] Neither of these cases is amenable to experimental study, and in each, the SCA clearly fails.

However, for the moment, let us focus on another application of the quadrupole moment: computation of the rate of emission of gravitational radiation. A freely rotating cow in vacuum will slow down as it loses energy and angular momentum to gravitational waves. We can compute the rate of spindown using the quadrupole formula for gravitational radiation. When written in terms of the components of the spherical quadrupole moment, this reads
$$
    \dot E_{\mathrm{quad}} =
    \frac{3G}{8\pi c^5}\sum_{m=-2}^2
    \left\langle\left|\overset{\cdots}{Q}_2^m\right|^2\right\rangle
    .
$$
In principle, one could take the spherical multipole coefficients from [ref:tab:mass-moments], transform them under rotations, and use this transformation behavior to evaluate $\overset{\cdots}{Q}_2^m$ for a rotation about a given axis. It is simpler, however, to work with Cartesian multipole coefficients for this purpose, since the Cartesian quadrupole coefficients take the form of a symmetric tensor with components $Q^{\mathrm{C}}_{ij} = \int\mathrm{d}^3\boldsymbol{\mathrm{x}}\,\rho(\boldsymbol{\mathrm{x}})\,(3x_ix_j - r^2\delta_{ij})$. The Cartesian quadrupole tensor of the benchmark cow in cow coordinates and benchmark units with $\rho(\boldsymbol{\mathrm{x}}) = 1$ is
$$
    Q^{\mathrm{C}} =
    \begin{pmatrix}
        4.23004 &  \phantom{-}0.83531 &  \phantom{-}0.00704 \\
        0.83531 & -1.64753 &  \phantom{-}0.00039 \\
        0.00704 &  \phantom{-}0.00039 & -2.58252
    \end{pmatrix}
    \times10^{-3}
    ,
$$
and the radiated power is given in terms of the $Q_{ij}^{\mathrm{C}}$ by
$$
    \dot E_{\mathrm{quad}} = \frac{G}{45c^5}
        \left\langle\overset{\cdots}{Q}^{\mathrm{C}}_{ij}\overset{\cdots}{Q}^{\mathrm{C}\,ij}
        \right\rangle
    .
$$
Now, suppose the cow is rotating about an axis $\boldsymbol{\mathrm{\hat a}}$ with angular frequency $\omega$. If $R_{ij}(\theta)$ denotes the rotation matrix about $\boldsymbol{\mathrm{\hat a}}$ by an angle $\theta$, then at any time $t$, we can write the quadrupole tensor as
$$
    Q_{ij}^{\mathrm{C}}(t) = R^{ik}(\omega t)\,R^{j\ell}(\omega t)\,
        Q_{k\ell}^{\mathrm{C}}(0)
    ,
$$
so the third derivative in time is given by
$$

    \overset{\cdots}{Q}_{ij}^{\mathrm{C}}(t) =
        Q_{k\ell}^{\mathrm{C}}(0)\,
        \partial_t^3\left[R^{ik}(\omega t)\,R^{j\ell}(\omega t)\right]
    .
$$
For a simple example, consider a constant rotation about the $y$ axis, as might be experienced by a cow wearing misaligned rollerskates. Here the rotation matrix is
$$
    R(\theta) =
    \begin{pmatrix}
        \cos\theta & 0 & \sin\theta \\
        0 & 1 & 0 \\
        -\sin\theta & 0 & \cos\theta
    \end{pmatrix}
    .
$$
Inserting this into [ref:eq:quadrupole-derivative] and averaging the result over a full period gives $\langle\overset{\cdots}{Q}^{\mathrm{C}}_{ij}\overset{\cdots}{Q}^{\mathrm{C}\,ij}\rangle \approx 0.00149 \omega^6$.

We can now restore physical dimensions to this result, which has units of $M^2L^4\omega^6$, where $M$ and $L$ are the units of mass and length for the benchmark cow. The average length of a cow is about 2.5\;\mathrm{m}, and the length in benchmark units is 1.044, so let us say that $L=2.39\;\mathrm{m}$. To determine the mass unit, recall that we set $\rho(\boldsymbol{\mathrm{x}}) = 1$, and the actual density of a cow is about 1\;\mathrm{g/cm^3}. Setting $M/L^3 = 1\;\mathrm{g/cm^3}$ gives $M = 13652\;\mathrm{kg}$. Replacing these quantities, in physical units, we obtain $\langle\overset{\cdots}{Q}^{\mathrm{C}}_{ij}\overset{\cdots}{Q}^{\mathrm{C}\,ij}\rangle \approx 9 \times 10^{12}\;\mathrm{kg^2m^4}\times\omega^6$, or
$$
    \dot E_{\mathrm{quad}} \approx 5.5 \times 10^{-41}\;\mathrm{erg/s} \times
        \left(\frac{\omega}{1\;\mathrm{Hz}}\right)^6
        .
$$
The inertia tensor in benchmark units is
$$
    I =
    \begin{pmatrix}
         \phantom{-}7.95079 & -2.78437 & -0.02348\\
        -2.78437 &  \phantom{-}27.5426 & -0.00130\\
        -0.02348 & -0.00130 &  \phantom{-}30.6593
    \end{pmatrix}
    \times10^{-4}
    ,
$$
so we can readily determine the kinetic energy associated with rotation about the $y$ axis, and thus obtain the timescale on which the rotation slows (neglecting the details of radiation of angular momentum):
$$
    \frac{E}{\dot E_{\mathrm{quad}}} \approx 1.9 \times 10^{49}\;\mathrm{s}
    \left(\frac{\omega}{1\;\mathrm{Hz}}\right)^{-4}
    .
$$
Since there is no gravitational radiation due to the monopole or dipole moments, this effect is completely lost when the mass distribution is regarded as purely monopolar. Thus, this is another case in which the SCA is wholly inadequate---under the SCA, one would incorrectly conclude that the cow never slows its spinning!

Again, I stress that while the spherical cow naturally corresponds to the spherically-symmetric monopole term, the monopole term actually represents the potential of a point mass, not a ball. The multipole expansion of the potential does not readily encode the full geometry of the source. In particular, the expansion is only well defined outside the source distribution---even including arbitrarily many terms, the potential of [ref:eq:spherical-multipole-potential] is sourceless everywhere away from the origin. We will study corrections to the spherical cow’s geometry in the next section.

![](sdf_0.png)

[^1]: A similar computation would apply for pigs in orbit, but these are only expected to be observed when pigs fly.