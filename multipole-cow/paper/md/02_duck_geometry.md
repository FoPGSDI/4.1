# Section 2: Duck Geometry and Multipole Characterization

## 2.1 The Duck Mesh

We adopt as our fiducial *Anas platyrhynchos* the DASSL rubber duck mesh [1], comprising 4732 vertices and 9460 triangular faces. The mesh is watertight and has been normalized to unit bounding-box scale. The original OBJ mesh was converted to OFF format via `trimesh` and centered at the geometric centroid.

The duck is oriented with the bill pointing along the negative $x$-axis, the body extending in the $+x$ direction, and the vertical (head-to-belly) along the $y$-axis. Unlike the cow of Lehmann (2025), which is roughly ellipsoidal, the duck possesses no continuous symmetry — the bill, neck, tail, and body break all axial and planar symmetries, making it a maximally general test case for non-spherical compact object physics.

## 2.2 Volume Integrals and Multipole Moments

We compute the duck's multipole moments using the validated pipeline of the cow reproduction study, which employs fan tetrahedralization from the center of mass and Hammer-Stroud quadrature (4-point rule, exact for degree-2 polynomials).

**Volume and center of mass.** The volume is computed via the divergence theorem:

$$V = \frac{1}{6}\sum_{\text{faces}} \mathbf{v}_0 \cdot (\mathbf{v}_1 \times \mathbf{v}_2)$$

The center of mass for uniform density $\rho = 1$:

$$\mathbf{x}_{\rm COM} = \frac{1}{V}\sum_{\text{tets}} V_{\rm tet}\,\mathbf{x}_{\rm centroid}$$

**Equivalent sphere radius:**

$$R_{\rm eq} = \left(\frac{3V}{4\pi}\right)^{1/3}$$

**Cartesian quadrupole tensor** (traceless, symmetric):

$$Q^C_{ij} = \int_{\mathcal{D}} \rho(\mathbf{x})(3x_i x_j - r^2 \delta_{ij})\,d^3\mathbf{x}$$

**Inertia tensor:**

$$I_{ij} = \int_{\mathcal{D}} \rho(\mathbf{x})(r^2 \delta_{ij} - x_i x_j)\,d^3\mathbf{x}$$

These satisfy the cross-check relation $Q^C_{ij} = \mathrm{Tr}(I)\,\delta_{ij} - 3I_{ij}$.

**Spherical multipole moments** up to $\ell = 10$:

$$Q_\ell^m = \int_{\mathcal{D}} \rho(\mathbf{x})\,\|\mathbf{x}\|^\ell\,Y_\ell^m(\hat{\mathbf{x}})\,d^3\mathbf{x}$$

computed via the surface-integral method (divergence theorem): $Q_\ell^m = \frac{1}{\ell+3}\oint_S (\mathbf{x}\cdot\hat{\mathbf{n}})\,r^\ell\,Y_\ell^m\,dA$.

**Shape-induced quadrupole parameter:**

$$\kappa_{\rm duck} = \frac{Q^C_{\rm eigenvalue}}{M\,R_{\rm eq}^2}$$

This is analogous to the spin-induced $\kappa$ for neutron stars (Laarakkers & Poisson 1999), where Kerr BHs have $\kappa = 1$ and NSs have $\kappa \sim 2$–$14$.

## 2.3 Surface Deformation Spectrum

We parametrize the duck boundary as a deformed sphere:

$$R_{\rm duck}(\theta,\phi) = R_0\left[1 + \sum_{\ell \geq 1}\sum_{m=-\ell}^{\ell} \varepsilon_{\ell m}\,Y_\ell^m(\theta,\phi)\right]$$

where $R_0 = R_{\rm eq}$ and the deformation coefficients $\varepsilon_{\ell m}$ are extracted via regularized least-squares fitting of the real spherical harmonic basis.

**Procedure:**
1. Center the duck mesh at the COM
2. Project each vertex radially onto $S^2$: $(\theta_i, \phi_i) = (\arccos(z_i/r_i),\, \arctan(y_i/x_i))$
3. Build the real SH basis matrix $Y_{\rm real}(\theta_i, \phi_i)$ up to $\ell_{\max} = 10$
4. Solve the Tikhonov-regularized least-squares problem:
$$\min_{\mathbf{c}} \|\mathbf{Y}\mathbf{c} - \mathbf{r}\|^2 + \alpha\sum_\ell \ell(\ell+1)\,|c_{\ell m}|^2$$
where $\alpha = 10^{-4}$ penalizes high-$\ell$ coefficients to suppress aliasing from the non-injective radial projection.

5. Convert to deformation coefficients: $\varepsilon_{\ell m} = c_{\ell m}/R_0 - \delta_{\ell 0}\delta_{m 0}$

**Deformation power spectrum:**

$$P_\ell = \sum_{m=-\ell}^{\ell} |\varepsilon_{\ell m}|^2$$

The dominant deformation is expected at $\ell = 2$ (overall elongation/flattening of the body) and $\ell = 4$ (neck and bill structure), with $P_\ell$ decreasing at higher $\ell$.

## 2.4 Verification

The following consistency checks are performed:

1. **Dipole vanishes:** $|Q_1^m| < 10^{-6}$ for all $m$ (origin at COM)
2. **Tracelessness:** $|\mathrm{Tr}(Q^C)| < 10^{-6}$
3. **Symmetry:** $Q^C_{ij} = Q^C_{ji}$, $I_{ij} = I_{ji}$
4. **Cross-check:** $Q^C_{ij} = \mathrm{Tr}(I)\,\delta_{ij} - 3I_{ij}$ within $< 1\%$
5. **Volume consistency:** $V = Q_0^0 \sqrt{4\pi}$ within $< 5\%$

## 2.5 Figures and Tables

- **Figure 1:** Duck mesh and SH reconstruction at $\ell_{\max} = 0, 2, 4, 8$, and full resolution (existing renders from `results/sphere_to_duck_transition.png`)
- **Figure 2:** Deformation power spectrum $P_\ell$ vs $\ell$ (bar chart)
- **Table 1:** Duck geometric and multipole parameters ($V$, $R_{\rm eq}$, $M$, principal moments $I_1, I_2, I_3$, $Q^C$ eigenvalues, $\kappa_{\rm duck}$ per axis)
- **Table 2:** Deformation coefficients $\varepsilon_{\ell m}$ for $\ell = 1$–$10$ (real and imaginary parts)

---

*Numerical values for Tables 1 and 2 are computed in `production-code/duck_baseline.json` and `production-code/duck_epsilon.dat`.*
