#!/usr/bin/env python3
"""
FEM Helmholtz Eigenvalue Solver for the Duck Interior
=====================================================

Solves  -nabla^2 psi = omega^2 psi  on the duck interior D
with Dirichlet BCs psi|_{dD} = 0, using P1 finite elements.

Steps:
  a) Load duck surface mesh
  b) Generate volumetric tet mesh (meshpy.tet or fan fallback)
  c) Assemble P1 stiffness K and mass M matrices
  d) Apply Dirichlet BCs
  e) Solve K psi = omega^2 M psi via shift-invert eigsh
  f) Classify modes by projecting onto Y_l^m
  g) Solve sphere (analytical) for comparison
  h) Print comparison table
  i) Generate eigenvalue spectrum plot
  j) Generate 3D mode shape visualization

Usage:
    python fem_helmholtz.py
"""

import os
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import spherical_jn, sph_harm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add code/ to path for mesh_io, render_sphere_to_cow
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "code"))

from mesh_io import load_off, compute_volume, compute_center_of_mass, tetrahedralize_fan


# ---------------------------------------------------------------------------
# Spherical Bessel zeros (reused from qnm_splitting.py)
# ---------------------------------------------------------------------------

def bessel_zero(ell, n, num_points=10000, r_max=50.0):
    """Find the n-th positive zero of j_ell(x) by bisection on a grid."""
    x = np.linspace(1e-8, r_max, num_points)
    vals = spherical_jn(ell, x)
    zeros = []
    for i in range(len(vals) - 1):
        if vals[i] * vals[i + 1] < 0:
            a, b = x[i], x[i + 1]
            for _ in range(60):
                mid = 0.5 * (a + b)
                if spherical_jn(ell, mid) * spherical_jn(ell, a) < 0:
                    b = mid
                else:
                    a = mid
            zeros.append(0.5 * (a + b))
            if len(zeros) == n:
                break
    if len(zeros) < n:
        raise RuntimeError(f"Could not find {n} zeros of j_{ell}")
    return zeros[n - 1]


# ---------------------------------------------------------------------------
# Tetrahedralization
# ---------------------------------------------------------------------------

def tetrahedralize_meshpy(V, F):
    """Generate a volumetric tet mesh using meshpy.tet (TetGen wrapper).

    Returns:
        tet_verts: (N, 3) array of vertex coordinates
        tets: (T, 4) array of tet vertex indices
        boundary_verts: set of vertex indices on the boundary
    """
    from meshpy.tet import MeshInfo, build

    mesh_info = MeshInfo()
    mesh_info.set_points(V.tolist())

    # Build facets from triangular faces
    facets = []
    for f in F:
        facets.append([list(f)])
    mesh_info.set_facets(facets)

    # Estimate a reasonable max volume from bounding box
    bbox_vol = np.prod(V.max(axis=0) - V.min(axis=0))
    max_vol = bbox_vol / 5000.0  # target ~5000 tets

    mesh = build(mesh_info, max_volume=max_vol,
                 options=MeshInfo.Options(switches="pq"))

    tet_verts = np.array(mesh.points, dtype=np.float64)
    tets = np.array(mesh.elements, dtype=np.int64)

    # Boundary vertices: those from the original surface mesh
    # (first len(V) vertices are preserved by TetGen)
    boundary_verts = set(range(len(V)))

    return tet_verts, tets, boundary_verts


def tetrahedralize_fallback(V, F, com):
    """Fan tetrahedralization fallback (each triangle + COM = one tet).

    Returns:
        tet_verts, tets, boundary_verts
    """
    tets, tet_verts = tetrahedralize_fan(V, F, com)
    # All original surface vertices are boundary vertices
    boundary_verts = set(range(len(V)))
    return tet_verts, tets, boundary_verts


# ---------------------------------------------------------------------------
# P1 FEM Assembly
# ---------------------------------------------------------------------------

def assemble_fem_matrices(tet_verts, tets):
    """Assemble global P1 FEM stiffness K and mass M matrices.

    For each tetrahedron with vertices v0, v1, v2, v3:
      - Volume = |det([v1-v0, v2-v0, v3-v0])| / 6
      - Basis function gradients from inverse Jacobian
      - K_local[i,j] = Volume * (grad phi_i . grad phi_j)
      - M_local[i,j] = Volume * (1/10 if i==j, 1/20 if i!=j)

    Returns:
        K: sparse stiffness matrix (N x N)
        M: sparse mass matrix (N x N)
    """
    n_verts = len(tet_verts)
    n_tets = len(tets)

    # Pre-allocate COO data
    # Each tet contributes 4x4 = 16 entries to both K and M
    rows = np.zeros(n_tets * 16, dtype=np.int64)
    cols = np.zeros(n_tets * 16, dtype=np.int64)
    K_data = np.zeros(n_tets * 16, dtype=np.float64)
    M_data = np.zeros(n_tets * 16, dtype=np.float64)

    # Reference gradients of P1 basis on reference tet
    # phi_0 = 1 - xi - eta - zeta
    # phi_1 = xi, phi_2 = eta, phi_3 = zeta
    # grad_ref = [[-1,-1,-1], [1,0,0], [0,1,0], [0,0,1]]
    grad_ref = np.array([[-1., -1., -1.],
                         [1., 0., 0.],
                         [0., 1., 0.],
                         [0., 0., 1.]])

    # Mass matrix template for P1 tet (consistent mass)
    # M_local[i,j] = Vol * (1/10 if i==j, 1/20 if i!=j)
    M_template = np.full((4, 4), 1.0 / 20.0)
    np.fill_diagonal(M_template, 1.0 / 10.0)

    for e in range(n_tets):
        idx = tets[e]  # 4 vertex indices
        v = tet_verts[idx]  # (4, 3)

        # Jacobian: J[k, :] = v[k+1] - v[0]  for k=0,1,2
        J = np.array([v[1] - v[0], v[2] - v[0], v[3] - v[0]])  # (3, 3)
        det_J = np.linalg.det(J)
        vol = abs(det_J) / 6.0

        if vol < 1e-20:
            continue  # degenerate tet

        # Inverse transpose of J for gradient transformation
        # grad phi_i (physical) = J^{-T} . grad_ref[i]
        J_inv_T = np.linalg.inv(J).T  # (3, 3)

        # Physical gradients: (4, 3)
        grad_phys = grad_ref @ J_inv_T  # (4, 3)

        # Local stiffness: K_local[i,j] = vol * (grad_phys[i] . grad_phys[j])
        K_local = vol * (grad_phys @ grad_phys.T)  # (4, 4)

        # Local mass
        M_local = vol * M_template

        # Scatter into global arrays
        offset = e * 16
        k = 0
        for i in range(4):
            for j in range(4):
                rows[offset + k] = idx[i]
                cols[offset + k] = idx[j]
                K_data[offset + k] = K_local[i, j]
                M_data[offset + k] = M_local[i, j]
                k += 1

    K = sparse.coo_matrix((K_data, (rows, cols)),
                          shape=(n_verts, n_verts)).tocsr()
    M = sparse.coo_matrix((M_data, (rows, cols)),
                          shape=(n_verts, n_verts)).tocsr()

    return K, M


def apply_dirichlet_bc(K, M, boundary_verts, n_verts):
    """Apply Dirichlet BCs by zeroing rows/cols and setting diagonal to 1.

    This shifts boundary eigenvalues to omega^2 = 1 and decouples them
    from interior modes.

    Returns modified K, M as CSR matrices.
    """
    K = K.tolil()
    M = M.tolil()

    for v in boundary_verts:
        K[v, :] = 0
        K[:, v] = 0
        K[v, v] = 1.0
        M[v, :] = 0
        M[:, v] = 0
        M[v, v] = 1.0

    return K.tocsr(), M.tocsr()


# ---------------------------------------------------------------------------
# Eigenvalue solver
# ---------------------------------------------------------------------------

def solve_eigenvalues(K, M, k=30, sigma=0.0):
    """Solve K psi = omega^2 M psi using shift-invert eigsh.

    Returns:
        eigenvalues: (k,) array of omega^2 values, sorted ascending
        eigenvectors: (N, k) array
    """
    try:
        eigenvalues, eigenvectors = eigsh(K, k=k, M=M, sigma=sigma,
                                         which="LM")
    except Exception as e:
        print(f"  eigsh with shift-invert failed: {e}")
        print("  Falling back to standard eigsh (smallest eigenvalues)...")
        eigenvalues, eigenvectors = eigsh(K, k=k, M=M, which="SM")

    # Sort by eigenvalue
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    return eigenvalues, eigenvectors


# ---------------------------------------------------------------------------
# Mode classification via SH projection
# ---------------------------------------------------------------------------

def classify_modes(eigenvectors, tet_verts, boundary_verts, com, l_max=6):
    """Project eigenvectors onto Y_l^m basis at surface vertices.

    Returns:
        dominant_ell: list of dominant ell for each mode
        power_spectra: list of P_l arrays for each mode
    """
    # Get boundary vertex positions relative to COM
    bv_list = sorted(boundary_verts)
    positions = tet_verts[bv_list] - com

    # Convert to spherical coordinates
    r = np.linalg.norm(positions, axis=1)
    r = np.maximum(r, 1e-15)
    theta = np.arccos(np.clip(positions[:, 2] / r, -1, 1))
    phi = np.arctan2(positions[:, 1], positions[:, 0])

    # Approximate solid angle weights (uniform for simplicity)
    d_omega = np.ones(len(bv_list)) * 4.0 * np.pi / len(bv_list)

    dominant_ell = []
    power_spectra = []

    n_modes = eigenvectors.shape[1]
    for alpha in range(n_modes):
        psi = eigenvectors[bv_list, alpha]

        P_l = np.zeros(l_max + 1)
        for ell in range(l_max + 1):
            for m in range(-ell, ell + 1):
                # sph_harm uses (m, ell, phi, theta) convention
                Ylm = sph_harm(m, ell, phi, theta)
                a_lm = np.sum(psi * np.conj(Ylm) * d_omega)
                P_l[ell] += abs(a_lm) ** 2

        dominant_ell.append(np.argmax(P_l))
        power_spectra.append(P_l)

    return dominant_ell, power_spectra


# ---------------------------------------------------------------------------
# Sphere analytical eigenvalues
# ---------------------------------------------------------------------------

def sphere_eigenvalues(R0, n_max=2, l_max=5):
    """Compute analytical Helmholtz eigenvalues for a sphere.

    Returns list of (omega, n, ell, degeneracy) sorted by omega.
    """
    results = []
    for n in range(1, n_max + 1):
        for ell in range(0, l_max + 1):
            try:
                z = bessel_zero(ell, n)
                omega = z / R0
                results.append((omega, n, ell, 2 * ell + 1))
            except RuntimeError:
                pass
    results.sort(key=lambda x: x[0])
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_eigenvalue_spectrum(duck_eigenvalues, duck_dominant_ell,
                            sphere_results, R0, outpath):
    """Plot eigenvalue spectrum: duck vs sphere as horizontal bars."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Color map for ell values
    ell_colors = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c",
                  3: "#d62728", 4: "#9467bd", 5: "#8c564b"}

    # Filter out boundary eigenvalues (omega^2 ~ 1 from Dirichlet BC)
    mask = duck_eigenvalues > 1e-3
    # Also filter out eigenvalues that are exactly 1.0 (BC artifacts)
    mask &= np.abs(duck_eigenvalues - 1.0) > 0.01
    duck_omega2 = duck_eigenvalues[mask]
    duck_omega = np.sqrt(np.maximum(duck_omega2, 0))
    duck_ell_filtered = [duck_dominant_ell[i] for i, m in enumerate(mask) if m]

    # Sphere levels
    x_sphere = 0.25
    x_duck = 0.75
    bar_width = 0.12

    # Find frequency range
    all_freqs = list(duck_omega)
    for omega, n, ell, deg in sphere_results:
        all_freqs.append(omega)
    if not all_freqs:
        print("  No eigenvalues to plot.")
        plt.close()
        return
    f_min = min(all_freqs) * 0.9
    f_max = max(all_freqs) * 1.05

    # Sphere levels (degenerate)
    for omega, n, ell, deg in sphere_results:
        if omega > f_max:
            continue
        color = ell_colors.get(ell, "#333333")
        ax.plot([x_sphere - bar_width, x_sphere + bar_width],
                [omega, omega], color=color, linewidth=2.5,
                solid_capstyle="round")
        ax.text(x_sphere - bar_width - 0.02, omega,
                f"$\\ell={ell}$ ({deg}$\\times$)",
                ha="right", va="center", fontsize=8, color=color)

    # Duck levels (split)
    for i, omega in enumerate(duck_omega):
        if omega > f_max:
            continue
        if i < len(duck_ell_filtered):
            ell = duck_ell_filtered[i]
        else:
            ell = -1
        color = ell_colors.get(ell, "#333333")
        ax.plot([x_duck - bar_width, x_duck + bar_width],
                [omega, omega], color=color, linewidth=1.8,
                solid_capstyle="round")

    # Connect sphere to duck
    sphere_by_ell = {}
    for omega, n, ell, deg in sphere_results:
        if ell not in sphere_by_ell:
            sphere_by_ell[ell] = []
        sphere_by_ell[ell].append(omega)

    for i, omega in enumerate(duck_omega):
        if omega > f_max or i >= len(duck_ell_filtered):
            continue
        ell = duck_ell_filtered[i]
        color = ell_colors.get(ell, "#333333")
        # Find nearest sphere level for this ell
        if ell in sphere_by_ell:
            nearest_sphere = min(sphere_by_ell[ell],
                                key=lambda w: abs(w - omega))
            ax.plot([x_sphere + bar_width + 0.01,
                     x_duck - bar_width - 0.01],
                    [nearest_sphere, omega],
                    color=color, linewidth=0.4, alpha=0.3, linestyle="--")

    # Labels
    ax.text(x_sphere, f_max * 1.02, "Sphere\n(degenerate)",
            ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.text(x_duck, f_max * 1.02, "Duck\n(split)",
            ha="center", va="bottom", fontsize=12, fontweight="bold")

    ax.set_xlim(0, 1)
    ax.set_ylim(f_min, f_max * 1.08)
    ax.set_ylabel(r"Frequency $\omega / R_0^{-1}$", fontsize=13)
    ax.set_title("FEM Helmholtz Eigenvalue Spectrum: Sphere vs. Duck",
                 fontsize=14, fontweight="bold")
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Legend
    for ell in sorted(ell_colors.keys()):
        ax.plot([], [], color=ell_colors[ell], linewidth=2,
                label=f"$\\ell = {ell}$")
    ax.legend(loc="lower right", fontsize=9)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"Eigenvalue spectrum saved to {outpath}")
    plt.close()


def plot_mode_shapes(eigenvectors, tet_verts, boundary_verts,
                     duck_eigenvalues, outpath):
    """Visualize the first few mode shapes on the duck surface."""
    try:
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError:
        print("  3D visualization requires mpl_toolkits; skipping mode shapes.")
        return

    bv_list = sorted(boundary_verts)
    positions = tet_verts[bv_list]

    # Filter valid modes
    mask = duck_eigenvalues > 1e-3
    mask &= np.abs(duck_eigenvalues - 1.0) > 0.01
    valid_indices = np.where(mask)[0]

    n_modes = min(6, len(valid_indices))
    if n_modes == 0:
        print("  No valid modes for 3D visualization.")
        return

    fig = plt.figure(figsize=(16, 9))

    for panel in range(n_modes):
        ax = fig.add_subplot(2, 3, panel + 1, projection="3d")
        mode_idx = valid_indices[panel]
        psi = eigenvectors[bv_list, mode_idx]
        omega2 = duck_eigenvalues[mode_idx]
        omega = np.sqrt(max(omega2, 0))

        # Normalize for color mapping
        psi_norm = psi / (np.max(np.abs(psi)) + 1e-15)

        scatter = ax.scatter(positions[:, 0], positions[:, 1],
                            positions[:, 2],
                            c=psi_norm, cmap="RdBu_r",
                            vmin=-1, vmax=1, s=1, alpha=0.8)
        ax.set_title(f"Mode {panel+1}\n$\\omega = {omega:.3f}$",
                     fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    fig.suptitle("Duck FEM Eigenmode Shapes (surface vertex values)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Mode shapes saved to {outpath}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    # =========================================================================
    # (a/b) Load duck mesh
    # =========================================================================
    duck_path = os.path.join(BASE_DIR, "data", "duck.off")
    print(f"Loading duck mesh from {duck_path}...")
    V, F = load_off(duck_path)
    n_verts = len(V)
    n_faces = len(F)
    print(f"  Vertices: {n_verts}, Faces: {n_faces}")

    volume = compute_volume(V, F)
    com = compute_center_of_mass(V, F)
    R_eq = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)
    print(f"  Volume: {volume:.8f}")
    print(f"  COM: ({com[0]:.6f}, {com[1]:.6f}, {com[2]:.6f})")
    print(f"  R_eq: {R_eq:.8f}")

    # =========================================================================
    # (c) Generate volumetric tet mesh
    # =========================================================================
    use_meshpy = False
    try:
        import meshpy.tet
        use_meshpy = True
    except ImportError:
        print("  meshpy not available. Attempting pip install...")
        try:
            import subprocess
            subprocess.check_call([sys.executable, "-m", "pip", "install",
                                   "meshpy", "--quiet"])
            import meshpy.tet
            use_meshpy = True
            print("  meshpy installed successfully.")
        except Exception as e:
            print(f"  Could not install meshpy: {e}")
            print("  Falling back to fan tetrahedralization.")

    if use_meshpy:
        print("Generating volumetric tet mesh with meshpy (TetGen)...")
        try:
            tet_verts, tets, boundary_verts = tetrahedralize_meshpy(V, F)
            print(f"  Tet vertices: {len(tet_verts)}, Tets: {len(tets)}")
        except Exception as e:
            print(f"  meshpy tetrahedralization failed: {e}")
            print("  Falling back to fan tetrahedralization.")
            use_meshpy = False

    if not use_meshpy:
        print("Using fan tetrahedralization (fallback)...")
        tet_verts, tets, boundary_verts = tetrahedralize_fallback(V, F, com)
        print(f"  Tet vertices: {len(tet_verts)}, Tets: {len(tets)}")
        print("  NOTE: Fan tetrahedralization is a crude approximation.")
        print("        Results are qualitative only. Install meshpy for")
        print("        accurate FEM eigenvalues.")

    n_total_verts = len(tet_verts)
    n_boundary = len(boundary_verts)
    n_interior = n_total_verts - n_boundary
    print(f"  Boundary vertices: {n_boundary}")
    print(f"  Interior vertices: {n_interior}")

    # =========================================================================
    # (d) Assemble FEM matrices
    # =========================================================================
    print("Assembling P1 FEM stiffness and mass matrices...")
    K, M = assemble_fem_matrices(tet_verts, tets)
    print(f"  K shape: {K.shape}, nnz: {K.nnz}")
    print(f"  M shape: {M.shape}, nnz: {M.nnz}")

    # =========================================================================
    # (e) Apply Dirichlet BCs
    # =========================================================================
    print("Applying Dirichlet boundary conditions...")
    K, M = apply_dirichlet_bc(K, M, boundary_verts, n_total_verts)

    # =========================================================================
    # (f) Solve eigenvalue problem
    # =========================================================================
    n_modes = min(30, n_interior - 1) if n_interior > 1 else 5
    n_modes = max(n_modes, 5)
    print(f"Solving generalized eigenvalue problem (k={n_modes})...")

    try:
        eigenvalues, eigenvectors = solve_eigenvalues(K, M, k=n_modes)
        print(f"  Found {len(eigenvalues)} eigenvalues.")
    except Exception as e:
        print(f"  Eigenvalue solve failed: {e}")
        print("  Attempting with fewer modes...")
        n_modes = min(10, n_interior - 1) if n_interior > 1 else 3
        n_modes = max(n_modes, 3)
        eigenvalues, eigenvectors = solve_eigenvalues(K, M, k=n_modes)
        print(f"  Found {len(eigenvalues)} eigenvalues.")

    # =========================================================================
    # (g) Classify modes by SH projection
    # =========================================================================
    print("Classifying modes via SH projection...")
    dominant_ell, power_spectra = classify_modes(
        eigenvectors, tet_verts, boundary_verts, com, l_max=6
    )

    # =========================================================================
    # (h) Sphere analytical eigenvalues
    # =========================================================================
    print(f"\nSphere analytical eigenvalues (R0 = {R_eq:.6f}):")
    sphere_results = sphere_eigenvalues(R_eq, n_max=2, l_max=5)

    # =========================================================================
    # (i) Print comparison table
    # =========================================================================
    print("\n" + "=" * 80)
    print("FEM Helmholtz Eigenvalue Results: Duck vs Sphere")
    print("=" * 80)

    # Filter duck eigenvalues (remove BC artifacts)
    mask = eigenvalues > 1e-3
    mask &= np.abs(eigenvalues - 1.0) > 0.01
    valid_idx = np.where(mask)[0]

    print(f"\n{'Mode':>5s}  {'omega^2_duck':>14s}  {'omega_duck':>12s}  "
          f"{'dom. ell':>8s}  {'omega_sphere':>13s}  {'ell_sphere':>10s}  "
          f"{'dw/w':>10s}")
    print("-" * 80)

    # Match duck modes to sphere levels
    sphere_freqs = [(w, n, ell, deg) for w, n, ell, deg in sphere_results]
    sphere_freqs.sort(key=lambda x: x[0])

    n_display = min(20, len(valid_idx))
    for i in range(n_display):
        idx = valid_idx[i]
        w2 = eigenvalues[idx]
        w = np.sqrt(max(w2, 0))
        ell_dom = dominant_ell[idx]

        # Find nearest sphere eigenvalue with same ell
        best_match = None
        best_dist = float("inf")
        for sw, sn, sell, sdeg in sphere_freqs:
            if sell == ell_dom:
                dist = abs(sw - w)
                if dist < best_dist:
                    best_dist = dist
                    best_match = (sw, sn, sell, sdeg)

        if best_match is not None:
            sw, sn, sell, sdeg = best_match
            dw_rel = (w - sw) / sw if sw > 0 else 0.0
            print(f"{i+1:5d}  {w2:14.6f}  {w:12.6f}  {ell_dom:8d}  "
                  f"{sw:13.6f}  {sell:10d}  {dw_rel:+10.4f}")
        else:
            print(f"{i+1:5d}  {w2:14.6f}  {w:12.6f}  {ell_dom:8d}  "
                  f"{'---':>13s}  {'---':>10s}  {'---':>10s}")

    print("-" * 80)

    # Also print sphere reference table
    print(f"\nSphere reference eigenvalues (R0 = {R_eq:.6f}):")
    print(f"{'n':>4s}  {'ell':>4s}  {'z_{n,ell}':>12s}  {'omega':>12s}  "
          f"{'degeneracy':>10s}")
    print("-" * 50)
    for omega, n, ell, deg in sphere_freqs[:15]:
        z = omega * R_eq
        print(f"{n:4d}  {ell:4d}  {z:12.6f}  {omega:12.6f}  {deg:10d}")

    # =========================================================================
    # (j) Eigenvalue spectrum plot
    # =========================================================================
    print("\nGenerating eigenvalue spectrum plot...")
    plot_path = os.path.join(results_dir, "fem_eigenvalue_spectrum.png")
    plot_eigenvalue_spectrum(eigenvalues, dominant_ell,
                            sphere_results, R_eq, plot_path)

    # =========================================================================
    # (k) 3D mode shape visualization
    # =========================================================================
    print("Generating 3D mode shape visualization...")
    mode_path = os.path.join(results_dir, "fem_mode_shapes.png")
    try:
        plot_mode_shapes(eigenvectors, tet_verts, boundary_verts,
                         eigenvalues, mode_path)
    except Exception as e:
        print(f"  Mode shape visualization failed: {e}")

    print("\n=== FEM Helmholtz analysis complete ===")


if __name__ == "__main__":
    main()
