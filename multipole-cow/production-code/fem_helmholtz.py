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
from scipy.special import spherical_jn
try:
    from scipy.special import sph_harm
except ImportError:
    from scipy.special import sph_harm_y as sph_harm

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

def tetrahedralize_tetgen(V, F):
    """Generate a volumetric tet mesh using the tetgen package + pyvista.

    Returns:
        tet_verts: (N, 3) array of vertex coordinates
        tets: (T, 4) array of tet vertex indices
        boundary_verts: set of vertex indices on the boundary
    """
    import tetgen
    import pyvista as pv

    faces_pv = np.column_stack([np.full(len(F), 3), F]).ravel()
    mesh = pv.PolyData(V, faces_pv)

    tg = tetgen.TetGen(mesh)
    tg.tetrahedralize(order=1, mindihedral=10, minratio=1.5)
    grid = tg.grid

    tet_verts = np.array(grid.points, dtype=np.float64)
    # Extract tet connectivity from unstructured grid
    cells = grid.cells.reshape(-1, 5)  # [4, v0, v1, v2, v3] per tet
    tets = cells[:, 1:].astype(np.int64)

    # Boundary vertices: original surface vertices (first len(V) are preserved)
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
    """Assemble global P1 FEM stiffness K and mass M matrices (vectorized).

    Returns:
        K: sparse stiffness matrix (N x N)
        M: sparse mass matrix (N x N)
    """
    n_verts = len(tet_verts)
    n_tets = len(tets)

    # Reference gradients: grad_ref[i] for basis i on reference tet
    grad_ref = np.array([[-1., -1., -1.],
                         [1., 0., 0.],
                         [0., 1., 0.],
                         [0., 0., 1.]])  # (4, 3)

    # Mass template
    M_template = np.full((4, 4), 1.0 / 20.0)
    np.fill_diagonal(M_template, 1.0 / 10.0)

    # Vectorized Jacobians: J[e, k, :] = v[k+1] - v[0]
    v = tet_verts[tets]  # (n_tets, 4, 3)
    J = v[:, 1:, :] - v[:, 0:1, :]  # (n_tets, 3, 3)

    det_J = np.linalg.det(J)  # (n_tets,)
    vol = np.abs(det_J) / 6.0  # (n_tets,)

    # Filter degenerate tets
    valid = vol > 1e-20
    if not np.all(valid):
        print(f"  Warning: {np.sum(~valid)} degenerate tets removed")

    # Inverse transpose for all tets: J_inv_T[e] = inv(J[e]).T
    J_inv_T = np.linalg.inv(J).transpose(0, 2, 1)  # (n_tets, 3, 3)

    # Physical gradients: grad_phys[e, i, :] = grad_ref[i] @ J_inv_T[e]
    # = (4,3) @ (n_tets,3,3) -> need einsum
    grad_phys = np.einsum("ij,ejk->eik", grad_ref, J_inv_T)  # (n_tets, 4, 3)

    # Local stiffness: K_local[e, i, j] = vol[e] * grad_phys[e,i] . grad_phys[e,j]
    K_local = np.einsum("eik,ejk->eij", grad_phys, grad_phys)  # (n_tets, 4, 4)
    K_local *= vol[:, None, None]

    # Local mass
    M_local = vol[:, None, None] * M_template[None, :, :]  # (n_tets, 4, 4)

    # Zero out degenerate tets
    K_local[~valid] = 0
    M_local[~valid] = 0

    # Build COO indices: (n_tets, 4, 4) -> flat
    ii = np.repeat(tets[:, :, None], 4, axis=2)  # (n_tets, 4, 4) row indices
    jj = np.repeat(tets[:, None, :], 4, axis=1)  # (n_tets, 4, 4) col indices

    K = sparse.coo_matrix((K_local.ravel(), (ii.ravel(), jj.ravel())),
                          shape=(n_verts, n_verts)).tocsr()
    M = sparse.coo_matrix((M_local.ravel(), (ii.ravel(), jj.ravel())),
                          shape=(n_verts, n_verts)).tocsr()

    return K, M


def apply_dirichlet_bc(K, M, boundary_verts, n_verts):
    """Apply Dirichlet BCs by eliminating boundary DOFs.

    Instead of zeroing rows/cols (slow for large matrices), extract
    the interior-interior block of K and M.

    Returns:
        K_int, M_int: sparse matrices for interior DOFs only
        interior_indices: array mapping interior DOF index -> global index
    """
    all_verts = set(range(n_verts))
    interior = sorted(all_verts - boundary_verts)
    interior_indices = np.array(interior, dtype=np.int64)

    K = K.tocsc()
    M = M.tocsc()

    # Extract interior-interior submatrix
    K_int = K[interior_indices][:, interior_indices]
    M_int = M[interior_indices][:, interior_indices]

    return K_int.tocsr(), M_int.tocsr(), interior_indices


# ---------------------------------------------------------------------------
# Eigenvalue solver
# ---------------------------------------------------------------------------

def solve_eigenvalues(K, M, k=30, sigma=1.0):
    """Solve K psi = omega^2 M psi using shift-invert eigsh.

    Returns:
        eigenvalues: (k,) array of omega^2 values, sorted ascending
        eigenvectors: (N, k) array
    """
    # Regularize M to handle zero-mass DOFs from degenerate tets
    diag_M = np.array(M.diagonal()).ravel()
    zero_mass = diag_M < 1e-20
    if np.any(zero_mass):
        n_zero = np.sum(zero_mass)
        print(f"  Regularizing {n_zero} zero-mass DOFs in M")
        eps = 1e-12 * np.max(diag_M)
        reg = sparse.diags(np.where(zero_mass, eps, 0.0))
        M = M + reg

    try:
        eigenvalues, eigenvectors = eigsh(K, k=k, M=M, sigma=sigma,
                                         which="LM")
    except Exception as e:
        print(f"  eigsh with sigma={sigma} failed: {e}")
        print("  Trying with larger sigma...")
        eigenvalues, eigenvectors = eigsh(K, k=k, M=M, sigma=100.0,
                                         which="LM")

    # Sort by eigenvalue
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    return eigenvalues, eigenvectors


# ---------------------------------------------------------------------------
# Mode classification via SH projection
# ---------------------------------------------------------------------------

def classify_modes(eigenvectors, tet_verts, boundary_verts, com, l_max=6):
    """Project eigenvectors onto Y_l^m basis at interior vertices.

    Uses interior vertices (where ψ ≠ 0) for angular decomposition.

    Returns:
        dominant_ell: list of dominant ell for each mode
        power_spectra: list of P_l arrays for each mode
    """
    all_verts = set(range(len(tet_verts)))
    interior = sorted(all_verts - boundary_verts)

    # Use a subsample of interior vertices for speed
    max_sample = 5000
    if len(interior) > max_sample:
        rng = np.random.RandomState(42)
        sample_idx = np.array(sorted(rng.choice(interior, max_sample, replace=False)))
    else:
        sample_idx = np.array(interior)

    positions = tet_verts[sample_idx] - com

    # Convert to spherical coordinates
    r = np.linalg.norm(positions, axis=1)
    r = np.maximum(r, 1e-15)
    theta = np.arccos(np.clip(positions[:, 2] / r, -1, 1))
    phi = np.arctan2(positions[:, 1], positions[:, 0])

    # Approximate volume weights (uniform for simplicity)
    d_omega = np.ones(len(sample_idx)) / len(sample_idx)

    dominant_ell = []
    power_spectra = []

    n_modes = eigenvectors.shape[1]
    for alpha in range(n_modes):
        psi = eigenvectors[sample_idx, alpha]

        P_l = np.zeros(l_max + 1)
        for ell in range(l_max + 1):
            for m in range(-ell, ell + 1):
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
    # No BC artifacts since we solve on interior DOFs only
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
    """Visualize mode shapes using a thin shell of interior vertices."""
    try:
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        print("  3D visualization requires mpl_toolkits; skipping mode shapes.")
        return

    bv_list = sorted(boundary_verts)
    bv_pos = tet_verts[bv_list]

    # Find interior vertices close to the surface (within 10% of R_eq)
    all_verts = set(range(len(tet_verts)))
    interior = sorted(all_verts - boundary_verts)

    from scipy.spatial import cKDTree
    tree = cKDTree(bv_pos)
    int_pos = tet_verts[interior]
    dists, _ = tree.query(int_pos)
    # Select near-surface interior vertices
    R_eq = (np.max(bv_pos, axis=0) - np.min(bv_pos, axis=0)).mean() / 2
    shell_mask = dists < 0.08 * R_eq
    shell_idx = np.array(interior)[shell_mask]

    if len(shell_idx) < 100:
        # Fallback: use all interior
        shell_idx = np.array(interior)
    positions = tet_verts[shell_idx]

    # Filter valid modes
    mask = duck_eigenvalues > 1e-3
    valid_indices = np.where(mask)[0]

    n_modes = min(6, len(valid_indices))
    if n_modes == 0:
        print("  No valid modes for 3D visualization.")
        return

    fig = plt.figure(figsize=(16, 9))

    for panel in range(n_modes):
        ax = fig.add_subplot(2, 3, panel + 1, projection="3d")
        mode_idx = valid_indices[panel]
        psi = eigenvectors[shell_idx, mode_idx]
        omega2 = duck_eigenvalues[mode_idx]
        omega = np.sqrt(max(omega2, 0))

        psi_norm = psi / (np.max(np.abs(psi)) + 1e-15)

        ax.scatter(positions[:, 0], positions[:, 1],
                   positions[:, 2],
                   c=psi_norm, cmap="RdBu_r",
                   vmin=-1, vmax=1, s=2, alpha=0.8)
        ax.set_title(f"Mode {panel+1}\n$\\omega = {omega:.3f}$",
                     fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    fig.suptitle("Duck FEM Eigenmode Shapes (near-surface interior vertices)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Mode shapes saved to {outpath}")
    plt.close()


# ---------------------------------------------------------------------------
# Deformed sphere for cross-validation with perturbation theory
# ---------------------------------------------------------------------------

def make_deformed_sphere(R0, epsilon_dict, n_subdiv=20):
    """Create a deformed sphere mesh R(theta,phi) = R0 * [1 + sum eps_lm Y_lm].

    Returns V, F as numpy arrays (OFF-compatible).
    """
    # Create icosphere by subdividing an icosahedron
    import trimesh
    sphere = trimesh.creation.icosphere(subdivisions=n_subdiv, radius=R0)
    V = np.array(sphere.vertices, dtype=np.float64)
    F = np.array(sphere.faces, dtype=np.int64)

    # Compute spherical coordinates
    r = np.linalg.norm(V, axis=1)
    r = np.maximum(r, 1e-15)
    theta = np.arccos(np.clip(V[:, 2] / r, -1, 1))
    phi = np.arctan2(V[:, 1], V[:, 0])

    # Compute deformation
    deform = np.zeros(len(V))
    for (ell, m), eps in epsilon_dict.items():
        Ylm = sph_harm(m, ell, phi, theta)
        deform += np.real(eps * Ylm)

    # Deform radially: R(theta,phi) = R0 * (1 + sum eps Y)
    scale = 1.0 + deform
    V *= scale[:, None]

    return V, F


def cross_validate_perturbative(R0, epsilon_dict, results_dir, label="0.1"):
    """Run FEM on a deformed sphere and compare with perturbative splitting.

    This is the key cross-validation: for small epsilon, the FEM eigenvalue
    shifts should match the Hadamard perturbation matrix eigenvalues.
    """
    print(f"\n{'=' * 72}")
    print(f"Cross-validation: FEM on deformed sphere (eps × {label})")
    print(f"{'=' * 72}")

    # Filter to even ell for the perturbation matrix
    eps_even = {k: v for k, v in epsilon_dict.items() if k[0] % 2 == 0}

    # Create deformed sphere mesh
    print("Creating deformed sphere mesh...")
    V, F = make_deformed_sphere(R0, epsilon_dict, n_subdiv=4)
    print(f"  Vertices: {len(V)}, Faces: {len(F)}")

    # Tetrahedralize
    print("Tetrahedralizing...")
    tet_verts, tets, boundary_verts = tetrahedralize_tetgen(V, F)
    n_interior = len(tet_verts) - len(boundary_verts)
    print(f"  Tet vertices: {len(tet_verts)}, Interior: {n_interior}")

    # Assemble and solve
    print("Assembling FEM matrices...")
    K, M = assemble_fem_matrices(tet_verts, tets)
    K_int, M_int, interior_indices = apply_dirichlet_bc(
        K, M, boundary_verts, len(tet_verts)
    )

    n_modes = min(30, len(interior_indices) - 2)
    print(f"Solving eigenvalue problem (k={n_modes})...")
    eigenvalues, eigvecs_int = solve_eigenvalues(K_int, M_int, k=n_modes)

    # Sphere analytical eigenvalues for reference
    sphere_eigs = sphere_eigenvalues(R0, n_max=2, l_max=4)

    # Group FEM eigenvalues near each sphere level
    print(f"\n{'Sphere level':>15s}  {'FEM eigenvalues (omega)':>40s}  {'Perturbative delta(w^2)':>30s}")
    print("-" * 90)

    from qnm_splitting import compute_splitting

    for omega_s, n, ell, deg in sphere_eigs[:6]:
        omega_s_sq = omega_s ** 2

        # Find FEM eigenvalues near this sphere level
        fem_omega = np.sqrt(np.maximum(eigenvalues, 0))
        nearby = np.where(np.abs(fem_omega - omega_s) / omega_s < 0.3)[0]

        # Perturbative splitting
        try:
            omega0_p, dw2_p, split_w_p = compute_splitting(
                ell, eps_even, n=n, R0=R0
            )
        except Exception:
            split_w_p = np.array([omega_s])
            dw2_p = np.array([0.0])

        fem_str = ", ".join(f"{fem_omega[i]:.4f}" for i in nearby[:deg+2])
        pert_str = ", ".join(f"{dw:.6f}" for dw in dw2_p[:5])

        print(f"  n={n} l={ell} ({deg}x) w={omega_s:.4f}  |  FEM: [{fem_str}]  |  dw2: [{pert_str}]")

    print("-" * 90)

    return eigenvalues


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
    use_tetgen = False
    try:
        import tetgen
        import pyvista
        use_tetgen = True
    except ImportError:
        print("  tetgen/pyvista not available.")

    if use_tetgen:
        print("Generating volumetric tet mesh with tetgen...")
        try:
            tet_verts, tets, boundary_verts = tetrahedralize_tetgen(V, F)
            print(f"  Tet vertices: {len(tet_verts)}, Tets: {len(tets)}")
        except Exception as e:
            print(f"  tetgen tetrahedralization failed: {e}")
            use_tetgen = False

    if not use_tetgen:
        print("Using fan tetrahedralization (fallback)...")
        tet_verts, tets, boundary_verts = tetrahedralize_fallback(V, F, com)
        print(f"  Tet vertices: {len(tet_verts)}, Tets: {len(tets)}")
        print("  WARNING: Fan fallback has ~1 interior vertex; results meaningless.")

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
    print("Applying Dirichlet boundary conditions (extracting interior block)...")
    K_int, M_int, interior_indices = apply_dirichlet_bc(
        K, M, boundary_verts, n_total_verts
    )
    print(f"  Interior DOF count: {len(interior_indices)}")
    print(f"  K_int shape: {K_int.shape}, nnz: {K_int.nnz}")

    # =========================================================================
    # (f) Solve eigenvalue problem
    # =========================================================================
    n_modes = min(50, len(interior_indices) - 2)
    n_modes = max(n_modes, 5)
    print(f"Solving generalized eigenvalue problem (k={n_modes})...")

    try:
        eigenvalues, eigvecs_int = solve_eigenvalues(K_int, M_int, k=n_modes)
        print(f"  Found {len(eigenvalues)} eigenvalues.")
    except Exception as e:
        print(f"  Eigenvalue solve failed: {e}")
        print("  Attempting with fewer modes...")
        n_modes = min(20, len(interior_indices) - 2)
        n_modes = max(n_modes, 5)
        eigenvalues, eigvecs_int = solve_eigenvalues(K_int, M_int, k=n_modes)
        print(f"  Found {len(eigenvalues)} eigenvalues.")

    # Map eigenvectors back to global indices (boundary = 0)
    eigenvectors = np.zeros((n_total_verts, eigvecs_int.shape[1]))
    eigenvectors[interior_indices, :] = eigvecs_int

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

    # Filter near-zero eigenvalues
    mask = eigenvalues > 1e-3
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

    # =========================================================================
    # (l) Cross-validation: FEM on deformed sphere vs perturbation theory
    # =========================================================================
    try:
        from qnm_splitting import load_epsilon_file
        eps_file_01 = os.path.join(script_dir, "duck_epsilon_0.1.dat")
        if os.path.exists(eps_file_01):
            eps_01, meta_01 = load_epsilon_file(eps_file_01)
            R0_val = meta_01.get("R0", R_eq)
            cross_validate_perturbative(R0_val, eps_01, results_dir, label="0.1")
    except Exception as e:
        print(f"\nCross-validation failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
