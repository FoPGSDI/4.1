"""
Render the gradual transition from a sphere to any mesh
at increasing spectral resolution using the mesh Laplacian eigenbasis.

This is the mesh analog of spherical harmonic expansion: the Laplacian
eigenvectors on the mesh play the role of Y_l^m, and truncating at
increasing numbers of eigenvectors gives progressively finer detail.

The eigenvalues lambda_k are analogous to l(l+1), so we label the
reconstructions by the equivalent l_max = sqrt(lambda_k).

Usage:
    python render_sphere_to_cow.py cow
    python render_sphere_to_cow.py bunny
"""

import sys
import os
import argparse
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.dirname(__file__))
from mesh_io import load_off, compute_center_of_mass


def cotangent_laplacian(V, F):
    """Compute the cotangent Laplacian matrix L and mass matrix M.

    L is the (N x N) cotangent-weight Laplacian (negative semi-definite).
    M is the (N x N) diagonal lumped mass matrix (vertex areas).

    Returns: L (sparse), M (sparse)
    """
    n = len(V)
    I, J, Wij = [], [], []

    for face in F:
        for k in range(3):
            i = face[k]
            j = face[(k + 1) % 3]
            opp = face[(k + 2) % 3]

            e1 = V[i] - V[opp]
            e2 = V[j] - V[opp]
            cos_angle = np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2) + 1e-30)
            cos_angle = np.clip(cos_angle, -1, 1)
            sin_angle = np.sqrt(1 - cos_angle**2) + 1e-30
            cot = cos_angle / sin_angle
            w = 0.5 * cot

            I.append(i); J.append(j); Wij.append(w)
            I.append(j); J.append(i); Wij.append(w)

    L = sparse.coo_matrix((Wij, (I, J)), shape=(n, n)).tocsc()
    # Diagonal: L_ii = -sum_j L_ij
    L = L - sparse.diags(np.array(L.sum(axis=1)).flatten())

    # Lumped mass matrix: 1/3 of incident face areas
    v0 = V[F[:, 0]]; v1 = V[F[:, 1]]; v2 = V[F[:, 2]]
    face_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    vertex_mass = np.zeros(n)
    for k in range(3):
        np.add.at(vertex_mass, F[:, k], face_areas / 3.0)
    M = sparse.diags(vertex_mass)

    return L, M


def spectral_decomposition(V_centered, L, M, n_modes):
    """Compute the spectral decomposition of vertex coordinates
    using the mesh Laplacian eigenbasis.

    Solves the generalized eigenvalue problem: L @ phi = -lambda * M @ phi

    Returns:
        eigenvalues: (n_modes,) sorted ascending (smallest first)
        eigenvectors: (N, n_modes) orthonormal w.r.t. M
        coefficients: (n_modes, 3) projection of V onto each eigenvector
    """
    # Solve for the smallest eigenvalues of -L (since L is negative semi-definite)
    # eigsh finds largest eigenvalues of M^-1 @ (-L), equivalent to smallest of L
    eigenvalues, eigenvectors = eigsh(-L, k=n_modes, M=M, sigma=0, which='LM')

    # Sort by eigenvalue (ascending)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Project coordinates onto eigenvectors: c_k = phi_k^T @ M @ V
    M_diag = np.array(M.diagonal())
    coefficients = np.zeros((n_modes, 3))
    for k in range(n_modes):
        for dim in range(3):
            coefficients[k, dim] = np.sum(eigenvectors[:, k] * M_diag * V_centered[:, dim])

    return eigenvalues, eigenvectors, coefficients


def reconstruct_at_n_modes(eigenvectors, coefficients, n):
    """Reconstruct vertex positions using the first n eigenmodes."""
    return eigenvectors[:, :n] @ coefficients[:n, :]


def render_mesh(ax, V, F, title, elev=25, azim=-60,
                global_range=None, global_rmin=None, global_rmax=None):
    """Render a triangulated surface on a 3D axis with coolwarm coloring."""
    r = np.linalg.norm(V, axis=1)
    face_r = (r[F[:, 0]] + r[F[:, 1]] + r[F[:, 2]]) / 3.0

    vmin = global_rmin if global_rmin is not None else r.min()
    vmax = global_rmax if global_rmax is not None else r.max()
    if vmax - vmin < 1e-10:
        vmax = vmin + 1.0
    cmap = plt.cm.coolwarm
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    face_colors = cmap(norm(face_r))

    triangles = V[F]
    poly = Poly3DCollection(triangles, facecolors=face_colors,
                            edgecolors='k', linewidths=0.05, alpha=1.0)
    ax.add_collection3d(poly)

    max_range = global_range if global_range is not None else np.abs(V).max() * 1.1
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)

    ax.set_title(title, fontsize=11, fontweight='bold', pad=2)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])


def render_transition(mesh_path, mesh_name, results_dir, n_modes_max=200,
                      mode_counts=None, elev=25, azim=-60):
    """Render sphere-to-mesh transition using mesh Laplacian spectral decomposition."""
    os.makedirs(results_dir, exist_ok=True)

    if mode_counts is None:
        # Number of eigenmodes for each panel (1=constant, 4≈dipole, 9≈quadrupole, ...)
        mode_counts = [1, 4, 9, 25, 50, 100]

    # Load and center mesh
    V, F = load_off(mesh_path)
    com = compute_center_of_mass(V, F)
    V_centered = V - com
    n_verts = len(V)

    print(f"Mesh '{mesh_name}': {n_verts} vertices, {len(F)} faces")

    # Compute cotangent Laplacian
    print(f"  Computing cotangent Laplacian...")
    L, M = cotangent_laplacian(V, F)

    # Clamp n_modes to available
    n_modes_max = min(n_modes_max, n_verts - 2)
    mode_counts = [min(n, n_modes_max) for n in mode_counts]

    # Spectral decomposition
    print(f"  Computing {n_modes_max} eigenmodes...")
    eigenvalues, eigenvectors, coefficients = spectral_decomposition(
        V_centered, L, M, n_modes_max)
    print(f"  Eigenvalue range: [{eigenvalues[0]:.4f}, {eigenvalues[-1]:.2f}]")

    # Equivalent ℓ for labeling: lambda_k ≈ l(l+1), so l ≈ sqrt(lambda_k)
    # For the first modes: mode 1 is constant (l=0), modes 2-4 are l=1, etc.

    # Reconstruct at each mode count
    labels = [f"{n} modes" for n in mode_counts] + ["Full mesh"]
    reconstructions = []
    for n in mode_counts:
        V_recon = reconstruct_at_n_modes(eigenvectors, coefficients, n)
        rmse = np.sqrt(np.mean((V_recon - V_centered)**2))
        reconstructions.append(V_recon)
        print(f"    {n:4d} modes: RMSE = {rmse:.6f}")
    reconstructions.append(V_centered)

    # Global bounds
    global_range = max(np.abs(Vr).max() for Vr in reconstructions) * 1.1
    all_radii = [np.linalg.norm(Vr, axis=1) for Vr in reconstructions]
    global_rmin = min(rad.min() for rad in all_radii)
    global_rmax = max(rad.max() for rad in all_radii)

    render_kw = dict(elev=elev, azim=azim, global_range=global_range,
                     global_rmin=global_rmin, global_rmax=global_rmax)

    # --- Multi-panel composite ---
    n_panels = len(reconstructions)
    fig = plt.figure(figsize=(4 * n_panels, 4.5))
    for i, (V_rec, label) in enumerate(zip(reconstructions, labels)):
        ax = fig.add_subplot(1, n_panels, i + 1, projection='3d')
        render_mesh(ax, V_rec, F, label, **render_kw)

    title = f"Sphere-to-{mesh_name.capitalize()} via Laplacian Spectral Decomposition"
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out_path = os.path.join(results_dir, f'sphere_to_{mesh_name}_transition.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {out_path}")

    # --- Individual renders ---
    for i, (V_rec, label) in enumerate(zip(reconstructions, labels)):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        render_mesh(ax, V_rec, F, label, **render_kw)
        plt.tight_layout()
        if i < len(mode_counts):
            fname = f'{mesh_name}_modes_{mode_counts[i]:04d}.png'
        else:
            fname = f'{mesh_name}_lmax_full.png'
        out_path = os.path.join(results_dir, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  Saved: {out_path}")

    print(f"Done rendering {mesh_name}!")


def main():
    parser = argparse.ArgumentParser(description="Render sphere-to-mesh spectral transition")
    parser.add_argument('mesh', nargs='?', default='cow',
                        help='Mesh name (cow, bunny) or path to .off file')
    parser.add_argument('--modes', type=int, default=200,
                        help='Maximum number of eigenmodes (default: 200)')
    parser.add_argument('--elev', type=float, default=25, help='Elevation angle')
    parser.add_argument('--azim', type=float, default=-60, help='Azimuth angle')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.dirname(__file__))
    results_dir = os.path.join(base_dir, 'results')
    data_dir = os.path.join(base_dir, 'data')

    if os.path.isfile(args.mesh):
        mesh_path = args.mesh
        mesh_name = os.path.splitext(os.path.basename(args.mesh))[0]
    else:
        mesh_name = args.mesh.lower()
        mesh_path = os.path.join(data_dir, f'{mesh_name}.off')
        if not os.path.isfile(mesh_path):
            print(f"Error: mesh file not found: {mesh_path}")
            sys.exit(1)

    render_transition(mesh_path, mesh_name, results_dir,
                      n_modes_max=args.modes, elev=args.elev, azim=args.azim)


if __name__ == '__main__':
    main()
