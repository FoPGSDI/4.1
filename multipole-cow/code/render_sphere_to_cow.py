"""
Render the gradual transition from a sphere to any mesh
at increasing spherical harmonic order l_max.

Approach: Project each vertex radially onto S^2 from the center of mass,
then fit SH coefficients for each Cartesian coordinate (x,y,z) using
regularized least-squares. The regularization prevents high-ℓ blowup
caused by the non-injective projection of non-star-shaped meshes.

Usage:
    python render_sphere_to_cow.py cow
    python render_sphere_to_cow.py bunny
"""

import sys
import os
import argparse
import numpy as np
from scipy.special import sph_harm_y
from scipy import sparse
from scipy.sparse.linalg import lsqr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.dirname(__file__))
from mesh_io import load_off, compute_center_of_mass


def cartesian_to_spherical(xyz):
    """Convert Cartesian to spherical. theta: polar [0,pi], phi: azimuthal [0,2pi]."""
    r = np.linalg.norm(xyz, axis=1)
    r_safe = np.maximum(r, 1e-15)
    theta = np.arccos(np.clip(xyz[:, 2] / r_safe, -1, 1))
    phi = np.arctan2(xyz[:, 1], xyz[:, 0]) % (2 * np.pi)
    return r, theta, phi


def build_sh_matrix(theta, phi, l_max):
    """Build the SH basis matrix Y[i, j] = Y_{l_j}^{m_j}(theta_i, phi_i).

    Returns: Y_real (N, n_coeffs) real matrix (stacking real and imag parts),
             coeff_index list of (l, m).
    """
    coeff_index = []
    for ell in range(l_max + 1):
        for m in range(-ell, ell + 1):
            coeff_index.append((ell, m))

    n_verts = len(theta)
    n_coeffs = len(coeff_index)

    # Build complex SH matrix, then convert to real system
    # For real-valued targets (x,y,z coords), we need:
    #   x_i = sum_j Re(c_j * Y_j(theta_i, phi_i))
    # This is a real linear system if we split c_j into real and imag parts.
    # But simpler: use real spherical harmonics directly.
    # Y_lm_real = Re(Y_lm) for m >= 0, Im(Y_l|m|) for m < 0 (with normalization)

    # Actually, use the real SH basis:
    #   S_{l,0} = Y_l^0 (real)
    #   S_{l,m,c} = (Y_l^m + (-1)^m Y_l^{-m}) / sqrt(2) for m > 0  (cos-like)
    #   S_{l,m,s} = (Y_l^m - (-1)^m Y_l^{-m}) / (i*sqrt(2)) for m > 0  (sin-like)
    real_index = []  # list of (l, m, 'c'/'s'/'0')
    for ell in range(l_max + 1):
        real_index.append((ell, 0, '0'))
        for m in range(1, ell + 1):
            real_index.append((ell, m, 'c'))
            real_index.append((ell, m, 's'))

    n_real = len(real_index)
    Y_real = np.zeros((n_verts, n_real))

    for j, (ell, m, kind) in enumerate(real_index):
        if kind == '0':
            Y_real[:, j] = np.real(sph_harm_y(ell, 0, theta, phi))
        elif kind == 'c':
            Ylm = sph_harm_y(ell, m, theta, phi)
            Ylnm = sph_harm_y(ell, -m, theta, phi)
            Y_real[:, j] = np.real((Ylm + (-1)**m * Ylnm) / np.sqrt(2))
        elif kind == 's':
            Ylm = sph_harm_y(ell, m, theta, phi)
            Ylnm = sph_harm_y(ell, -m, theta, phi)
            Y_real[:, j] = np.real((Ylm - (-1)**m * Ylnm) / (1j * np.sqrt(2)))

    return Y_real, real_index


def fit_sh_coefficients(Y_real, V_centered, real_index, l_max, alpha=1e-4):
    """Fit SH coefficients via Tikhonov-regularized least squares.

    Solves: min ||Y @ c - v||^2 + alpha * sum_l l(l+1) * |c_lm|^2
    for each coordinate independently.

    The regularization penalizes high-ℓ coefficients proportionally to l(l+1),
    preventing blowup from aliasing while preserving low-ℓ accuracy.

    Returns: coeffs (n_real, 3) array of real SH coefficients.
    """
    n_real = Y_real.shape[1]

    # Build regularization weights: penalty proportional to l(l+1)
    reg_weights = np.zeros(n_real)
    for j, (ell, m, kind) in enumerate(real_index):
        reg_weights[j] = alpha * ell * (ell + 1)

    # Augmented system: [Y; sqrt(reg) * I] @ c = [v; 0]
    D = sparse.diags(np.sqrt(reg_weights))

    coeffs = np.zeros((n_real, 3))
    for k in range(3):
        # Use scipy lsqr for efficiency
        result = lsqr(
            sparse.vstack([sparse.csc_matrix(Y_real), D]),
            np.concatenate([V_centered[:, k], np.zeros(n_real)]),
        )
        coeffs[:, k] = result[0]

    return coeffs


def reconstruct_at_lmax(Y_real, coeffs, real_index, l_max_cut):
    """Reconstruct vertex positions using only coefficients up to l_max_cut."""
    mask = np.array([ell <= l_max_cut for ell, m, kind in real_index])
    return Y_real[:, mask] @ coeffs[mask, :]


def render_mesh(ax, V, F, title, elev=25, azim=-60,
                global_range=None, global_rmin=None, global_rmax=None):
    """Render a triangulated surface on a 3D axis."""
    r = np.linalg.norm(V, axis=1)
    face_r = (r[F[:, 0]] + r[F[:, 1]] + r[F[:, 2]]) / 3.0

    vmin = global_rmin if global_rmin is not None else r.min()
    vmax = global_rmax if global_rmax is not None else r.max()
    if vmax - vmin < 1e-10:
        vmax = vmin + 1.0
    face_colors = plt.cm.coolwarm(plt.Normalize(vmin=vmin, vmax=vmax)(face_r))

    poly = Poly3DCollection(V[F], facecolors=face_colors,
                            edgecolors='k', linewidths=0.05, alpha=1.0)
    ax.add_collection3d(poly)

    mr = global_range if global_range is not None else np.abs(V).max() * 1.1
    ax.set_xlim(-mr, mr); ax.set_ylim(-mr, mr); ax.set_zlim(-mr, mr)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=2)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])


def render_transition(mesh_path, mesh_name, results_dir, l_max_full=20,
                      l_max_values=None, elev=25, azim=-60):
    """Render sphere-to-mesh transition using regularized SH fitting."""
    os.makedirs(results_dir, exist_ok=True)
    if l_max_values is None:
        l_max_values = [0, 1, 2, 4, 8, 16]

    V, F = load_off(mesh_path)
    com = compute_center_of_mass(V, F)
    V_centered = V - com
    print(f"Mesh '{mesh_name}': {len(V)} vertices, {len(F)} faces")

    _, theta, phi = cartesian_to_spherical(V_centered)

    # Build real SH basis matrix
    print(f"  Building SH basis matrix up to l_max = {l_max_full}...")
    Y_real, real_index = build_sh_matrix(theta, phi, l_max_full)
    n_coeffs = Y_real.shape[1]
    print(f"  Basis size: {n_coeffs} real SH functions")

    # Fit coefficients with regularization
    print(f"  Fitting regularized SH coefficients...")
    coeffs = fit_sh_coefficients(Y_real, V_centered, real_index, l_max_full, alpha=1e-4)

    # Reconstruct at each l_max
    labels = [f"$\\ell_{{\\max}} = {l}$" for l in l_max_values] + ["Full mesh"]
    reconstructions = []
    for lm in l_max_values:
        V_recon = reconstruct_at_lmax(Y_real, coeffs, real_index, lm)
        rmse = np.sqrt(np.mean((V_recon - V_centered)**2))
        reconstructions.append(V_recon)
        print(f"    l_max={lm:3d}: RMSE = {rmse:.6f}")
    reconstructions.append(V_centered)

    # Global bounds
    global_range = max(np.abs(Vr).max() for Vr in reconstructions) * 1.1
    all_r = [np.linalg.norm(Vr, axis=1) for Vr in reconstructions]
    global_rmin = min(r.min() for r in all_r)
    global_rmax = max(r.max() for r in all_r)

    kw = dict(elev=elev, azim=azim, global_range=global_range,
              global_rmin=global_rmin, global_rmax=global_rmax)

    # Multi-panel composite
    n_panels = len(reconstructions)
    fig = plt.figure(figsize=(4 * n_panels, 4.5))
    for i, (Vr, label) in enumerate(zip(reconstructions, labels)):
        ax = fig.add_subplot(1, n_panels, i + 1, projection='3d')
        render_mesh(ax, Vr, F, label, **kw)
    fig.suptitle(f"Sphere-to-{mesh_name.capitalize()} via Spherical Harmonic Expansion",
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(results_dir, f'sphere_to_{mesh_name}_transition.png')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {path}")

    # Individual renders
    for i, (Vr, label) in enumerate(zip(reconstructions, labels)):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        render_mesh(ax, Vr, F, label, **kw)
        plt.tight_layout()
        fname = f'{mesh_name}_lmax_{l_max_values[i]:03d}.png' if i < len(l_max_values) else f'{mesh_name}_lmax_full.png'
        path = os.path.join(results_dir, fname)
        fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  Saved: {path}")

    print(f"Done rendering {mesh_name}!")


def main():
    parser = argparse.ArgumentParser(description="Render sphere-to-mesh SH transition")
    parser.add_argument('mesh', nargs='?', default='cow',
                        help='Mesh name (cow, bunny) or path to .off file')
    parser.add_argument('--lmax', type=int, default=20,
                        help='Maximum SH order (default: 20)')
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
            print(f"Error: not found: {mesh_path}"); sys.exit(1)

    render_transition(mesh_path, mesh_name, results_dir,
                      l_max_full=args.lmax, elev=args.elev, azim=args.azim)


if __name__ == '__main__':
    main()
