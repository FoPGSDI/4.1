"""
Render duck with scaled deformation from sphere.

Correct approach: decompose the RADIAL function r(theta, phi) into SH,
then scale only the non-monopole (ell >= 1) part. This interpolates
between a perfect sphere (scale=0) and the full duck (scale=1).

The vertex positions are then reconstructed as:
    x_i = r_scaled(theta_i, phi_i) * n_hat_i
where n_hat_i is the unit direction of each vertex.
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
from render_sphere_to_cow import (
    cartesian_to_spherical, build_sh_matrix, render_mesh
)


def fit_radial_sh(r, Y_real, real_index, l_max, alpha=1e-4):
    """Fit SH coefficients for the radial function r(theta, phi).

    Returns: coeffs (n_real,) array of real SH coefficients for the scalar r.
    """
    n_real = Y_real.shape[1]
    reg_weights = np.zeros(n_real)
    for j, (ell, m, kind) in enumerate(real_index):
        reg_weights[j] = alpha * ell * (ell + 1)

    D = sparse.diags(np.sqrt(reg_weights))
    result = lsqr(
        sparse.vstack([sparse.csc_matrix(Y_real), D]),
        np.concatenate([r, np.zeros(n_real)]),
    )
    return result[0]


def render_scaled_deformation(mesh_path, mesh_name, results_dir,
                               scale=0.1, l_max_full=20,
                               l_max_values=None, elev=20, azim=-45):
    """Render sphere-to-mesh transition with scaled radial deformation."""
    os.makedirs(results_dir, exist_ok=True)
    if l_max_values is None:
        l_max_values = [0, 1, 2, 4, 8, 16, 32, 64, 128]

    V, F = load_off(mesh_path)
    com = compute_center_of_mass(V, F)
    V_centered = V - com
    print(f"Mesh '{mesh_name}': {len(V)} vertices, {len(F)} faces")

    # Compute radial distance and unit directions
    r_orig, theta, phi = cartesian_to_spherical(V_centered)
    n_hat = V_centered / np.maximum(r_orig, 1e-15)[:, None]  # unit direction vectors

    # Build SH basis for the radial function
    print(f"  Building SH basis matrix up to l_max = {l_max_full}...")
    Y_real, real_index = build_sh_matrix(theta, phi, l_max_full)
    n_coeffs = Y_real.shape[1]
    print(f"  Basis size: {n_coeffs} real SH functions")

    # Fit SH coefficients of the RADIAL function r(theta, phi)
    print(f"  Fitting radial SH coefficients...")
    r_coeffs = fit_radial_sh(r_orig, Y_real, real_index, l_max_full, alpha=1e-4)

    # Separate monopole (mean radius) from deformation
    r_monopole = np.zeros_like(r_coeffs)
    r_deformation = np.zeros_like(r_coeffs)
    for j, (ell, m, kind) in enumerate(real_index):
        if ell == 0:
            r_monopole[j] = r_coeffs[j]
        else:
            r_deformation[j] = r_coeffs[j]

    R0 = r_monopole[0] / Y_real[0, 0]  # mean radius (Y_00 = 1/sqrt(4pi))
    print(f"  Mean radius R_0 = {R0:.6f}")
    print(f"  Deformation RMS (original): {np.sqrt(np.mean((Y_real @ r_deformation)**2)):.6f}")
    print(f"  Deformation RMS (scaled):   {np.sqrt(np.mean((scale * Y_real @ r_deformation)**2)):.6f}")

    # Scaled radial coefficients: r_scaled = r_monopole + scale * r_deformation
    r_coeffs_scaled = r_monopole + scale * r_deformation

    def reconstruct_vertices(r_coeffs_subset, mask=None):
        """Reconstruct 3D vertices from radial SH coefficients."""
        if mask is not None:
            r_recon = Y_real[:, mask] @ r_coeffs_subset[mask]
        else:
            r_recon = Y_real @ r_coeffs_subset
        # Cartesian = r * n_hat
        return r_recon[:, None] * n_hat

    # Full scaled mesh
    V_full_scaled = reconstruct_vertices(r_coeffs_scaled)

    # Reconstruct at each l_max
    actual_lmax_values = [min(lm, l_max_full) for lm in l_max_values]
    # Deduplicate
    seen = set()
    deduped_lmax = []
    deduped_labels = []
    for lm, orig_lm in zip(actual_lmax_values, l_max_values):
        if lm not in seen:
            seen.add(lm)
            deduped_lmax.append(lm)
            deduped_labels.append(orig_lm)
    actual_lmax_values = deduped_lmax

    labels = [f"$\\ell_{{\\max}} = {l}$" for l in deduped_labels] + [f"Full ({scale:.0%} dev.)"]
    reconstructions = []
    for lm in actual_lmax_values:
        mask = np.array([ell <= lm for ell, m, kind in real_index])
        V_recon = reconstruct_vertices(r_coeffs_scaled, mask)
        rmse = np.sqrt(np.mean((V_recon - V_full_scaled)**2))
        reconstructions.append(V_recon)
        print(f"    l_max={lm:3d}: RMSE = {rmse:.6f}")
    reconstructions.append(V_full_scaled)

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
    fig.suptitle(f"Sphere-to-Duck ({scale:.0%} deformation) via Radial SH Expansion",
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(results_dir, f'sphere_to_duck_{scale:.0%}_transition.png')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {path}")

    # Individual renders
    for i, (Vr, label) in enumerate(zip(reconstructions, labels)):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        render_mesh(ax, Vr, F, label, **kw)
        plt.tight_layout()
        if i < len(actual_lmax_values):
            fname = f'duck_lmax_{deduped_labels[i]:03d}.png'
        else:
            fname = f'duck_lmax_full.png'
        path = os.path.join(results_dir, fname)
        fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  Saved: {path}")

    print(f"Done rendering {mesh_name} at {scale:.0%} deformation!")


if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.dirname(__file__))
    mesh_path = os.path.join(base_dir, 'data', 'duck.off')

    parser = argparse.ArgumentParser()
    parser.add_argument('--scale', type=float, default=0.1)
    parser.add_argument('--outdir', type=str, default=None)
    args = parser.parse_args()

    scale = args.scale
    if args.outdir:
        results_dir = args.outdir
    else:
        results_dir = os.path.join(base_dir, 'results', f'duck{scale}')

    render_scaled_deformation(
        mesh_path, 'duck',
        results_dir,
        scale=scale,
        l_max_full=20,
        l_max_values=[0, 1, 2, 4, 8, 16, 32, 64, 128],
        elev=20, azim=-45,
    )
