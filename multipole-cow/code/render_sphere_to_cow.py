"""
Render the gradual transition from a sphere to the cow mesh
at increasing multipole (spherical harmonic) order.
"""

import sys
import os
import numpy as np
from scipy.special import sph_harm_y
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Add code directory to path
sys.path.insert(0, os.path.dirname(__file__))
from mesh_io import load_off, compute_center_of_mass


def cartesian_to_spherical(xyz):
    """Convert Cartesian (x, y, z) to spherical (r, theta, phi).
    theta: polar angle [0, pi], phi: azimuthal [0, 2pi].
    """
    r = np.linalg.norm(xyz, axis=1)
    theta = np.arccos(np.clip(xyz[:, 2] / np.maximum(r, 1e-15), -1, 1))
    phi = np.arctan2(xyz[:, 1], xyz[:, 0]) % (2 * np.pi)
    return r, theta, phi


def compute_voronoi_areas(V, F):
    """Compute Voronoi area per vertex as 1/3 of the sum of incident face areas."""
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]
    face_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    vertex_areas = np.zeros(len(V))
    for k in range(3):
        np.add.at(vertex_areas, F[:, k], face_areas / 3.0)
    return vertex_areas


def compute_sh_coefficients(r, theta, phi, areas, l_max):
    """Compute spherical harmonic coefficients c_l^m via vertex-area weighted quadrature."""
    coeffs = {}
    for ell in range(l_max + 1):
        for m in range(-ell, ell + 1):
            # Y_l^m(theta, phi) — sph_harm_y returns complex values
            Ylm = sph_harm_y(ell, m, theta, phi)
            # c_l^m = sum_i f_r(theta_i, phi_i) * conj(Y_l^m(theta_i, phi_i)) * A_i
            coeffs[(ell, m)] = np.sum(r * np.conj(Ylm) * areas)
    return coeffs


def reconstruct_radii(coeffs, theta, phi, l_max):
    """Reconstruct radial function from SH coefficients up to l_max."""
    r_recon = np.zeros(len(theta), dtype=complex)
    for ell in range(l_max + 1):
        for m in range(-ell, ell + 1):
            if (ell, m) in coeffs:
                Ylm = sph_harm_y(ell, m, theta, phi)
                r_recon += coeffs[(ell, m)] * Ylm
    return np.real(r_recon)


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical to Cartesian."""
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.column_stack([x, y, z])


def render_mesh(ax, V, F, title, elev=25, azim=-60,
                global_range=None, global_rmin=None, global_rmax=None):
    """Render a triangulated surface on a 3D axis.

    Parameters
    ----------
    global_range : float, optional
        If given, use this as the axis limit for all three axes.
    global_rmin, global_rmax : float, optional
        If given, use these as the colormap normalization bounds.
    """
    # Compute radial distance for coloring
    r = np.linalg.norm(V, axis=1)

    # Face colors: average radius of face vertices
    face_r = (r[F[:, 0]] + r[F[:, 1]] + r[F[:, 2]]) / 3.0

    # Normalize for colormap
    vmin = global_rmin if global_rmin is not None else r.min()
    vmax = global_rmax if global_rmax is not None else r.max()
    if vmax - vmin < 1e-10:
        vmax = vmin + 1.0
    cmap = plt.cm.coolwarm
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    face_colors = cmap(norm(face_r))

    # Build polygon collection
    triangles = V[F]
    poly = Poly3DCollection(triangles, facecolors=face_colors,
                            edgecolors='k', linewidths=0.05, alpha=1.0)
    ax.add_collection3d(poly)

    # Set axis limits
    max_range = global_range if global_range is not None else np.abs(V).max() * 1.1
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)

    ax.set_title(title, fontsize=12, fontweight='bold', pad=2)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])


def main():
    base_dir = os.path.dirname(os.path.dirname(__file__))
    mesh_path = os.path.join(base_dir, 'data', 'cow.off')
    results_dir = os.path.join(base_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    # Load and center mesh
    V, F = load_off(mesh_path)
    com = compute_center_of_mass(V, F)
    V_centered = V - com

    print(f"Mesh: {len(V)} vertices, {len(F)} faces")
    print(f"Center of mass: {com}")

    # Spherical coordinates
    r, theta, phi = cartesian_to_spherical(V_centered)
    print(f"Radial range: [{r.min():.4f}, {r.max():.4f}]")

    # Compute solid angles on the unit sphere for each vertex.
    # Project vertices onto unit sphere, compute face areas there,
    # then use 1/3 of incident face solid angles as vertex weights.
    V_unit = V_centered / np.linalg.norm(V_centered, axis=1, keepdims=True)
    areas = compute_voronoi_areas(V_unit, F)
    # These should sum to ~4*pi for a closed mesh projected onto the sphere
    print(f"Total solid angle (vertex weights): {areas.sum():.4f} (4*pi = {4*np.pi:.4f})")
    # Normalize to exactly 4*pi
    areas = areas * (4 * np.pi / areas.sum())

    # Compute SH coefficients up to high l_max
    l_max_full = 30
    print(f"Computing SH coefficients up to l_max = {l_max_full}...")
    coeffs = compute_sh_coefficients(r, theta, phi, areas, l_max_full)
    print(f"  Number of coefficients: {len(coeffs)}")

    # Reconstruction orders
    l_max_values = [0, 1, 2, 4, 8, 16]
    labels = [f"$\\ell_{{\\max}} = {l}$" for l in l_max_values] + ["Full mesh"]

    # Reconstruct vertex positions for each l_max
    reconstructions = []
    for lm in l_max_values:
        r_recon = reconstruct_radii(coeffs, theta, phi, lm)
        # Clamp to avoid negative radii
        r_recon = np.maximum(r_recon, 0.01 * r.mean())
        V_recon = spherical_to_cartesian(r_recon, theta, phi)
        reconstructions.append(V_recon)
        print(f"  l_max={lm:3d}: r range [{r_recon.min():.4f}, {r_recon.max():.4f}]")

    # Add full mesh (original centered vertices)
    reconstructions.append(V_centered)

    # Compute global axis range and color bounds across all reconstructions
    global_range = max(np.abs(V_rec).max() for V_rec in reconstructions) * 1.1
    all_radii = [np.linalg.norm(V_rec, axis=1) for V_rec in reconstructions]
    global_rmin = min(rad.min() for rad in all_radii)
    global_rmax = max(rad.max() for rad in all_radii)
    print(f"Global axis range: {global_range:.4f}, color range: [{global_rmin:.4f}, {global_rmax:.4f}]")

    render_kw = dict(elev=25, azim=-60, global_range=global_range,
                     global_rmin=global_rmin, global_rmax=global_rmax)

    # --- Multi-panel figure ---
    n_panels = len(reconstructions)
    fig = plt.figure(figsize=(4 * n_panels, 4.5))
    for i, (V_rec, label) in enumerate(zip(reconstructions, labels)):
        ax = fig.add_subplot(1, n_panels, i + 1, projection='3d')
        render_mesh(ax, V_rec, F, label, **render_kw)

    fig.suptitle("Sphere-to-Cow Transition via Spherical Harmonic Expansion",
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out_path = os.path.join(results_dir, 'sphere_to_cow_transition.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {out_path}")

    # --- Individual high-res images ---
    for i, (V_rec, label) in enumerate(zip(reconstructions, labels)):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        lm_val = l_max_values[i] if i < len(l_max_values) else 'full'
        render_mesh(ax, V_rec, F, label, **render_kw)
        plt.tight_layout()
        if lm_val == 'full':
            fname = 'cow_lmax_full.png'
        else:
            fname = f'cow_lmax_{int(lm_val):03d}.png'
        out_path = os.path.join(results_dir, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"Saved: {out_path}")

    print("Done!")


if __name__ == '__main__':
    main()
