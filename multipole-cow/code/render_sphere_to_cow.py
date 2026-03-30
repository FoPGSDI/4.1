"""
Render the gradual transition from a sphere to any mesh
at increasing multipole (spherical harmonic) order.

Usage:
    python render_sphere_to_cow.py              # default: cow
    python render_sphere_to_cow.py cow
    python render_sphere_to_cow.py bunny
    python render_sphere_to_cow.py /path/to/mesh.off
"""

import sys
import os
import argparse
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
            Ylm = sph_harm_y(ell, m, theta, phi)
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
    """Render a triangulated surface on a 3D axis."""
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

    ax.set_title(title, fontsize=12, fontweight='bold', pad=2)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([1, 1, 1])


def render_transition(mesh_path, mesh_name, results_dir, l_max_full=30,
                      l_max_values=None, elev=25, azim=-60):
    """Render sphere-to-mesh transition for any mesh.

    Args:
        mesh_path: path to .off file
        mesh_name: short name (e.g. 'cow', 'bunny') for file naming
        results_dir: output directory for images
        l_max_full: maximum SH order for coefficient computation
        l_max_values: list of ℓ_max values for reconstruction panels
        elev, azim: viewing angles
    """
    os.makedirs(results_dir, exist_ok=True)

    if l_max_values is None:
        l_max_values = [0, 1, 2, 4, 8, 16]

    # Load and center mesh
    V, F = load_off(mesh_path)
    com = compute_center_of_mass(V, F)
    V_centered = V - com

    print(f"Mesh '{mesh_name}': {len(V)} vertices, {len(F)} faces")
    print(f"  Center of mass: {com}")

    # Spherical coordinates
    r, theta, phi = cartesian_to_spherical(V_centered)
    print(f"  Radial range: [{r.min():.4f}, {r.max():.4f}]")

    # Solid-angle weights on the unit sphere
    V_unit = V_centered / np.linalg.norm(V_centered, axis=1, keepdims=True)
    areas = compute_voronoi_areas(V_unit, F)
    areas = areas * (4 * np.pi / areas.sum())
    print(f"  Solid angle weights normalized to 4*pi")

    # Compute SH coefficients
    print(f"  Computing SH coefficients up to l_max = {l_max_full}...")
    coeffs = compute_sh_coefficients(r, theta, phi, areas, l_max_full)
    print(f"  Number of coefficients: {len(coeffs)}")

    # Reconstruct at each ℓ_max
    labels = [f"$\\ell_{{\\max}} = {l}$" for l in l_max_values] + ["Full mesh"]
    reconstructions = []
    for lm in l_max_values:
        r_recon = reconstruct_radii(coeffs, theta, phi, lm)
        r_recon = np.maximum(r_recon, 0.01 * r.mean())
        V_recon = spherical_to_cartesian(r_recon, theta, phi)
        reconstructions.append(V_recon)
        print(f"    l_max={lm:3d}: r range [{r_recon.min():.4f}, {r_recon.max():.4f}]")
    reconstructions.append(V_centered)

    # Global bounds
    global_range = max(np.abs(V_rec).max() for V_rec in reconstructions) * 1.1
    all_radii = [np.linalg.norm(V_rec, axis=1) for V_rec in reconstructions]
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

    title = f"Sphere-to-{mesh_name.capitalize()} Transition via Spherical Harmonic Expansion"
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
        if i < len(l_max_values):
            fname = f'{mesh_name}_lmax_{l_max_values[i]:03d}.png'
        else:
            fname = f'{mesh_name}_lmax_full.png'
        out_path = os.path.join(results_dir, fname)
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  Saved: {out_path}")

    print(f"Done rendering {mesh_name}!")


def main():
    parser = argparse.ArgumentParser(description="Render sphere-to-mesh SH transition")
    parser.add_argument('mesh', nargs='?', default='cow',
                        help='Mesh name (cow, bunny) or path to .off file')
    parser.add_argument('--lmax', type=int, default=30,
                        help='Maximum SH order for coefficients (default: 30)')
    parser.add_argument('--elev', type=float, default=25, help='Elevation angle')
    parser.add_argument('--azim', type=float, default=-60, help='Azimuth angle')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.dirname(__file__))
    results_dir = os.path.join(base_dir, 'results')
    data_dir = os.path.join(base_dir, 'data')

    # Resolve mesh path and name
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
                      l_max_full=args.lmax, elev=args.elev, azim=args.azim)


if __name__ == '__main__':
    main()
