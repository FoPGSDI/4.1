"""
Surface mapping methods for the cow mesh.

Implements the distance gradient flow and harmonic map methods for
mapping the cow surface to the sphere, plus SH decomposition of the
resulting map components (f_r, f_dtheta, f_dphi).

[FUTURE] - Full implementation deferred to Stage 2 for the gradient flow method.
The harmonic map framework and SH decomposition are provided.
"""

import numpy as np
from scipy.special import sph_harm_y
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import spsolve


def compute_sdf_trimesh(V, F, points):
    """Compute signed distance field using trimesh.

    Args:
        V: (N, 3) vertices
        F: (F, 3) faces
        points: (M, 3) query points

    Returns:
        sdf: (M,) signed distance values (negative inside)
    """
    import trimesh
    mesh = trimesh.Trimesh(vertices=V, faces=F, process=False)

    # Closest points on surface
    closest, distance, _ = trimesh.proximity.closest_point(mesh, points)

    # Determine sign: inside (negative) or outside (positive)
    # Use winding number or ray casting
    contains = mesh.contains(points)
    sign = np.where(contains, -1.0, 1.0)

    return sign * distance


def gradient_flow_rk4(V, F, surface_pts, t_max=2.0, dt=0.01):
    """Flow surface points along the gradient of the SDF using RK4.

    [FUTURE] - Full implementation requires careful SDF gradient computation.
    This is a placeholder that demonstrates the concept.

    Args:
        V, F: mesh vertices and faces
        surface_pts: (M, 3) starting points on the surface
        t_max: maximum flow time
        dt: time step

    Returns:
        flowed_pts: (M, 3) points after flow
    """
    import trimesh
    mesh = trimesh.Trimesh(vertices=V, faces=F, process=False)

    pts = surface_pts.copy()
    n_steps = int(t_max / dt)
    eps = 1e-4

    for _ in range(n_steps):
        # Numerical gradient via finite differences
        grad = np.zeros_like(pts)
        for ax in range(3):
            pts_plus = pts.copy()
            pts_minus = pts.copy()
            pts_plus[:, ax] += eps
            pts_minus[:, ax] -= eps
            d_plus = compute_sdf_trimesh(V, F, pts_plus)
            d_minus = compute_sdf_trimesh(V, F, pts_minus)
            grad[:, ax] = (d_plus - d_minus) / (2 * eps)

        # Normalize gradient
        gnorm = np.linalg.norm(grad, axis=1, keepdims=True)
        gnorm = np.maximum(gnorm, 1e-10)
        grad = grad / gnorm

        pts = pts + dt * grad

    return pts


def cotangent_laplacian(V, F):
    """Compute the cotangent Laplacian matrix for a triangle mesh.

    Returns:
        L: (N, N) sparse matrix (cotangent Laplacian)
        M: (N,) array of vertex areas (Voronoi/barycentric)
    """
    n = len(V)
    L = lil_matrix((n, n), dtype=np.float64)
    areas = np.zeros(n)

    for face in F:
        i, j, k = face
        vi, vj, vk = V[i], V[j], V[k]

        # Edge vectors
        eij = vj - vi
        eik = vk - vi
        ejk = vk - vj

        # Cotangents of angles at each vertex
        # Angle at i: between edges ij and ik
        cross_i = np.cross(eij, eik)
        area = 0.5 * np.linalg.norm(cross_i)
        if area < 1e-15:
            continue

        cot_i = np.dot(eij, eik) / (2 * area)
        cot_j = np.dot(-eij, ejk) / (2 * area)
        cot_k = np.dot(-eik, -ejk) / (2 * area)

        # Add cotangent weights
        L[j, k] += 0.5 * cot_i
        L[k, j] += 0.5 * cot_i
        L[i, k] += 0.5 * cot_j
        L[k, i] += 0.5 * cot_j
        L[i, j] += 0.5 * cot_k
        L[j, i] += 0.5 * cot_k

        # Barycentric area contribution
        a3 = area / 3.0
        areas[i] += a3
        areas[j] += a3
        areas[k] += a3

    # Set diagonal: L_ii = -sum_j L_ij
    L = csc_matrix(L)
    diag = -np.array(L.sum(axis=1)).flatten()
    from scipy.sparse import diags
    L = L + diags(diag, 0, shape=(n, n), format='csc')

    return L, areas


def harmonic_map_to_disk(V, F, delete_idx):
    """Compute a harmonic map from the mesh (with one vertex deleted) to the disk.

    [FUTURE] - This would use libigl's harmonic mapping.
    Placeholder for Stage 2 implementation.

    Args:
        V: (N, 3) vertices
        F: (F, 3) faces
        delete_idx: index of vertex to delete

    Returns:
        uv: (N, 2) UV coordinates (deleted vertex mapped to boundary)
    """
    # This is a placeholder - full implementation would use igl.harmonic
    print("[FUTURE] Harmonic map computation deferred to Stage 2")
    return None


def stereographic_inverse(uv):
    """Map points from the disk to the sphere via inverse stereographic projection.

    Projects from the plane (u, v) to the sphere, with the north pole
    as the projection point.

    Args:
        uv: (M, 2) array of UV coordinates

    Returns:
        sphere_pts: (M, 3) points on the unit sphere
    """
    u, v = uv[:, 0], uv[:, 1]
    r2 = u**2 + v**2
    denom = 1 + r2
    x = 2 * u / denom
    y = 2 * v / denom
    z = (r2 - 1) / denom  # South pole at (0,0,-1) for r2->inf
    return np.column_stack([x, y, z])


def extract_surface_components(cow_pts, sphere_coords, com):
    """Extract the (f_r, f_dtheta, f_dphi) components of the surface map.

    Args:
        cow_pts: (M, 3) points on the cow surface
        sphere_coords: (M, 2) (theta, phi) coordinates on the sphere
        com: (3,) center of mass

    Returns:
        f_r: (M,) radial distance from CoM
        f_dth: (M,) angular deviation in theta
        f_dph: (M,) angular deviation in phi
    """
    dx = cow_pts - com
    r = np.linalg.norm(dx, axis=1)
    # Spherical coords of cow points (y-up convention)
    r_safe = np.where(r > 0, r, 1.0)
    theta_cow = np.arccos(np.clip(dx[:, 1] / r_safe, -1, 1))
    phi_cow = np.arctan2(dx[:, 2], dx[:, 0])

    theta_sphere = sphere_coords[:, 0]
    phi_sphere = sphere_coords[:, 1]

    f_r = r
    f_dth = theta_cow - theta_sphere
    f_dph = phi_cow - phi_sphere

    # Wrap phi difference to [-pi, pi]
    f_dph = (f_dph + np.pi) % (2 * np.pi) - np.pi

    return f_r, f_dth, f_dph


def decompose_sh(f_vals, theta, phi, areas, ell_max=5):
    """Decompose a function on the sphere into spherical harmonics.

    f_lm = integral f(Omega) Y_l^m*(Omega) dOmega
         ~ sum_i f(Omega_i) Y_l^m*(Omega_i) * area_i

    Args:
        f_vals: (M,) function values at sample points
        theta: (M,) polar angles (colatitude)
        phi: (M,) azimuthal angles
        areas: (M,) area elements for each sample point
        ell_max: maximum multipole order

    Returns:
        dict[(l, m)] -> complex coefficient
    """
    coeffs = {}
    for ell in range(ell_max + 1):
        for m in range(-ell, ell + 1):
            Ylm = sph_harm_y(ell, m, theta, phi)
            # Inner product: <f | Ylm> = integral f * Ylm^* dOmega
            coeffs[(ell, m)] = np.sum(f_vals * np.conj(Ylm) * areas)
    return coeffs
