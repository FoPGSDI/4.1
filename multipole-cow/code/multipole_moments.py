"""
Multipole moment computation for the cow mesh.

Computes spherical multipole moments Q_l^m, Cartesian quadrupole tensor Q^C,
and the inertia tensor I, using tetrahedral quadrature over a fan
tetrahedralization from the center of mass.
"""

import numpy as np
from scipy.special import sph_harm_y

from mesh_io import tetrahedralize_fan


# =============================================================================
# Tetrahedral quadrature rules on the reference tetrahedron
# Reference tet: vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
# =============================================================================

def tet_quadrature_points(order=4):
    """Return (points, weights) for quadrature on the reference tetrahedron.

    order=4: 4-point rule (exact for degree 2 polynomials)
    order=14: 14-point rule (exact for degree 5 polynomials)

    Points are barycentric coordinates (4 values summing to 1).
    Weights sum to 1/6 (volume of reference tet).
    """
    if order == 4:
        # Hammer-Stroud 4-point rule, degree 2
        a = 0.1381966011250105
        b = 0.5854101966249685
        pts_bary = np.array([
            [a, a, a, b],
            [a, a, b, a],
            [a, b, a, a],
            [b, a, a, a],
        ])
        weights = np.array([1.0/4, 1.0/4, 1.0/4, 1.0/4]) / 6.0
    elif order == 14:
        # Shunn-Ham 14-point rule for tetrahedra, degree 5
        # Using Keast's 14-point rule (degree 4, but good for ell<=5 integrands)
        # Actually, let's use a well-known 14-point degree-5 rule.
        # Keast rule with 14 points:
        a1 = 0.0927352503108
        b1 = 0.7217942490674
        w1 = 0.0734930431163 / 6.0
        a2 = 0.3108859192633
        b2 = 0.0673422422101
        w2 = 0.0425460207771 / 6.0
        a3 = 0.0455397549780
        b3 = 0.4544602450220
        w3 = 0.0836167688816 / 6.0
        # 4 points of type (a1, a1, a1, b1) -- vertex-centered
        pts1 = np.array([
            [a1, a1, a1, b1],
            [a1, a1, b1, a1],
            [a1, b1, a1, a1],
            [b1, a1, a1, a1],
        ])
        w_1 = np.full(4, w1)
        # 4 points of type (a2, a2, a2, b2) -- vertex-centered
        pts2 = np.array([
            [a2, a2, a2, b2],
            [a2, a2, b2, a2],
            [a2, b2, a2, a2],
            [b2, a2, a2, a2],
        ])
        w_2 = np.full(4, w2)
        # 6 points of type (a3, a3, b3, b3) -- edge-centered
        pts3 = np.array([
            [a3, a3, b3, b3],
            [a3, b3, a3, b3],
            [a3, b3, b3, a3],
            [b3, a3, a3, b3],
            [b3, a3, b3, a3],
            [b3, b3, a3, a3],
        ])
        w_3 = np.full(6, w3)

        pts_bary = np.vstack([pts1, pts2, pts3])
        weights = np.concatenate([w_1, w_2, w_3])
    else:
        raise ValueError(f"Unsupported quadrature order: {order}")

    return pts_bary, weights


def integrate_over_tet(tet_verts, integrand, quad_order=4):
    """Integrate a function over a tetrahedron using quadrature.

    Args:
        tet_verts: (4, 3) array of tetrahedron vertices
        integrand: callable(points_array) -> array of complex values
            points_array has shape (n_quad, 3)
        quad_order: 4 or 14

    Returns:
        Complex scalar (integral value)
    """
    bary, weights = tet_quadrature_points(quad_order)

    # Map barycentric coords to physical coords: x = sum_i lambda_i * v_i
    physical_pts = bary @ tet_verts  # (n_quad, 3)

    # Signed volume of the tet = det(v1-v0, v2-v0, v3-v0) / 6
    # We use signed volume because the fan tetrahedralization from COM
    # has overlapping tets whose contributions must cancel via sign.
    d = tet_verts[1:] - tet_verts[0]
    signed_vol = np.linalg.det(d) / 6.0

    # Evaluate integrand
    vals = integrand(physical_pts)

    # Quadrature: weights already include 1/6 factor for reference tet
    # Scale by actual volume / reference volume
    # Reference tet volume = 1/6, so scale = signed_vol / (1/6) = 6*signed_vol
    result = np.sum(weights * vals) * 6.0 * signed_vol

    return result


def _tri_quadrature_points(order=7):
    """Quadrature points and weights for a reference triangle (0,0)-(1,0)-(0,1).

    Returns barycentric coords (n, 3) and weights (n,) summing to 0.5 (triangle area).
    """
    if order <= 3:
        # 4-point degree-3 rule
        pts_bary = np.array([
            [1/3, 1/3, 1/3],
            [0.6, 0.2, 0.2],
            [0.2, 0.6, 0.2],
            [0.2, 0.2, 0.6],
        ])
        weights = np.array([-27/96, 25/96, 25/96, 25/96])
    else:
        # 7-point degree-5 rule (Dunavant)
        a1 = 0.059715871789770
        b1 = 0.470142064105115
        a2 = 0.797426985353087
        b2 = 0.101286507323456
        pts_bary = np.array([
            [1/3, 1/3, 1/3],
            [a1, b1, b1],
            [b1, a1, b1],
            [b1, b1, a1],
            [a2, b2, b2],
            [b2, a2, b2],
            [b2, b2, a2],
        ])
        w0 = 0.1125
        w1 = 0.0661970763942530
        w2 = 0.0629695902724135
        weights = np.array([w0, w1, w1, w1, w2, w2, w2])

    return pts_bary, weights


def compute_Qlm(V, F, com, ell_max=5):
    """Compute spherical multipole moments Q_l^m.

    Q_l^m = integral rho(x) * |x|^l * Y_l^m(theta, phi) d^3x

    with rho = 1 (uniform density), and x measured from the center of mass.

    Uses the divergence theorem to convert to a surface integral:
    Since r^l Y_l^m is a homogeneous polynomial of degree l,
      Q_l^m = 1/(l+3) * integral_S (x . n) * r^l * Y_l^m dA

    where x is measured from COM and n is the outward face normal.
    """
    # Triangle quadrature for surface integration
    tri_bary, tri_weights = _tri_quadrature_points(order=7)
    n_quad = len(tri_weights)

    # Precompute face data
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]

    # Face normals (unnormalized, magnitude = 2 * face area)
    cross = np.cross(v1 - v0, v2 - v0)  # (n_faces, 3)
    face_areas = 0.5 * np.linalg.norm(cross, axis=1)  # (n_faces,)
    face_normals = cross / (2 * face_areas[:, None] + 1e-30)  # unit normals

    # Quadrature points on each face
    quad_pts = np.zeros((len(F), n_quad, 3))
    for q in range(n_quad):
        quad_pts[:, q, :] = (tri_bary[q, 0] * v0
                            + tri_bary[q, 1] * v1
                            + tri_bary[q, 2] * v2)

    # Displacement from COM
    dx = quad_pts - com  # (n_faces, n_quad, 3)

    # x . n for each face and quad point
    x_dot_n = np.sum(dx * face_normals[:, None, :], axis=2)  # (n_faces, n_quad)

    # Precompute r, theta, phi for all quad points
    # Standard spherical coordinates: z is polar axis
    # theta = arccos(z/r), phi = arctan2(y, x)
    # Cow coords: x=forward, y=up, z=right
    r = np.sqrt(np.sum(dx**2, axis=2))  # (n_faces, n_quad)
    r_safe = np.where(r > 0, r, 1.0)
    theta = np.arccos(np.clip(dx[:, :, 2] / r_safe, -1, 1))  # z is polar axis
    phi = np.arctan2(dx[:, :, 1], dx[:, :, 0])  # phi from x toward y

    Qlm = {}

    for ell in range(ell_max + 1):
        r_ell = r ** ell  # (n_faces, n_quad)

        for m in range(0, ell + 1):
            Ylm = sph_harm_y(ell, m, theta, phi)  # (n_faces, n_quad)

            # Integrand: (1/(l+3)) * (x.n) * r^l * Y_l^m
            integrand = x_dot_n * r_ell * Ylm / (ell + 3)  # (n_faces, n_quad)

            # Quadrature over each face, then sum over faces
            # integral_face f dA = 2*A * sum_q w_q * f(x_q)
            # (factor 2*A comes from Jacobian of reference triangle to physical triangle)
            face_integrals = 2.0 * face_areas * np.sum(
                tri_weights[None, :] * integrand, axis=1)  # (n_faces,)

            Qlm[(ell, m)] = np.sum(face_integrals)

        # Negative m via symmetry: Q_l^{-m} = (-1)^m * conj(Q_l^m)
        for m in range(1, ell + 1):
            Qlm[(ell, -m)] = (-1)**m * np.conj(Qlm[(ell, m)])

    return Qlm


def compute_cartesian_quadrupole(V, F, com):
    """Compute the Cartesian (traceless) quadrupole tensor.

    Q^C_ij = integral rho(x) (3 x_i x_j - r^2 delta_ij) d^3x

    with rho = 1, x relative to center of mass.
    """
    tets, tet_verts = tetrahedralize_fan(V, F, com)
    quad_order = 4

    QC = np.zeros((3, 3))

    for i_idx in range(3):
        for j_idx in range(i_idx, 3):
            total = 0.0
            for t in range(len(tets)):
                tv = tet_verts[tets[t]]

                def integrand(pts, _i=i_idx, _j=j_idx):
                    dx = pts - com
                    r2 = np.sum(dx**2, axis=1)
                    val = 3.0 * dx[:, _i] * dx[:, _j]
                    if _i == _j:
                        val -= r2
                    return val

                total += integrate_over_tet(tv, integrand, quad_order).real

            QC[i_idx, j_idx] = total
            QC[j_idx, i_idx] = total

    return QC


def compute_inertia_tensor(V, F, com):
    """Compute the inertia tensor.

    I_ij = integral rho(x) (r^2 delta_ij - x_i x_j) d^3x

    with rho = 1, x relative to center of mass.
    """
    tets, tet_verts = tetrahedralize_fan(V, F, com)
    quad_order = 4

    I = np.zeros((3, 3))

    for i_idx in range(3):
        for j_idx in range(i_idx, 3):
            total = 0.0
            for t in range(len(tets)):
                tv = tet_verts[tets[t]]

                def integrand(pts, _i=i_idx, _j=j_idx):
                    dx = pts - com
                    r2 = np.sum(dx**2, axis=1)
                    val = -dx[:, _i] * dx[:, _j]
                    if _i == _j:
                        val += r2
                    return val

                total += integrate_over_tet(tv, integrand, quad_order).real

            I[i_idx, j_idx] = total
            I[j_idx, i_idx] = total

    return I
