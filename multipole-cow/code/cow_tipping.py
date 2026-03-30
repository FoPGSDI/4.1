"""
Cow tipping analysis.

Determines the minimum force required to tip the cow, following the
analysis in Lehmann (2025). The cow rotates about the x-axis, with force
applied in the yz-plane.
"""

import numpy as np


def find_ground_contact(V, F, tol=1e-4):
    """Find vertices in contact with the ground.

    The ground plane is at the minimum y-coordinate. Vertices within
    `tol` of this minimum are considered in contact.

    Args:
        V: (N, 3) vertices
        F: (F, 3) faces
        tol: tolerance for ground contact

    Returns:
        contact_idx: array of vertex indices in contact with ground
    """
    y_min = V[:, 1].min()
    contact_idx = np.where(V[:, 1] < y_min + tol)[0]
    return contact_idx


def find_pivot_and_force_point(V, F, contact_idx, com):
    """Find the pivot point and force application point for tipping.

    The pivot is the ground contact point with the largest z-coordinate
    (furthest to the cow's right), as the tipping force pushes in +z.

    The force is applied at the intersection of the cow surface with the
    ray from the pivot through the barycenter.

    Args:
        V: (N, 3) vertices
        F: (F, 3) faces
        contact_idx: indices of ground contact vertices
        com: (3,) center of mass

    Returns:
        pivot: (3,) pivot point
        p_T: (3,) force application point on cow surface
        force_dir: (3,) unit direction of tipping force
    """
    # Pivot: ground contact with max z
    z_vals = V[contact_idx, 2]
    pivot_local_idx = np.argmax(z_vals)
    pivot = V[contact_idx[pivot_local_idx]].copy()

    # Ray from pivot through barycenter
    q_dir = com - pivot
    q_dir_norm = q_dir / np.linalg.norm(q_dir)

    # Find intersection with cow surface (approximate: find closest vertex
    # along the ray on the far side of the barycenter)
    # We look for the vertex that is closest to the ray, beyond the barycenter
    ray_origin = pivot
    # Project all vertices onto the ray
    dv = V - ray_origin
    t_proj = np.dot(dv, q_dir_norm)
    # Distance from ray
    proj_pts = ray_origin + t_proj[:, None] * q_dir_norm
    dist_from_ray = np.linalg.norm(V - proj_pts, axis=1)

    # Only consider vertices beyond the barycenter (t > |pivot - com|)
    t_com = np.linalg.norm(com - pivot)
    mask = t_proj > t_com * 0.5  # at least halfway
    if np.any(mask):
        # Among those, find closest to ray
        dist_masked = dist_from_ray.copy()
        dist_masked[~mask] = np.inf
        best_idx = np.argmin(dist_masked)
        p_T = V[best_idx].copy()
    else:
        # Fallback: use the point farthest from pivot along the ray direction
        best_idx = np.argmax(t_proj)
        p_T = V[best_idx].copy()

    # Force direction: perpendicular to both the ray q and x-hat
    x_hat = np.array([1.0, 0.0, 0.0])
    force_dir = np.cross(q_dir_norm, x_hat)
    force_norm = np.linalg.norm(force_dir)
    if force_norm > 1e-10:
        force_dir = force_dir / force_norm
    else:
        # Fallback
        force_dir = np.array([0.0, 0.0, 1.0])

    # Ensure force is in +z direction (tipping to the right)
    if force_dir[2] < 0:
        force_dir = -force_dir

    return pivot, p_T, force_dir


def compute_min_tipping_force(pivot, p_T, force_dir, com, weight):
    """Compute the minimum force required to start tipping.

    At the tipping threshold:
        F_T * |p_T - com| = F_N * |(pivot - com) x y_hat|

    where F_N = weight - F_T * (force_dir . z_hat)

    Solving for F_T:
        F_T = weight * |(p - b) x y_hat| / (|p_T - b| + |(p - b) x y_hat| * (F_dir . z_hat) / |force_dir|)

    Wait, more carefully:
    Torque from tipping: tau_T = F_T * |p_T - b|  (where b is barycenter)
    Torque from normal:  tau_N = F_N * |(p - b) x y_hat|
    F_N = w - F_T * (force_dir . y_hat)

    Actually the paper says: F_N = w - F_T . z_hat
    Let me re-read: "F_N = w - F_T . z_hat, where w is the weight"

    Torque balance: F_T * |p_T - b| = F_N * |(p - b) x y_hat|

    But the torques are about the pivot point p. Let me be more careful.

    The tipping torque about the x-axis through the pivot:
        tau_T = F_T * |arm_T|
    where arm_T is the moment arm of the force about the pivot.

    The restoring torque from normal force:
        tau_N = F_N * arm_N
    where arm_N = |(pivot - com) x y_hat| = |z_com - z_pivot|... no.

    Let me think in 2D (yz plane):
    - Gravity acts at CoM: force = -w * y_hat
    - Normal force at pivot: F_N * y_hat
    - Tipping force at p_T: F_T * force_dir

    Torques about pivot (taking x-axis moments):
    - Gravity: w * (z_com - z_pivot) [restoring, if z_com > z_pivot... no]
    Actually: gravity torque about pivot = w * (z_com - z_pivot) in the -x direction
    (tends to tip if z_com > z_pivot)

    Wait, the torque of gravity about the pivot in x-direction:
    tau_gravity = (com - pivot) x (-w * y_hat)
    = (-w) * (com - pivot) x y_hat
    The x-component: (-w) * [(y_com - y_pivot)*0 - (z_com - z_pivot)*1]...

    Let me just compute it with vectors.

    Args:
        pivot: (3,) pivot point on ground
        p_T: (3,) force application point
        force_dir: (3,) unit force direction
        com: (3,) center of mass
        weight: float, total weight (mg)

    Returns:
        F_min: minimum tipping force magnitude
    """
    y_hat = np.array([0.0, 1.0, 0.0])
    x_hat = np.array([1.0, 0.0, 0.0])

    # Moment arm of gravity about pivot (x-component of torque)
    r_com = com - pivot
    tau_gravity_x = np.cross(r_com, -weight * y_hat)[0]  # x-component

    # Moment arm of tipping force about pivot (x-component of torque)
    r_pT = p_T - pivot
    # Torque per unit force
    tau_tip_per_F = np.cross(r_pT, force_dir)[0]  # x-component

    # At threshold, total torque in x = 0:
    # tau_gravity_x + F_T * tau_tip_per_F = 0
    # But we also need to account for the change in normal force:
    # The normal force changes: F_N = weight - F_T * (force_dir . y_hat)
    # But the normal force acts at the pivot, so its torque about the pivot is zero!

    # Actually wait: F_N acts at pivot, so its torque about pivot = 0.
    # So we just need: tau_gravity + tau_tipping = 0

    # tau_gravity_x is the x-component of (com-pivot) x (-w y_hat)
    # This is negative (restoring) if the cow's CoM is above the pivot
    # tau_tip_per_F is the x-component of (p_T - pivot) x force_dir

    # For the cow to tip, we need the tipping torque to overcome gravity restoring torque
    # The restoring torque from gravity (keeping cow upright) has tau_gravity_x < 0
    # (because CoM is to the left of pivot in z, so gravity pulls it back)

    # Actually, let's think about signs:
    # If tipping to the right (+z), and pivot is at max z of contact,
    # then CoM is to the LEFT of pivot (z_com < z_pivot), so gravity restores.
    # The gravity torque about pivot in x:
    #   r_com x (-w y_hat) = (com - pivot) x (-w y_hat)

    # For tipping: we need F_T such that the net torque tips the cow
    if abs(tau_tip_per_F) < 1e-15:
        return np.inf

    F_min = -tau_gravity_x / tau_tip_per_F

    return abs(F_min)


def tipping_analysis(V, F, com, I):
    """Complete cow tipping analysis.

    Args:
        V: (N, 3) vertices
        F: (F, 3) faces
        com: (3,) center of mass
        I: (3,3) inertia tensor in benchmark units

    Returns:
        dict with all tipping results and physical units
    """
    # Physical parameters
    L_m = 2.39  # meters per benchmark unit
    L_cm = L_m * 100
    rho_cgs = 1.0  # g/cm^3
    M_g = rho_cgs * L_cm**3
    M_kg = M_g / 1000.0
    g = 9.81  # m/s^2

    weight_N = M_kg * g  # Newtons

    contact_idx = find_ground_contact(V, F)
    pivot, p_T, force_dir = find_pivot_and_force_point(V, F, contact_idx, com)
    F_min_benchmark = compute_min_tipping_force(pivot, p_T, force_dir, com, 1.0)

    # Convert force to physical units
    # In benchmark units, weight = 1 (dimensionless). F_min is in units of weight.
    F_min_N = F_min_benchmark * weight_N

    results = {
        'contact_vertices': len(contact_idx),
        'pivot': pivot.tolist(),
        'force_point': p_T.tolist(),
        'force_direction': force_dir.tolist(),
        'F_min_benchmark_units': float(F_min_benchmark),
        'F_min_N': float(F_min_N),
        'weight_N': float(weight_N),
        'F_min_over_weight': float(F_min_benchmark),
        'M_kg': float(M_kg),
        'com': com.tolist(),
        'human_sustainable_force_N': 500.0,
        'boxer_punch_N': 5000.0,
        'human_can_tip': F_min_N < 500.0,
        'boxer_can_tip': F_min_N < 5000.0,
    }

    return results
