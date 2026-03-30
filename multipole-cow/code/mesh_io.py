"""
Mesh I/O and basic geometric computations for the cow mesh.
"""

import numpy as np


def load_off(path):
    """Load an OFF mesh file.

    Returns:
        vertices: (N, 3) float array
        faces: (F, 3) int array
    """
    with open(path, 'r') as f:
        line = f.readline().strip()
        if line != 'OFF':
            raise ValueError(f"Not an OFF file: first line is '{line}'")
        # Read counts
        while True:
            line = f.readline().strip()
            if line and not line.startswith('#'):
                break
        parts = line.split()
        n_verts, n_faces, _ = int(parts[0]), int(parts[1]), int(parts[2])

        # Read vertices
        vertices = np.zeros((n_verts, 3), dtype=np.float64)
        for i in range(n_verts):
            line = f.readline().strip()
            vertices[i] = [float(x) for x in line.split()[:3]]

        # Read faces
        faces = np.zeros((n_faces, 3), dtype=np.int64)
        for i in range(n_faces):
            line = f.readline().strip()
            parts = line.split()
            # First number is vertex count per face (should be 3)
            nv = int(parts[0])
            if nv != 3:
                raise ValueError(f"Non-triangular face with {nv} vertices")
            faces[i] = [int(parts[j+1]) for j in range(3)]

    return vertices, faces


def compute_volume(V, F):
    """Compute volume of a closed triangular mesh using the divergence theorem.

    V = (1/6) * sum_faces v0 . (v1 x v2)
    """
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]
    # Signed volume of tetrahedra formed with origin
    vol = np.sum(v0 * np.cross(v1, v2)) / 6.0
    return abs(vol)


def compute_center_of_mass(V, F):
    """Compute center of mass of a uniform-density closed mesh.

    Uses the divergence theorem: each face contributes a tetrahedron with the origin.
    CoM_k = (1/(2*V)) * sum_faces integral of x_k over tet.

    For a tet with origin, the integral of x_k is:
    (signed_vol / 4) * (v0_k + v1_k + v2_k)... but more precisely we use
    the formula for centroid of uniform solid bounded by triangular mesh.
    """
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]

    # Signed volumes of tetrahedra with origin
    cross = np.cross(v1, v2)
    signed_vols = np.sum(v0 * cross, axis=1) / 6.0

    total_vol = np.sum(signed_vols)

    # Centroid: for each tet (origin, v0, v1, v2), centroid is (v0+v1+v2)/4
    # Weighted sum: CoM = sum(signed_vol_i * centroid_i) / total_vol
    centroids = (v0 + v1 + v2) / 4.0
    com = np.sum(signed_vols[:, None] * centroids, axis=0) / total_vol

    return com


def compute_bounding_box(V):
    """Compute axis-aligned bounding box.

    Returns:
        (min_corner, max_corner): each (3,) arrays
    """
    return V.min(axis=0), V.max(axis=0)


def tetrahedralize_fan(V, F, com):
    """Create a fan tetrahedralization from the center of mass.

    Each surface triangle, together with the CoM point, forms a tetrahedron.

    Returns:
        tets: (F_count, 4) int array - indices into tet_vertices
        tet_vertices: (N+1, 3) array - original vertices plus CoM appended
    """
    com_idx = len(V)
    tet_vertices = np.vstack([V, com.reshape(1, 3)])

    # Each face (i, j, k) -> tet (com_idx, i, j, k)
    n_faces = len(F)
    tets = np.zeros((n_faces, 4), dtype=np.int64)
    tets[:, 0] = com_idx
    tets[:, 1:] = F

    return tets, tet_vertices
