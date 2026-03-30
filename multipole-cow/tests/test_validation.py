"""
Validation test suite for "Higher multipoles of the cow" reproduction
(Lehmann 2025, arXiv:2504.00506).

Tests are grouped by functionality:
  - TestMeshIO: loading cow.off, geometry sanity checks
  - TestMultipoleMoments: spherical and Cartesian multipole moments
  - TestGWRadiation: gravitational-wave radiation quantities
  - TestSphHarmBasics: standalone spherical-harmonic identities (no code deps)
"""

import sys
import numpy as np
import pytest
from scipy.special import sph_harm_y

# Allow imports from the code/ directory (modules may not exist yet)
sys.path.insert(0, '/data/haiyangw/claude/4.1/multipole-cow/code')


# ---------------------------------------------------------------------------
# Mesh I/O and geometry
# ---------------------------------------------------------------------------
class TestMeshIO:
    """Tests that the cow mesh loads correctly and has expected dimensions."""

    def test_load_off(self):
        """cow.off loads and contains 2762 vertices and 5520 faces."""
        from mesh_io import load_off

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        assert verts.shape == (2762, 3), f"Expected 2762 vertices, got {verts.shape[0]}"
        assert faces.shape == (5520, 3), f"Expected 5520 faces, got {faces.shape[0]}"

    def test_bounding_box(self):
        """Bounding-box extents within 1 % of (1.044, 0.6397, 0.3403)."""
        from mesh_io import load_off

        verts, _ = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        extents = verts.max(axis=0) - verts.min(axis=0)
        expected = np.array([1.044, 0.6397, 0.3403])
        np.testing.assert_allclose(extents, expected, rtol=0.01)

    def test_volume_positive(self):
        """Mesh volume must be strictly positive."""
        from mesh_io import load_off
        from multipole import mesh_volume

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        vol = mesh_volume(verts, faces)
        assert vol > 0, f"Volume should be positive, got {vol}"

    def test_volume_from_monopole(self):
        r"""V = Q_0^0 * sqrt(4 pi) ~ 0.191, within 5 %."""
        from mesh_io import load_off
        from multipole import compute_Qlm

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        Q00 = compute_Qlm(verts, faces, l=0, m=0)
        volume = Q00.real * np.sqrt(4 * np.pi)
        np.testing.assert_allclose(volume, 0.191, rtol=0.05)


# ---------------------------------------------------------------------------
# Multipole moments
# ---------------------------------------------------------------------------
class TestMultipoleMoments:
    """Tests for spherical and Cartesian multipole moments of the cow."""

    def test_monopole(self):
        r"""Q_0^0 ~ 0.0539, within 5 %."""
        from mesh_io import load_off
        from multipole import compute_Qlm

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        Q00 = compute_Qlm(verts, faces, l=0, m=0)
        np.testing.assert_allclose(Q00.real, 0.0539, rtol=0.05)

    def test_dipole_vanishes(self):
        r"""|Q_1^m| < 1e-10 for m in {-1, 0, 1} (origin at centroid)."""
        from mesh_io import load_off
        from multipole import compute_Qlm

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        for m in [-1, 0, 1]:
            Q1m = compute_Qlm(verts, faces, l=1, m=m)
            assert abs(Q1m) < 1e-10, f"|Q_1^{m}| = {abs(Q1m)} is not negligible"

    def test_quadrupole_Q20(self):
        r"""Q_2^0 ~ -0.0029, within 10 %."""
        from mesh_io import load_off
        from multipole import compute_Qlm

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        Q20 = compute_Qlm(verts, faces, l=2, m=0)
        np.testing.assert_allclose(Q20.real, -0.0029, rtol=0.10)

    def test_Qlm_conjugation(self):
        r"""Q_l^{-m} = (-1)^m conj(Q_l^m) for l = 2, 3."""
        from mesh_io import load_off
        from multipole import compute_Qlm

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        for l in [2, 3]:
            for m in range(1, l + 1):
                Qlm = compute_Qlm(verts, faces, l=l, m=m)
                Ql_neg_m = compute_Qlm(verts, faces, l=l, m=-m)
                expected = ((-1) ** m) * np.conj(Qlm)
                np.testing.assert_allclose(
                    Ql_neg_m, expected, atol=1e-12,
                    err_msg=f"Conjugation failed for l={l}, m={m}",
                )

    def test_cartesian_traceless(self):
        r"""|Tr(Q^C)| < 1e-8 for the traceless Cartesian quadrupole."""
        from mesh_io import load_off
        from multipole import compute_Qlm, cartesian_quadrupole

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        QC = cartesian_quadrupole(verts, faces)
        assert abs(np.trace(QC)) < 1e-8, f"Tr(Q^C) = {np.trace(QC)}"

    def test_cartesian_symmetric(self):
        r"""Q^C_{ij} = Q^C_{ji}."""
        from mesh_io import load_off
        from multipole import cartesian_quadrupole

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        QC = cartesian_quadrupole(verts, faces)
        np.testing.assert_allclose(QC, QC.T, atol=1e-12)

    def test_inertia_symmetric(self):
        r"""I_{ij} = I_{ji}."""
        from mesh_io import load_off
        from multipole import inertia_tensor

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        I = inertia_tensor(verts, faces)
        np.testing.assert_allclose(I, I.T, atol=1e-12)

    def test_QC_inertia_crosscheck(self):
        r"""Q^C_{ij} = Tr(I) delta_{ij} - 3 I_{ij}, within 1 %."""
        from mesh_io import load_off
        from multipole import cartesian_quadrupole, inertia_tensor

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        QC = cartesian_quadrupole(verts, faces)
        I = inertia_tensor(verts, faces)
        expected = np.trace(I) * np.eye(3) - 3 * I
        np.testing.assert_allclose(QC, expected, rtol=0.01)


# ---------------------------------------------------------------------------
# GW radiation
# ---------------------------------------------------------------------------
class TestGWRadiation:
    """Gravitational-wave radiation tests."""

    def test_power_coefficient(self):
        r"""<Q'''_{ij} Q'''^{ij}> / omega^6 ~ 0.00149, within 10 %."""
        from mesh_io import load_off
        from multipole import cartesian_quadrupole
        from gw_radiation import power_coefficient

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        coeff = power_coefficient(verts, faces)
        np.testing.assert_allclose(coeff, 0.00149, rtol=0.10)

    def test_rotation_preserves_trace(self):
        r"""Tr(R Q^C R^T) = Tr(Q^C) for an arbitrary rotation R."""
        from mesh_io import load_off
        from multipole import cartesian_quadrupole

        verts, faces = load_off('/data/haiyangw/claude/4.1/multipole-cow/data/cow.off')
        QC = cartesian_quadrupole(verts, faces)

        # Build an arbitrary rotation matrix via QR decomposition
        rng = np.random.default_rng(42)
        A = rng.standard_normal((3, 3))
        R, _ = np.linalg.qr(A)
        # Ensure proper rotation (det = +1)
        if np.linalg.det(R) < 0:
            R[:, 0] *= -1

        rotated = R @ QC @ R.T
        np.testing.assert_allclose(np.trace(rotated), np.trace(QC), atol=1e-12)


# ---------------------------------------------------------------------------
# Spherical-harmonic basics (no project code dependencies)
# ---------------------------------------------------------------------------
class TestSphHarmBasics:
    """Standalone tests for spherical-harmonic identities using scipy."""

    def test_Y00(self):
        r"""Y_0^0 = 1 / sqrt(4 pi)."""
        val = sph_harm_y(0, 0, 0.5, 1.2)
        np.testing.assert_allclose(val, 1.0 / np.sqrt(4 * np.pi), atol=1e-14)

    def test_orthonormality(self):
        r"""Numerical check of int Y_l^m Y_{l'}^{m'}* d Omega = delta_{ll'} delta_{mm'}."""
        # Use Gauss-Legendre quadrature in theta and uniform grid in phi
        from numpy.polynomial.legendre import leggauss
        n_theta = 80
        n_phi = 160

        nodes, weights = leggauss(n_theta)
        theta = np.arccos(nodes)                  # theta in [0, pi]
        phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
        dphi = 2 * np.pi / n_phi

        pairs = [
            ((0, 0), (0, 0), 1.0),
            ((1, 0), (1, 0), 1.0),
            ((1, 1), (1, 1), 1.0),
            ((2, 1), (2, 1), 1.0),
            ((0, 0), (1, 0), 0.0),
            ((1, 1), (2, 1), 0.0),
            ((2, 0), (2, 1), 0.0),
        ]

        for (l1, m1), (l2, m2), expected in pairs:
            integral = 0.0 + 0.0j
            for i, th in enumerate(theta):
                Y1 = sph_harm_y(l1, m1, th, phi)
                Y2 = sph_harm_y(l2, m2, th, phi)
                integrand_phi = np.sum(Y1 * np.conj(Y2)) * dphi
                integral += integrand_phi * weights[i]  # sin(theta) absorbed in GL weight

            np.testing.assert_allclose(
                integral.real, expected, atol=1e-10,
                err_msg=f"Orthonormality failed for ({l1},{m1})-({l2},{m2})",
            )
            if expected == 0.0:
                assert abs(integral.imag) < 1e-10

    def test_conjugation(self):
        r"""Y_l^{-m} = (-1)^m (Y_l^m)* for several (l, m)."""
        theta_vals = [0.3, 1.0, 2.5]
        phi_vals = [0.0, 0.7, 4.1]

        for l in range(4):
            for m in range(1, l + 1):
                for th in theta_vals:
                    for ph in phi_vals:
                        Y_pos = sph_harm_y(l, m, th, ph)
                        Y_neg = sph_harm_y(l, -m, th, ph)
                        expected = ((-1) ** m) * np.conj(Y_pos)
                        np.testing.assert_allclose(
                            Y_neg, expected, atol=1e-14,
                            err_msg=f"Conjugation failed for l={l}, m={m}, "
                                    f"theta={th}, phi={ph}",
                        )
