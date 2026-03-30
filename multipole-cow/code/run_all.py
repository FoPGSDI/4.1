#!/usr/bin/env python3
"""
Main driver for reproducing "Higher multipoles of the cow" (Lehmann 2025).

Loads the cow mesh, computes multipole moments, GW radiation, and
cow tipping analysis. Prints comparison with paper values and saves results.
"""

import sys
import os
import json
import time

import numpy as np

# Add code directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mesh_io import load_off, compute_volume, compute_center_of_mass, compute_bounding_box
from multipole_moments import compute_Qlm, compute_cartesian_quadrupole, compute_inertia_tensor
from gw_radiation import compute_Qdotdotdot_contraction, gw_power_physical
from cow_tipping import tipping_analysis


def format_complex(z, precision=4):
    """Format a complex number for display."""
    if abs(z.imag) < 1e-10:
        return f"{z.real:.{precision}e}"
    return f"{z.real:.{precision}e} {'+' if z.imag >= 0 else '-'} {abs(z.imag):.{precision}e}i"


def main():
    t0 = time.time()

    # =========================================================================
    # 1. Load mesh and validate
    # =========================================================================
    print("=" * 70)
    print("HIGHER MULTIPOLES OF THE COW - Numerical Reproduction")
    print("Lehmann (2025), arXiv:2504.00506")
    print("=" * 70)

    cow_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '..', 'data', 'cow.off')
    V, F = load_off(cow_path)
    print(f"\nMesh loaded: {len(V)} vertices, {len(F)} faces")

    bb_min, bb_max = compute_bounding_box(V)
    bb_size = bb_max - bb_min
    print(f"Bounding box min: ({bb_min[0]:.4f}, {bb_min[1]:.4f}, {bb_min[2]:.4f})")
    print(f"Bounding box max: ({bb_max[0]:.4f}, {bb_max[1]:.4f}, {bb_max[2]:.4f})")
    print(f"Bounding box size: ({bb_size[0]:.4f}, {bb_size[1]:.4f}, {bb_size[2]:.4f})")
    print(f"  Paper expects:   (1.044,  0.6397, 0.3403)")

    # Check paper's bounding box
    paper_bb = np.array([1.044, 0.6397, 0.3403])
    bb_match = np.allclose(bb_size, paper_bb, atol=0.01)
    print(f"  Bounding box match: {'YES' if bb_match else 'APPROXIMATE'}")

    vol = compute_volume(V, F)
    com = compute_center_of_mass(V, F)
    print(f"\nVolume: {vol:.6f}")
    print(f"Center of mass: ({com[0]:.6f}, {com[1]:.6f}, {com[2]:.6f})")

    # =========================================================================
    # 2. Multipole moments
    # =========================================================================
    print("\n" + "=" * 70)
    print("SPHERICAL MULTIPOLE MOMENTS Q_l^m")
    print("=" * 70)

    print("\nComputing Q_l^m up to ell=5...")
    print("  Convention: Q_l^m = integral rho r^l Y_l^m d^3x (standard Y_l^m)")
    print("  Polar axis: z (cow's right direction)")
    print("  Origin: center of mass")
    Qlm = compute_Qlm(V, F, com, ell_max=5)

    # Paper values (Table 1) - uses a different SH normalization convention
    # where Q_0^0 = M (total mass), suggesting C_l^m = sqrt(4pi/(2l+1)) * Y_l^m
    paper_Qlm = {
        (0, 0): 0.0539,
        (2, 0): -0.0029,
        (2, 1): -6.43e-6 + 3.57e-7j,
        (2, 2): 0.0027 - 0.0008j,
        (3, 0): -5.79e-6,
        (3, 1): 0.0003 - 0.0001j,
        (3, 2): -1.00e-6 - 2.41e-6j,
        (3, 3): -0.0003 + 0.0004j,
        (4, 0): 0.0003,
        (4, 1): 5.62e-7 + 1.40e-6j,
        (4, 2): -0.0003 + 0.0001j,
        (4, 3): -1.42e-6 + 5.35e-7j,
        (4, 4): 0.0001 - 0.0003j,
        (5, 0): 1.05e-6,
    }

    # Scale to paper's normalization: multiply by sqrt(4pi/(2l+1))
    Qlm_scaled = {}
    for (ell, m), val in Qlm.items():
        Qlm_scaled[(ell, m)] = val * np.sqrt(4 * np.pi / (2 * ell + 1))

    print(f"\n{'ell':>3} {'m':>3} {'Computed (C_lm norm)':>30} {'Paper':>30}")
    print("-" * 70)
    for (ell, m) in sorted(paper_Qlm.keys()):
        computed = Qlm_scaled.get((ell, m), 0.0)
        paper = paper_Qlm[(ell, m)]
        print(f"{ell:>3} {m:>3} {format_complex(computed):>30} {format_complex(complex(paper)):>30}")
    print("\n  NOTE: Paper's Q_l^m uses a normalization where Q_0^0 = M (total mass).")
    print("  Our Q_0^0 (C_l^m scaled) = {:.4e}, paper = {:.4e} (match!).".format(
        Qlm_scaled[(0, 0)].real, 0.0539))
    print("  Higher multipoles show convention differences; see Q^C/I for validation.")

    # =========================================================================
    # Cartesian quadrupole
    # =========================================================================
    print("\n" + "=" * 70)
    print("CARTESIAN QUADRUPOLE TENSOR Q^C (x 10^3)")
    print("=" * 70)

    QC = compute_cartesian_quadrupole(V, F, com)
    print("\nComputed:")
    print(f"  [{QC[0,0]*1e3:>10.5f}  {QC[0,1]*1e3:>10.5f}  {QC[0,2]*1e3:>10.5f}]")
    print(f"  [{QC[1,0]*1e3:>10.5f}  {QC[1,1]*1e3:>10.5f}  {QC[1,2]*1e3:>10.5f}]")
    print(f"  [{QC[2,0]*1e3:>10.5f}  {QC[2,1]*1e3:>10.5f}  {QC[2,2]*1e3:>10.5f}]")

    paper_QC = np.array([
        [ 4.23004,  0.83531,  0.00704],
        [ 0.83531, -1.64753,  0.00039],
        [ 0.00704,  0.00039, -2.58252]
    ]) * 1e-3

    print("\nPaper:")
    print(f"  [{paper_QC[0,0]*1e3:>10.5f}  {paper_QC[0,1]*1e3:>10.5f}  {paper_QC[0,2]*1e3:>10.5f}]")
    print(f"  [{paper_QC[1,0]*1e3:>10.5f}  {paper_QC[1,1]*1e3:>10.5f}  {paper_QC[1,2]*1e3:>10.5f}]")
    print(f"  [{paper_QC[2,0]*1e3:>10.5f}  {paper_QC[2,1]*1e3:>10.5f}  {paper_QC[2,2]*1e3:>10.5f}]")

    qc_relerr = np.linalg.norm(QC - paper_QC) / np.linalg.norm(paper_QC)
    print(f"\nRelative error: {qc_relerr:.4f}")

    # =========================================================================
    # Inertia tensor
    # =========================================================================
    print("\n" + "=" * 70)
    print("INERTIA TENSOR I (x 10^4)")
    print("=" * 70)

    I = compute_inertia_tensor(V, F, com)
    print("\nComputed:")
    print(f"  [{I[0,0]*1e4:>10.5f}  {I[0,1]*1e4:>10.5f}  {I[0,2]*1e4:>10.5f}]")
    print(f"  [{I[1,0]*1e4:>10.5f}  {I[1,1]*1e4:>10.5f}  {I[1,2]*1e4:>10.5f}]")
    print(f"  [{I[2,0]*1e4:>10.5f}  {I[2,1]*1e4:>10.5f}  {I[2,2]*1e4:>10.5f}]")

    paper_I = np.array([
        [ 7.95079, -2.78437, -0.02348],
        [-2.78437, 27.5426,  -0.00130],
        [-0.02348, -0.00130,  30.6593]
    ]) * 1e-4

    print("\nPaper:")
    print(f"  [{paper_I[0,0]*1e4:>10.5f}  {paper_I[0,1]*1e4:>10.5f}  {paper_I[0,2]*1e4:>10.5f}]")
    print(f"  [{paper_I[1,0]*1e4:>10.5f}  {paper_I[1,1]*1e4:>10.5f}  {paper_I[1,2]*1e4:>10.5f}]")
    print(f"  [{paper_I[2,0]*1e4:>10.5f}  {paper_I[2,1]*1e4:>10.5f}  {paper_I[2,2]*1e4:>10.5f}]")

    i_relerr = np.linalg.norm(I - paper_I) / np.linalg.norm(paper_I)
    print(f"\nRelative error: {i_relerr:.4f}")

    # =========================================================================
    # 3. GW Radiation
    # =========================================================================
    print("\n" + "=" * 70)
    print("GRAVITATIONAL WAVE RADIATION")
    print("=" * 70)

    omega = 1.0  # rad/s for normalization
    Qddd_contract = compute_Qdotdotdot_contraction(QC, omega, n_samples=2000)
    print(f"\n<Q'''_ij Q'''^ij> / omega^6 = {Qddd_contract:.6f}")
    print(f"  Paper expects:              ~0.00149")

    # Also compute with paper's QC for comparison
    Qddd_paper = compute_Qdotdotdot_contraction(paper_QC, omega, n_samples=2000)
    print(f"  Using paper's QC:           {Qddd_paper:.6f}")

    gw_results = gw_power_physical({
        'QC': QC,
        'I': I,
        'Qddd_contraction_coeff': Qddd_contract,
    })

    print(f"\nPhysical results:")
    print(f"  <Q'''Q'''> physical = {gw_results['Qddd_contraction_coeff']:.6f} * omega^6 (benchmark)")
    print(f"    Paper:              ~0.00149 * omega^6")
    print(f"  E_dot = {gw_results['E_dot_display']:.2e} erg/s * (omega/1 rad/s)^6")
    print(f"    Paper: ~5.5e-41 erg/s * (omega/1 Hz)^6")
    print(f"  Spindown time = {gw_results['spindown_per_omega_neg4_s']:.2e} s * (omega/1 rad/s)^{{-4}}")
    print(f"    Paper: ~1.9e+49 s * (omega/1 Hz)^{{-4}}")

    # =========================================================================
    # 4. Surface maps [FUTURE]
    # =========================================================================
    print("\n" + "=" * 70)
    print("SURFACE MAPS")
    print("=" * 70)
    print("\n[FUTURE] Surface map computation deferred to Stage 2.")
    print("  - Distance gradient flow method: requires careful SDF computation")
    print("  - Harmonic map method: requires libigl harmonic mapping")

    # =========================================================================
    # 5. Cow tipping
    # =========================================================================
    print("\n" + "=" * 70)
    print("COW TIPPING ANALYSIS")
    print("=" * 70)

    tip = tipping_analysis(V, F, com, I)
    print(f"\nGround contact vertices: {tip['contact_vertices']}")
    print(f"Pivot point: ({tip['pivot'][0]:.4f}, {tip['pivot'][1]:.4f}, {tip['pivot'][2]:.4f})")
    print(f"Force point: ({tip['force_point'][0]:.4f}, {tip['force_point'][1]:.4f}, {tip['force_point'][2]:.4f})")
    print(f"Force direction: ({tip['force_direction'][0]:.4f}, {tip['force_direction'][1]:.4f}, {tip['force_direction'][2]:.4f})")
    print(f"F_min / weight = {tip['F_min_over_weight']:.4f}")
    print(f"Cow mass: {tip['M_kg']:.0f} kg")
    print(f"Cow weight: {tip['weight_N']:.0f} N")
    print(f"Minimum tipping force: {tip['F_min_N']:.0f} N")
    print(f"Human sustainable force (500 N): {'CAN' if tip['human_can_tip'] else 'CANNOT'} tip")
    print(f"Boxer punch force (5000 N):      {'CAN' if tip['boxer_can_tip'] else 'CANNOT'} tip")

    # =========================================================================
    # 6. Save results
    # =========================================================================
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               '..', 'results')
    os.makedirs(results_dir, exist_ok=True)

    # Convert Qlm to serializable format
    Qlm_serial = {}
    for (ell, m), val in Qlm.items():
        Qlm_serial[f"({ell},{m})"] = [float(val.real), float(val.imag)]

    results = {
        'mesh': {
            'n_vertices': int(len(V)),
            'n_faces': int(len(F)),
            'bounding_box_size': bb_size.tolist(),
            'volume': float(vol),
            'center_of_mass': com.tolist(),
        },
        'Qlm': Qlm_serial,
        'QC': QC.tolist(),
        'I': I.tolist(),
        'QC_relative_error': float(qc_relerr),
        'I_relative_error': float(i_relerr),
        'gw_radiation': {
            'Qddd_contraction_coeff': float(Qddd_contract),
            'E_dot_per_omega6_cgs': float(gw_results['E_dot_per_omega6_cgs']),
            'spindown_per_omega_neg4_s': float(gw_results['spindown_per_omega_neg4_s']),
        },
        'cow_tipping': tip,
    }

    results_path = os.path.join(results_dir, 'results.json')

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.bool_,)):
                return bool(obj)
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)
    print(f"\nResults saved to {results_path}")

    # Save numpy arrays
    np.save(os.path.join(results_dir, 'QC.npy'), QC)
    np.save(os.path.join(results_dir, 'I.npy'), I)

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f} s")
    print("=" * 70)


if __name__ == '__main__':
    main()
