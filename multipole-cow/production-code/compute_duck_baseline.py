"""
Phase 0: Compute all baseline multipole moments and geometric properties
for the duck mesh ("Quacky Normal Modes" paper).

Outputs:
  - production-code/duck_baseline.json
  - production-code/duck_epsilon.dat
  - results/duck_deformation_spectrum.png
  - results/duck_power_spectrum.png
  - production-code/progress/phase0_baseline.md
"""

import sys
import os
import json
import numpy as np

# Add code/ to path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, 'code'))

from mesh_io import load_off, compute_volume, compute_center_of_mass, compute_bounding_box
from multipole_moments import compute_Qlm, compute_cartesian_quadrupole, compute_inertia_tensor
from render_sphere_to_cow import build_sh_matrix, fit_sh_coefficients, cartesian_to_spherical


def main():
    # =========================================================================
    # 1. Load duck mesh and compute basic geometric properties
    # =========================================================================
    duck_path = os.path.join(BASE_DIR, 'data', 'duck.off')
    print(f"Loading duck mesh from {duck_path}...")
    V, F = load_off(duck_path)
    print(f"  Vertices: {len(V)}, Faces: {len(F)}")

    volume = compute_volume(V, F)
    com = compute_center_of_mass(V, F)
    bb_min, bb_max = compute_bounding_box(V)
    bb_size = bb_max - bb_min

    R_eq = (3 * volume / (4 * np.pi)) ** (1.0 / 3.0)

    print(f"  Volume: {volume:.10f}")
    print(f"  Center of mass: {com}")
    print(f"  Bounding box min: {bb_min}")
    print(f"  Bounding box max: {bb_max}")
    print(f"  Bounding box size: {bb_size}")
    print(f"  R_eq = (3V/4pi)^(1/3) = {R_eq:.10f}")

    # =========================================================================
    # 2. Compute multipole moments
    # =========================================================================
    ell_max = 10

    # Cartesian quadrupole Q^C
    print(f"\nComputing Cartesian quadrupole tensor Q^C...")
    QC = compute_cartesian_quadrupole(V, F, com)
    print(f"  Q^C =\n{QC}")
    QC_eigenvalues = np.sort(np.linalg.eigvalsh(QC))[::-1]
    print(f"  Q^C eigenvalues: {QC_eigenvalues}")

    # Inertia tensor I
    print(f"\nComputing inertia tensor I...")
    I = compute_inertia_tensor(V, F, com)
    print(f"  I =\n{I}")
    I_eigenvalues = np.sort(np.linalg.eigvalsh(I))[::-1]
    print(f"  I eigenvalues (principal moments): {I_eigenvalues}")

    # kappa_duck = Q^C_eigenvalues / (M * R_eq^2)
    # M = rho * V, and with rho=1, M = V
    M = volume
    kappa_duck = QC_eigenvalues / (M * R_eq**2)
    print(f"  kappa_duck = Q^C_eig / (M * R_eq^2) = {kappa_duck}")

    # Spherical multipole moments Q_lm
    print(f"\nComputing spherical multipole moments Q_lm up to ell={ell_max}...")
    Qlm = compute_Qlm(V, F, com, ell_max=ell_max)
    print(f"  Q_00 = {Qlm[(0,0)]:.10f}")
    print(f"  Q_10 = {Qlm[(1,0)]:.2e} (should be ~0 if centered at COM)")

    # =========================================================================
    # 3. Extract SH deformation coefficients epsilon_lm
    # =========================================================================
    print(f"\nExtracting SH deformation coefficients epsilon_lm...")
    l_max_sh = 10

    V_centered = V - com
    r_verts, theta, phi = cartesian_to_spherical(V_centered)

    # Build SH basis and fit radial function r(theta, phi)
    Y_real, real_index = build_sh_matrix(theta, phi, l_max_sh)
    print(f"  SH basis size: {Y_real.shape[1]} real functions")

    # Fit the radial function r(theta,phi) = sum c_lm S_lm(theta,phi)
    # We fit a single scalar (r) not 3 coordinates
    n_real = Y_real.shape[1]

    # Build regularization (same approach as render_sphere_to_cow)
    from scipy import sparse
    from scipy.sparse.linalg import lsqr

    alpha = 1e-4
    reg_weights = np.zeros(n_real)
    for j, (ell, m, kind) in enumerate(real_index):
        reg_weights[j] = alpha * ell * (ell + 1)

    D = sparse.diags(np.sqrt(reg_weights))

    result = lsqr(
        sparse.vstack([sparse.csc_matrix(Y_real), D]),
        np.concatenate([r_verts, np.zeros(n_real)]),
    )
    c_radial = result[0]

    # R_0 is the monopole coefficient divided by Y_00
    # Y_00 = 1/(2*sqrt(pi)), so R_0 = c_00 / Y_00... but actually
    # the fit gives: r(theta,phi) = sum c_j S_j(theta,phi)
    # The monopole term is c_00 * Y_00 (real SH), where Y_00 = 1/(2*sqrt(pi))
    # So R_0 = c_00 * Y_00 evaluated gives a constant R_0 across the sphere.
    # Actually: r = c_00 * Y_00 + higher terms
    # We want: r = R_0 * (1 + sum epsilon_lm Y_lm)
    # So R_0 = c_00 * Y_00... no, the (0,0) basis function IS Y_00.
    # r = c_00 Y_00 + ... = (c_00 / (2*sqrt(pi))) / (2*sqrt(pi)) ...
    # Let's think more carefully.
    #
    # The fit is: r(theta,phi) = sum_j c_j S_j(theta,phi)
    # where S_j are real SH. The (0,0) term is S_{0,0} = Y_0^0 = 1/(2*sqrt(pi)).
    # So the "mean radius" is R_0 = c_00 * S_{0,0} integrated...
    # Actually, if we average r over the sphere:
    # <r> = (1/4pi) int r dOmega = c_00 * (1/4pi) int Y_00 dOmega = c_00 * Y_00 = c_00/(2*sqrt(pi))
    # Wait, int Y_00 dOmega = sqrt(4pi) * delta_{l0 m0} ...
    # Actually for real SH orthonormality: int S_j S_k dOmega = delta_jk
    # So int r dOmega = int sum c_j S_j dOmega = c_00 * int S_00 dOmega
    # int S_00 dOmega = int 1/(2*sqrt(pi)) dOmega = 4pi/(2*sqrt(pi)) = 2*sqrt(pi)
    # So <r> = c_00 * 2*sqrt(pi) / (4*pi) = c_00 / (2*sqrt(pi))
    #
    # But we want R_0 such that r = R_0 (1 + sum eps_lm Y_lm)
    # This means: r = R_0 + R_0 * sum eps_lm Y_lm
    # Matching to: r = c_00 Y_00 + sum_{j>0} c_j S_j
    # R_0 = c_00 Y_00 = c_00 / (2*sqrt(pi))
    # Wait no. The l=0 term is c_00 * S_{0,0}(theta,phi) = c_00 * Y_0^0 = c_00/(2*sqrt(pi))
    # This IS the constant part. So R_0 = c_00 / (2*sqrt(pi)).
    # Then for l>0: R_0 * eps_j * S_j = c_j * S_j => eps_j = c_j / R_0

    Y00 = 1.0 / (2.0 * np.sqrt(np.pi))
    R_0 = c_radial[0] * Y00  # The l=0 real SH coefficient times Y_00 value
    print(f"  R_0 (mean radius) = {R_0:.10f}")
    print(f"  R_eq = {R_eq:.10f}")
    print(f"  R_0 / R_eq = {R_0/R_eq:.6f}")

    # epsilon_lm = c_j / R_0 for j > 0 (the deformation coefficients in real SH basis)
    # But we want to report in complex SH basis (ell, m) format
    # Let's convert: for real SH basis S_{l,0} = Y_l^0, S_{l,m,c}, S_{l,m,s}
    # to complex: Y_l^m = (S_{l,m,c} + i*S_{l,m,s}) * sqrt(2) * ...
    # Actually let's report both real basis and convert to complex.

    # Build epsilon in complex SH basis
    # Real SH to complex SH conversion:
    # Y_l^0 = S_{l,0}  (same)
    # Y_l^m = (S_{l,m,c} - i*S_{l,m,s}) / sqrt(2) * (-1)^m  for m > 0
    # Y_l^{-m} = (S_{l,m,c} + i*S_{l,m,s}) / sqrt(2)  for m > 0
    # So if r = R_0 + sum c_j S_j = R_0 (1 + sum eps_j S_j)
    # with eps_j = c_j / R_0
    # Then in complex basis:
    #   r = R_0 (1 + sum_{l,m} eps_lm^complex Y_l^m)

    # Collect real epsilon values
    eps_real = {}
    for j, (ell, m, kind) in enumerate(real_index):
        if ell == 0:
            continue  # skip monopole
        eps_real[(ell, m, kind)] = c_radial[j] / R_0

    # Convert to complex epsilon_lm
    epsilon_lm = {}
    for ell in range(1, l_max_sh + 1):
        # m = 0
        key0 = (ell, 0, '0')
        if key0 in eps_real:
            epsilon_lm[(ell, 0)] = complex(eps_real[key0], 0.0)

        for m in range(1, ell + 1):
            ec = eps_real.get((ell, m, 'c'), 0.0)
            es = eps_real.get((ell, m, 's'), 0.0)
            # Y_l^m = (-1)^m (S_c - i*S_s) / sqrt(2)
            # So eps_{l,m}^complex = (-1)^m (ec - i*es) / sqrt(2)
            epsilon_lm[(ell, m)] = (-1)**m * (ec - 1j*es) / np.sqrt(2)
            # Y_l^{-m} = (S_c + i*S_s) / sqrt(2)
            epsilon_lm[(ell, -m)] = (ec + 1j*es) / np.sqrt(2)

    # Print spectrum
    print(f"\n  Deformation spectrum |epsilon_lm|:")
    for ell in range(1, l_max_sh + 1):
        P_l = sum(abs(epsilon_lm[(ell, m)])**2 for m in range(-ell, ell+1))
        print(f"    l={ell:2d}: P_l = {P_l:.8e}, sqrt(P_l) = {np.sqrt(P_l):.6f}")

    # =========================================================================
    # 4. Save results
    # =========================================================================
    prod_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(BASE_DIR, 'results')
    progress_dir = os.path.join(prod_dir, 'progress')
    os.makedirs(progress_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # 4a. duck_baseline.json
    Qlm_serializable = {}
    for (ell, m), val in sorted(Qlm.items()):
        Qlm_serializable[f"({ell},{m})"] = [float(np.real(val)), float(np.imag(val))]

    epsilon_serializable = {}
    for (ell, m), val in sorted(epsilon_lm.items()):
        epsilon_serializable[f"({ell},{m})"] = [float(np.real(val)), float(np.imag(val))]

    baseline = {
        "mesh": {
            "file": "data/duck.off",
            "n_vertices": int(len(V)),
            "n_faces": int(len(F)),
        },
        "geometry": {
            "volume": float(volume),
            "center_of_mass": com.tolist(),
            "bounding_box_min": bb_min.tolist(),
            "bounding_box_max": bb_max.tolist(),
            "bounding_box_size": bb_size.tolist(),
            "R_eq": float(R_eq),
            "R_0_mean_radius": float(R_0),
        },
        "quadrupole_QC": {
            "tensor": QC.tolist(),
            "eigenvalues": QC_eigenvalues.tolist(),
        },
        "inertia_tensor_I": {
            "tensor": I.tolist(),
            "eigenvalues_principal_moments": I_eigenvalues.tolist(),
        },
        "kappa_duck": kappa_duck.tolist(),
        "Qlm": Qlm_serializable,
        "epsilon_lm": epsilon_serializable,
        "power_spectrum": {},
    }

    # Power spectrum
    for ell in range(1, l_max_sh + 1):
        P_l = sum(abs(epsilon_lm[(ell, m)])**2 for m in range(-ell, ell+1))
        baseline["power_spectrum"][f"l={ell}"] = float(P_l)

    json_path = os.path.join(prod_dir, 'duck_baseline.json')
    with open(json_path, 'w') as f:
        json.dump(baseline, f, indent=2)
    print(f"\nSaved: {json_path}")

    # 4b. duck_epsilon.dat
    dat_path = os.path.join(prod_dir, 'duck_epsilon.dat')
    with open(dat_path, 'w') as f:
        f.write("# Duck SH deformation coefficients epsilon_lm\n")
        f.write("# R(theta,phi) = R_0 [1 + sum epsilon_lm Y_lm(theta,phi)]\n")
        f.write(f"# R_0 = {R_0:.10f}\n")
        f.write(f"# R_eq = {R_eq:.10f}\n")
        f.write("# ell  m  epsilon_real  epsilon_imag\n")
        for ell in range(1, l_max_sh + 1):
            for m in range(-ell, ell + 1):
                val = epsilon_lm[(ell, m)]
                f.write(f"{ell:4d} {m:4d}  {np.real(val):+.12e}  {np.imag(val):+.12e}\n")
    print(f"Saved: {dat_path}")

    # 4c. Plots (use Agg backend)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Bar chart of |epsilon_lm| vs l
    fig, ax = plt.subplots(figsize=(12, 5))
    positions = []
    heights = []
    labels_bar = []
    colors = []
    cmap = plt.cm.tab10
    pos = 0
    for ell in range(1, l_max_sh + 1):
        for m in range(-ell, ell + 1):
            positions.append(pos)
            heights.append(abs(epsilon_lm[(ell, m)]))
            labels_bar.append(f"({ell},{m})")
            colors.append(cmap(ell % 10))
            pos += 1
        pos += 0.5  # gap between ell groups

    ax.bar(positions, heights, color=colors, width=0.8, edgecolor='none')
    ax.set_ylabel(r'$|\epsilon_{\ell m}|$', fontsize=13)
    ax.set_xlabel(r'$(\ell, m)$', fontsize=13)
    ax.set_title('Duck Deformation Spectrum: $|\\epsilon_{\\ell m}|$', fontsize=14)
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-6)

    # Add ell labels
    ell_centers = []
    idx = 0
    for ell in range(1, l_max_sh + 1):
        n_m = 2*ell + 1
        center = positions[idx] + (n_m - 1) / 2
        ell_centers.append((center, ell))
        idx += n_m

    ax.set_xticks([c for c, _ in ell_centers])
    ax.set_xticklabels([f'$\\ell={l}$' for _, l in ell_centers])
    plt.tight_layout()
    plot1_path = os.path.join(results_dir, 'duck_deformation_spectrum.png')
    fig.savefig(plot1_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {plot1_path}")

    # Power spectrum P_l vs l
    ells = list(range(1, l_max_sh + 1))
    P_l_values = []
    for ell in ells:
        P_l = sum(abs(epsilon_lm[(ell, m)])**2 for m in range(-ell, ell+1))
        P_l_values.append(P_l)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(ells, P_l_values, color='steelblue', edgecolor='black', width=0.7)
    ax.set_xlabel(r'$\ell$', fontsize=14)
    ax.set_ylabel(r'$P_\ell = \sum_m |\epsilon_{\ell m}|^2$', fontsize=13)
    ax.set_title('Duck Deformation Power Spectrum', fontsize=14)
    ax.set_yscale('log')
    ax.set_xticks(ells)
    plt.tight_layout()
    plot2_path = os.path.join(results_dir, 'duck_power_spectrum.png')
    fig.savefig(plot2_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {plot2_path}")

    # =========================================================================
    # 5. Summary report
    # =========================================================================
    summary_lines = [
        "# Phase 0: Duck Baseline — Quacky Normal Modes",
        "",
        "**Date:** 2026-03-30",
        "",
        "## Mesh Properties",
        "",
        f"- **File:** `data/duck.off`",
        f"- **Vertices:** {len(V)}",
        f"- **Faces:** {len(F)}",
        f"- **Volume:** V = {volume:.10e}",
        f"- **Center of mass:** ({com[0]:.6f}, {com[1]:.6f}, {com[2]:.6f})",
        f"- **Bounding box size:** ({bb_size[0]:.6f}, {bb_size[1]:.6f}, {bb_size[2]:.6f})",
        f"- **Equivalent radius:** R_eq = (3V/4pi)^(1/3) = {R_eq:.10f}",
        f"- **Mean radius (SH fit):** R_0 = {R_0:.10f}",
        "",
        "## Cartesian Quadrupole Tensor Q^C",
        "",
        "```",
        f"  [{QC[0,0]:+.10e}  {QC[0,1]:+.10e}  {QC[0,2]:+.10e}]",
        f"  [{QC[1,0]:+.10e}  {QC[1,1]:+.10e}  {QC[1,2]:+.10e}]",
        f"  [{QC[2,0]:+.10e}  {QC[2,1]:+.10e}  {QC[2,2]:+.10e}]",
        "```",
        "",
        f"Eigenvalues: {QC_eigenvalues[0]:+.10e}, {QC_eigenvalues[1]:+.10e}, {QC_eigenvalues[2]:+.10e}",
        "",
        "## Inertia Tensor I",
        "",
        "```",
        f"  [{I[0,0]:+.10e}  {I[0,1]:+.10e}  {I[0,2]:+.10e}]",
        f"  [{I[1,0]:+.10e}  {I[1,1]:+.10e}  {I[1,2]:+.10e}]",
        f"  [{I[2,0]:+.10e}  {I[2,1]:+.10e}  {I[2,2]:+.10e}]",
        "```",
        "",
        f"Principal moments (eigenvalues): {I_eigenvalues[0]:.10e}, {I_eigenvalues[1]:.10e}, {I_eigenvalues[2]:.10e}",
        "",
        "## Dimensionless Shape Parameters",
        "",
        f"- kappa_duck = Q^C_eigenvalues / (M * R_eq^2) = ({kappa_duck[0]:+.6f}, {kappa_duck[1]:+.6f}, {kappa_duck[2]:+.6f})",
        "",
        "## Spherical Multipole Moments Q_lm",
        "",
        f"- Q_00 = {Qlm[(0,0)]:.10e} (monopole, = V * Y_00 normalization)",
        f"- Q_10 = {Qlm[(1,0)]:.2e} (should vanish at COM)",
        "",
        "Full Q_lm up to l=10 stored in `duck_baseline.json`.",
        "",
        "## SH Deformation Spectrum",
        "",
        "Radial deformation: R(theta,phi) = R_0 [1 + sum epsilon_lm Y_lm]",
        "",
        "| l | P_l = sum_m |eps_lm|^2 | sqrt(P_l) |",
        "|---|---|---|",
    ]
    for ell in range(1, l_max_sh + 1):
        P_l = sum(abs(epsilon_lm[(ell, m)])**2 for m in range(-ell, ell+1))
        summary_lines.append(f"| {ell} | {P_l:.8e} | {np.sqrt(P_l):.6f} |")

    summary_lines += [
        "",
        "## Output Files",
        "",
        "- `production-code/duck_baseline.json` -- all numerical results",
        "- `production-code/duck_epsilon.dat` -- epsilon_lm in tabular format",
        "- `results/duck_deformation_spectrum.png` -- bar chart of |epsilon_lm|",
        "- `results/duck_power_spectrum.png` -- P_l vs l",
        "",
        "## Status",
        "",
        "Phase 0 baseline computation complete. All moments, eigenvalues, and",
        "deformation coefficients computed and saved. Ready for Phase 1 (normal mode analysis).",
    ]

    summary_path = os.path.join(progress_dir, 'phase0_baseline.md')
    with open(summary_path, 'w') as f:
        f.write('\n'.join(summary_lines) + '\n')
    print(f"Saved: {summary_path}")

    print("\n=== Phase 0 complete ===")


if __name__ == '__main__':
    main()
