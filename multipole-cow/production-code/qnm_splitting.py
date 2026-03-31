#!/usr/bin/env python3
"""
Quacky Normal Mode Splitting Calculator
========================================

Computes the perturbative QNM frequency splitting for a duck-shaped
compact object, given its spherical harmonic deformation coefficients
epsilon_{ell,m}.

The perturbation matrix for an (n, ell) multiplet is:

    V_{m1,m2} = -|R'_{n,ell}(R0)|^2 R0^3
                * sum_{ell',m'} eps_{ell',m'}
                * (-1)^m1 * sqrt((2*ell+1)^2 * (2*ell'+1) / (4*pi))
                * wigner_3j(ell, ell, ell', 0, 0, 0)
                * wigner_3j(ell, ell, ell', -m1, m2, m')

with the selection rule m' = m1 - m2.

Usage:
    python qnm_splitting.py

Outputs:
    - Printed table of split frequencies for ell=2,3,4
    - results/grotrian_diagram.png
"""

import os
import numpy as np
from scipy.special import spherical_jn
from sympy.physics.wigner import wigner_3j as sympy_wigner_3j
from functools import lru_cache
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Wigner 3j wrapper (cached for performance)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4096)
def w3j(j1, j2, j3, m1, m2, m3):
    """Evaluate the Wigner 3j symbol as a float."""
    return float(sympy_wigner_3j(j1, j2, j3, m1, m2, m3))


# ---------------------------------------------------------------------------
# Spherical Bessel zeros
# ---------------------------------------------------------------------------

def bessel_zero(ell, n, num_points=10000, r_max=50.0):
    """Find the n-th positive zero of j_ell(x) by bisection on a grid."""
    x = np.linspace(1e-8, r_max, num_points)
    vals = spherical_jn(ell, x)
    zeros = []
    for i in range(len(vals) - 1):
        if vals[i] * vals[i + 1] < 0:
            # Bisect
            a, b = x[i], x[i + 1]
            for _ in range(60):
                mid = 0.5 * (a + b)
                if spherical_jn(ell, mid) * spherical_jn(ell, a) < 0:
                    b = mid
                else:
                    a = mid
            zeros.append(0.5 * (a + b))
            if len(zeros) == n:
                break
    if len(zeros) < n:
        raise RuntimeError(f"Could not find {n} zeros of j_{ell}")
    return zeros[n - 1]


def bessel_deriv_at_zero(ell, n, R0=1.0):
    """Compute j'_ell(z_{n,ell}) / R0 where z is the n-th zero."""
    z = bessel_zero(ell, n)
    # j'_ell(x) via derivative=True
    jp = spherical_jn(ell, z, derivative=True)
    return jp / R0


# ---------------------------------------------------------------------------
# Perturbation matrix
# ---------------------------------------------------------------------------

def build_perturbation_matrix(ell, epsilon_dict, n=1, R0=1.0):
    """
    Build the (2*ell+1) x (2*ell+1) perturbation matrix V for the
    (n, ell) multiplet.

    Parameters
    ----------
    ell : int
        Angular quantum number of the multiplet being split.
    epsilon_dict : dict
        Deformation coefficients {(ell', m'): epsilon_value}.
        Only even ell' will contribute (selection rule).
    n : int
        Radial order (default 1 = f-mode).
    R0 : float
        Sphere radius.

    Returns
    -------
    V : ndarray, shape (2*ell+1, 2*ell+1)
        Hermitian perturbation matrix. Eigenvalues give delta(omega^2).
    """
    dim = 2 * ell + 1
    V = np.zeros((dim, dim), dtype=complex)

    # Radial prefactor: |R'_{n,ell}(R0)|^2 * R0^3
    Rp = bessel_deriv_at_zero(ell, n, R0)
    radial_prefactor = abs(Rp) ** 2 * R0 ** 3

    m_values = list(range(-ell, ell + 1))

    for i, m1 in enumerate(m_values):
        for j, m2 in enumerate(m_values):
            mp = m1 - m2  # selection rule

            val = 0.0 + 0.0j
            for (ellp, mp_key), eps in epsilon_dict.items():
                if mp_key != mp:
                    continue
                if ellp % 2 != 0:
                    continue  # parity selection
                if ellp > 2 * ell:
                    continue  # triangle inequality
                if ellp < 0:
                    continue

                # 3j symbols
                j1 = w3j(ell, ell, ellp, 0, 0, 0)
                if abs(j1) < 1e-15:
                    continue
                j2 = w3j(ell, ell, ellp, -m1, m2, mp)
                if abs(j2) < 1e-15:
                    continue

                angular = ((-1) ** m1
                           * np.sqrt((2 * ell + 1) ** 2 * (2 * ellp + 1)
                                     / (4 * np.pi))
                           * j1 * j2)
                val += eps * angular

            V[i, j] = -radial_prefactor * val

    # Enforce Hermiticity (should be by construction, but numerical safety)
    V = 0.5 * (V + V.conj().T)
    return V


# ---------------------------------------------------------------------------
# Splitting calculation
# ---------------------------------------------------------------------------

def compute_splitting(ell, epsilon_dict, n=1, R0=1.0):
    """
    Compute the split frequencies for the (n, ell) multiplet.

    Returns
    -------
    omega0 : float
        Unperturbed frequency.
    delta_omega2 : ndarray
        Eigenvalues of V (shifts in omega^2).
    split_omega : ndarray
        Perturbed frequencies omega_{n,ell,alpha}.
    """
    V = build_perturbation_matrix(ell, epsilon_dict, n, R0)
    delta_omega2 = np.sort(np.linalg.eigvalsh(V).real)

    z = bessel_zero(ell, n)
    omega0 = z / R0
    omega0_sq = omega0 ** 2

    split_omega = np.sqrt(np.maximum(omega0_sq + delta_omega2, 0.0))
    return omega0, delta_omega2, split_omega


# ---------------------------------------------------------------------------
# Main: table + Grotrian diagram
# ---------------------------------------------------------------------------

def print_results(epsilon_dict, ells=(2, 3, 4), n=1, R0=1.0):
    """Print a formatted table of QNM splitting results."""
    print("=" * 72)
    print("Quacky Normal Mode Splitting")
    print("=" * 72)
    print(f"\nDeformation coefficients:")
    for (lp, mp), val in sorted(epsilon_dict.items()):
        print(f"  epsilon_({lp},{mp:+d}) = {val: .4f}")
    print()

    all_results = {}
    for ell in ells:
        omega0, dw2, split_w = compute_splitting(ell, epsilon_dict, n, R0)
        all_results[ell] = (omega0, dw2, split_w)

        print(f"--- ell = {ell}  (degeneracy {2*ell+1}) ---")
        print(f"  Unperturbed:  omega_0 = {omega0:.6f}  (z = {omega0*R0:.6f})")
        print(f"  {'alpha':>6s}  {'delta(w^2)':>14s}  {'omega':>12s}  "
              f"{'delta_omega/omega':>18s}")
        for alpha in range(len(split_w)):
            dw = split_w[alpha] - omega0
            print(f"  {alpha:6d}  {dw2[alpha]:14.6f}  {split_w[alpha]:12.6f}  "
                  f"{dw/omega0:18.6f}")
        print()

    return all_results


def plot_grotrian(all_results, epsilon_dict, outpath="results/grotrian_diagram.png"):
    """
    Generate a Grotrian (energy-level) diagram showing the sphere-to-duck
    splitting of QNM multiplets.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    ells = sorted(all_results.keys())
    n_panels = len(ells)

    # Layout: sphere levels on left, duck levels on right
    x_sphere = 0.2
    x_duck = 0.8
    x_mid = 0.5
    level_width = 0.12

    colors = {2: "#2166ac", 3: "#b2182b", 4: "#1b7837"}

    y_offset = 0.0
    y_spacing = 1.0 / (n_panels + 1)

    # Normalize frequencies for plotting
    all_omega = []
    for ell in ells:
        omega0, dw2, split_w = all_results[ell]
        all_omega.append(omega0)
        all_omega.extend(split_w)

    omega_min = min(all_omega) * 0.95
    omega_max = max(all_omega) * 1.05

    def freq_to_y(omega):
        return (omega - omega_min) / (omega_max - omega_min)

    for ell in ells:
        omega0, dw2, split_w = all_results[ell]
        color = colors.get(ell, "#333333")
        deg = 2 * ell + 1

        # Sphere level (degenerate)
        y0 = freq_to_y(omega0)
        ax.plot([x_sphere - level_width, x_sphere + level_width],
                [y0, y0], color=color, linewidth=2.5, solid_capstyle="round")
        ax.text(x_sphere - level_width - 0.02, y0,
                f"$\\ell={ell}$\n$({deg}\\times)$",
                ha="right", va="center", fontsize=10, color=color)
        ax.text(x_sphere + level_width + 0.01, y0,
                f"$\\omega = {omega0:.3f}$",
                ha="left", va="center", fontsize=8, color="gray")

        # Duck levels (split)
        for alpha, w in enumerate(split_w):
            yd = freq_to_y(w)
            ax.plot([x_duck - level_width, x_duck + level_width],
                    [yd, yd], color=color, linewidth=2.0, solid_capstyle="round")
            ax.text(x_duck + level_width + 0.01, yd,
                    f"$\\omega = {w:.3f}$",
                    ha="left", va="center", fontsize=7, color="gray")
            # Connection line
            ax.plot([x_sphere + level_width + 0.01, x_duck - level_width - 0.01],
                    [y0, yd], color=color, linewidth=0.5, alpha=0.4,
                    linestyle="--")

    # Labels
    ax.text(x_sphere, 1.05, "Sphere\n(degenerate)", ha="center", va="bottom",
            fontsize=12, fontweight="bold", transform=ax.transAxes)
    ax.text(x_duck, 1.05, "Duck\n(split)", ha="center", va="bottom",
            fontsize=12, fontweight="bold", transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.1)
    ax.set_ylabel("Frequency $\\omega / R_0^{-1}$", fontsize=12)
    ax.set_title("Gravitational Zeeman Splitting: Sphere $\\to$ Duck QNMs",
                 fontsize=14, fontweight="bold")

    # Add epsilon annotation
    eps_str = ", ".join(
        [f"$\\varepsilon_{{{lp},{mp}}}={val}$"
         for (lp, mp), val in sorted(epsilon_dict.items())]
    )
    ax.text(0.5, -0.08, eps_str, ha="center", va="top",
            fontsize=9, transform=ax.transAxes, color="gray")

    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(left=True, labelleft=True)

    # Custom y-ticks from actual frequency range
    yticks_vals = np.linspace(0, 1, 6)
    yticks_labels = [f"{omega_min + v*(omega_max - omega_min):.2f}"
                     for v in yticks_vals]
    ax.set_yticks(yticks_vals)
    ax.set_yticklabels(yticks_labels)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"Grotrian diagram saved to {outpath}")
    plt.close()


def load_epsilon_file(filepath):
    """Load duck deformation coefficients from a .dat file.

    File format:
        # header lines starting with #
        ell  m  epsilon_real  epsilon_imag

    Returns:
        epsilon_dict: {(ell, m): complex} -- only even-ell entries
                      relevant for QNM splitting
        metadata: dict with R0, R_eq, scale_factor if present
    """
    epsilon_dict = {}
    metadata = {}

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Parse metadata from comments
                if "R_0 =" in line:
                    metadata["R0"] = float(line.split("=")[1].strip())
                elif "R_eq =" in line:
                    metadata["R_eq"] = float(line.split("=")[1].strip())
                elif "Scale factor:" in line:
                    metadata["scale_factor"] = float(line.split(":")[1].strip())
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            ell = int(parts[0])
            m = int(parts[1])
            eps_real = float(parts[2])
            eps_imag = float(parts[3])
            epsilon_dict[(ell, m)] = complex(eps_real, eps_imag)

    return epsilon_dict, metadata


def main():
    """Run the QNM splitting calculation with real duck coefficients."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Load real duck epsilon data at different scalings
    eps_files = {
        "0.1": os.path.join(script_dir, "duck_epsilon_0.1.dat"),
        "0.3": os.path.join(script_dir, "duck_epsilon_0.3.dat"),
        "full": os.path.join(script_dir, "duck_epsilon.dat"),
    }

    for label, filepath in eps_files.items():
        if not os.path.exists(filepath):
            print(f"Skipping {label}: {filepath} not found")
            continue

        print(f"\n{'#' * 72}")
        print(f"# Duck epsilon scale: {label}")
        print(f"# File: {os.path.basename(filepath)}")
        print(f"{'#' * 72}")

        epsilon_dict, metadata = load_epsilon_file(filepath)
        R0 = metadata.get("R0", 1.0)

        # Filter to even ell only (selection rule for splitting)
        eps_even = {k: v for k, v in epsilon_dict.items() if k[0] % 2 == 0}

        print(f"\nMetadata: R0={R0:.6f}, R_eq={metadata.get('R_eq', 'N/A')}")
        print(f"Total coefficients loaded: {len(epsilon_dict)}")
        print(f"Even-ell coefficients (active for splitting): {len(eps_even)}")

        # Compute and print results
        all_results = print_results(eps_even, ells=(2, 3, 4), R0=R0)

        # Generate Grotrian diagram
        outpath = os.path.join(script_dir, "results",
                               f"grotrian_diagram_{label}.png")
        plot_grotrian(all_results, eps_even, outpath=outpath)


if __name__ == "__main__":
    main()
