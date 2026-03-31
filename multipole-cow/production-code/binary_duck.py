"""
Section 5: Binary Duck on Post-Newtonian Orbits
================================================
Computes the quadrupole-monopole interaction, modified Kepler relation,
system quadrupole, GW power decomposition, and inspiral phase evolution
for a binary duck system.

Outputs:
  - results/binary_duck_power.png     : P_orb, P_body, P_cross vs r/R_eq
  - results/binary_duck_dephasing.png : Dephasing DN(f) for LIGO and LISA bands
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# Physical constants (CGS)
# ============================================================================
G_cgs = 6.67430e-8       # cm^3 g^-1 s^-2
c_cgs = 2.99792458e10    # cm s^-1
M_sun_cgs = 1.98892e33   # g
pc_cm = 3.08567758e18    # cm per parsec

# ============================================================================
# Load or define duck baseline parameters
# ============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

baseline_path = os.path.join(SCRIPT_DIR, 'duck_baseline.json')
if os.path.exists(baseline_path):
    with open(baseline_path) as f:
        baseline = json.load(f)
    QC_duck = np.array(baseline['quadrupole_QC']['tensor'])
    I_duck = np.array(baseline['inertia_tensor_I']['tensor'])
    R_eq = baseline['geometry']['R_eq']
    M_duck = baseline['geometry']['volume']  # benchmark: uniform density rho=1 => M = V
    print(f"Loaded duck baseline from {baseline_path}")
else:
    # Placeholder values from results.json eigenvalues (benchmark units)
    QC_duck = np.array([
        [ 4.183e-3,  8.511e-4,  2.353e-6],
        [ 8.511e-4, -1.615e-3,  1.696e-7],
        [ 2.353e-6,  1.696e-7, -2.569e-3],
    ])
    I_duck = np.array([
        [ 7.995e-4, -2.837e-4, -7.843e-7],
        [-2.837e-4,  2.732e-3, -5.653e-8],
        [-7.843e-7, -5.653e-8,  3.050e-3],
    ])
    M_duck = 0.347   # benchmark mass (from cow paper conventions)
    V_duck = 0.0535   # volume in benchmark units
    R_eq = (3.0 * V_duck / (4.0 * np.pi)) ** (1.0 / 3.0)
    print(f"Using placeholder duck parameters (Q^C, I from results.json)")

# Diagonalize Q^C for principal values
QC_eigenvalues = np.sort(np.linalg.eigvalsh(QC_duck))[::-1]
Q_max = QC_eigenvalues[0]
Q_hat = Q_max / (M_duck * R_eq**2)

print(f"Duck Q^C eigenvalues: {QC_eigenvalues}")
print(f"R_eq = {R_eq:.6f}")
print(f"Q_hat = {Q_hat:.6f}")


# ============================================================================
# 5.1  Quadrupole-Monopole energy U_QM
# ============================================================================
def U_QM(QC, M_companion, r, n_hat, R_orientation=None):
    """
    Quadrupole-monopole interaction energy (Eq. 29).

    Parameters
    ----------
    QC : (3,3) array
        Body-frame Cartesian quadrupole tensor (traceless, symmetric).
    M_companion : float
        Mass of the companion (monopole source).
    r : float
        Orbital separation.
    n_hat : (3,) array
        Unit separation vector in the inertial frame.
    R_orientation : (3,3) array or None
        Rotation matrix from body frame to inertial frame.
        If None, identity (body frame = inertial frame).

    Returns
    -------
    float : U_QM in the same unit system as the inputs.
    """
    if R_orientation is not None:
        QC_inertial = R_orientation @ QC @ R_orientation.T
    else:
        QC_inertial = QC

    Qnn = n_hat @ QC_inertial @ n_hat
    return -G_cgs * M_companion / (2.0 * r**3) * Qnn


# ============================================================================
# 5.2  Modified Kepler relation
# ============================================================================
def omega_kepler_modified(r, M_total, Q_nn_eff):
    """
    Modified Kepler angular frequency including quadrupole correction (Eq. 31).

    Parameters
    ----------
    r : float or array
        Orbital separation.
    M_total : float
        Total mass of the binary.
    Q_nn_eff : float
        Effective projected quadrupole: M_B * Q^A_nn + M_A * Q^B_nn.

    Returns
    -------
    omega : float or array
    """
    omega2 = G_cgs * M_total / r**3 * (1.0 + 15.0 * Q_nn_eff / (2.0 * M_total * r**2))
    return np.sqrt(np.maximum(omega2, 0.0))


# ============================================================================
# 5.3  System quadrupole I^sys(t) for tidally locked configuration
# ============================================================================
def rotation_matrix_z(phi):
    """Rotation by angle phi about z-axis."""
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def system_quadrupole_tidally_locked(t, omega, r, mu, QC_A, QC_B=None):
    """
    Compute the system quadrupole I^sys_ij(t) for a tidally locked binary
    on a circular orbit in the x-y plane (Eq. 37).

    Parameters
    ----------
    t : float or 1d array
        Time(s).
    omega : float
        Orbital angular frequency.
    r : float
        Orbital separation.
    mu : float
        Reduced mass.
    QC_A, QC_B : (3,3) arrays
        Body-frame quadrupole tensors. If QC_B is None, uses QC_A.

    Returns
    -------
    I_sys : (N, 3, 3) array
    """
    if QC_B is None:
        QC_B = QC_A.copy()

    t = np.atleast_1d(t)
    N = len(t)
    I_sys = np.zeros((N, 3, 3))

    for k, tk in enumerate(t):
        phi = omega * tk
        # Separation vector
        x = r * np.array([np.cos(phi), np.sin(phi), 0.0])

        # Orbital quadrupole (traceless)
        I_orb = mu * (np.outer(x, x) - (r**2 / 3.0) * np.eye(3))

        # Body quadrupoles rotated to inertial frame (tidally locked)
        Rz = rotation_matrix_z(phi)
        Q_A_inertial = Rz @ QC_A @ Rz.T
        Q_B_inertial = Rz @ QC_B @ Rz.T  # same rotation for equal-mass tidally locked

        I_body = Q_A_inertial + Q_B_inertial

        I_sys[k] = I_orb + I_body

    return I_sys


# ============================================================================
# GW power via numerical third time derivative
# ============================================================================
def numerical_third_derivative(I_sys, dt):
    """
    Compute the third time derivative of I_sys(t) using finite differences.

    Parameters
    ----------
    I_sys : (N, 3, 3) array
    dt : float

    Returns
    -------
    Iddd : (N-6, 3, 3) array  (trimmed by 3 on each side)
    """
    # Use central differences: f'''(x) ~ [-f(x-3h) + 3f(x-2h) - 3f(x-h) + ... ] / h^3
    # More robust: chain three first derivatives
    # First derivative (central, 2nd order)
    def d1(f, dt):
        return (f[2:] - f[:-2]) / (2.0 * dt)

    I1 = d1(I_sys, dt)
    I2 = d1(I1, dt)
    I3 = d1(I2, dt)
    return I3


def gw_power_from_Iddd(Iddd):
    """
    P_GW = G/(45 c^5) <I'''_ij I'''^ij>  (time-averaged).
    """
    contraction = np.einsum('nij,nij->n', Iddd, Iddd)
    return G_cgs / (45.0 * c_cgs**5) * np.mean(contraction)


def compute_power_decomposition(omega, r, mu, QC_A, QC_B=None, n_periods=4, n_samples_per_period=200):
    """
    Compute P_orb, P_body, P_cross by computing I^sys, I^orb, I^body
    separately and taking their third derivatives.

    Returns
    -------
    P_orb, P_body, P_cross : floats (CGS, erg/s)
    """
    if QC_B is None:
        QC_B = QC_A.copy()

    T = 2.0 * np.pi / omega
    N = n_periods * n_samples_per_period
    dt = n_periods * T / N
    t = np.arange(N) * dt

    # Build orbital-only quadrupole
    I_orb_arr = np.zeros((N, 3, 3))
    I_body_arr = np.zeros((N, 3, 3))

    for k in range(N):
        phi = omega * t[k]
        x = r * np.array([np.cos(phi), np.sin(phi), 0.0])
        I_orb_arr[k] = mu * (np.outer(x, x) - (r**2 / 3.0) * np.eye(3))

        Rz = rotation_matrix_z(phi)
        I_body_arr[k] = Rz @ QC_A @ Rz.T + Rz @ QC_B @ Rz.T

    I_sys_arr = I_orb_arr + I_body_arr

    Iddd_sys = numerical_third_derivative(I_sys_arr, dt)
    Iddd_orb = numerical_third_derivative(I_orb_arr, dt)
    Iddd_body = numerical_third_derivative(I_body_arr, dt)

    P_sys = gw_power_from_Iddd(Iddd_sys)
    P_orb = gw_power_from_Iddd(Iddd_orb)
    P_body = gw_power_from_Iddd(Iddd_body)

    # Cross term: P_sys = P_orb + P_body + P_cross
    P_cross = P_sys - P_orb - P_body

    return P_orb, P_body, P_cross


# ============================================================================
# 5.5  Phase evolution: dephasing delta_Psi_Q(f)
# ============================================================================
def delta_psi_Q(f, M_total_cgs, eta, Q_hat_eff):
    """
    Quadrupole correction to the SPA phase (Eq. 51).

    delta_Psi_Q = -(75 / 64 eta) * Q_hat_eff * (pi * M_chirp * f)^{1/3}

    Parameters
    ----------
    f : float or array
        GW frequency in Hz.
    M_total_cgs : float
        Total mass in grams.
    eta : float
        Symmetric mass ratio.
    Q_hat_eff : float
        Dimensionless effective quadrupole parameter.

    Returns
    -------
    delta_psi : float or array (radians)
    """
    M_chirp = eta**(3.0 / 5.0) * M_total_cgs
    # Convert M_chirp to seconds: M_s = G M / c^3
    M_chirp_s = G_cgs * M_chirp / c_cgs**3
    v_param = (np.pi * M_chirp_s * f) ** (1.0 / 3.0)
    return -75.0 / (64.0 * eta) * Q_hat_eff * v_param


def dephasing_cycles(f, M_total_cgs, eta, Q_hat_eff):
    """
    Number of dephasing cycles: DN = |delta_Psi_Q(f) - delta_Psi_Q(f_low)| / (2 pi).
    Returns DN(f) with respect to the lowest frequency in the array.
    """
    psi = delta_psi_Q(f, M_total_cgs, eta, Q_hat_eff)
    return np.abs(psi - psi[0]) / (2.0 * np.pi)


# ============================================================================
# Main: produce plots
# ============================================================================
def main():
    print("\n" + "=" * 70)
    print("Binary Duck on Post-Newtonian Orbits")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Setup: equal-mass binary with duck quadrupole, in CGS
    # -----------------------------------------------------------------------
    # Use 1.4 M_sun per duck (neutron-star-like)
    M_A = 1.4 * M_sun_cgs
    M_B = 1.4 * M_sun_cgs
    M_total = M_A + M_B
    mu = M_A * M_B / M_total
    eta = mu / M_total   # = 0.25 for equal mass

    # Physical duck radius: use 12 km (NS-like)
    R_phys = 1.2e6  # 12 km in cm

    # Scale the benchmark Q^C to physical units
    # Q^C_phys = Q^C_bench * (M_phys / M_bench) * (R_phys / R_bench)^2
    # But we need Q^C in physical (CGS) units: [g cm^2]
    # Q_hat is dimensionless, so Q^C_phys = Q_hat * M_A * R_phys^2
    QC_phys = Q_hat * M_A * R_phys**2 * np.diag([1.0, -QC_eigenvalues[1]/QC_eigenvalues[0],
                                                    -QC_eigenvalues[2]/QC_eigenvalues[0]])
    # For tidally locked: Q_nn = Q_xx = largest eigenvalue direction
    Q_nn = Q_hat * M_A * R_phys**2
    Q_nn_eff = M_B * Q_nn + M_A * Q_nn  # equal-mass, both ducks

    print(f"\nBinary parameters:")
    print(f"  M_A = M_B = {M_A/M_sun_cgs:.1f} M_sun")
    print(f"  R_phys = {R_phys/1e5:.0f} km")
    print(f"  Q_hat = {Q_hat:.6f}")
    print(f"  Q_nn (phys) = {Q_nn:.4e} g cm^2")
    print(f"  eta = {eta:.4f}")

    # -----------------------------------------------------------------------
    # Plot 1: GW power decomposition vs r/R_phys
    # -----------------------------------------------------------------------
    print("\nComputing GW power decomposition...")
    r_over_R = np.logspace(np.log10(3.0), np.log10(100.0), 30)
    r_arr = r_over_R * R_phys

    P_orb_arr = np.zeros_like(r_arr)
    P_body_arr = np.zeros_like(r_arr)
    P_cross_arr = np.zeros_like(r_arr)

    for i, r in enumerate(r_arr):
        omega = omega_kepler_modified(r, M_total, Q_nn_eff)
        Po, Pb, Pc = compute_power_decomposition(omega, r, mu, QC_phys, QC_phys,
                                                  n_periods=4, n_samples_per_period=200)
        P_orb_arr[i] = Po
        P_body_arr[i] = Pb
        P_cross_arr[i] = abs(Pc)  # can be negative; plot abs on log scale
        if i % 10 == 0:
            print(f"  r/R = {r_over_R[i]:.1f}: P_orb={Po:.3e}, P_body={Pb:.3e}, P_cross={Pc:.3e}")

    # Identify sign of cross term for annotation
    omega_test = omega_kepler_modified(r_arr[0], M_total, Q_nn_eff)
    _, _, Pc_test = compute_power_decomposition(omega_test, r_arr[0], mu, QC_phys, QC_phys)
    cross_sign = "negative (destructive)" if Pc_test < 0 else "positive (constructive)"
    print(f"  Cross-term sign at r/R=3: {cross_sign}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(r_over_R, P_orb_arr, 'b-', lw=2, label=r'$P_{\mathrm{orb}}$ (Peters-Mathews)')
    ax.loglog(r_over_R, P_body_arr, 'r--', lw=2, label=r'$P_{\mathrm{body}}$')
    ax.loglog(r_over_R, P_cross_arr, 'g-.', lw=2,
              label=r'$|P_{\mathrm{cross}}|$ (' + cross_sign.split('(')[1].rstrip(')') + ')')
    ax.set_xlabel(r'$r / R_{\mathrm{eq}}$', fontsize=14)
    ax.set_ylabel(r'$P_{\mathrm{GW}}$ [erg/s]', fontsize=14)
    ax.set_title('Binary Duck: GW Power Decomposition (Tidally Locked)', fontsize=14)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    # Add scaling annotations
    ax.text(0.05, 0.15, r'$P_{\mathrm{orb}} \propto r^{-5}$', transform=ax.transAxes, fontsize=11, color='b')
    ax.text(0.05, 0.08, r'$P_{\mathrm{body}} \propto r^{-9}$', transform=ax.transAxes, fontsize=11, color='r')
    ax.text(0.05, 0.01, r'$P_{\mathrm{cross}} \propto r^{-7}$', transform=ax.transAxes, fontsize=11, color='g')

    power_path = os.path.join(RESULTS_DIR, 'binary_duck_power.png')
    fig.tight_layout()
    fig.savefig(power_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {power_path}")

    # -----------------------------------------------------------------------
    # Plot 2: Dephasing DN(f) for LIGO and LISA bands
    # -----------------------------------------------------------------------
    print("\nComputing dephasing...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- LIGO band: stellar-mass ducks ---
    M_ligo = 2.8 * M_sun_cgs
    eta_ligo = 0.25
    f_ligo = np.logspace(np.log10(10.0), np.log10(1000.0), 500)  # Hz

    DN_ligo = dephasing_cycles(f_ligo, M_ligo, eta_ligo, Q_hat)

    ax = axes[0]
    ax.loglog(f_ligo, DN_ligo, 'b-', lw=2)
    ax.axhline(1.0, color='k', ls=':', lw=1, alpha=0.5)
    ax.text(15, 1.3, r'$\Delta N = 1$ cycle', fontsize=10, alpha=0.6)
    ax.set_xlabel('$f$ [Hz]', fontsize=13)
    ax.set_ylabel(r'$\Delta N$ [cycles]', fontsize=13)
    ax.set_title(r'LIGO band ($M = 2.8\, M_\odot$)', fontsize=13)
    ax.grid(True, which='both', alpha=0.3)
    ax.text(0.05, 0.9, rf'$\hat{{Q}} = {Q_hat:.4f}$', transform=ax.transAxes, fontsize=11)
    ax.text(0.05, 0.82, rf'$\Delta N_{{max}} = {DN_ligo[-1]:.4f}$', transform=ax.transAxes, fontsize=11)

    # --- LISA band: supermassive ducks ---
    M_lisa = 1.0e7 * M_sun_cgs
    eta_lisa = 0.25
    f_lisa = np.logspace(np.log10(1e-4), np.log10(0.1), 500)  # Hz

    DN_lisa = dephasing_cycles(f_lisa, M_lisa, eta_lisa, Q_hat)

    ax = axes[1]
    ax.loglog(f_lisa, DN_lisa, 'r-', lw=2)
    ax.axhline(1.0, color='k', ls=':', lw=1, alpha=0.5)
    ax.text(2e-4, 1.3, r'$\Delta N = 1$ cycle', fontsize=10, alpha=0.6)
    ax.set_xlabel('$f$ [Hz]', fontsize=13)
    ax.set_ylabel(r'$\Delta N$ [cycles]', fontsize=13)
    ax.set_title(r'LISA band ($M = 10^7\, M_\odot$)', fontsize=13)
    ax.grid(True, which='both', alpha=0.3)
    ax.text(0.05, 0.9, rf'$\hat{{Q}} = {Q_hat:.4f}$', transform=ax.transAxes, fontsize=11)
    ax.text(0.05, 0.82, rf'$\Delta N_{{max}} = {DN_lisa[-1]:.4f}$', transform=ax.transAxes, fontsize=11)

    dephasing_path = os.path.join(RESULTS_DIR, 'binary_duck_dephasing.png')
    fig.tight_layout()
    fig.savefig(dephasing_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {dephasing_path}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Duck Q_hat = {Q_hat:.6f}")
    print(f"  Power decomposition at r/R = 10:")
    idx10 = np.argmin(np.abs(r_over_R - 10.0))
    print(f"    P_orb   = {P_orb_arr[idx10]:.4e} erg/s")
    print(f"    P_body  = {P_body_arr[idx10]:.4e} erg/s")
    print(f"    P_cross = {P_cross_arr[idx10]:.4e} erg/s  (|P_cross|)")
    print(f"    P_cross/P_orb = {P_cross_arr[idx10]/P_orb_arr[idx10]:.4e}")
    print(f"    P_body/P_orb  = {P_body_arr[idx10]/P_orb_arr[idx10]:.4e}")
    print(f"  LIGO dephasing (10-1000 Hz): DN = {DN_ligo[-1]:.6f} cycles")
    print(f"  LISA dephasing (0.1 mHz - 0.1 Hz): DN = {DN_lisa[-1]:.6f} cycles")
    print("=" * 70)


if __name__ == '__main__':
    main()
