"""
Gravitational wave radiation computation for a rotating cow.

Computes the time-averaged third-derivative contraction of the Cartesian
quadrupole tensor and the resulting GW power, for rotation about the y-axis.
"""

import numpy as np


def rotation_matrix_y(theta):
    """Rotation matrix for rotation about the y-axis by angle theta.

    R(theta) = [[cos(theta), 0, sin(theta)],
                [0,          1, 0          ],
                [-sin(theta), 0, cos(theta)]]
    """
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])


def compute_Qdotdotdot_contraction(QC, omega, n_samples=1000):
    """Compute time-averaged <Q'''_ij Q'''^ij> for rotation about y-axis.

    Uses numerical differentiation of Q_ij(t) = R(wt) Q R(wt)^T
    via dense sampling over one period, computing the third time derivative
    numerically and averaging the contraction.

    Args:
        QC: (3,3) Cartesian quadrupole tensor at t=0
        omega: angular frequency
        n_samples: number of time samples for averaging

    Returns:
        float: time-averaged <Q'''_ij Q'''^ij> (includes omega^6 dependence)
    """
    T = 2 * np.pi / omega
    dt = T / n_samples

    # Compute Q(t) at many time points
    times = np.linspace(0, T, n_samples + 4, endpoint=False)
    Q_series = np.zeros((len(times), 3, 3))

    for k, t in enumerate(times):
        R = rotation_matrix_y(omega * t)
        Q_series[k] = R @ QC @ R.T

    # Third derivative via finite differences (4th-order central)
    # d^3f/dt^3 ~ (-f(t-2h) + 2f(t-h) - 2f(t+h) + f(t+2h)) / (2h^3)
    h = dt
    Qddd_contraction = []

    for k in range(2, len(times) - 2):
        Qddd = (-Q_series[k-2] + 2*Q_series[k-1] - 2*Q_series[k+1] + Q_series[k+2]) / (2 * h**3)
        contraction = np.sum(Qddd * Qddd)
        Qddd_contraction.append(contraction)

    return np.mean(Qddd_contraction)


def gw_power_physical(benchmark_result):
    """Compute GW power in physical (CGS) units.

    Args:
        benchmark_result: dict with keys:
            'QC': (3,3) Cartesian quadrupole tensor in benchmark units
            'I': (3,3) inertia tensor in benchmark units
            'Qddd_contraction_coeff': <Q'''Q'''> / omega^6 coefficient

    Returns:
        dict with:
            'Qddd_contraction_coeff': coefficient C such that <Q'''Q'''> = C * omega^6
            'E_dot_cgs': power in erg/s per (omega/Hz)^6
            'spindown_time_s': E/Edot in seconds per (omega/Hz)^{-4}
    """
    # Physical constants
    G = 6.674e-8       # cm^3 g^-1 s^-2
    c = 2.998e10       # cm/s

    # Benchmark -> physical unit conversion
    # Cow length ~ 2.5 m, benchmark x-extent = 1.044
    L = 2.39  # meters = 239 cm
    L_cgs = L * 100  # cm

    # rho = 1 g/cm^3, M = rho * L^3
    M_g = 1.0 * L_cgs**3  # grams

    # <Q'''Q'''> has units M^2 L^4 omega^6 in benchmark
    coeff_benchmark = benchmark_result['Qddd_contraction_coeff']
    coeff_phys = coeff_benchmark * M_g**2 * L_cgs**4  # g^2 cm^4

    # E_dot = G / (45 c^5) * <Q'''Q'''>
    E_dot_per_omega6 = G / (45 * c**5) * coeff_phys  # erg/s per (rad/s)^6

    # Convert to per Hz^6: omega = 2*pi*f, so omega^6 = (2pi)^6 * f^6
    # But paper uses omega in Hz directly (angular frequency)
    # Actually paper says (omega / 1 Hz)^6 where omega is angular frequency
    E_dot_cgs = E_dot_per_omega6  # erg/s per (omega / (rad/s))^6

    # Spindown time: E / Edot where E = (1/2) I_yy omega^2
    I_yy_benchmark = benchmark_result['I'][1, 1]
    I_yy_phys = I_yy_benchmark * M_g * L_cgs**2  # g cm^2

    E_per_omega2 = 0.5 * I_yy_phys  # erg per (rad/s)^2
    spindown_per_omega_neg4 = E_per_omega2 / E_dot_per_omega6  # s per (rad/s)^{-4}

    return {
        'Qddd_contraction_coeff': coeff_benchmark,
        'E_dot_per_omega6_cgs': E_dot_cgs,
        'E_dot_display': E_dot_cgs,  # erg/s at omega=1 rad/s
        'spindown_per_omega_neg4_s': spindown_per_omega_neg4,
        'I_yy_phys_g_cm2': I_yy_phys,
        'M_g': M_g,
        'L_cm': L_cgs,
    }
