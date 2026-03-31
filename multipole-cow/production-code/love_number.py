#!/usr/bin/env python3
"""
Tidal Love Numbers of the Duck
===============================

Computes:
1. TOV solution for a Gamma=2 polytrope (K=100, rho_c = 1.28e-3)
2. Coupled TOV + y(r) Riccati integration for the tidal perturbation
3. y = R h_2'/h_2 at the surface
4. k_2 via the Hinderer (2008) formula
5. Lambda = (2/3) k_2 C^{-5}
6. Newtonian k_2 for n=0, 0.5, 1, 1.5, 2, 3 polytropes (Brooker & Olle 1955)
7. Elastic k_2 for a rubber duck
8. Plot of k_2 vs compactness C

The Riccati variable y(r) = r h_2'(r)/h_2(r) satisfies a first-order ODE
derived from the Hinderer (2008) master equation.  This avoids numerical
issues with the second-order h_2 ODE.

References:
  - Hinderer, ApJ 677, 1216 (2008)
  - Hinderer, Lackey, Lang & Read, PRD 81, 123016 (2010)
  - Flanagan & Hinderer, PRD 77, 021502 (2008)
  - Damour & Nagar, PRD 80, 084035 (2009)
  - Postnikov, Prakash & Lattimer, PRD 82, 024016 (2010)
  - Brooker & Olle, MNRAS 115, 101 (1955)
  - Love, Proc. R. Soc. A 82, 73 (1909)

All quantities in geometric units (G = c = 1) unless noted.
"""

import numpy as np
from scipy.integrate import solve_ivp
import os


# ---------------------------------------------------------------------------
# TOV + y-Riccati coupled system
# ---------------------------------------------------------------------------

def tov_y_rhs(r, state, K, Gamma):
    """
    Right-hand side for the coupled TOV + y-Riccati system.

    State: [p, m, nu, y]
      p   = pressure
      m   = enclosed gravitational mass M(r)
      nu  = metric function (g_tt = -e^nu)
      y   = r * h_2'(r) / h_2(r)  (Riccati variable)

    The y-Riccati ODE is obtained by substituting y = r h'/h into the
    Hinderer (2008) master equation for h_2.  See Postnikov et al. (2010)
    eq. (12) and Hinderer et al. (2010) eq. (14).

    For a polytrope p = K * rho_rest^Gamma, total energy density
    epsilon = rho_rest + p/(Gamma - 1).
    """
    p, m, nu_val, y_val = state

    if p <= 0 or r <= 1e-15:
        return [0, 0, 0, 0]

    # Rest-mass density from EoS
    rho_rest = (p / K) ** (1.0 / Gamma)
    eps = rho_rest + p / (Gamma - 1.0)  # total energy density

    if r <= 2.0 * m:
        return [0, 0, 0, 0]

    e_lambda = 1.0 / (1.0 - 2.0 * m / r)

    # TOV equations
    dp_dr = -(eps + p) * (m + 4.0 * np.pi * r**3 * p) / (r * (r - 2.0 * m))
    dm_dr = 4.0 * np.pi * r**2 * eps
    dnu_dr = 2.0 * (m + 4.0 * np.pi * r**3 * p) / (r * (r - 2.0 * m))

    # d(epsilon)/dp for the EoS
    drho_rest_dp = 1.0 / (K * Gamma * rho_rest ** (Gamma - 1.0))
    deps_dp = drho_rest_dp + 1.0 / (Gamma - 1.0)

    # Coefficients of the h_2 master ODE: h'' + A(r) h' + B(r) h = 0
    # where A and B are from Hinderer (2008) eq. (6)
    A = 2.0 / r + e_lambda * (2.0 * m / r**2 + 4.0 * np.pi * r * (p - eps))

    B = -(
        6.0 * e_lambda / r**2
        - 4.0 * np.pi * e_lambda * (5.0 * eps + 9.0 * p + (eps + p) * deps_dp)
        + dnu_dr**2
    )

    # y-Riccati equation: y' = -(y^2 + y + r^2 B) / r - A * y
    # Derivation:
    #   h = f, h' = y f/r, h'' = (y' f/r + y^2 f/r^2 - y f/r^2)
    #   [since h'' = d/dr(yh/r) = y'h/r + y h'/r - yh/r^2 = y'h/r + y^2 h/r^2 - yh/r^2]
    #   Substituting into h'' + A h' + B h = 0:
    #   (y'/r + y^2/r^2 - y/r^2) + A*y/r + B = 0
    #   y' = -y^2/r + y/r - A*r*y/r - B*r = -(y^2 - y)/r - A*y - B*r
    #   WRONG SIGN: h'' + A h' + B h = 0 means B is as defined above (negative of the
    #   "potential" term).  Let me redo:
    #
    #   From h'' + A h' + B h = 0 with y = r h'/h:
    #   h'/h = y/r, h''/h = y'/r + y^2/r^2 - y/r^2
    #   => y'/r + y^2/r^2 - y/r^2 + A*y/r + B = 0
    #   => y' = -y^2/r + y/r - A*y - B*r
    #
    # That is: dy/dr = -y^2/r + y/r - A*y - B*r
    #        = -(y^2 - y)/r - A*y - B*r

    dy_dr = -(y_val**2 - y_val) / r - A * y_val - B * r

    return [dp_dr, dm_dr, dnu_dr, dy_dr]


def surface_event(r, state, K, Gamma):
    """Event: p = 0 (stellar surface)."""
    return state[0] - 1e-15


surface_event.terminal = True
surface_event.direction = -1


def integrate_tov_y(K, Gamma, rho_c, r_start=1e-5, r_max=50.0):
    """
    Integrate the coupled TOV + y-Riccati system from center to surface.

    Returns: (R_star, M_star, C, y_surface, k2, Lambda)
    """
    # Central values
    p_c = K * rho_c**Gamma
    eps_c = rho_c + p_c / (Gamma - 1.0)

    # Near-origin initial conditions
    r0 = r_start
    m0 = (4.0 * np.pi / 3.0) * eps_c * r0**3
    p0 = p_c - (2.0 * np.pi / 3.0) * (eps_c + p_c) * (eps_c + 3.0 * p_c) * r0**2
    nu0 = 0.0  # arbitrary

    # y near origin: h_2 ~ r^2, so h_2' ~ 2r, y = r * 2r / r^2 = 2
    y0 = 2.0

    state0 = [p0, m0, nu0, y0]

    sol = solve_ivp(
        tov_y_rhs,
        [r0, r_max],
        state0,
        args=(K, Gamma),
        method="RK45",
        events=surface_event,
        max_step=0.005,
        rtol=1e-10,
        atol=1e-13,
    )

    if sol.t_events[0].size == 0:
        raise RuntimeError("Surface not reached; increase r_max or check parameters.")

    R_star = sol.t_events[0][0]
    p_s, M_star, nu_s, y_surf = sol.y_events[0][0]

    C = M_star / R_star

    # Hinderer formula for k_2
    k2 = hinderer_k2(C, y_surf)
    Lambda = (2.0 / 3.0) * k2 * C ** (-5)

    return R_star, M_star, C, y_surf, k2, Lambda


def hinderer_k2(C, y):
    """
    Tidal apsidal constant k_2 from Hinderer (2008), eq. (11).

    Parameters:
        C : compactness M/R
        y : logarithmic derivative R h_2'/h_2 at the surface
    """
    num = (8.0 / 5.0) * C**5 * (1.0 - 2.0 * C)**2 * (2.0 + 2.0 * C * (y - 1.0) - y)

    den = (
        2.0 * C * (6.0 - 3.0 * y + 3.0 * C * (5.0 * y - 8.0))
        + 4.0 * C**3 * (
            13.0 - 11.0 * y + C * (3.0 * y - 2.0) + 2.0 * C**2 * (1.0 + y)
        )
        + 3.0 * (1.0 - 2.0 * C)**2
        * (2.0 - y + 2.0 * C * (y - 1.0))
        * np.log(1.0 - 2.0 * C)
    )

    return num / den


# ---------------------------------------------------------------------------
# Newtonian Love numbers (Brooker & Olle 1955)
# ---------------------------------------------------------------------------

def newtonian_k2_table():
    """
    Return Newtonian k_2 values from Brooker & Olle (1955).
    Computed from the Clairaut-Radau equation.
    """
    return {
        0.0: 0.7500,   # uniform density (exact: 3/4)
        0.5: 0.4489,   # Brooker & Olle
        1.0: 0.2599,   # Brooker & Olle
        1.5: 0.1433,   # Brooker & Olle
        2.0: 0.0728,   # Brooker & Olle
        3.0: 0.01157,  # Brooker & Olle
    }


# ---------------------------------------------------------------------------
# Elastic Love number for rubber duck
# ---------------------------------------------------------------------------

def elastic_k2(mu_shear, rho, g, R):
    """
    Elastic Love number for a solid body (Love 1909).

    k_2 = (3/2) / (1 + 19*mu_shear / (2*rho*g*R))
    """
    rigidity_ratio = 19.0 * mu_shear / (2.0 * rho * g * R)
    return 1.5 / (1.0 + rigidity_ratio)


# ---------------------------------------------------------------------------
# k_2 vs compactness scan
# ---------------------------------------------------------------------------

def scan_k2_vs_C(K=100.0, Gamma=2.0, rho_c_min=1e-4, rho_c_max=5e-3, N=50):
    """Compute k_2 and Lambda for a range of central densities."""
    rho_c_arr = np.geomspace(rho_c_min, rho_c_max, N)
    results = {"rho_c": rho_c_arr}
    for key in ["C", "k2", "Lambda", "M", "R"]:
        results[key] = np.full(N, np.nan)

    for i, rho_c in enumerate(rho_c_arr):
        try:
            R_s, M_s, C, y_s, k2, Lam = integrate_tov_y(K, Gamma, rho_c)
            results["C"][i] = C
            results["k2"][i] = k2
            results["Lambda"][i] = Lam
            results["M"][i] = M_s
            results["R"][i] = R_s
        except Exception:
            pass

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_k2_vs_C(results, save_path="results/love_number_vs_C.png"):
    """Plot k_2 vs compactness C."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C = results["C"]
    k2 = results["k2"]
    mask = np.isfinite(C) & np.isfinite(k2) & (k2 > 0)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(C[mask], k2[mask], "b-", linewidth=2, label=r"$\Gamma=2$ polytrope (TOV)")

    # Mark the fiducial duck point
    idx = np.argmin(np.abs(C[mask] - 0.145))
    C_m, k2_m = C[mask], k2[mask]
    ax.plot(C_m[idx], k2_m[idx], "r*", markersize=15,
            label=rf"Nuclear duck ($C={C_m[idx]:.3f}$, $k_2={k2_m[idx]:.4f}$)")

    # BH limit
    ax.axhline(y=0, color="k", linestyle="--", alpha=0.5, label="BH limit ($k_2=0$)")

    # Newtonian limit
    ax.axhline(y=0.260, color="gray", linestyle=":", alpha=0.5,
               label=r"Newtonian $k_2$ ($n=1$, $\Gamma=2$)")

    ax.set_xlabel(r"Compactness $C = M/R$", fontsize=14)
    ax.set_ylabel(r"Tidal apsidal constant $k_2$", fontsize=14)
    ax.set_title(r"Tidal Love Number vs Compactness ($\Gamma=2$ Polytrope)", fontsize=14)
    ax.legend(fontsize=11, loc="upper right")
    ax.set_xlim(0, 0.35)
    ax.set_ylim(-0.01, 0.30)
    ax.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {save_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("  Tidal Love Numbers of the Duck")
    print("=" * 72)
    print()

    # -----------------------------------------------------------------------
    # 1. Fiducial nuclear duck: TOV + y-Riccati integration
    # -----------------------------------------------------------------------
    K, Gamma, rho_c = 100.0, 2.0, 1.28e-3
    print(f"EoS: Gamma={Gamma}, K={K}, rho_c={rho_c:.4e}")
    print()

    R_s, M_s, C, y_s, k2, Lambda = integrate_tov_y(K, Gamma, rho_c)

    print("--- Fiducial Nuclear Duck (TOV + y-Riccati) ---")
    print(f"  R_star    = {R_s:.4f}  (geometric units)")
    print(f"  M_star    = {M_s:.4f}  (geometric units)")
    print(f"  C = M/R   = {C:.5f}")
    print(f"  y_surface = {y_s:.5f}")
    print(f"  k_2       = {k2:.6f}")
    print(f"  Lambda    = {Lambda:.1f}")
    print()

    # -----------------------------------------------------------------------
    # 2. Newtonian k_2 for polytropes (Brooker & Olle 1955)
    # -----------------------------------------------------------------------
    print("--- Newtonian k_2 for Polytropes (Brooker & Olle 1955) ---")
    print(f"{'n':>5s}  {'Gamma':>10s}  {'k_2^N':>10s}")
    print("-" * 32)

    lookup = newtonian_k2_table()
    for n_val in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]:
        k2_lu = lookup[n_val]
        if n_val == 0:
            gamma_str = "inf"
        else:
            gamma_str = f"{1.0 + 1.0 / n_val:.4f}"
        print(f"{n_val:5.1f}  {gamma_str:>10s}  {k2_lu:10.4f}")
    print()

    # -----------------------------------------------------------------------
    # 3. Elastic k_2 for rubber duck
    # -----------------------------------------------------------------------
    # PVC (polyvinyl chloride), the actual material of bath ducks
    mu_rubber = 3.0e9       # Pa (shear modulus of rigid PVC)
    rho_rubber = 1.3e3      # kg/m^3
    g_earth = 9.8           # m/s^2
    R_duck = 0.04           # m

    k2_rubber = elastic_k2(mu_rubber, rho_rubber, g_earth, R_duck)
    rigidity = 19.0 * mu_rubber / (2.0 * rho_rubber * g_earth * R_duck)

    print("--- Rubber Duck (Elastic Love Number) ---")
    print(f"  mu_shear         = {mu_rubber:.1e} Pa")
    print(f"  rho              = {rho_rubber:.1f} kg/m^3")
    print(f"  g                = {g_earth:.1f} m/s^2")
    print(f"  R                = {R_duck:.3f} m")
    print(f"  19mu/(2 rho g R) = {rigidity:.2e}")
    print(f"  k_2              = {k2_rubber:.2e}")
    print()

    # -----------------------------------------------------------------------
    # 4. Comparison table
    # -----------------------------------------------------------------------
    print("=" * 72)
    print("  Comparison Table: k_2 and Lambda")
    print("=" * 72)
    print(f"{'Object':30s}  {'C':>10s}  {'k_2':>10s}  {'Lambda':>12s}")
    print("-" * 72)
    print(f"{'Nuclear duck (fiducial)':30s}  {C:10.4f}  {k2:10.5f}  {Lambda:12.1f}")
    print(f"{'Nuclear duck (range)':30s}  {'0.10-0.20':>10s}  {'0.06-0.10':>10s}  {'290-880':>10s}")
    print(f"{'Rubber duck (bath toy)':30s}  {'~1e-25':>10s}  {k2_rubber:10.2e}  {'---':>12s}")
    print(f"{'NS - SLy (1.4 Msun)':30s}  {'0.176':>10s}  {'0.091':>10s}  {'306':>12s}")
    print(f"{'NS - APR (1.4 Msun)':30s}  {'0.170':>10s}  {'0.085':>10s}  {'261':>12s}")
    print(f"{'NS - GW170817':30s}  {'~0.16':>10s}  {'---':>10s}  {'190+390-120':>12s}")
    print(f"{'Black hole':30s}  {'0.500':>10s}  {'0':>10s}  {'0':>12s}")
    print()

    # -----------------------------------------------------------------------
    # 5. k_2 vs compactness scan + plot
    # -----------------------------------------------------------------------
    print("Computing k_2 vs C scan...")
    results = scan_k2_vs_C()

    save_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results", "love_number_vs_C.png"
    )
    plot_k2_vs_C(results, save_path=save_path)

    print()
    print("--- k_2 vs C Scan ---")
    print(f"{'rho_c':>10s}  {'C':>8s}  {'M':>8s}  {'R':>8s}  {'k_2':>10s}  {'Lambda':>12s}")
    print("-" * 65)
    for i in range(len(results["rho_c"])):
        if np.isfinite(results["C"][i]):
            print(
                f"{results['rho_c'][i]:10.4e}  "
                f"{results['C'][i]:8.4f}  "
                f"{results['M'][i]:8.4f}  "
                f"{results['R'][i]:8.4f}  "
                f"{results['k2'][i]:10.5f}  "
                f"{results['Lambda'][i]:12.1f}"
            )

    print()
    print("Done.")


if __name__ == "__main__":
    main()
