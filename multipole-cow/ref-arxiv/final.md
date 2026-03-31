# Reference Papers for "Quacky Normal Modes" -- Consolidated Summary

This directory contains arXiv source and markdown conversions of three key reference papers for the duck QNM project.

## Papers

### 1. Kokkotas & Schmidt (1999) -- `9909058/`

**"Quasi-Normal Modes of Stars and Black Holes"**, Living Rev. Rel. 2, 2
[arXiv:gr-qc/9909058](https://arxiv.org/abs/gr-qc/9909058)

**What it provides:**
- Complete mode classification: f, p, g, w, r-modes with frequency/damping scalings
- Perturbation equations for relativistic stars (polar and axial)
- The f-mode frequency formula: omega^2 = 2*l*(l-1)/(2l+1) * M/R^3 (Newtonian uniform density)
- Cowling and inverse Cowling approximations
- Table of typical QNM frequencies and damping times

**Key equation for us:** The f-mode frequency ~ sqrt(M/R^3) scaling is the "quack frequency" baseline.

---

### 2. Andersson & Kokkotas (1998) -- `9711088/`

**"Towards gravitational-wave asteroseismology"**, MNRAS 299, 1059
[arXiv:gr-qc/9711088](https://arxiv.org/abs/gr-qc/9711088)

**What it provides:**
- Universal (EOS-independent) fitting formulas for f-mode, p-mode, and w-mode frequencies and damping times as functions of (M, R)
- Parameter estimation accuracy: combining f + w mode data yields (M, R) to (2%, 5%)
- Framework for the inverse problem: observed QNMs -> stellar parameters

**Key equations for us:**
- f-mode: omega_f (kHz) = 0.78 + 1.635 * sqrt(M_bar / R_bar^3)
- f-mode damping: 1/tau_f (s) = (M_bar^3/R_bar^4) * [22.85 - 14.65*(M_bar/R_bar)]
- w-mode: omega_w (kHz) = (1/R_bar) * [20.92 - 9.14*(M_bar/R_bar)]

These are the spherical baselines. Duck shape splits and shifts these frequencies.

---

### 3. Hinderer (2008) -- `0711.2420/`

**"Tidal Love numbers of neutron stars"**, ApJ 677, 1216
[arXiv:0711.2420](https://arxiv.org/abs/0711.2420)

**What it provides:**
- Definition of tidal Love number k_2: Q_ij = -lambda * E_ij, with k_2 = (3/2) G lambda R^{-5}
- The h_2 ODE for static perturbations inside a relativistic star
- Closed-form k_2 expression in terms of compactness C = M/R and boundary value y = R*H'(R)/H(R)
- Newtonian limit: k_2^N = (2-y)/(2(y+3))
- Numerical tables and fitting formula for polytropic stars

**Key equation for us:** The k_2 formula (Eq. 20) generalizes to a tensor for non-spherical stars. The duck's anisotropic shape creates direction-dependent tidal response.

---

## How These Connect to the Duck Paper

The three papers form a hierarchy:

1. **Kokkotas & Schmidt** provides the theoretical framework (perturbation theory, mode classification)
2. **Andersson & Kokkotas** provides the quantitative baseline (universal relations for spherical stars)
3. **Hinderer** provides the tidal deformability framework (Love numbers)

For the duck-shaped star:
- The spherical QNM frequencies from (1) and (2) get **split** by the duck multipole moments, analogous to Zeeman splitting
- The scalar Love number from (3) becomes a **tensor**, with the duck's shape introducing anisotropic tidal response
- The "melting" transition (duck -> sphere) corresponds to the splitting going to zero and the Love tensor becoming isotropic

## Notes on arXiv IDs

The user-provided arXiv IDs required correction:
- "9904024" -> actual paper is **gr-qc/9909058** (Kokkotas & Schmidt)
- "9706075" -> actual paper is **gr-qc/9711088** (Andersson & Kokkotas)
- "0805.3235" -> actual paper is **0711.2420** (Hinderer)
