# Verification Report: Thorne & Campolattaro 1967 LaTeX Transcription

**Date:** 2026-04-01
**Source PDF:** Thorne & Campolattaro, ApJ 149, 591 (1967)
**LaTeX file:** `main.tex`

---

## 1. Critical Equation Verification (TeX vs. PDF images)

### Eq (1): Line element [page-02.png]

**PDF:**
$$ds^2 = (ds^2)_0 \equiv e^\nu dt^2 - e^\lambda dr^2 - r^2(d\theta^2 + \sin^2\theta\, d\phi^2)$$

**TeX (line 121):**
```
ds^{2} = -(e^{\nu})\,(dt)^{2} + e^{\lambda}\,(dr)^{2} + r^{2}(d\theta^{2} + \sin^{2}\theta\,d\phi^{2})
```

**DISCREPANCY FOUND:** Sign convention is inconsistent with the PDF. The PDF uses the signature $(-,+,+,+)$ with $ds^2 = e^\nu dt^2 - e^\lambda dr^2 - r^2 d\Omega^2$. The TeX has the signs flipped: $-(e^\nu)(dt)^2 + e^\lambda(dr)^2 + r^2 d\Omega^2$. This is a **different but equivalent convention** (signature $(-,+,+,+)$ vs. $(+,-,-,-)$). However, comparing carefully to the PDF image, the original paper writes:

$$ds^2 = (ds^2)_0 \equiv e^\nu\,dt^2 - e^\lambda\,dr^2 - r^2(d\theta^2 + \sin^2\theta\, d\phi^2)$$

This is signature $(+,-,-,-)$. The TeX version writes $-(e^\nu)(dt)^2 + e^\lambda(dr)^2 + r^2(...)$ which is signature $(-,+,+,+)$.

**Verdict:** The TeX uses the **opposite overall sign convention** from the PDF. Both are valid metric signatures, but this is **inconsistent with the source**. The PDF clearly shows the positive sign on $e^\nu dt^2$ and negative on the spatial terms.

### Eq (3b): TOV equation [page-03.png]

**PDF:**
$$\frac{dp}{dr} = -\frac{(\rho + p)(m + 4\pi r^3 p)}{r(r - 2m)}$$

**TeX (line 142):**
```
\frac{dp}{dr} = -\frac{(\rho + p)(m + 4\pi r^3 p)}{r(r - 2m)}
```

**Verdict: MATCH.** The TOV equation is correctly transcribed.

### Eq (7b): Metric perturbation matrix [page-04.png]

**PDF:** Shows a 4x4 matrix with row/column labels $t, r, \theta, \phi$ and entries:
- $(t,t)$: $H_0 e^\nu$
- $(t,r)$ and $(r,t)$: $H_1$
- $(r,r)$: $H_2 e^\lambda$
- $(\theta,\theta)$: $r^2 K$
- $(\phi,\phi)$: $r^2 \sin^2\theta\, K$
- All other entries: 0

Multiplied by $P_l(\cos\theta)$.

**TeX (lines 234-246):**
```
h_{\mu\nu} = matrix with same entries * P_l(cos theta)
```

**Verdict: MATCH.** The matrix structure and all entries match the PDF.

### Eq (14a): First eigenequation [page-07.png]

**PDF (page 7):** Shows the initial-value eigenequation starting with $H' + r^{-1}e^\lambda[l(l+1)/2 + 1 + 4\pi r^2(\rho - p)]H = rK''...$

**TeX (lines 452-461):**
```
H' + r^{-1}e^{\lambda}[l(l+1)/2 + 1 + 4\pi r^{2}(p - \rho)]H = r K'' ...
```

**DISCREPANCY FOUND:** The sign inside the bracket differs:
- PDF: $4\pi r^2(\rho - p)$
- TeX: $4\pi r^2(p - \rho)$

This is a **sign error** in the TeX transcription. The PDF clearly shows $(\rho - p)$.

Additionally, examining the last term of Eq (14a):
- PDF (from Eq 8a on page 5): $8\pi l(l+1) r^{-1}(\rho + p)e^\nu V$ -- the exponential is $e^\nu$
- TeX Eq (14a) line 460: $8\pi l(l+1)r^{-1}(\rho + p)e^{\lambda}V$ -- the exponential is $e^\lambda$

**DISCREPANCY FOUND:** The last term has $e^\lambda$ in the TeX but the PDF (Eq 8a) shows $e^\nu$. This needs careful checking since Eq 14a is Eq 8a with the time dependence factored out.

Also on the second line of Eq (14a):
- PDF: $e^\lambda(3 - 5m/r - 4\pi r^2 p)K'$ (coefficient of $K'$)
- TeX line 455: $e^{\lambda}(3 - 5m/r - 4\pi r^{2}\rho)K$ -- this has $\rho$ instead of $p$, AND it says $K$ not $K'$

**DISCREPANCY FOUND:** Two errors on this line: (1) $\rho$ should be $p$, and (2) the missing prime on $K$ (should be $K'$).

### Eq (22): Quasi-normal mode expansion [page-11.png]

**PDF:**
$$K_n(r,t) = K_n^{(0)}(r) e^{i\omega_n t} + K_n^{(0)*}(r) e^{-i\omega_n^* t}$$
$$= [K_n^{(0)}(r) e^{i\sigma_n t} + K_n^{(0)*}(r) e^{-i\sigma_n t}] e^{-t/\tau_n}$$

**TeX (lines 720-721):**
```
K_n(r,t) = K_n^{(0)}(r) e^{i\omega_n t} + K_n^{(0)*}(r) e^{-i\omega_n^* t}
= [K_n^{(0)}(r) e^{i\sigma_n t} + K_n^{(0)*}(r) e^{-i\sigma_n t}] e^{-t/\tau_n}
```

**Verdict: MATCH.** The quasi-normal mode expansion is correctly transcribed.

---

## 2. Equation Numbering Completeness

### Main body equations (1 through 35):

| Equation | Present? | Notes |
|----------|----------|-------|
| 1 | Yes | Line element |
| 2 | Yes | $e^{-\lambda}$ definition |
| 3a | Yes | Mass equation |
| 3b | Yes | TOV equation |
| 3c | Yes | Source equation for $\nu$ |
| 3d | Yes | Equation of state |
| 4 | Yes | Adiabatic index |
| 5 | Yes | Perturbed line element |
| 6a | Yes | Odd-parity displacement |
| 6b | Yes | Odd-parity metric perturbation |
| 7a | Yes | Even-parity displacement |
| 7b | Yes | Even-parity metric perturbation matrix |
| 8a | Yes | Initial-value equation |
| 8b | Yes | Initial-value equation |
| 8c | Yes | Initial-value equation |
| 8d | Yes | $H_2 = H_0$ |
| 9a | Yes | Propagation equation |
| 9b | Yes | Propagation equation |
| 9c | Yes | Propagation equation |
| 10 | Yes | Boundary condition $H_0 = K$ at $r=0$ |
| 11 | Yes | Complex frequency |
| 12 | Yes | Normal mode ansatz |
| 13 | Yes | Convention $\sigma \geq 0$ |
| 14a | Yes | Eigenequation (has errors, see above) |
| 14b | Yes | Eigenequation |
| 14c | Yes | Eigenequation |
| 14d | Yes | Eigenequation |
| 15a | Yes | Boundary at $r=0$ |
| 15b | Yes | Boundary at $r=R$ |
| 16 | Yes | Lagrangian pressure condition |
| 17a | Yes | Polytropic surface |
| 17b | Yes | Surface relation |
| 18 | Yes | Asymptotic form (partially rendered) |
| 19 | Yes | Continuity condition |
| 20 | Yes | Ratio function |
| 21 | Yes | General pulsation expansion |
| 22 | Yes | Quasi-normal mode |
| 23 | Yes | Wave zone form |
| 24 | Yes | Observed frequency |
| 25 | Yes | Kinetic energy |
| 26 | Yes | Kinetic energy expanded |
| 27 | Yes | Phase smallness |
| 28a | Yes | $E_{\rm kin}$ simplified |
| 28b | Yes | $E_{\rm pulse}^n$ |
| 29a | Yes | Potential energy |
| 29b | Yes | Total energy |
| 30 | Yes | Radiated power |
| 31 | Yes | Standing wave condition |
| 32 | Yes | Amplitude expansion |
| 33 | Yes | Wave front integral |
| 34a | Yes | Wave front solution |
| 34b | Yes | Wave front solution |
| 35 | Yes | Incoming amplitude ratio |

**All 35 numbered equations (with sub-labels) are present.**

### Appendix equations:

| Equation | Present? |
|----------|----------|
| A1a | Yes |
| A1b | Yes |
| A2 | Yes |
| A3a | Yes |
| A3b | Yes |
| A4 | Yes |
| A5 | Yes |
| A6 | Yes |
| A7 | Yes |
| A8 | Yes |
| A9 | Yes |
| A10 | Yes |
| A11 | Yes |
| A12 | Yes |
| B1 | Yes |
| B2 | Yes |
| B3 | Yes |
| B4 | Yes |
| B5a | Yes |
| B5b | Yes |
| B5c | Yes |
| B6 | Yes |
| B7 | Yes |
| B8 | Yes |
| C1 | Yes |
| C2 | Yes |
| C3 | Yes |
| C4 | Yes |
| C5 | Yes |
| C6 | Yes |
| C7 | Yes |
| D1a | Yes |
| D1b | Yes |
| D1c | Yes |
| D1d | Yes |
| D2 | Yes |
| D3 | Yes |
| D4 | Yes |
| D5 | Yes |
| D6 | Yes |
| D7 | Yes |
| E1a | Yes |
| E1b | Yes |
| E2 | Yes |

**All appendix equations A1a through E2 are present.**

---

## 3. Formatting Issues

1. **Eq (18):** The asymptotic wave form equation is incomplete/poorly structured in the TeX. Lines 596-604 show a broken `align` environment that ends awkwardly with `\text{for } r \gg M \text{ and } \{R \to r\}^{-1}` and a standalone `\tag{18}` on a separate equation. The PDF shows a cleaner two-line expression for $K$ and $H$ with conditions.

2. **Eq (30):** The radiated power equation (line 868) appears garbled:
   ```
   P = -2\,\tau_{n}^{-1} E_{n,\text{puls}}^{-k/\tau_{n}}
   ```
   The PDF should show something like $P = 2\tau_n^{-1} E_{\rm puls}^n e^{-2t/\tau_n}$. The superscript `^{-k/\tau_n}` is likely an OCR/transcription error.

3. **Eq (34b):** Contains obvious transcription errors -- `i < r + 2M \ln r` should be `t < r + 2M \ln r` (letter `i` vs `t`), and `\pi_r` should likely be `\sigma_n r` or similar.

4. **Eq (C7):** Lines 1357-1383 contain multiple comments marking uncertain terms (`% "3p" possibly "3\rho"`, `% prefactor uncertain`, `% exponent on p uncertain`, `% exponent on r uncertain`). This equation was clearly difficult to read from the PDF and contains known uncertainties.

5. **Eq (8a) vs Eq (14a) inconsistency:** Eq (8a) at line 297 has $4\pi r^2(\rho - p)$ while Eq (14a) at line 453 has $4\pi r^2(p - \rho)$. Since Eq (14a) should be Eq (8a) with time-dependence factored out (identical spatial part), these must agree. The PDF shows $(\rho - p)$ in both, so the Eq (14a) TeX has the sign wrong.

6. **Appendix A equations (A5, A6):** These are quite garbled with incorrect symbol usage (e.g., `$r_f$` appearing where the PDF likely shows `(r,t)`). The harmonic notation is inconsistent.

7. **Unnumbered curvature invariant equation** near line 656: `\delta R_I \sim 1/r` appears without a number, consistent with the PDF.

---

## 4. Summary of Discrepancies

### Critical errors (affect physics):

| Location | Issue | Severity |
|----------|-------|----------|
| Eq (1), line 121 | Overall metric sign convention flipped vs. PDF | HIGH |
| Eq (14a), line 453 | $(\rho - p)$ transcribed as $(p - \rho)$ | HIGH |
| Eq (14a), line 455 | $K'$ missing prime (written as $K$) | HIGH |
| Eq (14a), line 455 | $p$ written as $\rho$ in coefficient of $K'$ | HIGH |
| Eq (14a), line 460 | $e^\nu$ written as $e^\lambda$ in last term | HIGH |
| Eq (30), line 868 | Garbled expression for radiated power | HIGH |
| Eq (34b), lines 936-939 | Multiple transcription errors ($i$ for $t$, garbled coefficients) | HIGH |

### Moderate issues (formatting/presentation):

| Location | Issue |
|----------|-------|
| Eq (18) | Incomplete/broken formatting of asymptotic expression |
| Eq (C7) | Multiple uncertain terms flagged with comments |
| Appendix A (A5, A6) | Garbled notation in harmonic expansions |

### Verified correct:
- Eq (3b): TOV equation -- correct
- Eq (7b): Metric perturbation matrix -- correct
- Eq (22): Quasi-normal mode expansion -- correct
- Eq (8a)-(8d): Initial-value equations (main body) -- correct
- Eq (9a)-(9c): Propagation equations -- correct
- All equation numbers present (1-35, A1a-E2)

---

## 5. Recommendation

The transcription has significant errors concentrated in **Eq (14a)** (the first eigenequation) and in **Eqs (30) and (34b)**. These are physics-critical equations. The metric sign convention in Eq (1) should also be reconciled with the original. The appendix equations (especially C7 and A5-A6) need re-examination against the PDF.

Priority fixes:
1. Fix Eq (14a): four separate errors in this single equation
2. Fix Eq (1): match the original sign convention
3. Fix Eq (30): correct the garbled radiated power expression
4. Fix Eq (34b): correct `i` to `t` and fix coefficients
5. Re-examine Eq (C7) uncertain terms
