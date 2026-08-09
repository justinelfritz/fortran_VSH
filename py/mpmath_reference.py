#!/usr/bin/env python3
"""
Independent, arbitrary-precision reference implementation of the five
routine families plotted in the convergence_accuracy figure (Legendre, SSH,
VSH toroidal, VSH poloidal up/down), plus (added later, for the rotation
infrastructure's own stability investigation -- see STABILITY_FINDINGS.md
Part 3) the Wigner small-d rotation matrix.

This exists because src/vsh.f90's "single-mode" functions (SSH, VSH_TOR,
VSH_POL_UP, VSH_POL_DN) all normalize the *unnormalized* ASSOC_LEGENDRE
Bonnet recurrence after the fact, while the batch "_ALL" functions use the
separate Holmes-Featherstone ASSOC_LEGENDRE_NORM_ALL recurrence that folds
normalization into every step. Diffing those two Fortran code paths against
each other (the old convergence_accuracy figure) only shows where the
single-mode path degrades -- it says nothing about whether either path is
actually close to the true value.

Importantly, `legendre_ref` (unnormalized) is NOT the right reference for
ASSOC_LEGENDRE_NORM_ALL specifically -- normalization is the actual
stabilizing feature of that routine (it's why it stays accurate to
l~2700 instead of ~1400), so validating it against an unnormalized
reference would defeat the point. `legendre_norm_ref` applies the same
normalization to the reference side for exactly this comparison.

The unnormalized P_l^m is computed here by `_legendre_recurrence`: the same
three-term Bonnet recurrence as src/vsh.f90's ASSOC_LEGENDRE, but evaluated
at mpmath's arbitrary precision (50 decimal digits, unbounded exponent)
instead of double precision. mpmath.legenp (its built-in associated
Legendre function, hypergeometric-series based) was tried first but fails
to converge for extreme (l, m) near the diagonal at l~300 -- a limitation
of that specific series-evaluation method, not of the underlying math. The
recurrence itself is mathematically exact; the Fortran double-precision
version only breaks down from floating-point effects (overflow past
~1.8e308, or accumulated rounding) that don't exist at 50-digit precision
with unbounded exponent, so re-running the *identical* formula at high
precision isolates exactly that effect and remains a legitimate
independent ground truth for "how far can double precision be trusted."
The theta-derivative for the SSH/VSH families uses mpmath.diff (numerical
differentiation) rather than a second closed-form recurrence.

The five public functions mirror the combination formulas documented in
src/vsh.f90 (VSH_TOR = PVSH_TOR, VSH_POL_UP, VSH_POL_DN) exactly, since
those are definitions of the VSH basis, not something to re-derive
independently -- only the underlying P_l^m/Y_l^m/dY_l^m evaluation is
independent of the Fortran source.

Run standalone to self-check against known closed forms before trusting it
as ground truth elsewhere:
    python3 py/mpmath_reference.py
"""

import mpmath as mp

mp.mp.dps = 50


def _legendre_recurrence(l, m, x):
    """
    Unnormalized P_l^m(x) via the same Bonnet recurrence as
    src/vsh.f90's ASSOC_LEGENDRE_DP (src/vsh.f90:309-354), evaluated at
    mpmath's working precision. m is assumed >= 0 (all sample points used
    by py/plot_convergence.py sample non-negative m, matching how
    BATCH_ALM_CONS/BATCH_SSH_CONS exercise the Fortran side).
    """
    if m > l or abs(x) > 1:
        return mp.mpf(0)
    somx2 = mp.sqrt((1 - x) * (1 + x))
    pmm = mp.mpf(1)
    for i in range(1, m + 1):
        pmm = -pmm * (2 * i - 1) * somx2
    if l == m:
        return pmm
    pmm1 = x * (2 * m + 1) * pmm
    if l == m + 1:
        return pmm1
    plm = pmm1
    for i in range(m + 2, l + 1):
        plm = (x * (2 * i - 1) * pmm1 - (i + m - 1) * pmm) / (i - m)
        pmm = pmm1
        pmm1 = plm
    return plm


def _y_and_grad(l, m, theta, phi):
    """
    (Y_l^m, dY/dtheta, (1/sin(theta))*dY/dphi) at arbitrary precision.
    Y_l^m via `_legendre_recurrence` (Condon-Shortley, matches
    src/vsh.f90); the theta-derivative via mpmath.diff (numerical,
    independent of any closed-form derivative recurrence).
    """
    def Y(th):
        x = mp.cos(th)
        norm = mp.sqrt((2 * l + 1) / (4 * mp.pi)) * \
            mp.sqrt(mp.factorial(l - m) / mp.factorial(l + m))
        return norm * _legendre_recurrence(l, m, x) * mp.exp(1j * m * phi)

    y = Y(theta)
    gth = mp.diff(Y, theta)
    gph = 1j * m * y / mp.sin(theta)
    return y, gth, gph


def legendre_ref(l, m, x):
    """Unnormalized P_l^m(x), independent of ASSOC_LEGENDRE(_ALL)."""
    return complex(_legendre_recurrence(l, abs(m), mp.mpf(x)))


def legendre_norm_ref(l, m, x):
    """
    N_l^m * P_l^m(x), independent of ASSOC_LEGENDRE_NORM_ALL -- the
    normalization is the actual stabilizing feature of that routine (folds
    into every recurrence step, keeps values O(1) instead of letting the
    unnormalized P_l^m grow/shrink combinatorially), so it must be applied
    on the reference side too rather than comparing unnormalized values
    only. Same N_l^m = sqrt((2l+1)/4pi * (l-m)!/(l+m)!) as ASSOC_LEGENDRE_NORM_ALL
    documents (src/vsh.f90:1220-1230) and as SSH's own normalization.
    """
    m = abs(m)
    x = mp.mpf(x)
    norm = mp.sqrt((2 * l + 1) / (4 * mp.pi)) * \
        mp.sqrt(mp.factorial(l - m) / mp.factorial(l + m))
    return float(norm * _legendre_recurrence(l, m, x))


def ssh_ref(l, m, theta, phi):
    """Y_l^m(theta, phi), independent of SSH/SSH_ALL."""
    y, _, _ = _y_and_grad(l, m, mp.mpf(theta), mp.mpf(phi))
    return complex(y)


def vsh_tor_ref(l, m, theta, phi):
    """VSH_TOR(l,m,theta,phi) = PVSH_TOR; zero at l=0 (see src/vsh.f90)."""
    if l == 0:
        return (0j, 0j, 0j)
    _, gth, gph = _y_and_grad(l, m, mp.mpf(theta), mp.mpf(phi))
    scale = 1 / mp.sqrt(l * (l + 1))
    return (0j, complex(1j * gph * scale), complex(-1j * gth * scale))


def vsh_pol_dn_ref(l, m, theta, phi):
    """VSH_POL_DN(l,m,theta,phi); zero at l=0 (see src/vsh.f90)."""
    if l == 0:
        return (0j, 0j, 0j)
    y, gth, gph = _y_and_grad(l, m, mp.mpf(theta), mp.mpf(phi))
    sc_r = mp.sqrt(mp.mpf(l) / (2 * l + 1))
    sc_th = 1 / mp.sqrt((2 * l + 1) * l)
    return (complex(sc_r * y), complex(sc_th * gth), complex(sc_th * gph))


def vsh_pol_up_ref(l, m, theta, phi):
    """VSH_POL_UP(l,m,theta,phi); nonzero radial part at l=0 (see src/vsh.f90)."""
    y, gth, gph = _y_and_grad(l, m, mp.mpf(theta), mp.mpf(phi))
    sc_r = -mp.sqrt((l + 1) / mp.mpf(2 * l + 1))
    if l == 0:
        return (complex(sc_r * y), 0j, 0j)
    sc_th = 1 / mp.sqrt((2 * l + 1) * (l + 1))
    return (complex(sc_r * y), complex(sc_th * gth), complex(sc_th * gph))


def wigner_d_small_ref(l, mp_, m, beta):
    """
    Wigner small-d matrix element d^l_{m'm}(beta), independent of
    src/vsh.f90's WIGNER_D_SMALL_DP. Same closed-form finite-sum formula
    (Wigner's explicit sum -- see that function's docstring), but
    evaluated via mpmath's *direct* arbitrary-precision factorials
    (mp.factorial) rather than WIGNER_D_SMALL_DP's log-factorial/
    incremental-term-update technique -- that technique exists purely to
    dodge double-precision overflow, which mpmath's unbounded exponent
    range doesn't need, so this is a genuinely different code path, not
    a higher-precision rerun of the same arithmetic sequence (contrast
    `_legendre_recurrence`, which deliberately *does* replay the same
    recurrence at higher precision -- see that function's docstring for
    why that distinction matters there).

    Used by py/plot_wigner_d_error.py as ground truth for characterizing
    WIGNER_D_SMALL_DP's catastrophic-cancellation accuracy limit at
    large l (see STABILITY_FINDINGS.md, Part 3).
    """
    beta = mp.mpf(beta)
    coshb = mp.cos(beta / 2)
    sinhb = mp.sin(beta / 2)
    prefactor = mp.sqrt(mp.factorial(l + m) * mp.factorial(l - m) *
                         mp.factorial(l + mp_) * mp.factorial(l - mp_))
    kmin = max(0, m - mp_)
    kmax = min(l + m, l - mp_)
    total = mp.mpf(0)
    for k in range(kmin, kmax + 1):
        denom = (mp.factorial(l + m - k) * mp.factorial(k) *
                  mp.factorial(l - mp_ - k) * mp.factorial(mp_ - m + k))
        term = prefactor / denom * coshb**(2 * l + m - mp_ - 2 * k) * \
            sinhb**(mp_ - m + 2 * k)
        total += term if k % 2 == 0 else -term
    return float(total)


if __name__ == '__main__':
    checks = []

    checks.append(('P_1^0(0.6)', legendre_ref(1, 0, 0.6), 0.6))
    checks.append(('P_1^1(0.6)', legendre_ref(1, 1, 0.6),
                   -mp.sqrt(1 - 0.6**2)))
    checks.append(('Y_0^0', ssh_ref(0, 0, 0.9, 0.3),
                   1 / mp.sqrt(4 * mp.pi)))
    checks.append(('Y_1^0(theta=0.9)', ssh_ref(1, 0, 0.9, 0.0),
                   mp.sqrt(3 / (4 * mp.pi)) * mp.cos(0.9)))
    checks.append(('N_1^0*P_1^0(0.6)', legendre_norm_ref(1, 0, 0.6),
                   mp.sqrt(3 / (4 * mp.pi)) * 0.6))
    checks.append(('d^0_00(0.7)', wigner_d_small_ref(0, 0, 0, 0.7), 1.0))
    checks.append(('d^1_00(0.7)', wigner_d_small_ref(1, 0, 0, 0.7),
                   mp.cos(0.7)))
    checks.append(('d^1_10(0.7)', wigner_d_small_ref(1, 1, 0, 0.7),
                   mp.sin(0.7) / mp.sqrt(2)))

    ok = True
    for name, got, expect in checks:
        err = abs(complex(got) - complex(expect))
        status = 'OK' if err < 1e-30 else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f'{status}  {name:20s} got={got}  expect={expect}  err={err:.2e}')

    print('\nAll checks passed.' if ok else '\nSOME CHECKS FAILED.')
