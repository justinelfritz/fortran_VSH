#!/usr/bin/env python3
"""
Compute validation statistics from VSH test suite output files and write
LaTeX newcommand definitions to validation/validation_values.tex.

Run from the project root:
    python3 py/compute_validation.py

Or from the py/ directory:
    python3 compute_validation.py

After the test suite has been run to regenerate the validation data:
    cmake --build build && ctest --test-dir build
"""

import re
import os
import numpy as np

# ── path setup ────────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)
VDIR       = os.path.join(ROOT_DIR, 'validation')


# ── parsing helpers ───────────────────────────────────────────────────────────

CPAT = re.compile(r'\(([^,)]+),([^)]+)\)')


def load_plain(fname, vdir=VDIR):
    """
    Load a whitespace-delimited file into a list of string-token rows,
    skipping lines that start with '#'. `vdir` defaults to the ctest-gated
    validation/ directory, but callers (e.g. the convergence sweep, which
    lives entirely outside ctest) can point it elsewhere.
    """
    rows = []
    with open(os.path.join(vdir, fname)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append(line.split())
    return rows


def parse_test_file(fname):
    """
    Parse TEST1.dat / TEST2.dat.
    Format per line: theta  (re,im)_analytic  (re,im)_vsh  (re,im)_gw
    Returns four numpy arrays: theta, analytic, vsh_form, gw_form.
    """
    thetas, analytic, vsh_vals, gw_vals = [], [], [], []
    with open(os.path.join(VDIR, fname)) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Everything before the first '(' is theta
            first = line.index('(')
            theta = float(line[:first].strip())
            pairs = CPAT.findall(line)
            if len(pairs) < 3:
                continue
            cvals = [complex(float(r), float(i)) for r, i in pairs]
            thetas.append(theta)
            analytic.append(cvals[0])
            vsh_vals.append(cvals[1])
            gw_vals.append(cvals[2])
    return (np.array(thetas), np.array(analytic),
            np.array(vsh_vals), np.array(gw_vals))


# ── per-file statistics ───────────────────────────────────────────────────────

def cross_validation_stats(fname):
    """
    Max relative error for Tests 1 and 2.
    Interior points only (avoids pole regions where the function is zero).
    """
    theta, analytic, vsh_form, gw_form = parse_test_file(fname)
    interior = (theta > 0.05) & (theta < np.pi - 0.05)
    scale = np.abs(analytic[interior]).max()
    if scale == 0.0:
        return 0.0, 0.0
    vsh_rel = np.abs(analytic[interior] - vsh_form[interior]).max() / scale
    gw_rel  = np.abs(analytic[interior] - gw_form[interior]).max()  / scale
    return vsh_rel, gw_rel


def alp_consistency_pairs(fname, vdir=VDIR):
    """
    (L, abs_diff) per row from an ASSOC_LEGENDRE_ALL or
    DDX_ASSOC_LEGENDRE_ALL batch file.
    Format: L  M  X  batch_val  single_val  abs_diff
    """
    return [(int(row[0]), float(row[-1])) for row in load_plain(fname, vdir)]


def alp_consistency_stats(fname):
    return max(err for _, err in alp_consistency_pairs(fname))


def ssh_consistency_pairs(fname, vdir=VDIR):
    """
    (L, |batch - single|) per row from batch_ssh_cons.dat.
    Format: L  M  theta  phi  re_batch  im_batch  re_single  im_single
    """
    pairs = []
    for row in load_plain(fname, vdir):
        batch  = complex(float(row[4]), float(row[5]))
        single = complex(float(row[6]), float(row[7]))
        pairs.append((int(row[0]), abs(batch - single)))
    return pairs


def ssh_consistency_stats(fname):
    return max(err for _, err in ssh_consistency_pairs(fname))


def vsh_consistency_pairs(fname, vdir=VDIR):
    """
    (L, max component abs-difference) per row from any batch_*_cons.dat
    file. Format: L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph
    """
    return [(int(row[0]), max(float(row[-3]), float(row[-2]), float(row[-1])))
            for row in load_plain(fname, vdir)]


def vsh_consistency_stats(fname):
    return max(err for _, err in vsh_consistency_pairs(fname))


def max_by_degree(pairs, overflow_threshold=1e10):
    """
    Reduce a list of (L, err) pairs to (sorted L values, max err per L).
    A degree where any mode's error is non-finite (NaN/Inf) OR exceeds
    `overflow_threshold` is reported as NaN rather than an arbitrary max()
    result. The threshold matters because the unnormalized single-mode
    recurrence doesn't jump straight from "correct" to Inf as l, m near the
    diagonal grow -- it passes through a few degrees of astronomically large
    but still-finite values (e.g. ~1e291) on its way to overflow, and those
    values are just as meaningless as Inf while also being large enough to
    wreck a shared log-scale axis. Plotted with a plain line, NaN reads as a
    gap -- an honest "breaks down here" marker rather than a misleading
    number.
    """
    degrees = sorted(set(l for l, _ in pairs))
    errs = []
    for d in degrees:
        vals = [e for l, e in pairs if l == d]
        broken = any(not np.isfinite(v) or v > overflow_threshold for v in vals)
        errs.append(np.nan if broken else max(vals))
    return degrees, errs


def ssh_ortho_stats():
    """
    Max diagonal error |<Y_lm, Y_lm> - 1| and
    max off-diagonal error |<Y_lm, Y_l'm'>| from ssh_ortho.dat.
    Format: L1  M  L2  inner_product  expected
    """
    rows = load_plain('ssh_ortho.dat')
    diag, offdiag = [], []
    for row in rows:
        ip  = float(row[3])
        exp = int(row[4])
        (diag if exp == 1 else offdiag).append(abs(ip - exp))
    return max(diag), max(offdiag)


def pvsh_ortho_stats():
    """
    Max |DOT(PVSH_POL, PVSH_TOR)| from pvsh_pol_tor_ortho.dat.
    Format: L  M  theta  phi  |dot|
    """
    rows = load_plain('pvsh_pol_tor_ortho.dat')
    vals = [float(row[-1]) for row in rows]
    return max(vals) if vals else 0.0


def vsh_inversion_stats():
    """
    Max reconstruction errors from vsh_pol_inversion.dat.
    Format: L  M  theta  phi  max_absdiff_pol  max_absdiff_rad
    """
    rows = load_plain('vsh_pol_inversion.dat')
    pol = [float(row[-2]) for row in rows]
    rad = [float(row[-1]) for row in rows]
    return max(pol), max(rad)


def wigner_d_spotcheck_stats():
    """
    Max |Fortran WIGNER_D - sympy Rotation.D| over wigner_d_spotcheck.dat,
    an independent cross-check against a genuinely different
    library/algorithm (see py/sympy_reference.py's module docstring for
    why this is kept separate from py/mpmath_reference.py's role).
    Format: L  MP  M  alpha  beta  gamma  re(D)  im(D)

    Manual/investigative, like py/mpmath_reference.py's own self-check --
    not wired into main()'s validation_values.tex output, since this
    rotation infrastructure isn't part of the submitted manuscript.
    Requires sympy (see py/requirements.txt); import is local to this
    function so the rest of this module has no hard sympy dependency.
    """
    from sympy_reference import wigner_d_ref

    rows = load_plain('wigner_d_spotcheck.dat')
    max_err = 0.0
    for row in rows:
        l, mp, m = int(row[0]), int(row[1]), int(row[2])
        alpha, beta, gamma = float(row[3]), float(row[4]), float(row[5])
        fortran_val = complex(float(row[6]), float(row[7]))
        sympy_val = wigner_d_ref(l, mp, m, alpha, beta, gamma)
        max_err = max(max_err, abs(fortran_val - sympy_val))
    return max_err


# ── LaTeX formatting ──────────────────────────────────────────────────────────

def latex_sci(val):
    """
    Format a non-negative float as compact LaTeX scientific notation.
    Values below 1e-18 are reported as the double-precision floor.
    """
    if val == 0.0:
        return r'$0$'
    if val < 1e-18:
        return r'$<\!10^{-18}$'
    exp  = int(np.floor(np.log10(val)))
    mant = val / 10.0**exp
    if abs(mant - 1.0) < 0.06:
        return f'$10^{{{exp}}}$'
    return f'${mant:.1f}\\times10^{{{exp}}}$'


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    stats = {}

    # Cross-validation (Tests 1 and 2)
    t1_vsh, t1_gw = cross_validation_stats('TEST1.dat')
    t2_vsh, t2_gw = cross_validation_stats('TEST2.dat')
    stats['TestOneVSH'] = t1_vsh
    stats['TestOneGW']  = t1_gw
    stats['TestTwoVSH'] = t2_vsh
    stats['TestTwoGW']  = t2_gw

    # ALP batch consistency
    stats['ALMmax']  = alp_consistency_stats('batch_alm_cons.dat')
    stats['DALMmax'] = alp_consistency_stats('batch_dalm_cons.dat')

    # SSH batch consistency
    stats['SSHmax'] = ssh_consistency_stats('batch_ssh_cons.dat')

    # VSH batch consistency (max over r, theta, phi components)
    vsh_files = [
        ('GradSSHmax',  'batch_grad_ssh_cons.dat'),
        ('LSHmax',      'batch_l_ssh_cons.dat'),
        ('PVSHRadmax',  'batch_pvsh_rad_cons.dat'),
        ('PVSHPolmax',  'batch_pvsh_pol_cons.dat'),
        ('PVSHTormax',  'batch_pvsh_tor_cons.dat'),
        ('VSHTormax',   'batch_vsh_tor_cons.dat'),
        ('VSHUpmax',    'batch_vsh_pol_up_cons.dat'),
        ('VSHDnmax',    'batch_vsh_pol_dn_cons.dat'),
    ]
    for key, fname in vsh_files:
        stats[key] = vsh_consistency_stats(fname)

    # Orthogonality tests
    diag_err, offdiag_err = ssh_ortho_stats()
    stats['SSHOrthoDiag']    = diag_err
    stats['SSHOrthoOffdiag'] = offdiag_err
    stats['PVSHOrtho']       = pvsh_ortho_stats()

    # VSH_POL inversion
    pol_err, rad_err = vsh_inversion_stats()
    stats['VSHInvPol'] = pol_err
    stats['VSHInvRad'] = rad_err

    # Print summary
    print('=== Validation Statistics ===\n')
    width = max(len(k) for k in stats)
    for k, v in stats.items():
        print(f'  {k:{width}s} : {v:.3e}')

    # Wigner-D vs. sympy cross-check -- informational only (see
    # wigner_d_spotcheck_stats' docstring for why this stays out of
    # validation_values.tex/`stats` above rather than joining them).
    try:
        wigd_err = wigner_d_spotcheck_stats()
        print(f'\nWigner-D vs. sympy cross-check: max abs err = {wigd_err:.3e}')
    except ImportError:
        print('\n(sympy not installed -- skipping Wigner-D vs. sympy cross-check)')

    # Write validation_values.tex
    outpath = os.path.join(VDIR, 'validation_values.tex')
    with open(outpath, 'w') as fh:
        fh.write('% Auto-generated by py/compute_validation.py\n')
        fh.write('% Re-run after: cmake --build build && ctest --test-dir build\n\n')
        for k, v in stats.items():
            fh.write(f'\\newcommand{{\\val{k}}}{{{latex_sci(v)}}}\n')
    print(f'\nWrote {outpath}')


if __name__ == '__main__':
    main()
