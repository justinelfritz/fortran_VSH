#!/usr/bin/env python3
"""
Plot the convergence_accuracy figure: max absolute error against an
INDEPENDENT arbitrary-precision reference (py/mpmath_reference.py) for the
NORMALIZED batch ("_ALL") routines -- ASSOC_LEGENDRE_NORM_ALL, SSH_ALL,
VSH_TOR_ALL, VSH_POL_UP_ALL, VSH_POL_DN_ALL -- restricted to their verified
safe range l<=2000 (see STABILITY_FINDINGS.md, Part 1).

This replaces an earlier version of this figure that mixed the normalized
batch routines together with the *unnormalized* single-mode/batch code
paths in the same panels. That mixing caused two problems, both resolved
by splitting into two figures:
  1. The unnormalized series needed relative error (their raw magnitude
     grows combinatorially even when accurate -- see
     plot_convergence_unnormalized.py), while the normalized series are
     fine with absolute error (bounded O(1) always) -- one metric can't
     honestly serve both.
  2. The unnormalized series' much lower safe ceiling (l<=150) forced an
     awkward y-axis cap and, before that, an outright crash from the raw
     magnitude reaching ~1e300+ (see git history) -- none of which has
     anything to do with whether the normalized routines (the ones this
     figure is actually about) are trustworthy.
See plot_convergence_unnormalized.py for the unnormalized companion figure.

Only ONE series per panel now (no single-mode comparison): there is no
single-mode counterpart to ASSOC_LEGENDRE_NORM_ALL, and single-mode
SSH/VSH_TOR/VSH_POL_UP/VSH_POL_DN are built on the *unnormalized*
recurrence (see their docstrings in src/vsh.f90) -- plotting them here
would misleadingly imply they share this figure's l<=2000 safe range, when
they're actually limited to l<=150 (that's exactly what the unnormalized
companion figure shows instead).

Reads convergence/reference_*.dat (produced by the opt-in vsh_convergence
executable), filtered to l<=2000 -- the sparse sample grid already extends
to l=4000 (for the investigation that led to STABILITY_FINDINGS.md), but
l=2000-4000 is where ASSOC_LEGENDRE_NORM_ALL's theta-dependent instability
can appear (Part 1), so this figure deliberately stops at the verified-safe
ceiling rather than showing that cliff -- STABILITY_FINDINGS.md and its own
heatmap (py/plot_stability.py) already document it thoroughly.

Run from the project root, after generating the sweep data:
    cmake -B build -DVSH_BUILD_CONVERGENCE=ON && cmake --build build
    ./build/vsh_convergence
    python3 py/plot_convergence.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import mpmath_reference as mpref
from compute_validation import max_by_degree
from plotstyle import BLUE, INK_PRIMARY, INK_SECONDARY, SURFACE, \
    style_axes, savefig_pair, ROOT_DIR

CONVERGENCE_DIR = os.path.join(ROOT_DIR, 'convergence')
SAFE_LMAX = 2000


def _safe_abs(diff):
    """
    abs() that treats OverflowError as +inf instead of crashing. Observed
    empirically at extreme (l, m) in earlier, wider sweeps -- kept here
    defensively even though this figure's l<=2000 filter should stay well
    clear of it (see plot_convergence_unnormalized.py's docstring for
    where this was first hit).
    """
    try:
        return abs(diff)
    except OverflowError:
        return float('inf')


def _load_rows(fname):
    rows = []
    with open(os.path.join(CONVERGENCE_DIR, fname)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append(line.split())
    return rows


def _alm_pairs():
    """(L, err) pairs for batch ASSOC_LEGENDRE_NORM_ALL vs. the equally-normalized reference."""
    pairs = []
    for row in _load_rows('reference_alm.dat'):
        l, m, x = int(row[0]), int(row[1]), float(row[2])
        if l > SAFE_LMAX:
            continue
        p_batch_norm = float(row[5])
        ref_norm = mpref.legendre_norm_ref(l, m, x)
        pairs.append((l, _safe_abs(p_batch_norm - ref_norm)))
    return pairs


def _ssh_pairs():
    """(L, err) pairs for batch SSH_ALL."""
    pairs = []
    for row in _load_rows('reference_ssh.dat'):
        l, m, theta, phi = int(row[0]), int(row[1]), float(row[2]), float(row[3])
        if l > SAFE_LMAX:
            continue
        b_re, b_im = float(row[6]), float(row[7])
        ref = mpref.ssh_ref(l, m, theta, phi)
        pairs.append((l, _safe_abs(complex(b_re, b_im) - ref)))
    return pairs


def _vsh_pairs(fname, ref_fn):
    """(L, err) pairs for a batch VSH_*_ALL family (max over 3 vector components)."""
    pairs = []
    for row in _load_rows(fname):
        l, m, theta, phi = int(row[0]), int(row[1]), float(row[2]), float(row[3])
        if l > SAFE_LMAX:
            continue
        vals = [float(v) for v in row[4:16]]
        b = [complex(vals[6], vals[7]), complex(vals[8], vals[9]),
             complex(vals[10], vals[11])]
        ref = ref_fn(l, m, theta, phi)
        pairs.append((l, max(_safe_abs(b[i] - ref[i]) for i in range(3))))
    return pairs


# Same five routine families as py/plot_benchmark.py.
FAMILIES = [
    ('Legendre', _alm_pairs),
    ('SSH', _ssh_pairs),
    ('VSH toroidal',
     lambda: _vsh_pairs('reference_vsh_tor.dat', mpref.vsh_tor_ref)),
    ('VSH poloidal (up)',
     lambda: _vsh_pairs('reference_vsh_pol_up.dat', mpref.vsh_pol_up_ref)),
    ('VSH poloidal (dn)',
     lambda: _vsh_pairs('reference_vsh_pol_dn.dat', mpref.vsh_pol_dn_ref)),
]


def plot_convergence():
    """
    Max abs error vs. the independent mpmath reference, per family, for
    the normalized batch routine only (one line per panel). max_by_degree
    (compute_validation.py) reduces the sparse M/theta samples at each
    Lmax to a single worst-case point per degree.
    """
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.4), sharey=True,
                              facecolor=SURFACE)
    for ax, (title, pairs_fn) in zip(axes, FAMILIES):
        pairs = pairs_fn()
        l_vals, err_vals = max_by_degree(pairs)
        style_axes(ax)
        ax.plot(l_vals, np.clip(err_vals, 1e-18, None), color=BLUE,
                marker='o', markersize=4, linewidth=1.5,
                solid_capstyle='round', zorder=3)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.9, SAFE_LMAX * 1.15)
        ax.set_xlabel(r'$\ell$', fontsize=9, color=INK_SECONDARY)
        ax.set_title(title, color=INK_PRIMARY, fontsize=10)
    axes[0].set_ylabel('max abs error vs.\nindependent (mpmath) reference',
                        fontsize=9, color=INK_SECONDARY)
    fig.suptitle('Normalized batch routines (ASSOC_LEGENDRE_NORM_ALL, '
                 'SSH_ALL, VSH_*_ALL), verified-safe range',
                 color=INK_PRIMARY, fontsize=11, y=1.02)
    fig.tight_layout()
    savefig_pair(fig, 'convergence_accuracy')


if __name__ == '__main__':
    plot_convergence()
