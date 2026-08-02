#!/usr/bin/env python3
"""
Plot the convergence_accuracy_unnormalized figure: max absolute error
against an independent arbitrary-precision reference (py/mpmath_reference.py)
for the unnormalized Legendre recurrence (ASSOC_LEGENDRE, ASSOC_LEGENDRE_ALL),
restricted to its verified safe range l<=150 (see STABILITY_FINDINGS.md,
Part 2).

This is a single-panel companion to convergence_accuracy_normalized.py's
5-panel figure. Only Legendre is shown here (not SSH/VSH) because every
single-mode SSH/VSH function's stability is inherited entirely from this
same unnormalized recurrence (see the @warning notes added to their
docstrings in src/vsh.f90) -- Legendre is the root cause, not one case
among five independent ones.

Deliberately scoped to l<=150 rather than sweeping until it breaks: that
breakdown is already thoroughly documented in STABILITY_FINDINGS.md, and a
manuscript accuracy figure showing the *usable* range is more informative
than one dominated by a cliff. Both single-mode (ASSOC_LEGENDRE) and batch
(ASSOC_LEGENDRE_ALL) are plotted -- they overlap almost exactly, since it's
the same recurrence -- to keep the color convention consistent with the
normalized figure and the benchmark figures.

Reads stability/unnorm_diagonal.dat (produced by the opt-in
vsh_unnorm_stability executable), diagonal (m=l, the worst case -- see
STABILITY_FINDINGS.md Part 2) samples only.

Run from the project root, after generating the scan data:
    cmake -B build -DVSH_BUILD_STABILITY=ON && cmake --build build
    ./build/vsh_unnorm_stability
    python3 py/plot_convergence_unnormalized.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import mpmath_reference as mpref
from compute_validation import max_by_degree
from plotstyle import BLUE, ORANGE, INK_PRIMARY, INK_SECONDARY, SURFACE, \
    style_axes, savefig_pair, ROOT_DIR

STABILITY_DIR = os.path.join(ROOT_DIR, 'stability')
SAFE_LMAX = 150


def _load_rows(fname):
    rows = []
    with open(os.path.join(STABILITY_DIR, fname)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append(line.split())
    return rows


def _diagonal_pairs():
    """
    (L, relative_err) pairs for single-mode ASSOC_LEGENDRE and batch
    ASSOC_LEGENDRE_ALL, diagonal (m=l) samples only, restricted to the
    verified-safe l<=150. Relative error, not absolute -- the unnormalized
    P_l^l(x) grows combinatorially with l even within the safe range (that's
    the whole reason ASSOC_LEGENDRE_NORM_ALL exists), so absolute error
    would just track that growth and say nothing about accuracy, exactly
    the pitfall STABILITY_FINDINGS.md Part 2 worked around when actually
    characterizing this range.
    """
    single, batch = [], []
    for row in _load_rows('unnorm_diagonal.dat'):
        theta_deg, l, m, mfrac = float(row[0]), int(row[1]), int(row[2]), float(row[3])
        if mfrac != 1.0 or l > SAFE_LMAX:
            continue
        p_single, p_batch = float(row[4]), float(row[5])
        x = np.cos(np.radians(theta_deg))
        ref = mpref.legendre_ref(l, m, x).real
        scale = max(abs(ref), 1e-300)
        single.append((l, abs(p_single - ref) / scale))
        batch.append((l, abs(p_batch - ref) / scale))
    return single, batch


def plot_convergence_unnormalized():
    single_pairs, batch_pairs = _diagonal_pairs()
    l_single, err_single = max_by_degree(single_pairs)
    l_batch, err_batch = max_by_degree(batch_pairs)

    fig, ax = plt.subplots(figsize=(5.5, 4.2), facecolor=SURFACE)
    style_axes(ax)
    ax.plot(l_single, np.clip(err_single, 1e-18, None), color=ORANGE,
            marker='o', markersize=4, linewidth=1.5, solid_capstyle='round',
            label='Single-mode', zorder=3)
    ax.plot(l_batch, np.clip(err_batch, 1e-18, None), color=BLUE,
            marker='o', markersize=4, linewidth=1.5, solid_capstyle='round',
            label='Batch (_ALL)', zorder=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.9, SAFE_LMAX * 1.1)
    ax.set_xlabel(r'$\ell$', fontsize=10, color=INK_SECONDARY)
    ax.set_ylabel('max relative error vs.\nindependent (mpmath) reference',
                  fontsize=9, color=INK_SECONDARY)
    ax.set_title(f'Legendre (unnormalized), worst-case $m=\\ell$, '
                 f'verified-safe range', color=INK_PRIMARY, fontsize=10)
    ax.legend(loc='upper left', frameon=False, fontsize=9,
              labelcolor=INK_SECONDARY)
    fig.tight_layout()
    savefig_pair(fig, 'convergence_accuracy_unnormalized')


if __name__ == '__main__':
    plot_convergence_unnormalized()
