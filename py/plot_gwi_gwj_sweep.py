#!/usr/bin/env python3
"""
Plot the Geppert-Wiebicke coefficient sweep worked example (2e): GWI
(vs. the brute-force Gaunt-coefficient integral) and GWJ (vs. the
DOT(PVSH_POL,PVSH_TOR) pointwise identity), each swept across mode
combinations and aggregated by J1+J2.

Reads examples/gwi_gwj_sweep/*.dat (produced by the gwi_gwj_sweep
executable) and writes gwi_gwj_sweep_accuracy.{pdf,png} to
tex/Copernicus-EGU/figures/.

Run from the project root, after generating the example data:
    cmake -B build -DVSH_BUILD_EXAMPLES=ON && cmake --build build
    ./build/gwi_gwj_sweep
    python3 py/plot_gwi_gwj_sweep.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import BLUE, ORANGE, INK_PRIMARY, INK_SECONDARY, SURFACE, \
    AXIS_LINE, style_axes, savefig_pair, ROOT_DIR

EX_DIR = os.path.join(ROOT_DIR, 'examples', 'gwi_gwj_sweep')
MACHINE_EPS = np.finfo(np.float64).eps


def max_by_sum(j1, j2, absdiff):
    """Reduce per-row (J1,J2,absdiff) to (sorted J1+J2, max absdiff per sum)."""
    total = (j1 + j2).astype(int)
    sums = sorted(set(total))
    return sums, [absdiff[total == s].max() for s in sums]


if __name__ == '__main__':
    gwi_j1, gwi_m1, gwi_j2, gwi_m2, gwi_l, gwi_m, _, _, _, _, gwi_err = \
        np.loadtxt(os.path.join(EX_DIR, 'gwi_check.dat'), comments='#',
                   unpack=True)
    gwj_j1, gwj_m1, gwj_j2, gwj_m2, _, _, _, _, _, _, gwj_err = \
        np.loadtxt(os.path.join(EX_DIR, 'gwj_check.dat'), comments='#',
                   unpack=True)

    gwi_sums, gwi_max = max_by_sum(gwi_j1, gwi_j2, gwi_err)
    gwj_sums, gwj_max = max_by_sum(gwj_j1, gwj_j2, gwj_err)

    fig, ax = plt.subplots(figsize=(6.5, 4.5), facecolor=SURFACE)
    style_axes(ax)
    ax.axhline(MACHINE_EPS, color=AXIS_LINE, linewidth=1.2, linestyle=':',
               zorder=2, label='double-precision epsilon')
    ax.plot(gwi_sums, np.clip(gwi_max, 1e-18, None), 'o-', color=BLUE,
            markersize=6, linewidth=2, solid_capstyle='round',
            label='GWI vs. Gaunt-coefficient integral', zorder=3)
    ax.plot(gwj_sums, np.clip(gwj_max, 1e-18, None), 'o-', color=ORANGE,
            markersize=6, linewidth=2, solid_capstyle='round',
            label='GWJ vs. DOT(PVSH_POL,PVSH_TOR) identity', zorder=3)
    ax.set_yscale('log')
    ax.set_xlabel(r'$J_1+J_2$', fontsize=12, color=INK_SECONDARY)
    ax.set_ylabel('max |closed form - independent check|', fontsize=12,
                  color=INK_SECONDARY)
    ax.set_xticks([0,2,4,6,8,10,12])
    ax.set_xticklabels([0,2,4,6,8,10,12],fontsize=11,color=INK_SECONDARY)
    ax.set_yticks([1.e-17, 1.e-16, 1.e-15, 1.e-14])
    ax.set_yticklabels([r'$10^{-17}$', r'$10^{-16}$',r'$10^{-15}$', r'$10^{-14}$'],fontsize=11,color=INK_SECONDARY)
    ax.set_title('Geppert-Wiebicke coefficients across parameter ranges',
                 color=INK_PRIMARY, fontsize=12)
    ax.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    savefig_pair(fig, 'gwi_gwj_sweep_accuracy')
