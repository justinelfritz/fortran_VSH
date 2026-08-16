#!/usr/bin/env python3
"""
Plot the pure-toroidal VSH spectral decomposition worked example (2f):
known-mode recovery (Part A) and spectral truncation convergence (Part B).
Companion to py/plot_vsh_decomposition.py (example 2d, radial family) --
same two-panel layout and styling, applied to the toroidal-only field
built in src/examples/vsh_decomposition_tor.f90 (eight known toroidal
modes for Part A; a quadrupole (l=2,m=0) background plus a Gaussian bump
shaped in both theta and phi for Part B).

Reads examples/vsh_decomposition_tor/*.dat (produced by the
vsh_decomposition_tor executable) and writes
vsh_decomposition_tor_convergence.{pdf,png} to tex/Copernicus-EGU/figures/.

Run from the project root, after generating the example data:
    cmake -B build -DVSH_BUILD_EXAMPLES=ON && cmake --build build
    ./build/vsh_decomposition_tor
    python3 py/plot_vsh_decomposition_tor.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import BLUE, ORANGE, AQUA, INK_PRIMARY, INK_SECONDARY, \
    SURFACE, AXIS_LINE, style_axes, savefig_pair, ROOT_DIR

EX_DIR = os.path.join(ROOT_DIR, 'examples', 'vsh_decomposition_tor')


if __name__ == '__main__':
    L, M, ar_re, ar_im, ap_re, ap_im, at_re, at_im, expected = np.loadtxt(
        os.path.join(EX_DIR, 'known_modes.dat'), comments='#', unpack=True)
    lmax_list, max_err = np.loadtxt(
        os.path.join(EX_DIR, 'convergence.dat'), comments='#', unpack=True)

    mag_rad = np.abs(ar_re + 1j*ar_im)
    mag_pol = np.abs(ap_re + 1j*ap_im)
    mag_tor = np.abs(at_re + 1j*at_im)
    mode_idx = np.arange(1, len(L)+1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)

    style_axes(ax1)
    ax1.axhline(1.0, color=AXIS_LINE, linewidth=1.2, linestyle=':', zorder=2,
                label='expected a_tor (1.0 or 0)')
    # Small horizontal jitter so the three families (which coincide
    # exactly at 0 -- or, for toroidal, at 1.0 -- by construction) don't
    # fully hide each other. Unlike example 2d, only a_tor is ever
    # expected nonzero here: a_rad/a_pol pinned at the floor throughout is
    # the direct visual confirmation that this field really is purely
    # toroidal, not just that the known coefficients come back right.
    jitter = 0.22
    ax1.plot(mode_idx-jitter, np.clip(mag_rad, 1e-18, None), 'o', color=BLUE,
             markersize=4, alpha=0.8, linestyle='-', label='radial', zorder=3)
    ax1.plot(mode_idx, np.clip(mag_pol, 1e-18, None), 'o', color=ORANGE,
             markersize=4, alpha=0.8, linestyle='-',label='poloidal', zorder=3)
    ax1.plot(mode_idx+jitter, np.clip(mag_tor, 1e-18, None), 'o', color=AQUA,
             markersize=4, alpha=0.8, linestyle='--',label='toroidal', zorder=3)
    ax1.set_yscale('log')
    ax1.set_ylim(top=6.0)
    ax1.set_xticks([5*i for i in range(0,6)])
    ax1.set_xticklabels([5*i for i in range(0,6)],fontsize=11,color=INK_SECONDARY)
    ax1.set_yticks([1.e-16, 1.e-12, 1.e-8, 1.e-4, 1.e0])
    ax1.set_yticklabels([r'$10^{-16}$', r'$10^{-12}$', r'$10^{-8}$',r'$10^{-4}$', r'$10^{0}$'],fontsize=11,color=INK_SECONDARY)
    
    ax1.set_xlabel('mode index (ordered by l, m)', fontsize=11,
                    color=INK_SECONDARY)
    ax1.set_ylabel('|recovered coefficient|', fontsize=11,
                    color=INK_SECONDARY)
    ax1.set_title('Part A: known-mode recovery (8 toroidal modes, l<=2)',
                  color=INK_PRIMARY, fontsize=11)
    ax1.legend(frameon=True, fontsize=10, labelcolor=INK_SECONDARY, ncol=2,
               loc='center right')

    style_axes(ax2)
    ax2.plot(lmax_list, np.clip(max_err, 1e-18, None), color=BLUE,
              marker='o', markersize=6, linewidth=2, solid_capstyle='round',
              zorder=3)
    ax2.set_yscale('log')
    ax2.set_xticks([5*i for i in range(1,7)])
    ax2.set_xticklabels([5*i for i in range(1,7)],fontsize=11,color=INK_SECONDARY)
    ax2.set_yticks([1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e-0])
    ax2.set_yticklabels([r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$',r'$10^{-2}$', r'$10^{-1}$', r'$10^{-0}$'],fontsize=11,color=INK_SECONDARY)
    
    ax2.set_xlabel(r'$L_{max}$ retained', fontsize=11, color=INK_SECONDARY)
    ax2.set_ylabel('max reconstruction error', fontsize=11,
                    color=INK_SECONDARY)
    ax2.set_title('Part B: spectral truncation convergence\n'
                  '(quadrupole + toroidal Gaussian bump)', color=INK_PRIMARY,
                  fontsize=11)

    fig.tight_layout()
    savefig_pair(fig, 'vsh_decomposition_tor_convergence')
