#!/usr/bin/env python3
"""
Heatmap of ASSOC_LEGENDRE_NORM_ALL's numerical stability across (l, theta):
peak |PNORM| observed at each (theta, l) checkpoint, log-scaled and clipped.
See STABILITY_FINDINGS.md at the repo root for the full writeup -- this
figure is the visual companion to that document's onset-l table.

Reads stability/stability_map.dat (produced by the opt-in vsh_stability
executable).

Run from the project root, after generating the scan data:
    cmake -B build -DVSH_BUILD_STABILITY=ON && cmake --build build
    ./build/vsh_stability
    python3 py/plot_stability.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import INK_PRIMARY, INK_SECONDARY, INK_MUTED, SURFACE, \
    savefig_pair, ROOT_DIR

STABILITY_DIR = os.path.join(ROOT_DIR, 'stability')

# Values above this are indistinguishable "broken" for display purposes
# (some are literal double-precision Infinity) -- clipping keeps the
# colormap's dynamic range meaningful instead of one pixel dominating it.
LOG_PEAK_CLIP = 20.0


def plot_stability():
    data = np.loadtxt(os.path.join(STABILITY_DIR, 'stability_map.dat'),
                       comments='#')
    theta_deg, l, broken, mmin, mmax, mmin_over_l, mmax_over_l, peak = data.T

    thetas = np.array(sorted(set(theta_deg)))
    ls = np.array(sorted(set(l)))
    grid = np.full((len(thetas), len(ls)), np.nan)
    theta_idx = {t: i for i, t in enumerate(thetas)}
    l_idx = {v: i for i, v in enumerate(ls)}
    for th, lv, pk in zip(theta_deg, l, peak):
        grid[theta_idx[th], l_idx[lv]] = np.log10(pk) if pk > 0 else -18

    grid_clipped = np.clip(grid, None, LOG_PEAK_CLIP)

    fig, ax = plt.subplots(figsize=(8, 5.5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    mesh = ax.pcolormesh(ls, thetas, grid_clipped, shading='nearest',
                          cmap='magma', vmin=0, vmax=LOG_PEAK_CLIP)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r'$\log_{10}|N_\ell^m P_\ell^m(\cos\theta)|$ (clipped at '
                   f'{LOG_PEAK_CLIP:.0f})', fontsize=9, color=INK_SECONDARY)
    cbar.ax.tick_params(colors=INK_MUTED, labelsize=8)

    ax.set_xlabel(r'$\ell$', fontsize=10, color=INK_SECONDARY)
    ax.set_ylabel(r'$\theta$ (degrees)', fontsize=10, color=INK_SECONDARY)
    ax.set_title('ASSOC_LEGENDRE_NORM_ALL stability: peak normalized value '
                 'per (l, theta)', color=INK_PRIMARY, fontsize=11)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8)
    ax.axhline(90, color='white', linewidth=0.5, alpha=0.3, zorder=2)

    fig.tight_layout()
    savefig_pair(fig, 'stability_map')


if __name__ == '__main__':
    plot_stability()
