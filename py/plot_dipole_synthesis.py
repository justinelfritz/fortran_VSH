#!/usr/bin/env python3
"""
Plot the dipole-synthesis worked example (2a): the exterior vacuum dipole
field built from FORTVSH's PVSH_RAD/PVSH_POL at l=1,m=0, and its overlap
with l=0,2,3 (purity check).

Reads examples/dipole_synthesis/*.dat (produced by the dipole_synthesis
executable; see the plan's "Real-world usage demonstration" section) and
writes dipole_synthesis_accuracy.{pdf,png} to tex/Copernicus-EGU/figures/.

Run from the project root, after generating the example data:
    cmake -B build -DVSH_BUILD_EXAMPLES=ON && cmake --build build
    ./build/dipole_synthesis
    python3 py/plot_dipole_synthesis.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import BLUE, ORANGE, INK_PRIMARY, INK_SECONDARY, \
    SURFACE, AXIS_LINE, style_axes, savefig_pair, ROOT_DIR

EX_DIR = os.path.join(ROOT_DIR, 'examples', 'dipole_synthesis')
MACHINE_EPS = np.finfo(np.float64).eps


def load(fname):
    return np.loadtxt(os.path.join(EX_DIR, fname), comments='#', unpack=True)


if __name__ == '__main__':
    theta, br_lib, bth_lib, br_ana, bth_ana = load('field_profile.dat')
    l_vals, overlap_rad, overlap_pol = load('mode_overlap.dat')

    err_br  = np.clip(np.abs(br_lib - br_ana), 1e-18, None)
    err_bth = np.clip(np.abs(bth_lib - bth_ana), 1e-18, None)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4), facecolor=SURFACE)

    style_axes(ax1)
    ax1.plot(theta, br_lib, color=BLUE, linewidth=2, solid_capstyle='round',
              label=r'$B_r$', zorder=3)
    ax1.plot(theta, bth_lib, color=ORANGE, linewidth=2, solid_capstyle='round',
              label=r'$B_\theta$', zorder=3)
    ax1.set_xlabel(r'$\theta$ (rad)', fontsize=11, color=INK_SECONDARY)
    ax1.set_ylabel(r'field (units of $\mu/R_0^3$)', fontsize=11,
                    color=INK_SECONDARY)
    ax1.set_xticks([0.0, np.pi/4.0, np.pi/2.0, 3.0*np.pi/4.0, np.pi])
    ax1.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'],fontsize=11,color=INK_SECONDARY)
    ax1.set_yticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    ax1.set_yticklabels([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0],fontsize=11,color=INK_SECONDARY)
    ax1.set_title('Synthesized exterior dipole', color=INK_PRIMARY, fontsize=11)
    ax1.legend(frameon=True, fontsize=12, labelcolor=INK_SECONDARY)

    style_axes(ax2)
    ax2.axhline(MACHINE_EPS, color=AXIS_LINE, linewidth=1.2, linestyle=':',
                zorder=2, label='double-precision epsilon')
    ax2.plot(theta, err_br, 'o', color=BLUE, markersize=3.5, alpha=0.75,
              label=r'$|B_r^{lib}-B_r^{analytic}|$', zorder=3)
    ax2.plot(theta, err_bth, 'o', color=ORANGE, markersize=3.5, alpha=0.75,
              label=r'$|B_\theta^{lib}-B_\theta^{analytic}|$', zorder=3)
    ax2.set_yscale('log')
    ax2.set_ylim(top=1e-14)
    ax2.set_xlabel(r'$\theta$ (rad)', fontsize=11, color=INK_SECONDARY)
    ax2.set_ylabel('absolute error', fontsize=11, color=INK_SECONDARY)
    ax2.set_xticks([0.0, np.pi/4.0, np.pi/2.0, 3.0*np.pi/4.0, np.pi])
    ax2.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'],fontsize=11,color=INK_SECONDARY)
    ax2.set_yticks([1.e-18, 1.e-17, 1.e-16, 1.e-15, 1.e-14])
    ax2.set_yticklabels([r'$10^{-18}$', r'$10^{-17}$', r'$10^{-16}$',r'$10^{-15}$', r'$10^{-14}$'],fontsize=11,color=INK_SECONDARY)      
    ax2.set_title('Library vs. closed-form error', color=INK_PRIMARY, fontsize=11)
    ax2.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY)

    style_axes(ax3)
    width = 0.35
    x = np.arange(len(l_vals))
    floor = 1e-18
    ax3.bar(x - width/2, np.clip(overlap_rad, floor, None), width,
            color=BLUE, label='radial overlap', zorder=3)
    ax3.bar(x + width/2, np.clip(overlap_pol, floor, None), width,
            color=ORANGE, label='horizontal overlap', zorder=3)
    ax3.set_yscale('log')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f'{int(l)}' for l in l_vals],fontsize=11,color=INK_SECONDARY)
    ax3.set_yticks([1.e-16, 1.e-13, 1.e-10, 1.e-7, 1.e-4, 1.e-1])
    ax3.set_yticklabels([r'$10^{-16}$', r'$10^{-13}$', r'$10^{-10}$',r'$10^{-7}$', r'$10^{-4}$', r'$10^{-1}$'],fontsize=11,color=INK_SECONDARY)   
    ax3.set_xlabel(r'$\ell$', fontsize=12, color=INK_SECONDARY)
    ax3.set_ylabel('|overlap with target field|', fontsize=11,
                    color=INK_SECONDARY)
    ax3.set_title('Mode purity: only $\\ell=1$ survives', color=INK_PRIMARY,
                  fontsize=11)
    ax3.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    savefig_pair(fig, 'dipole_synthesis_accuracy')
