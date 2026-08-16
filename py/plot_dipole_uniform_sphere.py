#!/usr/bin/env python3
"""
Plot the uniformly-magnetized-sphere worked example (2b): a uniform
interior field matched to an exterior dipole at r=R, showing Br continuous
and Btheta discontinuous (bound surface current) across the boundary.

Reads examples/dipole_uniform_sphere/field_profile.dat (produced by the
dipole_uniform_sphere executable) and writes
dipole_uniform_sphere_match.{pdf,png} to figures/.

Run from the project root, after generating the example data:
    cmake -B build -DVSH_BUILD_EXAMPLES=ON && cmake --build build
    ./build/dipole_uniform_sphere
    python3 py/plot_dipole_uniform_sphere.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import BLUE, ORANGE, AQUA, INK_PRIMARY, INK_SECONDARY, \
    SURFACE, AXIS_LINE, style_axes, savefig_pair, ROOT_DIR

EX_DIR = os.path.join(ROOT_DIR, 'examples', 'dipole_uniform_sphere')
R_SPHERE = 1.0
THETA_LABELS = [r'$\theta=0$', r'$\theta=\pi/4$', r'$\theta=\pi/2$']
COLORS = [BLUE, ORANGE, AQUA]
MACHINE_EPS = np.finfo(np.float64).eps


if __name__ == '__main__':
    r, theta, br_lib, bth_lib, br_ana, bth_ana = np.loadtxt(
        os.path.join(EX_DIR, 'field_profile.dat'), comments='#', unpack=True)

    n_theta = len(THETA_LABELS)
    n_r = len(r) // n_theta
    r_grid = r[:n_r]
    br_blocks = br_lib.reshape(n_theta, n_r)
    bth_blocks = bth_lib.reshape(n_theta, n_r)
    err_br_blocks = np.clip(np.abs(br_lib - br_ana), 1e-18, None).reshape(n_theta, n_r)
    err_bth_blocks = np.clip(np.abs(bth_lib - bth_ana), 1e-18, None).reshape(n_theta, n_r)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, figsize=(10, 7.5), facecolor=SURFACE)

    for ax, blocks, ylabel, title in (
        (ax1, br_blocks, r'$B_r$', r'$B_r$: continuous at $r=R$'),
        (ax2, bth_blocks, r'$B_\theta$', r'$B_\theta$: jump at $r=R$ (surface current)'),
    ):
        style_axes(ax)
        for block, color, label in zip(blocks, COLORS, THETA_LABELS):
            ax.plot(r_grid, block, color=color, linewidth=2,
                    solid_capstyle='round', label=label, zorder=3)
        ax.axvline(R_SPHERE, color=AXIS_LINE, linewidth=1.2, linestyle='--',
                   zorder=2)
        ax.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        ax.set_xticklabels([0.5, 1.0, 1.5, 2.0, 2.5, 3.0],fontsize=11,color=INK_SECONDARY)
        if ax == ax1:
            ax.set_yticks([0.5*i for i in range(0,5)])
            ax.set_yticklabels([0.5*i for i in range(0,5)],fontsize=11,color=INK_SECONDARY)
        if ax == ax2:
            ax.set_yticks([-2.0+0.5*i for i in range(0,7)])
            ax.set_yticklabels([-2.0+0.5*i for i in range(0,7)],fontsize=11,color=INK_SECONDARY)
        ax.set_xlabel(r'$r/R$', fontsize=11, color=INK_SECONDARY)
        ax.set_ylabel(ylabel, fontsize=12, color=INK_SECONDARY)
        ax.set_title(title, color=INK_PRIMARY, fontsize=11)
        ax.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY)

    for ax, blocks, ylabel, title in (
        (ax3, err_br_blocks, r'$|B_r^{lib}-B_r^{analytic}|$',
         'Library vs. closed-form error'),
        (ax4, err_bth_blocks, r'$|B_\theta^{lib}-B_\theta^{analytic}|$',
         'Library vs. closed-form error'),
    ):
        style_axes(ax)
        ax.axhline(MACHINE_EPS, color=AXIS_LINE, linewidth=1.2, linestyle=':',
                   zorder=2, label='double-precision epsilon')
        for block, color, label in zip(blocks, COLORS, THETA_LABELS):
            ax.plot(r_grid, block, 'o', color=color, markersize=3.5,
                    alpha=0.75, label=label, zorder=3)
        ax.axvline(R_SPHERE, color=AXIS_LINE, linewidth=1.2, linestyle='--',
                   zorder=2)
        ax.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        ax.set_xticklabels([0.5, 1.0, 1.5, 2.0, 2.5, 3.0],fontsize=11,color=INK_SECONDARY)
        ax.set_yticks([1.e-18, 1.e-17, 1.e-16, 1.e-15, 1.e-14])
        ax.set_yticklabels([r'$10^{-18}$', r'$10^{-17}$', r'$10^{-16}$',r'$10^{-15}$', r'$10^{-14}$'],fontsize=11,color=INK_SECONDARY)
        ax.set_yscale('log')
        ax.set_ylim(top=1e-14)
        ax.set_xlabel(r'$r/R$', fontsize=12, color=INK_SECONDARY)
        ax.set_ylabel(ylabel, fontsize=12, color=INK_SECONDARY)
        ax.set_title(title, color=INK_PRIMARY, fontsize=11)
        ax.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY)

    fig.tight_layout()
    savefig_pair(fig, 'dipole_uniform_sphere_match')
