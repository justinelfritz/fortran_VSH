#!/usr/bin/env python3
"""
Plot cross-validation accuracy: percent difference between the three
independent derivations behind Tests 1 and 2 (closed-form analytic,
VSH-inner-product, and Geppert & Wiebicke forms) across the full
colatitude domain. This is the visual companion to the scalar max-error
values already reported via validation_values.tex -- it shows the three
forms agree to numerical noise everywhere, not just at that one max.

Run from the project root, after the test suite has produced
validation/TEST1.dat and validation/TEST2.dat:
    cmake --build build && ctest --test-dir build
    python3 py/plotError.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from compute_validation import parse_test_file
from plotstyle import BLUE, ORANGE, INK_PRIMARY, INK_SECONDARY, \
    style_axes, savefig_pair

EPS = 1.0e-15   # reference values below this are treated as exactly zero


def pct_diff(reference, comparison):
    """
    |comparison - reference| / |reference| * 100, elementwise.
    Points where |reference| < EPS are set to 0 (avoids blow-up at the
    isolated zeros of the analytic form); the result is floored at a
    small positive value so it still plots on a log axis.
    """
    ref = np.asarray(reference)
    cmp = np.asarray(comparison)
    denom = np.abs(ref)
    diff = np.abs(cmp - ref)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = np.where(denom < EPS, 0.0, 100.0 * diff / denom)
    out = np.nan_to_num(out, nan=0.0)
    return np.clip(out, 1e-18, None)


def plot_test(ax, fname, title):
    theta, analytic, vsh_form, gw_form = parse_test_file(fname)
    vsh_err = pct_diff(analytic, vsh_form)
    gw_err  = pct_diff(analytic, gw_form)
    style_axes(ax)
    ax.plot(theta, vsh_err, color=BLUE, linewidth=2, solid_capstyle='round',
            label='VSH inner-product form', zorder=3)
    ax.plot(theta, gw_err, color=ORANGE, linewidth=2, solid_capstyle='round',
            label='Geppert-Wiebicke form', zorder=3)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\theta$ (rad)', fontsize=9, color=INK_SECONDARY)
    ax.set_title(title, color=INK_PRIMARY, fontsize=10)


if __name__ == '__main__':
    from plotstyle import SURFACE

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), facecolor=SURFACE)
    plot_test(axes[0], 'TEST1.dat', 'Test 1')
    plot_test(axes[1], 'TEST2.dat', 'Test 2')
    axes[0].set_ylabel('difference from analytic form (%)', fontsize=9,
                        color=INK_SECONDARY)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.06), fontsize=9,
               labelcolor=INK_SECONDARY)
    fig.tight_layout(rect=(0, 0, 1, 0.90))

    savefig_pair(fig, 'cross_validation_accuracy')
