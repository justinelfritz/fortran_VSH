#!/usr/bin/env python3
"""
Plot FORTVSH batch-vs-naive-loop timing benchmarks.

Reads benchmark/bench_*.dat (produced by vsh_benchmark; see the
"Benchmarking" section of README.md) and writes two figures to
figures/:

  benchmark_scaling.pdf  -- wall time per evaluation vs Lmax, batch vs
                            naive per-mode loop, one panel per routine
  benchmark_speedup.pdf  -- batch/loop speedup vs Lmax, overlaid across
                            routines

bench_main.f90's LMAX_LIST extends to 4000 (from the investigation that
produced STABILITY_FINDINGS.md), but each series here is truncated to its
OWN verified-safe range rather than plotted all the way out: "Naive loop"
(single-mode, built on the unnormalized recurrence) stops at
NAIVE_SAFE_LMAX=150, "Batch (_ALL)" (normalized recurrence) stops at
BATCH_SAFE_LMAX=2000 (STABILITY_FINDINGS.md Parts 2 and 1 respectively).
Wall-clock timing itself doesn't care whether the computed value is
correct -- a NaN takes just as long to produce as a valid float -- so nothing
would visibly break if these were left unfiltered, but showing timing data
for a regime where the manuscript elsewhere states the routine gives wrong
answers would be a strange thing to leave in. Truncating each line at its
own limit rather than capping the whole figure at 150 also makes a real
point: batch isn't just faster, it's the only one still valid past l=150.

Run from the project root, after generating benchmark data:
    cmake -B build -DVSH_BUILD_BENCHMARK=ON && cmake --build build
    ./build/vsh_benchmark
    python3 py/plot_benchmark.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plotstyle import BLUE, ORANGE, AQUA, YELLOW, INK_PRIMARY, \
    INK_SECONDARY, SURFACE, style_axes, savefig_pair

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)
BENCH_DIR  = os.path.join(ROOT_DIR, 'benchmark')

# Verified-safe ranges from STABILITY_FINDINGS.md (Parts 2 and 1).
NAIVE_SAFE_LMAX = 150
BATCH_SAFE_LMAX = 2000

# All five benchmarked routines, for the per-routine scaling small multiples.
ROUTINES = [
    ('bench_legendre.dat',   'Legendre'),
    ('bench_ssh.dat',        'SSH'),
    ('bench_vsh_tor.dat',    'VSH toroidal'),
    ('bench_vsh_pol_up.dat', 'VSH poloidal (up)'),
    ('bench_vsh_pol_dn.dat', 'VSH poloidal (dn)'),
]

# VSH_POL_UP and VSH_POL_DN are mirror-image computations with near-identical
# cost, so the combined speedup overlay drops POL_DN to avoid two
# indistinguishable overlapping lines; its data is still in the scaling
# small multiples above.
SPEEDUP_ROUTINES = [
    ('bench_legendre.dat',   'Legendre',      BLUE),
    ('bench_ssh.dat',        'SSH',           ORANGE),
    ('bench_vsh_tor.dat',    'VSH toroidal',  AQUA),
    ('bench_vsh_pol_up.dat', 'VSH poloidal',  YELLOW),
]


def load(fname):
    lmax, n_modes, t_batch, t_loop, speedup = np.loadtxt(
        os.path.join(BENCH_DIR, fname), comments='#', unpack=True)
    return lmax, n_modes, t_batch, t_loop, speedup


def plot_scaling():
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.4), sharey=True,
                              facecolor=SURFACE)
    for ax, (fname, label) in zip(axes, ROUTINES):
        lmax, _, t_batch, t_loop, _ = load(fname)
        naive_mask = lmax <= NAIVE_SAFE_LMAX
        batch_mask = lmax <= BATCH_SAFE_LMAX
        style_axes(ax)
        ax.plot(lmax[naive_mask], t_loop[naive_mask], color=ORANGE,
                marker='o', markersize=5, linewidth=2,
                solid_capstyle='round', label='Naive loop', zorder=3)
        ax.plot(lmax[batch_mask], t_batch[batch_mask], color=BLUE,
                marker='o', markersize=5, linewidth=2,
                solid_capstyle='round', label='Batch (_ALL)', zorder=3)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(label, color=INK_PRIMARY, fontsize=12)
        ax.set_xlabel(r'$L_{max}$', fontsize=12, color=INK_SECONDARY)
        ax.set_xticks([1.0, 10.0, 100.0, 1000.0])
        ax.set_xticklabels([r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$'],fontsize=11,color=INK_SECONDARY)
        ax.set_yticks([1.e-8, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1])
        ax.set_yticklabels([r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$',r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$',r'$10^{-2}$', r'$10^{-1}$'],fontsize=11,color=INK_SECONDARY)
    axes[0].set_ylabel('wall time per evaluation (s)', fontsize=12,
                        color=INK_SECONDARY)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.06), fontsize=12,
               labelcolor=INK_SECONDARY)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def plot_speedup():
    """
    Speedup = T_loop/T_batch requires both series to mean something at the
    same Lmax, so this is truncated to NAIVE_SAFE_LMAX (the tighter of the
    two ranges) rather than each line independently -- unlike plot_scaling,
    there's no valid "naive" baseline left to speed up from past l=150.
    """
    fig, ax = plt.subplots(figsize=(6.2, 4.2), facecolor=SURFACE)
    style_axes(ax)
    last_lmax = None
    for fname, label, color in SPEEDUP_ROUTINES:
        lmax, _, _, _, speedup = load(fname)
        mask = lmax <= NAIVE_SAFE_LMAX
        lmax, speedup = lmax[mask], speedup[mask]
        last_lmax = lmax
        ax.plot(lmax, speedup, color=color, marker='o', markersize=6,
                linewidth=2, solid_capstyle='round', zorder=3)
        ax.annotate(label, xy=(lmax[-1], speedup[-1]),
                    xytext=(8, 0), textcoords='offset points',
                    va='center', fontsize=9, color=color, clip_on=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(right=last_lmax[-1] * 3.2)
    ax.set_xlabel(r'$L_{max}$', fontsize=10, color=INK_SECONDARY)
    ax.set_ylabel('speedup (naive loop / batch)', fontsize=10,
                  color=INK_SECONDARY)
    ax.set_title('Batch vs. naive per-mode loop speedup '
                 '(naive-loop verified-safe range)',
                 color=INK_PRIMARY, fontsize=11)
    fig.tight_layout()
    return fig


if __name__ == '__main__':
    savefig_pair(plot_scaling(), 'benchmark_scaling')
    savefig_pair(plot_speedup(), 'benchmark_speedup')
