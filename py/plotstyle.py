"""
Shared matplotlib styling for FORTVSH validation/benchmark figures, so
the convergence, accuracy, and performance plots read as one consistent
set in the manuscript. Palette is the validated default from the
dataviz-skill reference palette (colorblind-safe fixed hue order).
"""

import os
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)
FIG_DIR    = os.path.join(ROOT_DIR, 'tex', 'Copernicus-EGU', 'figures')

# ── categorical palette, fixed order ─────────────────────────────────────────
BLUE    = '#2a78d6'
ORANGE  = '#eb6834'
AQUA    = '#1baf7a'
YELLOW  = '#eda100'
MAGENTA = '#e87ba4'

# ── chrome & ink ──────────────────────────────────────────────────────────────
INK_PRIMARY   = '#0b0b0b'
INK_SECONDARY = "#0b0b0b"
INK_MUTED     = '#898781'
GRIDLINE      = '#e1e0d9'
AXIS_LINE     = "#3a3a3a"
SURFACE       = '#fcfcfb'


def style_axes(ax):
    """Apply the shared chart chrome (surface, gridlines, spines, ticks)."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, which='major', color=GRIDLINE, linewidth=0.7, zorder=0)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_color(AXIS_LINE)
    ax.tick_params(colors=INK_MUTED, labelsize=8)


def savefig_pair(fig, name):
    """Save fig as both PDF (for LaTeX) and PNG (for quick viewing) into FIG_DIR."""
    os.makedirs(FIG_DIR, exist_ok=True)
    pdf_path = os.path.join(FIG_DIR, f'{name}.pdf')
    png_path = os.path.join(FIG_DIR, f'{name}.png')
    fig.savefig(pdf_path, dpi=300, facecolor=fig.get_facecolor(),
                bbox_inches='tight')
    fig.savefig(png_path, dpi=150, facecolor=fig.get_facecolor(),
                bbox_inches='tight')
    print(f'Wrote {pdf_path}')
    print(f'Wrote {png_path}')
