#!/usr/bin/env python3
"""
Line-plot companion to plot_vsh_decomposition_tor_field.py's 3D sphere:
Re(B_theta), Im(B_theta), Re(B_phi), Im(B_phi) vs. theta (at a fixed phi)
and vs. phi (at a fixed theta), for the same Part B target field
(src/examples/vsh_decomposition_tor.f90's FIELD_QUAD_BUMP).

The field is complex (see plot_vsh_decomposition_tor_field.py's module
docstring for why: the axisymmetric quadrupole term is purely imaginary,
the bump term purely real, verified against mpmath_reference.py). Splitting
into real/imaginary parts makes that phase split directly visible: Im(B_theta)
is exactly zero everywhere (the quadrupole's theta-component is identically
zero, and the bump is purely real), while Im(B_phi) traces the quadrupole's
sin(theta)*cos(theta) band on its own, separate from the bump's real-valued
contribution.

Left panel (theta-dependence): phi=pi/3. phi=pi/2 was deliberately avoided
here -- it's a special value in these equations (it maximizes sin(phi) and
zeroes cos(phi), so it isolates the bump's sin(phi) term in B_theta from
its cos(phi) term in B_phi rather than showing them mixed together as they
generically are). pi/3 gives both sin(phi) and cos(phi) nonzero, a more
representative slice. The theta axis runs the full [0,pi] range,
deliberately including the poles: the field was redesigned (see
FIELD_QUAD_BUMP's docstring) so that B_theta and B_phi vanish *exactly*
there, by construction -- every term in the closed form carries an
explicit sin(theta), sin(theta)**2, or sin(2*theta) factor. Unlike an
earlier, pole-centered version of this bump, there is no
coordinate-singularity ambiguity left to route around: since the field
genuinely goes to zero at theta=0,pi for *every* phi, "the value at the
pole" is unambiguous. The curves should visibly taper to zero at both ends.

Right panel (phi-dependence): theta=pi/3, off the equator (theta=pi/2) for
the same reason -- the equator is a node of the axisymmetric quadrupole's
Im(B_phi) (proportional to sin(theta)*cos(theta), zero at theta=pi/2) and
of the bump's B_phi generating-function derivative, so a theta=pi/2 slice
would flatten out structure that's actually there elsewhere. theta=pi/3
still sits well inside the bump's core (one bump-width, sigma=0.4, from
the equator) while showing the quadrupole's genuine phi-independent
contribution to B_phi.

Reads no data files -- evaluates the same closed-form field directly, like
plot_vsh_decomposition_tor_field.py.

Run from the project root:
    python3 py/plot_vsh_decomposition_tor_profiles.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plot_vsh_decomposition_tor_field import field_quad_bump
from plotstyle import BLUE, ORANGE, INK_PRIMARY, INK_SECONDARY, \
    SURFACE, style_axes, savefig_pair

PHI_FIXED = np.pi / 3
THETA_FIXED = np.pi / 3


def _plot_re_im(ax, x, bt, bp):
    ax.plot(x, bt.real, color=BLUE, linewidth=2, linestyle='-',
             solid_capstyle='round', label=r'$\mathrm{Re}(B_\theta)$', zorder=3)
    ax.plot(x, bt.imag, color=BLUE, linewidth=2, linestyle='--',
             dash_capstyle='round', label=r'$\mathrm{Im}(B_\theta)$', zorder=3)
    ax.plot(x, bp.real, color=ORANGE, linewidth=2, linestyle='-',
             solid_capstyle='round', label=r'$\mathrm{Re}(B_\phi)$', zorder=3)
    ax.plot(x, bp.imag, color=ORANGE, linewidth=2, linestyle='--',
             dash_capstyle='round', label=r'$\mathrm{Im}(B_\phi)$', zorder=3)


def plot_profiles():
    theta = np.linspace(0.0, np.pi, 400)
    phi = np.linspace(0.0, 2.0*np.pi, 400)

    bt_th, bp_th = field_quad_bump(theta, PHI_FIXED)
    bt_ph, bp_ph = field_quad_bump(THETA_FIXED, phi)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)

    style_axes(ax1)
    _plot_re_im(ax1, theta, bt_th, bp_th)
    ax1.set_xlabel(r'$\theta$ (rad)', fontsize=11, color=INK_SECONDARY)
    ax1.set_ylabel('field component', fontsize=11, color=INK_SECONDARY)
    ax1.set_xticks([0.0, np.pi/4.0, np.pi/2.0, 3.0*np.pi/4.0, np.pi])
    ax1.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'],fontsize=11,color=INK_SECONDARY)
    ax1.set_yticks([-2 + i for i in range(0,5)])
    ax1.set_yticklabels([-2 + i for i in range(0,5)],fontsize=11,color=INK_SECONDARY)
    ax1.set_title(r'$\theta$-dependence at $\phi=\pi/3$',
                  color=INK_PRIMARY, fontsize=12)
    ax1.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY, ncol=2)

    style_axes(ax2)
    _plot_re_im(ax2, phi, bt_ph, bp_ph)
    ax2.set_xlabel(r'$\phi$ (rad)', fontsize=11, color=INK_SECONDARY)
    ax2.set_ylabel('field component', fontsize=11, color=INK_SECONDARY)
    ax2.set_xticks([0.0, np.pi/2.0, np.pi, 3.0*np.pi/2.0, 2.*np.pi])
    ax2.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'],fontsize=11,color=INK_SECONDARY)
    ax2.set_yticks([0.5*i for i in range(0,5)])
    ax2.set_yticklabels([0.5*i for i in range(0,5)],fontsize=11,color=INK_SECONDARY)   
    ax2.set_title(r'$\phi$-dependence at $\theta=\pi/3$',
                  color=INK_PRIMARY, fontsize=12)
    ax2.legend(frameon=True, fontsize=11, labelcolor=INK_SECONDARY, ncol=2)

    fig.tight_layout()
    savefig_pair(fig, 'vsh_decomposition_tor_profiles')


if __name__ == '__main__':
    plot_profiles()
