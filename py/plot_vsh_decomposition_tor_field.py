#!/usr/bin/env python3
"""
3D visualization of the Part B target field from
src/examples/vsh_decomposition_tor.f90 (FIELD_QUAD_BUMP): a sphere colored
by field magnitude |B|, with vector arrows overplotted at sample points to
show the local tangent-field direction. Companion, qualitative-structure
counterpart to py/plot_vsh_decomposition_tor.py's quantitative
accuracy/convergence figure -- and to py/plot_gaussian_bump_sphere.py,
which visualizes example 2d's (purely radial) target field the same way.

The field is purely tangential (B_r=0 always, by construction -- see
src/examples/vsh_decomposition_tor.f90 and its companion .md). The bump
is an equator-centered Gaussian, generated from T(theta,phi) =
exp(-(theta-pi/2)**2/2sigma**2) * sin(theta)**2 * (1+0.5*cos(phi)) via
B = r_hat x grad_perp(T), so that (unlike an earlier, pole-centered
version of this field) B_theta and B_phi vanish *exactly* at both poles
-- every term below carries an explicit sin(theta), sin(theta)**2, or
sin(2*theta) factor:
    B_theta(theta,phi) = 0.5*exp(-(theta-pi/2)**2/2sigma**2)*sin(theta)*sin(phi)
    B_phi(theta,phi)   = A_QUAD * PVSH_TOR(2,0,theta,phi)_phi
                          + exp(-(theta-pi/2)**2/2sigma**2) *
                            [-(theta-pi/2)/sigma**2*sin(theta)**2 + sin(2*theta)]
                            * (1+0.5*cos(phi))
sigma=0.4, A_QUAD=1.0, matching the Fortran source exactly. The
quadrupole background needed no changes: PVSH_TOR(2,0,...)'s phi-component
is already proportional to sin(theta)*cos(theta), zero at both poles for
any axisymmetric toroidal mode.

A subtlety worth being explicit about: this field is complex-valued (the
same complex VSH convention the whole library and Part A/B's tests use),
and its two pieces are systematically out of phase with each other --
verified directly against the library's own PVSH_TOR via
py/mpmath_reference.py's vsh_tor_ref: the axisymmetric (l=2,m=0)
quadrupole term is purely IMAGINARY (PVSH_TOR's built-in i-phase
convention has no counterpart to cancel it at m=0), while the Gaussian
bump term was constructed directly as a real closed form. There is no
single global phase that makes both pieces simultaneously real, so this
plot draws arrows from Re(B) (the theta=0-phase snapshot -- a standard,
defensible convention for visualizing one instant of a complex mode
expansion) while coloring the sphere by the full complex magnitude |B|
(phase-independent, so it reflects both pieces properly). The practical
upshot: the arrows are dominated by the bump's structure (the
quadrupole's real-part contribution is exactly zero at this phase), but
the quadrupole's presence is still visible in the background color as a
band peaking at mid-latitudes (~sin(theta)cos(theta)), vanishing at both
the poles and the equator.

PNG only (no companion PDF), matching plot_gaussian_bump_sphere.py's
one-off-illustrative-figure convention.

Run from the project root:
    python3 py/plot_vsh_decomposition_tor_field.py
"""

import os
import sys
import types
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Same mpl_toolkits fix as plot_gaussian_bump_sphere.py -- see that
# script's comments for the full explanation (a system/pip matplotlib
# version conflict in this environment, not something to work around
# differently here).
_pip_mpl_toolkits = os.path.join(
    os.path.dirname(os.path.dirname(matplotlib.__file__)), 'mpl_toolkits')
_mpl_toolkits_pkg = types.ModuleType('mpl_toolkits')
_mpl_toolkits_pkg.__path__ = [_pip_mpl_toolkits]
sys.modules['mpl_toolkits'] = _mpl_toolkits_pkg

from mpl_toolkits.mplot3d import Axes3D  # noqa: E402
import matplotlib.projections as _projections  # noqa: E402
_projections.register_projection(Axes3D)

from plotstyle import SURFACE, INK_PRIMARY, INK_MUTED, GRIDLINE, MAGENTA  # noqa: E402

SIGMA = 0.4
A_QUAD = 1.0


def field_quad_bump(theta, phi):
    """
    (B_theta, B_phi) of FIELD_QUAD_BUMP, complex-valued, matching
    src/examples/vsh_decomposition_tor.f90 exactly. Closed form for the
    quadrupole's phi-component verified against mpmath_reference.py's
    vsh_tor_ref(2,0,...) (exact match, both magnitude and phase). The
    bump's closed form was verified against a finite-difference check on
    its generating function T, including exact vanishing at theta=0,pi.
    """
    b_theta_quad = 0j  # PVSH_TOR(2,0,...) theta-component is exactly 0 (m=0)
    b_phi_quad = 1j * A_QUAD * 3*np.sin(theta)*np.cos(theta) \
        * np.sqrt(5/(4*np.pi)) / np.sqrt(6)

    env = np.exp(-(theta-np.pi/2)**2 / (2.0*SIGMA**2))
    sin2th = np.sin(theta)**2
    b_theta_bump = 0.5 * env * np.sin(theta) * np.sin(phi)
    b_phi_bump = env * (-(theta-np.pi/2)/SIGMA**2 * sin2th + np.sin(2*theta)) \
        * (1 + 0.5*np.cos(phi))

    return b_theta_quad + b_theta_bump, b_phi_quad + b_phi_bump


def plot_field():
    n_theta, n_phi = 140, 200
    theta = np.linspace(0.0, np.pi, n_theta)
    phi = np.linspace(0.0, 2.0*np.pi, n_phi)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing='ij')

    b_theta, b_phi = field_quad_bump(theta_grid, phi_grid)
    magnitude = np.sqrt(np.abs(b_theta)**2 + np.abs(b_phi)**2)

    x = np.sin(theta_grid) * np.cos(phi_grid)
    y = np.sin(theta_grid) * np.sin(phi_grid)
    z = np.cos(theta_grid)

    norm = plt.Normalize(vmin=0.0, vmax=magnitude.max())
    colors = plt.cm.viridis(norm(magnitude))

    fig = plt.figure(figsize=(8.5, 7.2), facecolor=SURFACE)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(SURFACE)

    ax.plot_surface(x, y, z, facecolors=colors, rstride=1, cstride=1,
                     linewidth=0, edgecolor='none', antialiased=False,
                     shade=False, zorder=1)

    # ── Tangent-direction arrows at sample points ───────────────────────
    # Re(B) is used for arrow direction -- see module docstring for why
    # (the field is complex; this is the theta=0-phase snapshot).
    #
    # Restricted to the camera-facing hemisphere: mplot3d does its own
    # internal depth-sorting between artist types (a quiver collection vs.
    # a plot_surface) and doesn't reliably occlude one against the other
    # the way a true 3D renderer would -- arrows on the far side ended up
    # rendering as a tangled cluster in front of the sphere instead of
    # being hidden behind it. Filtering to points whose outward normal
    # faces the camera sidesteps the ambiguity rather than fighting it
    # (the same lesson learned with the wall-projection curves in
    # plot_gaussian_bump_sphere.py, applied here at the point level).
    elev, azim = 20, 40
    elev_r, azim_r = np.radians(elev), np.radians(azim)
    view_dir = np.array([np.cos(elev_r)*np.cos(azim_r),
                          np.cos(elev_r)*np.sin(azim_r),
                          np.sin(elev_r)])

    theta_arrows = np.linspace(0.2, np.pi-0.2, 14)
    phi_arrows = np.linspace(0.0, 2.0*np.pi, 24, endpoint=False)
    TH_A, PH_A = np.meshgrid(theta_arrows, phi_arrows, indexing='ij')

    normal = np.stack([np.sin(TH_A)*np.cos(PH_A),
                        np.sin(TH_A)*np.sin(PH_A),
                        np.cos(TH_A)], axis=-1)
    facing = normal @ view_dir > 0.15  # small margin back from the limb
    TH_A, PH_A = TH_A[facing], PH_A[facing]

    bt_a, bp_a = field_quad_bump(TH_A, PH_A)
    bt_a, bp_a = bt_a.real, bp_a.real

    # theta_hat and phi_hat in Cartesian, at each arrow anchor point.
    theta_hat = np.stack([np.cos(TH_A)*np.cos(PH_A),
                           np.cos(TH_A)*np.sin(PH_A),
                           -np.sin(TH_A)], axis=-1)
    phi_hat = np.stack([-np.sin(PH_A), np.cos(PH_A), np.zeros_like(PH_A)],
                        axis=-1)
    vec = bt_a[..., None]*theta_hat + bp_a[..., None]*phi_hat

    r_anchor = 1.02  # slightly proud of the surface to avoid z-fighting
    xa = r_anchor * np.sin(TH_A)*np.cos(PH_A)
    ya = r_anchor * np.sin(TH_A)*np.sin(PH_A)
    za = r_anchor * np.cos(TH_A)

    arrow_scale = 0.35
    ax.quiver(xa, ya, za, vec[..., 0], vec[..., 1], vec[..., 2],
              length=arrow_scale, normalize=False, color=MAGENTA,
              linewidth=1.3, zorder=5)

    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    mappable.set_array(magnitude)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.55, pad=0.06)
    cbar.set_label('|B| (field magnitude)', color=INK_PRIMARY, fontsize=10)
    cbar.ax.tick_params(colors=INK_PRIMARY, labelsize=8)

    ax.set_title(
        'Part B target field: quadrupole + toroidal Gaussian bump\n'
        'color = |B|, arrows = Re(B) local direction',
        color=INK_PRIMARY, fontsize=11)
    # Coordinate box (panes + ticks + labels) instead of a bare floating
    # sphere -- gives a visual sense of scale/orientation. Panes and grid
    # lines styled to match the project's light chrome rather than
    # matplotlib's default mid-gray.
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRIDLINE)
        axis._axinfo['grid']['color'] = GRIDLINE
        axis._axinfo['grid']['linewidth'] = 0.6
    ax.set_xlabel('x', color=INK_PRIMARY, fontsize=9)
    ax.set_ylabel('y', color=INK_PRIMARY, fontsize=9)
    ax.set_zlabel('z', color=INK_PRIMARY, fontsize=9)
    ax.tick_params(colors=INK_MUTED, labelsize=7)
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=20, azim=40)

    fig.tight_layout()
    outpath = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'tex', 'Copernicus-EGU', 'figures', 'vsh_decomposition_tor_field.png')
    fig.savefig(outpath, dpi=200, facecolor=SURFACE, bbox_inches='tight')
    print(f'Wrote {outpath}')


if __name__ == '__main__':
    plot_field()
