#!/usr/bin/env python3
"""
3D visualization of the Gaussian bump test field from
src/examples/vsh_decomposition.f90's Part B (FIELD_GAUSSIAN, the smooth,
effectively-infinite-spectral-content field used to demonstrate spectral
convergence):

    BUMP(theta, phi) = exp(-theta**2 / (2*sigma**2))
                        * (1 + 0.5*sin(theta)*cos(phi)),   sigma = 0.4

(purely radial: F = (BUMP, 0, 0), matching FIELD_GAUSSIAN exactly). The
sphere's radius is modestly displaced by the field value so the bump
reads as real 3D geometry rather than a flat colored sphere, and colored
by that same value on a red-blue scale, so shape and color reinforce the
same magnitude information.

Two wall-projected profile curves (matplotlib's zdir/offset "flatten onto
a wall" technique, same idea as the matplotlib gallery's 3D contour-walls
example) show the theta- and phi-dependence separately:
  - back wall (constant y): <B>_phi(theta), the phi-AVERAGED field -- this
    is exactly the pure Gaussian envelope exp(-theta**2/2sigma**2), since
    cos(phi) integrates to zero over a full period, so it isolates the
    m=0 (phi-independent) content.
  - side wall (constant x): B(theta_peak, phi), a raw slice through the
    field's true peak colatitude (theta_peak ~ 0.0768 rad, found earlier
    by direct optimization -- NOT exactly the pole, since the cos(phi)
    term's linear growth briefly outpaces the Gaussian's quadratic falloff
    near theta=0), showing the cos(phi) modulation directly.
Unlike the gallery example, these are ax.plot curves rather than
ax.contour: a sphere's surface isn't a single-valued function of two
Cartesian axes the way the gallery's z=f(x,y) test data is, so a literal
ax.contour(X,Y,Z,zdir=...) call doesn't apply here -- contour needs a
genuine 2D height field, and these are 1D profiles. The zdir/offset
wall-placement mechanism is identical either way.

PNG only (no companion PDF) -- a one-off illustrative figure, not part of
the matched PDF/PNG figure-pair convention the other plot_*.py scripts
follow.

Run from the project root:
    python3 py/plot_gaussian_bump_sphere.py
"""

import os
import sys
import types
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# This environment has two mpl_toolkits installs on sys.path: a pip one
# matching the installed matplotlib, and an incompatible apt/dist-packages
# one that Python's import system prefers regardless of sys.path order
# (it has a real __init__.py, so it wins over the pip install's namespace
# package, per Python's package-resolution rules). Point mpl_toolkits at
# the correct directory explicitly -- the one next to matplotlib itself,
# which we know resolved correctly -- scoped to this process only; doesn't
# touch sys.path or anything else dist-packages legitimately provides
# (e.g. cycler).
_pip_mpl_toolkits = os.path.join(
    os.path.dirname(os.path.dirname(matplotlib.__file__)), 'mpl_toolkits')
_mpl_toolkits_pkg = types.ModuleType('mpl_toolkits')
_mpl_toolkits_pkg.__path__ = [_pip_mpl_toolkits]
sys.modules['mpl_toolkits'] = _mpl_toolkits_pkg

from mpl_toolkits.mplot3d import Axes3D  # noqa: E402
import matplotlib.projections as _projections  # noqa: E402
_projections.register_projection(Axes3D)  # matplotlib.pyplot's own import
# attempt (above) already ran and silently failed before our patch, which
# leaves the '3d' projection unregistered even though the class import
# now succeeds -- register it explicitly.

from plotstyle import SURFACE, INK_PRIMARY, INK_SECONDARY, AQUA, MAGENTA, \
    FIG_DIR  # noqa: E402

SIGMA = 0.4
L = 2.0  # half-extent of the bounding box the walls sit on


def bump(theta, phi):
    return np.exp(-theta**2 / (2.0 * SIGMA**2)) * \
        (1.0 + 0.5 * np.sin(theta) * np.cos(phi))


def find_theta_peak():
    """theta of the field's true maximum at phi=0 (see prior analysis:
    NOT exactly the pole -- cos(phi)'s linear growth briefly outpaces the
    Gaussian's quadratic falloff near theta=0)."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda t: -bump(t, 0.0), bounds=(0.0, 1.0),
                           method='bounded', options={'xatol': 1e-10})
    return res.x


def plot_gaussian_bump_sphere():
    n_theta, n_phi = 140, 200
    theta = np.linspace(0.0, np.pi, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing='ij')

    field = bump(theta_grid, phi_grid)

    # Subtle radial exaggeration so the bump reads as real geometry without
    # overwhelming the sphere's own shape -- the field itself is purely
    # radial in the source example, so displacing the radius is a faithful
    # reading of "a bump", not an invented effect. Kept small deliberately:
    # a large exaggeration turns the polar Gaussian cap into a dominant
    # egg/teardrop shape rather than a bump on an otherwise round sphere.
    exaggeration = 0.10
    r = 1.0 + exaggeration * (field - field.min())

    x = r * np.sin(theta_grid) * np.cos(phi_grid)
    y = r * np.sin(theta_grid) * np.sin(phi_grid)
    z = r * np.cos(theta_grid)

    norm = plt.Normalize(vmin=field.min(), vmax=field.max())
    colors = plt.cm.RdBu_r(norm(field))

    fig = plt.figure(figsize=(8.5, 7.2), facecolor=SURFACE)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(SURFACE)

    # ── Flat wall panels, drawn first (furthest back) ──────────────────
    # Explicit rectangular surfaces at x=-L and y=-L give the profile
    # curves below a visible plane to sit on -- matplotlib's default 3D
    # axis panes are tied to set_axis_off()/tick behavior in ways that
    # made them unreliable here, so these are drawn directly instead.
    wall_span = np.array([-L, L])
    WY, WZ = np.meshgrid(wall_span, wall_span)
    ax.plot_surface(np.full_like(WY, -L), WY, WZ, color=INK_SECONDARY,
                     alpha=0.06, linewidth=0, shade=False, zorder=0)
    WX, WZ2 = np.meshgrid(wall_span, wall_span)
    ax.plot_surface(WX, np.full_like(WX, -L), WZ2, color=INK_SECONDARY,
                     alpha=0.06, linewidth=0, shade=False, zorder=0)

    ax.plot_surface(x, y, z, facecolors=colors, rstride=1, cstride=1,
                     linewidth=0, edgecolor='none', antialiased=False,
                     shade=False, zorder=1)

    # ── Wall-projected profile curves ──────────────────────────────────
    # Map each profile's (angle, value) onto the wall's two in-plane
    # coordinates: angle spans the wall's full width [-L, L], value rises
    # from a shared floor at z=-L. Explicit fmax normalizes both curves to
    # the same height scale for fair visual comparison.
    #
    # Two things were needed to keep the curves from being cut through by
    # the opaque sphere: (1) height is capped so the curves' Z-range stays
    # entirely below the sphere's own Z-extent (~+-1.1) -- matplotlib's
    # mplot3d does its own internal depth-sorting between artists and
    # largely ignores manually set zorder for a plot_surface vs. a line,
    # so relying on zorder alone doesn't work; (2) that alone wasn't
    # sufficient either -- the (x=-L, y=-L) wall corner line sits directly
    # behind the sphere (at the origin) for some camera azimuths regardless
    # of Z-separation, since the sphere's near-camera silhouette is
    # magnified by perspective and can still cover a farther-back curve
    # even at a different height. view_init's azim=125 below was chosen
    # specifically because it avoids that alignment -- the two are only
    # verified compatible together, not independently.
    fmax = field.max()
    height = 0.7

    theta_profile = np.array([np.trapz(bump(t * np.ones_like(phi), phi),
                                        phi) / (2 * np.pi) for t in theta])
    x_wall_theta = -L + (theta / np.pi) * (2 * L)
    z_wall_theta = -L + (theta_profile / fmax) * height
    ax.plot(x_wall_theta, z_wall_theta, zs=-L, zdir='y',
            color=AQUA, linewidth=2.2, zorder=3,
            label=r'$\langle B\rangle_\phi(\theta)$ (= pure Gaussian envelope)')

    theta_peak = find_theta_peak()
    phi_profile = bump(theta_peak, phi)
    y_wall_phi = -L + (phi / (2.0 * np.pi)) * (2 * L)
    z_wall_phi = -L + (phi_profile / fmax) * height
    ax.plot(y_wall_phi, z_wall_phi, zs=-L, zdir='x',
            color=MAGENTA, linewidth=2.2, zorder=3,
            label=r'$B(\theta_{peak},\phi)$, $\theta_{peak}\approx0.077$')

    mappable = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    mappable.set_array(field)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.55, pad=0.02)
    cbar.set_label('field amplitude', color=INK_PRIMARY, fontsize=10)
    cbar.ax.tick_params(colors=INK_PRIMARY, labelsize=8)

    ax.set_title(
        'Gaussian bump test field on the sphere\n'
        r'$\exp(-\theta^2/2\sigma^2)\,(1+0.5\sin\theta\cos\phi)$, '
        r'$\sigma=0.4$',
        color=INK_PRIMARY, fontsize=11)
    ax.legend(loc='upper left', fontsize=8, frameon=False,
              labelcolor=INK_SECONDARY, bbox_to_anchor=(0.02, 0.98))
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-L, L)
    ax.set_axis_off()
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=16, azim=125)

    fig.tight_layout()
    outpath = os.path.join(FIG_DIR, 'gaussian_bump_sphere.png')
    fig.savefig(outpath, dpi=200, facecolor=SURFACE, bbox_inches='tight')
    print(f'Wrote {outpath}')


if __name__ == '__main__':
    plot_gaussian_bump_sphere()
