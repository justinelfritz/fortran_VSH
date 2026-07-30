---
project: FORTVSH
summary: A Fortran library for scalar and vector spherical harmonics (VSH), Legendre functions, and angular-momentum coupling coefficients, built for spectral MHD applications in geophysics and astrophysics.
author: Justin G. Elfritz
email: 45849387+justinelfritz@users.noreply.github.com
project_github: https://github.com/justinelfritz/FORTVSH
version: 0.3.0
year: 2026
src_dir: src
output_dir: docs
display: public
         protected
         private
source: true
graph: true
coloured_edges: true
search: true
sort: permission-alpha
---

FORTVSH is a lightweight, physics-first Fortran library for computing
scalar spherical harmonics (SSH), two complementary vector spherical
harmonic (VSH) bases, and the Legendre functions and angular-momentum
coupling coefficients (Clebsch-Gordan, Wigner 3-j, and Geppert-Wiebicke)
that underlie them.

## Two vector bases, one scalar harmonic

Both VSH bases are built from the same underlying scalar spherical
harmonic \( Y_\ell^m \):

- The **polar basis** ([[PVSH_RAD]]/[[PVSH_POL]]/[[PVSH_TOR]]), indexed
  directly by \( (\ell,m) \) -- the natural choice for the
  poloidal/toroidal decomposition of a divergence-free vector field in
  MHD.
- The **standard basis** ([[VSH_TOR]]/[[VSH_POL_DN]]/[[VSH_POL_UP]]),
  indexed by the total-angular-momentum-coupled \( J=\ell,\ell\mp1 \) --
  the form familiar from the quantum-mechanical vector-spherical-harmonic
  literature.

Every routine has a single-point form (one mode at one point) and a
batch `_ALL` form (every mode at one point, in a single call).

## Where to start

- [[VSH]] -- the library itself: every spherical harmonic, Legendre
  function, and coupling coefficient routine.
- [[TESTS]] -- the numerical validation suite backing the accuracy
  claims in the accompanying manuscript.
- [[BENCHMARK]] -- timing comparisons between the batch `_ALL` routines
  and naive per-mode loops.
- `src/examples/` -- five worked examples, from a basic magnetic-dipole
  synthesis through a full arbitrary-field spectral decomposition and a
  Geppert-Wiebicke coupling-coefficient sweep.

## Building

```bash
cmake -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

See `README.md` in the repository root for the full build, installation,
and usage guide, including the optional benchmark (`VSH_BUILD_BENCHMARK`)
and example (`VSH_BUILD_EXAMPLES`) targets.
