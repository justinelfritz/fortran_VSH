# Changelog

All notable changes to FORTVSH are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**Release convention:**
Bump `VERSION` in `CMakeLists.txt` → update this file → commit → merge to
`master` → `git tag -a vX.Y.Z -m "Release vX.Y.Z: <summary>"`.
PATCH for bug fixes · MINOR for new API additions · MAJOR for breaking changes.

---

## [0.3.1] — 2026-08-02

### Added
- GitHub Actions CI (`.github/workflows/ci.yml`): `build-and-test` job on every
  push/PR — one combined configure with every opt-in CMake option enabled,
  `ctest`, running all worked examples, a benchmark smoke-test, install +
  pkg-config + a `find_package(FORTVSH)` downstream-consumer check.
- GitHub Actions docs publishing (`.github/workflows/docs.yml`): `build-docs`
  (FORD-warning regression gate on push to `master`/`develop` and on every PR)
  and `deploy-docs` (publishes HTML API docs to GitHub Pages on push to
  `master` only).
- `VSH_BUILD_CONVERGENCE` CMake option (default `OFF`): `vsh_convergence`
  executable (`src/convergence_main.f90`), an accuracy sweep feeding the
  manuscript's convergence figures. Not part of `ctest`.
- `VSH_BUILD_STABILITY` CMake option (default `OFF`): `vsh_stability` /
  `vsh_unnorm_stability` executables (`src/stability_main.f90`,
  `src/unnorm_stability_main.f90`) mapping the numerically safe `(l, m,
  theta)` range of the normalized and unnormalized Legendre recurrences,
  respectively. Not part of `ctest`.
- `STABILITY_FINDINGS.md`: full writeup of both stability investigations —
  methodology, per-degree onset tables, and the practical safe-range
  guidance now reflected in the affected routines' docstrings.
- `py/mpmath_reference.py`: independent arbitrary-precision reference
  (Legendre/SSH/VSH, via a high-precision rerun of the library's own Bonnet
  recurrence rather than a separate ground-truth implementation) backing
  both the stability investigation and the rebuilt convergence-accuracy
  figures.
- `py/plot_convergence_unnormalized.py`, `py/plot_stability.py`: new
  figures (`convergence_accuracy_unnormalized`, `stability_map`).
- `py/requirements.txt`: first explicit Python dependency manifest for the
  `py/` scripts (`numpy`, `matplotlib`, `mpmath`).

### Fixed
- `ASSOC_LEGENDRE_NORM_ALL` and the 10 batch routines built on it
  (`DDX_ASSOC_LEGENDRE_NORM_ALL`, `SSH_ALL`, `GRAD_SSH_ALL`, `L_SSH_ALL`,
  `PVSH_RAD_ALL`, `PVSH_POL_ALL`, `PVSH_TOR_ALL`, `VSH_TOR_ALL`,
  `VSH_POL_UP_ALL`, `VSH_POL_DN_ALL`): docstrings corrected from an
  unverified "stable to l~2700" claim to a verified range — unconditionally
  safe to `l<=2000`; beyond that, safety depends on `theta` (see
  `STABILITY_FINDINGS.md`).
- `ASSOC_LEGENDRE`, `ASSOC_LEGENDRE_ALL`, and the 12 single-mode/derivative
  routines built on the unnormalized recurrence: docstrings corrected from
  an unverified "stable to l~1400" claim to a verified `l<=150` safe range
  (worst case at the equator).
- `vsh_benchmark`'s naive-loop timing measurement: each naive-loop pass
  costs O(NX·Lmax³) (every single-mode call is itself an O(l) recurrence),
  so extending the sweep to `Lmax=4000` made a single measurement run for
  hours. Added `LOOP_LMAX_CAP=200` (`src/benchmark.f90`) so naive-loop
  timing stops there; batch timing is unaffected and still runs to
  `Lmax=2000`.

### Changed
- `convergence_accuracy` figure rebuilt: now diffs the normalized batch
  routines against the independent arbitrary-precision reference (rather
  than diffing batch against single-mode, which only showed where the
  *less* stable single-mode path degraded, not whether either was close to
  the true value), scoped to the verified-safe `l<=2000` range. The
  unnormalized family got its own new figure,
  `convergence_accuracy_unnormalized`, scoped to `l<=150`.
- `benchmark_scaling` / `benchmark_speedup` figures: each series now
  truncated to its own verified-safe range (naive loop `l<=150`, batch
  `l<=2000`) instead of both plotted uniformly out to `Lmax=4000`.
- `README.md`: corrected the `ASSOC_LEGENDRE_NORM_ALL` API table's stale
  "stable to l ≈ 2700" claim and the Benchmarking section's `Lmax` sweep
  description to match the current sweep and the naive-loop cap; added a
  Numerical stability section pointing to `STABILITY_FINDINGS.md`.
- `.gitignore`: added `/benchmark/`, `/convergence/`, `/stability/`
  (regenerable output from the three opt-in CMake targets above, same
  category as the existing `/validation/*.dat` entry) and `__pycache__/`.

---

## [0.3.0] — 2026-07-30

### Added
- FORD ([forddocs.readthedocs.io](https://forddocs.readthedocs.io/)) HTML API
  documentation: `project.md` project file plus full `!>`-style docstrings
  (LaTeX-rendered governing equations, per-argument descriptions, and
  cross-references) across every public routine in `vsh.f90`, `tests.f90`,
  `benchmark.f90`, and the five example programs. Build via
  `pip install ford && ford project.md`; the generated `docs/` output is not
  committed.
- `VSH_BUILD_BENCHMARK` CMake option (default `OFF`): builds `vsh_benchmark`
  (`src/benchmark.f90` + `src/bench_main.f90`), timing every batch (`_ALL`)
  routine against a naive per-mode loop across a sweep of `Lmax`, with each
  measurement auto-calibrated against wall-clock time. Not part of `ctest`
  (timings are machine/compiler-dependent) — a manually-run target writing to
  `benchmark/`.
- `VSH_BUILD_EXAMPLES` CMake option (default `OFF`): builds five worked-example
  executables under `src/examples/`: `dipole_synthesis`,
  `dipole_uniform_sphere`, and `dipole_flux_expulsion` (a static poloidal-dipole
  magnetostatics progression — field synthesis and mode-purity check,
  uniformly-magnetized-sphere boundary matching, and idealized
  superconducting-core flux expulsion), `gwi_gwj_sweep` (Geppert-Wiebicke
  coupling-coefficient validation across swept parameter ranges), and
  `vsh_decomposition` (arbitrary-field VSH spectral decomposition,
  reconstruction, and convergence — the first routine in the repo to exercise
  `SHGLQ`). Each prints its own `PASS`/`FAIL` self-check and writes to
  `examples/<name>/`.
- `py/plotstyle.py`: shared matplotlib styling (validated colorblind-safe
  palette) used by all new plotting scripts, plus `py/plot_benchmark.py` and
  one `py/plot_*.py` script per new example, feeding figures into the
  accompanying manuscript.

### Fixed
- `GWJ_DP`: guarded a previously-unguarded `sqrt()` term that could go
  negative (producing `NaN`, or crashing under `-ffpe-trap` builds) for
  out-of-range `(J1,J2,L)` mode triples, matching the triangle-inequality
  guard `CGCOEFF` already had. Purely additive — identical results for every
  previously-valid input.
- Split ~62 grouped dummy-argument declarations (e.g.
  `INTEGER, INTENT(IN) :: L, K`) into one variable per line throughout
  `vsh.f90`, `tests.f90`, `benchmark.f90`, and `vsh_decomposition.f90`. Works
  around a FORD bug where grouped declarations caused each variable's
  generated documentation to accumulate the descriptions of all preceding
  variables in the same group. Pure declaration-syntax change — no behavior
  change.

---

## [0.2.0] — 2026-06-06

### Added
- Formal pass/fail exit code from `vsh_test`: all 14 test subroutines in
  `tests.f90` accept a new `INTEGER(KIND=i4), INTENT(OUT) :: STATUS` argument
  and evaluate numerical thresholds internally (relative error for
  `BATCH_ALM_CONS` / `BATCH_DALM_CONS`; absolute error for all others;
  quadrature-aware tolerance for `SSH_ORTHO`).
- `main.f90` rewritten as a proper test driver: collects `STATUS` from each
  subroutine, checks `TEST1` / `TEST2` agreement with `MAXVAL(ABS(...))`,
  prints a per-test `PASS` / `FAIL` line, and exits via `STOP 1` if any test
  fails — making ctest a meaningful go/no-go signal.
- CMake fixture `create_validation_dir` runs before `vsh_validation` and
  creates `validation/` if absent, so `ctest` succeeds in a fresh clone.
- `VSH_REGRESSION_TEST` CMake option (default `OFF`): when enabled, adds a
  `vsh_regression` ctest that compares every generated `.dat` file against its
  committed reference counterpart using `cmake/check_regression.py`.
- `cmake/check_regression.py`: Python script for column-by-column numeric
  comparison of two directories of `.dat` files to a relative tolerance of
  10⁻¹⁰.
- `validation/reference/`: 16 committed reference `.dat` files serving as the
  regression baseline.
- `py/compute_validation.py`: Python script that reads all validation output
  files, computes per-test maximum errors, and writes
  `validation/validation_values.tex` — a set of `\newcommand` definitions for
  direct inclusion in the whitepaper.
- `py/plotError.py`: matplotlib script for plotting `TEST1` percentage
  differences between the analytic and VSH-inner-product results.

### Changed
- `validation/*.dat` added to `.gitignore`; generated outputs are no longer
  tracked (reference copies live in `validation/reference/`).
- README Testing section updated: automated pass/fail output described,
  `batch_ssh_cons.dat` column layout corrected, optional regression comparison
  documented.

---

## [0.1.0] — 2026-06-06

First versioned release. The VSH mathematics are unchanged from the original
library; this release modernises the Fortran source and build infrastructure
to produce an installable, standards-conforming package.

### Added
- `VSH_VERSION` Fortran module exposing `VERSION_MAJOR`, `VERSION_MINOR`,
  `VERSION_PATCH`, and `VERSION_STRING` constants (generated by CMake at
  configure time from the `project(VERSION ...)` declaration).
- CMake build system (`CMakeLists.txt`) replacing the legacy makefile:
  static library `libvsh.a`, optional shared `libvsh.so` (`-DVSH_BUILD_SHARED=ON`),
  `ctest` integration, and `cmake --install` support.
- pkg-config file `fortvsh.pc` installed to `lib/pkgconfig/`.
- CMake package config enabling `find_package(FORTVSH REQUIRED)` /
  `target_link_libraries(myapp FORTVSH::vsh)` in downstream projects.
- Precision-generic `INTERFACE` blocks for all 33 public routines; existing
  dp implementations renamed to `X_DP` specifics. Future sp variants require
  only a new `X_SP` routine and one `MODULE PROCEDURE X_SP` line per block.
- Library/executable separation: KINDS + GLOBALS + VSH compile into
  `libvsh.a`; TESTS + MAIN link against the library and are not installed.

### Changed
- Source converted from fixed-form (`.f`) to free-form Fortran 90 (`.f90`).
- `IMPLICIT NONE` added throughout all modules and subroutines.
- `INTENT(IN)`, `INTENT(OUT)`, `INTENT(INOUT)` added to all dummy arguments
  across `vsh.f90` and `tests.f90`.
- `PRIVATE` default + explicit `PUBLIC` declarations in every module
  (KINDS, GLOBALS, VSH, TESTS) to document and enforce the public API surface.
- All `USE` statements updated to `USE … ONLY:` to prevent namespace pollution.
