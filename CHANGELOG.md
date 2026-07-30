# Changelog

All notable changes to FORTVSH are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**Release convention:**
Bump `VERSION` in `CMakeLists.txt` → update this file → commit → merge to
`master` → `git tag -a vX.Y.Z -m "Release vX.Y.Z: <summary>"`.
PATCH for bug fixes · MINOR for new API additions · MAJOR for breaking changes.

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
