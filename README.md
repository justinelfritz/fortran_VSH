# FORTVSH

Fortran library for Legendre polynomials and derivatives, scalar spherical
harmonics (SSH) and angular derivatives, vector spherical harmonics (VSH)
and angular derivatives, Clebsch-Gordan and Wigner 3j coefficients, and
Gauss-Legendre quadrature support.

Supporting modules `KINDS` and `GLOBALS` provide datatype definitions and
shared physical constants (π, imaginary unit j).

Please review the LaTeX whitepaper in `/tex/` for normalization conventions
and sign choices that propagate through all numerical routines (not yet available).

---

## Requirements

| Tool | Minimum version |
|------|----------------|
| gfortran | 9 |
| CMake | 3.14 |

No external Fortran libraries are required. CMake is the only build system
dependency. pkg-config is optional (only needed to consume the installed library
via `pkg-config --cflags --libs FORTVSH`).

---

## Quick start

```bash
cmake -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

---

## Building from source

An out-of-source build is recommended so that generated files do not mix with
source files.

```bash
# Configure (Release mode by default)
cmake -B build

# Optional: debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug

# Optional: also build shared library (libvsh.so)
cmake -B build -DVSH_BUILD_SHARED=ON

# Compile
cmake --build build

# Run the validation suite
ctest --test-dir build --output-on-failure
```

Build artifacts:

| Path | Contents |
|------|----------|
| `build/libvsh.a` | Static library |
| `build/libvsh.so` | Shared library (if `-DVSH_BUILD_SHARED=ON`) |
| `build/mod/` | Fortran module files (`.mod`) for the installed library |
| `build/obj/` | Module files for test-only code (not installed) |
| `build/vsh_test` | Test executable |

---

## Installing

```bash
cmake --install build --prefix /usr/local
```

This copies:

| Destination | Contents |
|-------------|----------|
| `<prefix>/lib/libvsh.a` | Static library |
| `<prefix>/include/vsh/` | Module files: `kinds.mod`, `globals.mod`, `vsh.mod`, `vsh_version.mod` |
| `<prefix>/lib/pkgconfig/fortvsh.pc` | pkg-config file |
| `<prefix>/lib/cmake/FORTVSH/` | CMake package config files |

---

## Using the library

### From a CMake project

```cmake
find_package(FORTVSH REQUIRED)
target_link_libraries(myapp PRIVATE FORTVSH::vsh)
```

Set `FORTVSH_DIR` (or `CMAKE_PREFIX_PATH`) to the install prefix if
`find_package` cannot locate the library automatically:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/usr/local
```

### From pkg-config / manual linking

```bash
# Compiler and linker flags
pkg-config --cflags --libs FORTVSH
# → -I/usr/local/include/vsh  -L/usr/local/lib  -lvsh

# Example manual compile
gfortran $(pkg-config --cflags FORTVSH) mycode.f90 $(pkg-config --libs FORTVSH) -o myapp
```

### Minimal usage example

```fortran
PROGRAM EXAMPLE
  USE VSH,         ONLY: SSH, GRAD_SSH
  USE VSH_VERSION, ONLY: VERSION_STRING
  USE KINDS,       ONLY: dp
  IMPLICIT NONE
  COMPLEX(KIND=dp) :: Y, GY(3)
  PRINT *, 'FORTVSH version: ', VERSION_STRING
  Y  = SSH(2, 1, 1.0_dp, 0.5_dp)       ! Y_2^1(θ=1, φ=0.5)
  GY = GRAD_SSH(2, 1, 1.0_dp, 0.5_dp)  ! ∇_Ω Y_2^1
  PRINT *, 'Y_2^1 = ', Y
END PROGRAM EXAMPLE
```

---

## Source layout

| File | Module | Purpose |
|------|--------|---------|
| `src/kinds.f90` | `KINDS` | Kind parameters (`sp`, `dp`, `i4`, …) |
| `src/globals.f90` | `GLOBALS` | Shared constants (π, j) |
| `src/vsh.f90` | `VSH` | All spherical harmonic and VSH routines |
| `src/vsh_version.f90.in` | `VSH_VERSION` | Version constants; filled by CMake at configure time |
| `src/tests.f90` | `TESTS` | Validation and consistency test routines (not installed) |
| `src/main.f90` | — | Driver: calls test routines, writes output (not installed) |

`TESTS` and `main.f90` are compiled only into the `vsh_test` executable and
are not part of `libvsh.a` or any installed target.

---

## API reference

All public names are generic interfaces: the current implementations are
double-precision (`dp`) specifics named `X_DP`. Adding a single-precision
variant later requires only writing `X_SP` and inserting
`MODULE PROCEDURE X_SP` into the existing `INTERFACE` block — no call-site
changes.

### Version

| Name | Type | Description |
|------|------|-------------|
| `VSH_VERSION % VERSION_STRING` | `CHARACTER(LEN=*)` | Full version string, e.g. `"0.1.0"` |
| `VSH_VERSION % VERSION_MAJOR` | `INTEGER(i4)` | Major version |
| `VSH_VERSION % VERSION_MINOR` | `INTEGER(i4)` | Minor version |
| `VSH_VERSION % VERSION_PATCH` | `INTEGER(i4)` | Patch version |

### Legendre polynomials

| Name | Type | Description |
|------|------|-------------|
| `LEGENDRE(L, X)` | Function → `dp` | Legendre polynomial P_l(x) |
| `DDX_LEGENDRE(L, X)` | Function → `dp` | d/dx P_l(x) |
| `ASSOC_LEGENDRE(L, K, X)` | Function → `dp` | Associated Legendre P_l^k(x), Condon-Shortley |
| `DDX_ASSOC_LEGENDRE(L, K, X)` | Function → `dp` | d/dx P_l^k(x) |

### Batch Legendre (0 ≤ l ≤ lmax, 0 ≤ m ≤ l)

Output arrays use `PLM_INDEX(l,m) = l*(l+1)/2 + m + 1`.

| Name | Description |
|------|-------------|
| `ASSOC_LEGENDRE_ALL(P, LMAX, X)` | All P_l^m(x), unnormalized, Bonnet recurrence |
| `DDX_ASSOC_LEGENDRE_ALL(DP, P, LMAX, X)` | All d/dx P_l^m(x), requires precomputed `P` |
| `ASSOC_LEGENDRE_NORM_ALL(PNORM, LMAX, X)` | All 4π-normalized P̄_l^m(x) via Holmes & Featherstone (2002) recurrence; stable to l ≈ 2700 |
| `DDX_ASSOC_LEGENDRE_NORM_ALL(DPNORM, PNORM, LMAX, X)` | All d/dx P̄_l^m(x), requires precomputed `PNORM` |

### Scalar spherical harmonics Y_l^m(θ, φ)

| Name | Type | Description |
|------|------|-------------|
| `SSH(L, K, THETA, PHI)` | Function → `complex dp` | Single Y_l^k |
| `GRAD_SSH(L, K, THETA, PHI)` | Function → `complex dp (3)` | Angular gradient ∇_Ω Y_l^k |
| `L_SSH(L, K, THETA, PHI)` | Function → `complex dp (3)` | r̂ × ∇_Ω Y_l^k |

### Batch SSH (0 ≤ l ≤ lmax, −l ≤ m ≤ l)

Output arrays use `YLM_INDEX(l,m) = l² + l + m + 1`.
VSH batch output arrays are shaped `(3, (lmax+1)²)` — first index is
spherical component (r, θ, φ).

| Name | Description |
|------|-------------|
| `SSH_ALL(YLM, LMAX, THETA, PHI)` | All Y_l^m; uses H&F normalized ALPs internally |
| `GRAD_SSH_ALL(OUT, LMAX, THETA, PHI)` | All ∇_Ω Y_l^m |
| `L_SSH_ALL(OUT, LMAX, THETA, PHI)` | All r̂ × ∇_Ω Y_l^m |

### Polar vector spherical harmonics

| Name | Type | Description |
|------|------|-------------|
| `PVSH_RAD(L, K, THETA, PHI)` | Function → `complex dp (3)` | Radial PVSH Ỹ_l^{k(−1)} = r̂ Y_l^k |
| `PVSH_POL(L, K, THETA, PHI)` | Function → `complex dp (3)` | Poloidal PVSH Ỹ_l^{k(+1)} ∝ ∇_Ω Y_l^k |
| `PVSH_TOR(L, K, THETA, PHI)` | Function → `complex dp (3)` | Toroidal PVSH Ỹ_l^{k(0)} ∝ r̂ × ∇_Ω Y_l^k |
| `PVSH_RAD_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch radial PVSH |
| `PVSH_POL_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch poloidal PVSH |
| `PVSH_TOR_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch toroidal PVSH |

### Standard vector spherical harmonics

| Name | Type | Description |
|------|------|-------------|
| `VSH_TOR(L, K, THETA, PHI)` | Function → `complex dp (3)` | Toroidal VSH (= PVSH_TOR) |
| `VSH_POL_UP(L, K, THETA, PHI)` | Function → `complex dp (3)` | J = L+1 poloidal VSH |
| `VSH_POL_DN(L, K, THETA, PHI)` | Function → `complex dp (3)` | J = L−1 poloidal VSH |
| `VSH_TOR_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch toroidal VSH |
| `VSH_POL_UP_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch J = L+1 poloidal VSH |
| `VSH_POL_DN_ALL(OUT, LMAX, THETA, PHI)` | Subroutine | Batch J = L−1 poloidal VSH |

### Angular momentum algebra and utilities

| Name | Type | Description |
|------|------|-------------|
| `CGCOEFF(J1,M1,J2,M2,J3,M3)` | Function → `dp` | Clebsch-Gordan coefficient |
| `SYMBOL3J(J1,M1,J2,M2,J3,M3)` | Function → `dp` | Wigner 3j symbol |
| `GWI(J1,M1,J2,M2,L,M)` | Function → `complex dp` | Geppert-Wiebicke integral I |
| `GWJ(J1,M1,J2,M2,L,M)` | Function → `complex dp` | Geppert-Wiebicke integral J |
| `DOT(VSH1, VSH2)` | Function → `complex dp` | Bilinear dot product of two 3-vectors |
| `LOG_FACT(N)` | Function → `dp` | ln(N!) via log-gamma; avoids factorial overflow |
| `PLM_INDEX(L, M)` | Function → `i4` | 1D index for (l,m) in P batch arrays |
| `YLM_INDEX(L, M)` | Function → `i4` | 1D index for (l,m) in Y/VSH batch arrays |
| `SHGLQ(ZERO, W, LMAX)` | Subroutine | Gauss-Legendre nodes (cos θ) and weights for exact integration of polynomials up to degree 2·lmax+1 |

---

## Testing

The validation suite runs automatically via `ctest`:

```bash
ctest --test-dir build --output-on-failure
```

The test executable writes all output to `validation/`. Every output file
begins with a `#` header line followed by column labels.

### Cross-validation tests (Tests 1 and 2)

`TEST1.DAT` and `TEST2.DAT` compare three independent computations of the
same VSH inner-product expression:

| Column | Contents |
|--------|----------|
| 1 | θ (radians) |
| 2–3 | Analytic result (re, im) |
| 4–5 | VSH inner product (re, im) |
| 6–7 | Geppert-Wiebicke form (re, im) |

All three columns should agree to near machine precision across all θ.

### Batch consistency tests

Each `BATCH_*_CONS` subroutine compares the batch `*_ALL` output against
the corresponding single-mode function for every `(l, m, θ)` in the grid.

| Output file | Compares | Columns |
|-------------|----------|---------|
| `validation/batch_alm_cons.dat` | `ASSOC_LEGENDRE_ALL` vs `ASSOC_LEGENDRE` | `L  M  X  P_batch  P_single  abs_diff` |
| `validation/batch_dalm_cons.dat` | `DDX_ASSOC_LEGENDRE_ALL` vs `DDX_ASSOC_LEGENDRE` | `L  M  X  dP_batch  dP_single  abs_diff` |
| `validation/batch_ssh_cons.dat` | `SSH_ALL` vs `SSH` | `L  M  theta  phi  re(batch)  im(batch)  re(single)  im(single)` |
| `validation/batch_grad_ssh_cons.dat` | `GRAD_SSH_ALL` vs `GRAD_SSH` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_l_ssh_cons.dat` | `L_SSH_ALL` vs `L_SSH` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_pvsh_rad_cons.dat` | `PVSH_RAD_ALL` vs `PVSH_RAD` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_pvsh_pol_cons.dat` | `PVSH_POL_ALL` vs `PVSH_POL` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_pvsh_tor_cons.dat` | `PVSH_TOR_ALL` vs `PVSH_TOR` (l ≥ 1) | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_vsh_tor_cons.dat` | `VSH_TOR_ALL` vs `VSH_TOR` (l ≥ 1) | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_vsh_pol_up_cons.dat` | `VSH_POL_UP_ALL` vs `VSH_POL_UP` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |
| `validation/batch_vsh_pol_dn_cons.dat` | `VSH_POL_DN_ALL` vs `VSH_POL_DN` | `L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph` |

> `PVSH_TOR` and `VSH_TOR` single-mode functions produce 0/0 at l = 0
> (the l = 0 toroidal VSH is identically zero). The batch routines correctly
> return zero; the consistency tests begin at l = 1 to avoid NaN comparisons.

### Mathematical property tests

| Output file | Property verified |
|-------------|-------------------|
| `validation/ssh_ortho.dat` | SSH orthonormality via midpoint quadrature; inner product should be 0 (off-diagonal) or 1 (diagonal) |
| `validation/pvsh_pol_tor_ortho.dat` | `DOT(PVSH_POL, PVSH_TOR) = 0` pointwise for all (l, m, θ); expected value ≤ ~10⁻¹⁵ |
| `validation/vsh_pol_inversion.dat` | VSH_POL_UP/DN invert to PVSH_POL and PVSH_RAD via the orthogonal rotation; columns `max_absdiff_pol` and `max_absdiff_rad` should be ≤ ~10⁻¹⁴ |

### Interpreting results

`vsh_test` checks every test against a numerical threshold and prints a
per-test `PASS` / `FAIL` line followed by an overall `RESULT:` summary.
ctest marks the run failed if the executable exits nonzero, so no manual
inspection of the `.dat` files is needed for a go/no-go decision. The
output files remain available for plotting or further analysis.

### Optional regression comparison

Committed reference outputs live in `validation/reference/`. Enable the
comparison test with:

```bash
cmake -B build -DVSH_REGRESSION_TEST=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

This adds a `vsh_regression` ctest that runs `cmake/check_regression.py`,
which compares every `validation/*.dat` file against its reference
counterpart to a relative tolerance of 10⁻¹⁰.

## Benchmarking

The batch (`_ALL`) routines exist to avoid redundantly recomputing shared
recurrence terms when evaluating every mode at a point — `vsh_benchmark`
measures that payoff directly, timing each batch routine against a naive
loop over the corresponding single-mode routine, swept across increasing
`Lmax`. It is an opt-in build target, not a `ctest`, since wall-clock
timings are machine- and compiler-dependent and shouldn't gate a pass/fail
build:

```bash
cmake -B build -DVSH_BUILD_BENCHMARK=ON
cmake --build build
./build/vsh_benchmark
```

Run it from the repository root (as shown) so it can create `benchmark/`
alongside `validation/`. Each measurement auto-calibrates its repeat count
against wall-clock time, so results stay resolvable from `Lmax=5` up
through `Lmax=320` without a hardcoded repetition count.

| File | Comparison |
|---|---|
| `benchmark/bench_legendre.dat`     | `ASSOC_LEGENDRE_NORM_ALL` vs. looped `ASSOC_LEGENDRE` |
| `benchmark/bench_ssh.dat`          | `SSH_ALL` vs. looped `SSH` |
| `benchmark/bench_vsh_tor.dat`      | `VSH_TOR_ALL` vs. looped `VSH_TOR` |
| `benchmark/bench_vsh_pol_up.dat`   | `VSH_POL_UP_ALL` vs. looped `VSH_POL_UP` |
| `benchmark/bench_vsh_pol_dn.dat`   | `VSH_POL_DN_ALL` vs. looped `VSH_POL_DN` |

Each file has columns `LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP`,
one row per `Lmax` in `{5, 10, 20, 40, 80, 160, 320}`, timed per grid-point
evaluation over a 50-point sweep in the routine's natural argument
(colatitude, or `x=cos(theta)` for the Legendre case).
