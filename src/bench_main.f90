!> Driver for the [[BENCHMARK]] suite (opt-in `VSH_BUILD_BENCHMARK` CMake
!> target, not part of `ctest`): runs every `BENCH_*` routine across a
!> fixed sweep of \( \ell_{max} \) values, writing one row per
!> \( \ell_{max} \) to `./benchmark/bench_*.dat` (consumed by
!> `py/plot_benchmark.py` for the manuscript's performance figure).
PROGRAM BENCH_MAIN
USE KINDS,     ONLY: i4
USE BENCHMARK, ONLY: &
  BENCH_LEGENDRE, BENCH_SSH, &
  BENCH_VSH_TOR, BENCH_VSH_POL_UP, BENCH_VSH_POL_DN
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: NX = 50   !- grid points swept per Lmax
!- Capped at 2000 (BATCH_SAFE_LMAX in plot_benchmark.py, matching
!- STABILITY_FINDINGS.md Part 1) -- no need to spend time computing batch
!- timings past the range the manuscript figure actually plots. The naive
!- loop is separately capped inside BENCHMARK (LOOP_LMAX_CAP=200; see
!- src/benchmark.f90) since its per-measurement cost is O(NX*Lmax^3), not
!- O(NX*Lmax^2) like batch's, and became a multi-hour runaway at Lmax=4000
!- before that cap was added.
INTEGER(KIND=i4), PARAMETER :: NLMAX = 12
INTEGER(KIND=i4), DIMENSION(NLMAX), PARAMETER :: &
  LMAX_LIST = (/1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 2000/)
INTEGER(KIND=i4) :: I

OPEN(UNIT=51, FILE="./benchmark/bench_legendre.dat")
OPEN(UNIT=52, FILE="./benchmark/bench_ssh.dat")
OPEN(UNIT=53, FILE="./benchmark/bench_vsh_tor.dat")
OPEN(UNIT=54, FILE="./benchmark/bench_vsh_pol_up.dat")
OPEN(UNIT=55, FILE="./benchmark/bench_vsh_pol_dn.dat")

WRITE(51,'(A)') '# LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP'
WRITE(52,'(A)') '# LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP'
WRITE(53,'(A)') '# LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP'
WRITE(54,'(A)') '# LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP'
WRITE(55,'(A)') '# LMAX  N_MODES  T_BATCH_sec  T_LOOP_sec  SPEEDUP'

DO I = 1, NLMAX
  WRITE(*,'(A,I0)') "Benchmarking Lmax=", LMAX_LIST(I)
  CALL BENCH_LEGENDRE  (LMAX_LIST(I), NX, 51)
  CALL BENCH_SSH       (LMAX_LIST(I), NX, 52)
  CALL BENCH_VSH_TOR   (LMAX_LIST(I), NX, 53)
  CALL BENCH_VSH_POL_UP(LMAX_LIST(I), NX, 54)
  CALL BENCH_VSH_POL_DN(LMAX_LIST(I), NX, 55)
ENDDO

CLOSE(51)
CLOSE(52)
CLOSE(53)
CLOSE(54)
CLOSE(55)

WRITE(*,'(A)') "RESULT: benchmark complete - see ./benchmark/*.dat"

END PROGRAM BENCH_MAIN
