!> Driver for mapping the UNNORMALIZED Bonnet recurrence's
!> (`ASSOC_LEGENDRE`/`ASSOC_LEGENDRE_ALL`, and everything single-mode built
!> on it) usable accuracy range across \( (\ell, m, \theta) \) (opt-in
!> `VSH_BUILD_STABILITY` CMake target, same as [[STABILITY_MAIN]] -- not
!> part of `ctest`). See `STABILITY_FINDINGS.md` at the repo root.
!>
!> Unlike the normalized-recurrence scan, this one doesn't hunt for a
!> literal overflow -- the seed \( P_m^m(x) = (-1)^m(2m-1)!!\sin^m\theta \)
!> grows/shrinks combinatorially with \( m \), so raw magnitude alone is a
!> poor proxy for accuracy: catastrophic cancellation in the later
!> three-term upward recurrence in \( \ell \) can destroy relative
!> precision well before the value is anywhere near double-precision
!> overflow. So this driver writes RAW single-mode and batch values (not a
!> broken/not-broken flag) for `py/plot_stability_unnorm.py` to diff
!> against `py/mpmath_reference.py`'s `legendre_ref` -- the same
!> arbitrary-precision rerun of this exact recurrence used to validate
!> ASSOC_LEGENDRE_NORM_ALL's stable range.
!>
!> For each `theta`, one call to `ASSOC_LEGENDRE_ALL` at `LMAX_SCAN`
!> yields the whole `l<=LMAX_SCAN, m<=l` table; three `m` fractions of
!> `l` (1.0 = diagonal/worst-case, 0.75, 0.5) are sampled per `l` to
!> confirm the diagonal really is the worst case rather than assuming it.
PROGRAM UNNORM_STABILITY_MAIN
USE KINDS,   ONLY: i4, dp
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: ASSOC_LEGENDRE, ASSOC_LEGENDRE_ALL, PLM_INDEX
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: LMAX_SCAN = 200
INTEGER(KIND=i4), PARAMETER :: NTHETA = 179
INTEGER(KIND=i4), PARAMETER :: NFRAC = 3
REAL(KIND=dp), DIMENSION(NFRAC), PARAMETER :: MFRAC = (/1.0_dp, 0.75_dp, 0.5_dp/)

REAL(KIND=dp), ALLOCATABLE :: PBATCH(:)
REAL(KIND=dp) :: THETA_DEG, THETA, X, PSING
INTEGER(KIND=i4) :: ITH, L, M, IF_, MPREV

ALLOCATE(PBATCH((LMAX_SCAN+1)*(LMAX_SCAN+2)/2))

OPEN(UNIT=41, FILE="./stability/unnorm_diagonal.dat")
WRITE(41,'(A)') '# theta_deg  l  m  m_frac  P_single  P_batch'

DO ITH = 1, NTHETA
  THETA_DEG = REAL(ITH, dp)
  THETA = THETA_DEG * pi / 180.0_dp
  X = DCOS(THETA)
  CALL ASSOC_LEGENDRE_ALL(PBATCH, LMAX_SCAN, X)

  DO L = 1, LMAX_SCAN
    MPREV = -1
    DO IF_ = 1, NFRAC
      M = NINT(MFRAC(IF_) * L)
      IF (M == MPREV) CYCLE
      MPREV = M
      PSING = ASSOC_LEGENDRE(L, M, X)
      WRITE(41,*) THETA_DEG, L, M, MFRAC(IF_), PSING, PBATCH(PLM_INDEX(L,M))
    END DO
  END DO
  IF (MOD(ITH, 20) == 0) WRITE(*,'(A,F6.1,A)') "  theta=", THETA_DEG, " done"
END DO

CLOSE(41)
DEALLOCATE(PBATCH)
WRITE(*,'(A)') "RESULT: unnormalized scan complete - see ./stability/unnorm_diagonal.dat"

END PROGRAM UNNORM_STABILITY_MAIN
