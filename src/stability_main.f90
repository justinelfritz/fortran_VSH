!> Driver for mapping ASSOC_LEGENDRE_NORM_ALL's numerical stability across
!> the \( (\ell, m, \theta) \) parameter space (opt-in `VSH_BUILD_STABILITY`
!> CMake target, not part of `ctest`). See `STABILITY_FINDINGS.md` at the
!> repo root for the full writeup -- this exists because the documented
!> "stable to \( \ell\sim2700 \)" claim for [[ASSOC_LEGENDRE_NORM_ALL]]
!> turned out to depend on \( \theta \) as well as \( \ell \): for \( \theta
!> \) close to the poles and \( m \) in a band tied to the classical
!> turning point \( \sin\theta \approx m/\ell \) of the associated Legendre
!> ODE, the forward-column recurrence can blow up by 100+ orders of
!> magnitude even at moderate \( \ell \).
!>
!> For each \( \theta \) sample, one call to `ASSOC_LEGENDRE_NORM_ALL` at
!> `LMAX_SCAN` yields the entire \( 0\le\ell\le\text{LMAX\_SCAN},\,0\le
!> m\le\ell \) table in one pass (the recurrence is forward in \( \ell \)
!> for fixed \( m \), so there is no need to re-call per \( \ell \)
!> checkpoint). `BREAK_THRESHOLD` is set far above any legitimate
!> normalized value (which should be at most \( O(\sqrt{\ell}) \), ~70 at
!> \( \ell=5000 \)) and far below the observed blowups (\( >10^{100} \)),
!> so it cleanly separates "legitimate" from "broken" without false
!> positives.
PROGRAM STABILITY_MAIN
USE KINDS,   ONLY: i4, dp
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: ASSOC_LEGENDRE_NORM_ALL, PLM_INDEX
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: LMAX_SCAN = 5000
REAL(KIND=dp),    PARAMETER :: BREAK_THRESHOLD = 1000.0_dp
INTEGER(KIND=i4), PARAMETER :: NTHETA = 179
INTEGER(KIND=i4), PARAMETER :: LCHECK_STEP = 50

REAL(KIND=dp), ALLOCATABLE :: PNORM(:)
REAL(KIND=dp) :: THETA_DEG, THETA, X, VAL, PEAK
INTEGER(KIND=i4) :: ITH, L, M, MMIN, MMAX
LOGICAL :: BROKEN

ALLOCATE(PNORM((LMAX_SCAN+1)*(LMAX_SCAN+2)/2))

OPEN(UNIT=31, FILE="./stability/stability_map.dat")
WRITE(31,'(A)') &
    '# theta_deg  l  broken  m_min_broken  m_max_broken'// &
    '  m_min_over_l  m_max_over_l  peak_abs_value'

DO ITH = 1, NTHETA
  THETA_DEG = REAL(ITH, dp)
  THETA = THETA_DEG * pi / 180.0_dp
  X = DCOS(THETA)
  CALL ASSOC_LEGENDRE_NORM_ALL(PNORM, LMAX_SCAN, X)

  DO L = LCHECK_STEP, LMAX_SCAN, LCHECK_STEP
    BROKEN = .FALSE.
    MMIN = -1
    MMAX = -1
    PEAK = 0.0_dp
    DO M = 0, L
      VAL = ABS(PNORM(PLM_INDEX(L, M)))
      IF (VAL > PEAK) PEAK = VAL
      IF (VAL > BREAK_THRESHOLD) THEN
        BROKEN = .TRUE.
        IF (MMIN < 0) MMIN = M
        MMAX = M
      END IF
    END DO
    IF (BROKEN) THEN
      WRITE(31,*) THETA_DEG, L, 1, MMIN, MMAX, &
          REAL(MMIN,dp)/REAL(L,dp), REAL(MMAX,dp)/REAL(L,dp), PEAK
    ELSE
      WRITE(31,*) THETA_DEG, L, 0, -1, -1, -1.0_dp, -1.0_dp, PEAK
    END IF
  END DO
  IF (MOD(ITH, 20) == 0) WRITE(*,'(A,F6.1,A)') "  theta=", THETA_DEG, " done"
END DO

CLOSE(31)
DEALLOCATE(PNORM)
WRITE(*,'(A)') "RESULT: stability map complete - see ./stability/stability_map.dat"

END PROGRAM STABILITY_MAIN
