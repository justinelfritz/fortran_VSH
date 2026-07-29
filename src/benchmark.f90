!     Timing benchmarks comparing the batch (_ALL) Legendre/SSH/VSH
!     routines against naive per-mode loops over the corresponding
!     single-mode routines, evaluated across a grid of NX colatitude/
!     argument points -- mirroring how a spectral MHD transform would
!     call these routines once per grid point every timestep. Each
!     measurement is auto-calibrated: the full grid sweep is repeated
!     until elapsed wall time crosses MIN_TIME, so timings stay
!     resolvable across the full Lmax sweep used for the manuscript's
!     performance figure, without a hardcoded repeat count.

MODULE BENCHMARK
USE KINDS,   ONLY: dp, i4, i8
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: &
  ASSOC_LEGENDRE, ASSOC_LEGENDRE_NORM_ALL, &
  SSH, SSH_ALL, &
  VSH_TOR, VSH_TOR_ALL, &
  VSH_POL_UP, VSH_POL_UP_ALL, &
  VSH_POL_DN, VSH_POL_DN_ALL
IMPLICIT NONE
PRIVATE
PUBLIC :: &
  BENCH_LEGENDRE, BENCH_SSH, &
  BENCH_VSH_TOR, BENCH_VSH_POL_UP, BENCH_VSH_POL_DN

REAL(KIND=dp),    PARAMETER :: MIN_TIME = 0.2_dp
INTEGER(KIND=i4), PARAMETER :: NREP_CAP = 1000000

CONTAINS

!     Wall-clock seconds elapsed, via SYSTEM_CLOCK with an INTEGER(i8)
!     count so the calibration loops below cannot roll the counter over.
  FUNCTION WALL_TIME() RESULT(T)
  IMPLICIT NONE
  REAL(KIND=dp) :: T
  INTEGER(KIND=i8) :: TCOUNT, TRATE
  CALL SYSTEM_CLOCK(TCOUNT, TRATE)
  T = REAL(TCOUNT, dp) / REAL(TRATE, dp)
  RETURN
  END FUNCTION WALL_TIME

!     Benchmark ASSOC_LEGENDRE_NORM_ALL (batch) vs. looped ASSOC_LEGENDRE
!     (single-mode) for 0<=l<=lmax, 0<=m<=l, swept over NX points
!     x in [-0.9,0.9]. Appends one row to OUTUNIT:
!     LMAX  N_MODES  T_BATCH  T_LOOP  SPEEDUP
  SUBROUTINE BENCH_LEGENDRE(LMAX, NX, OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4) :: L, M, IX, R, NREP, PSIZE
  REAL(KIND=dp) :: X, DX, T0, T1, T_BATCH, T_LOOP, ACC
  REAL(KIND=dp), ALLOCATABLE :: PNORM(:)

  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(PNORM(PSIZE))
  DX = 1.8_dp / (NX - 1)

  ! ── Batch: one ASSOC_LEGENDRE_NORM_ALL call per grid point ──────────────
  NREP = 1
  DO
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        X = -0.9_dp + IX*DX
        CALL ASSOC_LEGENDRE_NORM_ALL(PNORM, LMAX, X)
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_BATCH = (T1-T0) / (NREP*NX)
  ACC = SUM(PNORM)
  IF (ACC == -1.0d300) WRITE(*,'(A)') ""  !- defeat dead-code elimination

  ! ── Loop: ASSOC_LEGENDRE called once per (l,m) at each grid point ───────
  NREP = 1
  DO
    ACC = 0.d0
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        X = -0.9_dp + IX*DX
        DO L = 0, LMAX
          DO M = 0, L
            ACC = ACC + ASSOC_LEGENDRE(L, M, X)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_LOOP = (T1-T0) / (NREP*NX)
  IF (ACC == -1.0d300) WRITE(*,'(A)') ""

  WRITE(OUTUNIT,'(I6,1X,I8,1X,ES14.6,1X,ES14.6,1X,F12.3)') &
    LMAX, PSIZE, T_BATCH, T_LOOP, T_LOOP/T_BATCH

  DEALLOCATE(PNORM)
  RETURN
  END SUBROUTINE BENCH_LEGENDRE

!     Benchmark SSH_ALL (batch) vs. looped SSH (single-mode) for
!     0<=l<=lmax, -l<=m<=l, swept over NX colatitude points in
!     (0,pi) at a fixed longitude. Appends one row to OUTUNIT:
!     LMAX  N_MODES  T_BATCH  T_LOOP  SPEEDUP
  SUBROUTINE BENCH_SSH(LMAX, NX, OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4) :: L, M, IX, R, NREP, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH, T0, T1, T_BATCH, T_LOOP
  COMPLEX(KIND=dp) :: ACC
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:)

  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM))
  PHI = 0.123d0   !- arbitrary fixed longitude
  DTH = pi / (NX - 1)

  NREP = 1
  DO
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        CALL SSH_ALL(YLM, LMAX, THETA, PHI)
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_BATCH = (T1-T0) / (NREP*NX)
  ACC = SUM(YLM)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  NREP = 1
  DO
    ACC = DCMPLX(0.d0,0.d0)
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        DO L = 0, LMAX
          DO M = -L, L
            ACC = ACC + SSH(L, M, THETA, PHI)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_LOOP = (T1-T0) / (NREP*NX)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  WRITE(OUTUNIT,'(I6,1X,I8,1X,ES14.6,1X,ES14.6,1X,F12.3)') &
    LMAX, NYLM, T_BATCH, T_LOOP, T_LOOP/T_BATCH

  DEALLOCATE(YLM)
  RETURN
  END SUBROUTINE BENCH_SSH

!     Benchmark VSH_TOR_ALL (batch) vs. looped VSH_TOR (single-mode)
!     for 0<=l<=lmax, -l<=m<=l, swept over NX colatitude points.
!     Appends one row to OUTUNIT: LMAX  N_MODES  T_BATCH  T_LOOP  SPEEDUP
  SUBROUTINE BENCH_VSH_TOR(LMAX, NX, OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4) :: L, M, IX, R, NREP, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH, T0, T1, T_BATCH, T_LOOP
  COMPLEX(KIND=dp) :: ACC
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)

  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3,NYLM))
  PHI = 0.123d0
  DTH = pi / (NX - 1)

  NREP = 1
  DO
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        CALL VSH_TOR_ALL(OUT, LMAX, THETA, PHI)
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_BATCH = (T1-T0) / (NREP*NX)
  ACC = SUM(OUT)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  NREP = 1
  DO
    ACC = DCMPLX(0.d0,0.d0)
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        DO L = 0, LMAX
          DO M = -L, L
            ACC = ACC + SUM(VSH_TOR(L, M, THETA, PHI))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_LOOP = (T1-T0) / (NREP*NX)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  WRITE(OUTUNIT,'(I6,1X,I8,1X,ES14.6,1X,ES14.6,1X,F12.3)') &
    LMAX, NYLM, T_BATCH, T_LOOP, T_LOOP/T_BATCH

  DEALLOCATE(OUT)
  RETURN
  END SUBROUTINE BENCH_VSH_TOR

!     Benchmark VSH_POL_UP_ALL (batch) vs. looped VSH_POL_UP
!     (single-mode) for 0<=l<=lmax, -l<=m<=l, swept over NX colatitude
!     points. Appends one row to OUTUNIT:
!     LMAX  N_MODES  T_BATCH  T_LOOP  SPEEDUP
  SUBROUTINE BENCH_VSH_POL_UP(LMAX, NX, OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4) :: L, M, IX, R, NREP, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH, T0, T1, T_BATCH, T_LOOP
  COMPLEX(KIND=dp) :: ACC
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)

  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3,NYLM))
  PHI = 0.123d0
  DTH = pi / (NX - 1)

  NREP = 1
  DO
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        CALL VSH_POL_UP_ALL(OUT, LMAX, THETA, PHI)
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_BATCH = (T1-T0) / (NREP*NX)
  ACC = SUM(OUT)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  NREP = 1
  DO
    ACC = DCMPLX(0.d0,0.d0)
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        DO L = 0, LMAX
          DO M = -L, L
            ACC = ACC + SUM(VSH_POL_UP(L, M, THETA, PHI))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_LOOP = (T1-T0) / (NREP*NX)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  WRITE(OUTUNIT,'(I6,1X,I8,1X,ES14.6,1X,ES14.6,1X,F12.3)') &
    LMAX, NYLM, T_BATCH, T_LOOP, T_LOOP/T_BATCH

  DEALLOCATE(OUT)
  RETURN
  END SUBROUTINE BENCH_VSH_POL_UP

!     Benchmark VSH_POL_DN_ALL (batch) vs. looped VSH_POL_DN
!     (single-mode) for 0<=l<=lmax, -l<=m<=l, swept over NX colatitude
!     points. Appends one row to OUTUNIT:
!     LMAX  N_MODES  T_BATCH  T_LOOP  SPEEDUP
  SUBROUTINE BENCH_VSH_POL_DN(LMAX, NX, OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4) :: L, M, IX, R, NREP, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH, T0, T1, T_BATCH, T_LOOP
  COMPLEX(KIND=dp) :: ACC
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)

  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3,NYLM))
  PHI = 0.123d0
  DTH = pi / (NX - 1)

  NREP = 1
  DO
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        CALL VSH_POL_DN_ALL(OUT, LMAX, THETA, PHI)
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_BATCH = (T1-T0) / (NREP*NX)
  ACC = SUM(OUT)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  NREP = 1
  DO
    ACC = DCMPLX(0.d0,0.d0)
    T0 = WALL_TIME()
    DO R = 1, NREP
      DO IX = 0, NX-1
        THETA = IX*DTH
        DO L = 0, LMAX
          DO M = -L, L
            ACC = ACC + SUM(VSH_POL_DN(L, M, THETA, PHI))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    T1 = WALL_TIME()
    IF (T1-T0 >= MIN_TIME .OR. NREP >= NREP_CAP) EXIT
    NREP = NREP*2
  ENDDO
  T_LOOP = (T1-T0) / (NREP*NX)
  IF (ACC == (-1.0d300,0.d0)) WRITE(*,'(A)') ""

  WRITE(OUTUNIT,'(I6,1X,I8,1X,ES14.6,1X,ES14.6,1X,F12.3)') &
    LMAX, NYLM, T_BATCH, T_LOOP, T_LOOP/T_BATCH

  DEALLOCATE(OUT)
  RETURN
  END SUBROUTINE BENCH_VSH_POL_DN

END MODULE BENCHMARK
