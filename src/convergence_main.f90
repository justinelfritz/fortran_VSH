!> Driver for the convergence-vs-degree accuracy sweep (opt-in
!> `VSH_BUILD_CONVERGENCE` CMake target, not part of `ctest`): calls the
!> same [[TESTS]] batch-consistency routines the core `ctest` gate uses
!> (`BATCH_ALM_CONS`, `BATCH_SSH_CONS`, `BATCH_VSH_TOR_CONS`,
!> `BATCH_VSH_POL_UP_CONS`, `BATCH_VSH_POL_DN_CONS`), but at
!> \( \ell_{max}=300 \) instead of the small, fixed values `main.f90` uses
!> to gate correctness. Each routine already loops internally over every
!> degree \( 0\le\ell\le\ell_{max} \) in one call, so this produces dense
!> per-degree error data up to 300 without needing an outer sweep here.
!>
!> Deliberately does **not** check the `STATUS` each routine returns:
!> batch-vs-single error legitimately grows with \( \ell \) (that's the
!> point of the resulting figure), so failing a `1e-12`-style tolerance at
!> high \( \ell \) is expected, not a regression. Writes to
!> `./convergence/*.dat` -- a directory kept separate from
!> `./validation/*.dat` specifically so this sweep can never interact with
!> `VSH_REGRESSION_TEST`'s comparison against `validation/reference/`.
!> Consumed by `py/plot_convergence.py` for the manuscript's convergence
!> figure.
PROGRAM CONVERGENCE_MAIN
USE KINDS,  ONLY: i4, dp
USE TESTS,  ONLY: &
  BATCH_ALM_CONS, BATCH_SSH_CONS, &
  BATCH_VSH_TOR_CONS, BATCH_VSH_POL_UP_CONS, BATCH_VSH_POL_DN_CONS
USE GLOBALS, ONLY: pi
USE VSH,    ONLY: &
  PLM_INDEX, YLM_INDEX, &
  ASSOC_LEGENDRE, ASSOC_LEGENDRE_ALL, ASSOC_LEGENDRE_NORM_ALL, &
  SSH, SSH_ALL, &
  VSH_TOR, VSH_TOR_ALL, VSH_POL_UP, VSH_POL_UP_ALL, &
  VSH_POL_DN, VSH_POL_DN_ALL
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: LMAX = 300
INTEGER(KIND=i4), PARAMETER :: NX_ALM = 20   !- x-grid density (matches main.f90)
INTEGER(KIND=i4), PARAMETER :: NTH    = 15   !- theta-grid density (matches main.f90)
INTEGER(KIND=i4) :: STATUS

!- Sparse sample grid for the independent-reference accuracy check (see
!- py/mpmath_reference.py): matches bench_main.f90's LMAX_LIST so this
!- figure and benchmark_scaling share the same x-axis sample points. Writes
!- RAW single-mode and batch values here (not diffs -- unlike the dense
!- sweep above) since py/plot_convergence.py diffs each independently
!- against an arbitrary-precision mpmath reference rather than against
!- each other.
INTEGER(KIND=i4), PARAMETER :: NREFL = 16
INTEGER(KIND=i4), DIMENSION(NREFL), PARAMETER :: &
  REFL_LIST = (/1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 2000, &
                2500, 3000, 3500, 4000/)
INTEGER(KIND=i4), PARAMETER :: NX_REF = 3
REAL(KIND=dp), DIMENSION(NX_REF), PARAMETER :: &
  X_REF_LIST = (/-0.6_dp, 0.0_dp, 0.6_dp/)
INTEGER(KIND=i4), PARAMETER :: NTH_REF = 5
REAL(KIND=dp), DIMENSION(NTH_REF), PARAMETER :: &
  TH_REF_LIST = (/pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6/)
REAL(KIND=dp), PARAMETER :: PHI_REF = pi/4
INTEGER(KIND=i4) :: L, M, IL, IX, ITH, IM, NMSAMP, MSAMP(3)
REAL(KIND=dp) :: XX, TH, PSING
COMPLEX(KIND=dp) :: YSING, VSING(3)
REAL(KIND=dp), ALLOCATABLE :: PBATCH(:), PNORMBATCH(:)
COMPLEX(KIND=dp), ALLOCATABLE :: YBATCH(:), TORBATCH(:,:), UPBATCH(:,:), DNBATCH(:,:)

WRITE(*,'(A,I0,A)') "Sweeping batch-consistency error up to Lmax=", LMAX, "..."

OPEN(UNIT=11, FILE="./convergence/batch_alm_cons.dat")
CALL BATCH_ALM_CONS(LMAX, NX_ALM, 11, STATUS)
CLOSE(11)
WRITE(*,'(A)') "  Legendre done"

OPEN(UNIT=12, FILE="./convergence/batch_ssh_cons.dat")
CALL BATCH_SSH_CONS(LMAX, NTH, 12, STATUS)
CLOSE(12)
WRITE(*,'(A)') "  SSH done"

OPEN(UNIT=13, FILE="./convergence/batch_vsh_tor_cons.dat")
CALL BATCH_VSH_TOR_CONS(LMAX, NTH, 13, STATUS)
CLOSE(13)
WRITE(*,'(A)') "  VSH toroidal done"

OPEN(UNIT=14, FILE="./convergence/batch_vsh_pol_up_cons.dat")
CALL BATCH_VSH_POL_UP_CONS(LMAX, NTH, 14, STATUS)
CLOSE(14)
WRITE(*,'(A)') "  VSH poloidal (up) done"

OPEN(UNIT=15, FILE="./convergence/batch_vsh_pol_dn_cons.dat")
CALL BATCH_VSH_POL_DN_CONS(LMAX, NTH, 15, STATUS)
CLOSE(15)
WRITE(*,'(A)') "  VSH poloidal (dn) done"

WRITE(*,'(A)') "RESULT: convergence sweep complete - see ./convergence/*.dat"

WRITE(*,'(A)') "Sampling single-mode/batch values for independent-reference check..."

OPEN(UNIT=21, FILE="./convergence/reference_alm.dat")
WRITE(21,'(A)') &
    '# L  M  X  P_single  P_batch  P_batch_norm'// &
    '  -- P_single/P_batch are unnormalized (ASSOC_LEGENDRE/_ALL);'// &
    ' P_batch_norm is ASSOC_LEGENDRE_NORM_ALL (no single-mode normalized'// &
    ' counterpart exists) -- see py/mpmath_reference.py legendre_norm_ref'
OPEN(UNIT=22, FILE="./convergence/reference_ssh.dat")
WRITE(22,'(A)') '# L  M  THETA  PHI  single_re  single_im  batch_re  batch_im'
OPEN(UNIT=23, FILE="./convergence/reference_vsh_tor.dat")
WRITE(23,'(A)') '# L  M  THETA  PHI'// &
    '  single_r_re single_r_im single_th_re single_th_im single_ph_re single_ph_im'// &
    '  batch_r_re batch_r_im batch_th_re batch_th_im batch_ph_re batch_ph_im'
OPEN(UNIT=24, FILE="./convergence/reference_vsh_pol_up.dat")
WRITE(24,'(A)') '# L  M  THETA  PHI'// &
    '  single_r_re single_r_im single_th_re single_th_im single_ph_re single_ph_im'// &
    '  batch_r_re batch_r_im batch_th_re batch_th_im batch_ph_re batch_ph_im'
OPEN(UNIT=25, FILE="./convergence/reference_vsh_pol_dn.dat")
WRITE(25,'(A)') '# L  M  THETA  PHI'// &
    '  single_r_re single_r_im single_th_re single_th_im single_ph_re single_ph_im'// &
    '  batch_r_re batch_r_im batch_th_re batch_th_im batch_ph_re batch_ph_im'

DO IL = 1, NREFL
  L = REFL_LIST(IL)

  !- dedup M samples {0, L/2, L}
  MSAMP(1) = 0
  MSAMP(2) = L/2
  MSAMP(3) = L
  NMSAMP = 1
  DO IM = 2, 3
    IF (MSAMP(IM) /= MSAMP(NMSAMP)) THEN
      NMSAMP = NMSAMP + 1
      MSAMP(NMSAMP) = MSAMP(IM)
    END IF
  END DO

  !- Legendre
  DO IX = 1, NX_REF
    XX = X_REF_LIST(IX)
    ALLOCATE(PBATCH((L+1)*(L+2)/2))
    ALLOCATE(PNORMBATCH((L+1)*(L+2)/2))
    CALL ASSOC_LEGENDRE_ALL(PBATCH, L, XX)
    CALL ASSOC_LEGENDRE_NORM_ALL(PNORMBATCH, L, XX)
    DO IM = 1, NMSAMP
      M = MSAMP(IM)
      PSING = ASSOC_LEGENDRE(L, M, XX)
      WRITE(21, *) L, M, XX, PSING, PBATCH(PLM_INDEX(L, M)), &
          PNORMBATCH(PLM_INDEX(L, M))
    END DO
    DEALLOCATE(PBATCH)
    DEALLOCATE(PNORMBATCH)
  END DO

  !- SSH, VSH toroidal, VSH poloidal (up/dn) share one THETA loop -- their
  !- batch calls are all keyed off the same (L, THETA, PHI_REF)
  DO ITH = 1, NTH_REF
    TH = TH_REF_LIST(ITH)
    ALLOCATE(YBATCH((L+1)**2))
    ALLOCATE(TORBATCH(3, (L+1)**2))
    ALLOCATE(UPBATCH(3, (L+1)**2))
    ALLOCATE(DNBATCH(3, (L+1)**2))
    CALL SSH_ALL(YBATCH, L, TH, PHI_REF)
    CALL VSH_TOR_ALL(TORBATCH, L, TH, PHI_REF)
    CALL VSH_POL_UP_ALL(UPBATCH, L, TH, PHI_REF)
    CALL VSH_POL_DN_ALL(DNBATCH, L, TH, PHI_REF)

    DO IM = 1, NMSAMP
      M = MSAMP(IM)

      YSING = SSH(L, M, TH, PHI_REF)
      WRITE(22, *) L, M, TH, PHI_REF, &
          DREAL(YSING), DIMAG(YSING), &
          DREAL(YBATCH(YLM_INDEX(L,M))), DIMAG(YBATCH(YLM_INDEX(L,M)))

      VSING = VSH_TOR(L, M, TH, PHI_REF)
      WRITE(23, *) L, M, TH, PHI_REF, &
          DREAL(VSING(1)), DIMAG(VSING(1)), &
          DREAL(VSING(2)), DIMAG(VSING(2)), &
          DREAL(VSING(3)), DIMAG(VSING(3)), &
          DREAL(TORBATCH(1,YLM_INDEX(L,M))), DIMAG(TORBATCH(1,YLM_INDEX(L,M))), &
          DREAL(TORBATCH(2,YLM_INDEX(L,M))), DIMAG(TORBATCH(2,YLM_INDEX(L,M))), &
          DREAL(TORBATCH(3,YLM_INDEX(L,M))), DIMAG(TORBATCH(3,YLM_INDEX(L,M)))

      VSING = VSH_POL_UP(L, M, TH, PHI_REF)
      WRITE(24, *) L, M, TH, PHI_REF, &
          DREAL(VSING(1)), DIMAG(VSING(1)), &
          DREAL(VSING(2)), DIMAG(VSING(2)), &
          DREAL(VSING(3)), DIMAG(VSING(3)), &
          DREAL(UPBATCH(1,YLM_INDEX(L,M))), DIMAG(UPBATCH(1,YLM_INDEX(L,M))), &
          DREAL(UPBATCH(2,YLM_INDEX(L,M))), DIMAG(UPBATCH(2,YLM_INDEX(L,M))), &
          DREAL(UPBATCH(3,YLM_INDEX(L,M))), DIMAG(UPBATCH(3,YLM_INDEX(L,M)))

      VSING = VSH_POL_DN(L, M, TH, PHI_REF)
      WRITE(25, *) L, M, TH, PHI_REF, &
          DREAL(VSING(1)), DIMAG(VSING(1)), &
          DREAL(VSING(2)), DIMAG(VSING(2)), &
          DREAL(VSING(3)), DIMAG(VSING(3)), &
          DREAL(DNBATCH(1,YLM_INDEX(L,M))), DIMAG(DNBATCH(1,YLM_INDEX(L,M))), &
          DREAL(DNBATCH(2,YLM_INDEX(L,M))), DIMAG(DNBATCH(2,YLM_INDEX(L,M))), &
          DREAL(DNBATCH(3,YLM_INDEX(L,M))), DIMAG(DNBATCH(3,YLM_INDEX(L,M)))
    END DO
    DEALLOCATE(YBATCH, TORBATCH, UPBATCH, DNBATCH)
  END DO

  WRITE(*,'(A,I0,A)') "  L=", L, " done"
END DO

CLOSE(21)
CLOSE(22)
CLOSE(23)
CLOSE(24)
CLOSE(25)

WRITE(*,'(A)') &
    "RESULT: independent-reference sample data complete - see ./convergence/reference_*.dat"

END PROGRAM CONVERGENCE_MAIN
