!     Worked example (2c): idealized superconducting-core flux expulsion.
!     Ties to Elfritz (2016), cited in the manuscript, on field expulsion
!     from superconducting neutron-star cores. Simplified to the standard
!     "perfectly conducting sphere in a uniform ambient field" textbook
!     problem: a sphere of radius R totally excludes flux (idealized
!     Meissner/Type-I limit -- real neutron-star cores are believed
!     Type-II, with flux confined to quantized fluxoids rather than fully
!     excluded; that fuller picture is not computed here). The boundary
!     condition Br(R)=0 for all theta (no field may cross into a perfect
!     diamagnet) fixes the induced dipole moment in closed form:
!         MU_IND = -B0_AMB * R^3 / 2
!     which superposes on the uniform ambient field to give the exterior
!     field. Interior field is exactly zero. Purely magnetostatic: no
!     induction equation, no time evolution.

PROGRAM DIPOLE_FLUX_EXPULSION
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: PVSH_RAD, PVSH_POL
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: NR = 61
INTEGER(KIND=i4), PARAMETER :: NTHETA_LIST = 3
REAL(KIND=dp),    PARAMETER :: R_SPHERE = 1.0_dp
REAL(KIND=dp),    PARAMETER :: B0_AMB = 1.0_dp     !- ambient uniform field
REAL(KIND=dp),    PARAMETER :: MU_IND = -B0_AMB*R_SPHERE**3/2.0_dp   !- induced dipole
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.2_dp, R_MAX = 3.0_dp
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-12_dp

REAL(KIND=dp) :: THETA_LIST(NTHETA_LIST)
REAL(KIND=dp) :: A_AMB, C_IND, R, THETA, PHI, DR, S, DSDR
REAL(KIND=dp) :: BR_LIB, BTH_LIB, BR_ANALYTIC, BTH_ANALYTIC, MAX_ERR
REAL(KIND=dp) :: BR_OUT_AT_R, BTH_OUT_AT_R, JUMP_EXPECTED
COMPLEX(KIND=dp), DIMENSION(3) :: RADVEC, POLVEC
INTEGER(KIND=i4) :: IR, ITH, N_FAIL

THETA_LIST = (/ 0.0_dp, 0.25_dp*pi, 0.5_dp*pi /)
PHI = 0.0_dp
N_FAIL = 0

A_AMB = B0_AMB / (2.0_dp*DSQRT(3.0_dp/(4.0_dp*pi)))   !- ambient term: S = A_AMB * r^2
C_IND = MU_IND / DSQRT(3.0_dp/(4.0_dp*pi))             !- induced term: S = C_IND / r

! ── Field profile across the boundary, library vs. closed form ──────────────
OPEN(UNIT=65, FILE="./examples/dipole_flux_expulsion/field_profile.dat")
WRITE(65,'(A)') '# r  theta  Br_lib  Btheta_lib  Br_analytic  Btheta_analytic'

DR = (R_MAX - R_MIN) / (NR - 1)
MAX_ERR = 0.0_dp
DO ITH = 1, NTHETA_LIST
  THETA = THETA_LIST(ITH)
  DO IR = 0, NR-1
    R = R_MIN + IR*DR
    IF (R < R_SPHERE) THEN
      ! Perfect flux exclusion: field is exactly zero inside.
      BR_LIB = 0.0_dp
      BTH_LIB = 0.0_dp
      BR_ANALYTIC = 0.0_dp
      BTH_ANALYTIC = 0.0_dp
    ELSE
      S    = A_AMB*R**2 + C_IND/R
      DSDR = 2.0_dp*A_AMB*R - C_IND/R**2
      BR_ANALYTIC  = B0_AMB*DCOS(THETA) + 2.0_dp*MU_IND*DCOS(THETA)/R**3
      BTH_ANALYTIC = -B0_AMB*DSIN(THETA) +        MU_IND*DSIN(THETA)/R**3

      RADVEC = PVSH_RAD(1, 0, THETA, PHI)
      POLVEC = PVSH_POL(1, 0, THETA, PHI)
      BR_LIB  = (2.0_dp*S/R**2) * DBLE(RADVEC(1))
      BTH_LIB = (DSQRT(2.0_dp)/R * DSDR) * DBLE(POLVEC(2))
    END IF

    MAX_ERR = MAX(MAX_ERR, ABS(BR_LIB-BR_ANALYTIC), ABS(BTH_LIB-BTH_ANALYTIC))
    WRITE(65,*) R, THETA, BR_LIB, BTH_LIB, BR_ANALYTIC, BTH_ANALYTIC
  ENDDO
ENDDO
CLOSE(65)

IF (MAX_ERR > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  flux_expulsion_accuracy  max_err=", MAX_ERR
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  flux_expulsion_accuracy  max_err=", MAX_ERR
END IF

! ── Boundary check at r=R: Br(R)=0 (defining condition), Btheta jump ────────
OPEN(UNIT=66, FILE="./examples/dipole_flux_expulsion/boundary_check.dat")
WRITE(66,'(A)') '# theta  Br_out(R)  Btheta_out(R)  jump_Btheta  jump_expected'
DO ITH = 1, NTHETA_LIST
  THETA = THETA_LIST(ITH)

  S    = A_AMB*R_SPHERE**2 + C_IND/R_SPHERE
  DSDR = 2.0_dp*A_AMB*R_SPHERE - C_IND/R_SPHERE**2
  RADVEC = PVSH_RAD(1, 0, THETA, PHI)
  POLVEC = PVSH_POL(1, 0, THETA, PHI)
  BR_OUT_AT_R  = (2.0_dp*S/R_SPHERE**2) * DBLE(RADVEC(1))
  BTH_OUT_AT_R = (DSQRT(2.0_dp)/R_SPHERE * DSDR) * DBLE(POLVEC(2))

  JUMP_EXPECTED = -1.5_dp*B0_AMB*DSIN(THETA)   !- Btheta_out(R) - 0 (interior)

  WRITE(66,*) THETA, BR_OUT_AT_R, BTH_OUT_AT_R, BTH_OUT_AT_R, JUMP_EXPECTED

  IF (ABS(BR_OUT_AT_R) > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  Br_zero_at_surface  Br=", BR_OUT_AT_R
  END IF
  IF (ABS(BTH_OUT_AT_R-JUMP_EXPECTED) > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  Btheta_jump  absdiff=", ABS(BTH_OUT_AT_R-JUMP_EXPECTED)
  END IF
ENDDO
CLOSE(66)
IF (N_FAIL == 0) WRITE(*,'(A)') "PASS  Br_zero_and_Btheta_jump (all theta)"

! ── Summary ──────────────────────────────────────────────────────────────────
IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - dipole_flux_expulsion"
END IF

END PROGRAM DIPOLE_FLUX_EXPULSION
