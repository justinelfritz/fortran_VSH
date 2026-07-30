!     Worked example (2b): uniformly magnetized sphere <-> exterior dipole.
!     Classic textbook magnetostatic boundary-value problem (e.g. Griffiths
!     Ex. 6.1 / Jackson): a sphere of radius R with uniform magnetization
!     has a genuinely uniform interior field and a pure exterior dipole
!     field. In the normalized units used here (mu_0/4*pi = 1, i.e.
!     Br_exterior = 2*MU*cos(theta)/r^3), the magnetic-scalar-potential
!     solution gives the closed-form relation
!         B0 (interior uniform field) = 2*MU/R^3
!     Because there is no free surface current in this relation's
!     derivation but the interior/exterior fields are NOT continuous in
!     Btheta (physically: continuity is broken by the sphere's own bound
!     magnetization surface current, K = M x n-hat) -- Br IS continuous at
!     r=R, but Btheta has a jump. Both are checked below using FORTVSH's
!     PVSH_RAD/PVSH_POL evaluated with the interior (S_in=A*r^2, uniform
!     field) and exterior (S_out=C/r, dipole) stream functions.
!     Purely static/magnetostatic: no induction equation, no time
!     evolution.

PROGRAM DIPOLE_UNIFORM_SPHERE
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: PVSH_RAD, PVSH_POL
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: NR = 61
INTEGER(KIND=i4), PARAMETER :: NTHETA_LIST = 3
REAL(KIND=dp),    PARAMETER :: R_SPHERE = 1.0_dp
REAL(KIND=dp),    PARAMETER :: MU = 1.0_dp        !- exterior dipole moment
REAL(KIND=dp),    PARAMETER :: B0 = 2.0_dp*MU/R_SPHERE**3   !- interior uniform field
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.2_dp, R_MAX = 3.0_dp
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-12_dp

REAL(KIND=dp) :: THETA_LIST(NTHETA_LIST)
REAL(KIND=dp) :: A_STREAM, C_STREAM, R, THETA, PHI, DR, S, DSDR
REAL(KIND=dp) :: BR_LIB, BTH_LIB, BR_ANALYTIC, BTH_ANALYTIC, MAX_ERR
REAL(KIND=dp) :: BR_IN, BTH_IN, BR_OUT, BTH_OUT, JUMP_BTH, JUMP_EXPECTED
COMPLEX(KIND=dp), DIMENSION(3) :: RADVEC, POLVEC
INTEGER(KIND=i4) :: IR, ITH, N_FAIL

THETA_LIST = (/ 0.0_dp, 0.25_dp*pi, 0.5_dp*pi /)
PHI = 0.0_dp
N_FAIL = 0

A_STREAM = B0 / (2.0_dp*DSQRT(3.0_dp/(4.0_dp*pi)))   !- S_in(r) = A_STREAM * r^2
C_STREAM = MU / DSQRT(3.0_dp/(4.0_dp*pi))             !- S_out(r) = C_STREAM / r

! ── Field profile across the boundary, library vs. closed form ──────────────
OPEN(UNIT=63, FILE="./examples/dipole_uniform_sphere/field_profile.dat")
WRITE(63,'(A)') '# r  theta  Br_lib  Btheta_lib  Br_analytic  Btheta_analytic'

DR = (R_MAX - R_MIN) / (NR - 1)
MAX_ERR = 0.0_dp
DO ITH = 1, NTHETA_LIST
  THETA = THETA_LIST(ITH)
  DO IR = 0, NR-1
    R = R_MIN + IR*DR
    IF (R < R_SPHERE) THEN
      S = A_STREAM * R**2
      DSDR = 2.0_dp * A_STREAM * R
      BR_ANALYTIC  =  B0*DCOS(THETA)
      BTH_ANALYTIC = -B0*DSIN(THETA)
    ELSE
      S = C_STREAM / R
      DSDR = -C_STREAM / R**2
      BR_ANALYTIC  = 2.0_dp*MU*DCOS(THETA)/R**3
      BTH_ANALYTIC =         MU*DSIN(THETA)/R**3
    END IF

    RADVEC = PVSH_RAD(1, 0, THETA, PHI)
    POLVEC = PVSH_POL(1, 0, THETA, PHI)
    BR_LIB  = (2.0_dp*S/R**2) * DBLE(RADVEC(1))
    BTH_LIB = (DSQRT(2.0_dp)/R * DSDR) * DBLE(POLVEC(2))

    MAX_ERR = MAX(MAX_ERR, ABS(BR_LIB-BR_ANALYTIC), ABS(BTH_LIB-BTH_ANALYTIC))
    WRITE(63,*) R, THETA, BR_LIB, BTH_LIB, BR_ANALYTIC, BTH_ANALYTIC
  ENDDO
ENDDO
CLOSE(63)

IF (MAX_ERR > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  uniform_sphere_accuracy  max_err=", MAX_ERR
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  uniform_sphere_accuracy  max_err=", MAX_ERR
END IF

! ── Boundary check at r=R: Br continuous, Btheta jumps by 3*MU*sin(theta)/R^3
OPEN(UNIT=64, FILE="./examples/dipole_uniform_sphere/boundary_check.dat")
WRITE(64,'(A)') '# theta  Br_in  Br_out  absdiff_Br  Btheta_in  Btheta_out  jump_Btheta  jump_expected'
DO ITH = 1, NTHETA_LIST
  THETA = THETA_LIST(ITH)

  RADVEC = PVSH_RAD(1, 0, THETA, PHI)
  POLVEC = PVSH_POL(1, 0, THETA, PHI)

  S    = A_STREAM * R_SPHERE**2
  DSDR = 2.0_dp * A_STREAM * R_SPHERE
  BR_IN  = (2.0_dp*S/R_SPHERE**2) * DBLE(RADVEC(1))
  BTH_IN = (DSQRT(2.0_dp)/R_SPHERE * DSDR) * DBLE(POLVEC(2))

  S    = C_STREAM / R_SPHERE
  DSDR = -C_STREAM / R_SPHERE**2
  BR_OUT  = (2.0_dp*S/R_SPHERE**2) * DBLE(RADVEC(1))
  BTH_OUT = (DSQRT(2.0_dp)/R_SPHERE * DSDR) * DBLE(POLVEC(2))

  JUMP_BTH = BTH_OUT - BTH_IN
  JUMP_EXPECTED = 3.0_dp*MU*DSIN(THETA)/R_SPHERE**3

  WRITE(64,*) THETA, BR_IN, BR_OUT, ABS(BR_OUT-BR_IN), &
              BTH_IN, BTH_OUT, JUMP_BTH, JUMP_EXPECTED

  IF (ABS(BR_OUT-BR_IN) > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  Br_continuity  absdiff=", ABS(BR_OUT-BR_IN)
  END IF
  IF (ABS(JUMP_BTH-JUMP_EXPECTED) > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  Btheta_jump  absdiff=", ABS(JUMP_BTH-JUMP_EXPECTED)
  END IF
ENDDO
CLOSE(64)
IF (N_FAIL == 0) WRITE(*,'(A)') "PASS  Br_continuity_and_Btheta_jump (all theta)"

! ── Summary ──────────────────────────────────────────────────────────────────
IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - dipole_uniform_sphere"
END IF

END PROGRAM DIPOLE_UNIFORM_SPHERE
