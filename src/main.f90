!> Driver for the full [[TESTS]] validation suite (run via `ctest`, target
!> `vsh_validation`): runs every cross-validation, batch-consistency,
!> identity, and Wigner-D/rotation check in [[TESTS]], writes each one's
!> full per-point data to `./validation/*.dat`, prints a `PASS`/`FAIL`
!> line per check, and exits non-zero if any check failed (the signal
!> `ctest` uses for a pass/fail build gate).
PROGRAM MAIN
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE TESTS,   ONLY: &
  TEST1_ANALYTIC, TEST1_VSH, TEST1_GW, &
  TEST2_ANALYTIC, TEST2_VSH, TEST2_GW, &
  BATCH_ALM_CONS, BATCH_DALM_CONS, BATCH_SSH_CONS, &
  SSH_ORTHO, BATCH_GRAD_SSH_CONS, BATCH_L_SSH_CONS, &
  BATCH_PVSH_RAD_CONS, BATCH_PVSH_POL_CONS, BATCH_PVSH_TOR_CONS, &
  BATCH_VSH_TOR_CONS, BATCH_VSH_POL_UP_CONS, BATCH_VSH_POL_DN_CONS, &
  PVSH_POL_TOR_ORTHO, VSH_POL_INVERSION, &
  WIGNER_D_SMALL_IDENTITY, WIGNER_D_SMALL_SYMMETRY, &
  WIGNER_D_IDENTITY, WIGNER_D_UNITARITY, WIGNER_D_COMPOSITION, &
  ROTATE_SSH_INVERSION, ROTATE_PVSH_INVERSION, ROTATE_VSH_STD_INVERSION, &
  ROTATE_PVSH_SPECTRUM_GT, ROTATE_VSH_STD_SPECTRUM_GT, &
  WIGNER_D_SPOTCHECK
IMPLICIT NONE

INTEGER(KIND=i4) :: NTHETA, I, STATUS, N_FAIL
PARAMETER (NTHETA=100)
COMPLEX(KIND=dp), DIMENSION(NTHETA) :: TEST1A, TEST1B, TEST1C
COMPLEX(KIND=dp), DIMENSION(NTHETA) :: TEST2A, TEST2B, TEST2C
REAL(KIND=dp)    :: MAX_ERR

N_FAIL = 0

! ── Cross-validation tests ──────────────────────────────────────────────────
TEST1A = TEST1_ANALYTIC(NTHETA)
TEST1B = TEST1_VSH(NTHETA)
TEST1C = TEST1_GW(NTHETA)

TEST2A = TEST2_ANALYTIC(NTHETA)
TEST2B = TEST2_VSH(NTHETA)
TEST2C = TEST2_GW(NTHETA)

OPEN(UNIT=11, FILE="./validation/TEST1.dat")
OPEN(UNIT=12, FILE="./validation/TEST2.dat")
DO I=1,NTHETA
  WRITE(11,*) (I-1)*pi/(NTHETA-1),TEST1A(I),TEST1B(I),TEST1C(I)
ENDDO
DO I=1,NTHETA
  WRITE(12,*) (I-1)*pi/(NTHETA-1),TEST2A(I),TEST2B(I),TEST2C(I)
ENDDO
CLOSE(11)
CLOSE(12)

MAX_ERR = MAXVAL(ABS(TEST1A - TEST1B))
IF (MAX_ERR > 1.0E-12_dp) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  test1_analytic_vs_vsh       max_err=", MAX_ERR
ELSE
  WRITE(*,'(A)') "PASS  test1_analytic_vs_vsh"
END IF

MAX_ERR = MAXVAL(ABS(TEST1A - TEST1C))
IF (MAX_ERR > 1.0E-12_dp) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  test1_analytic_vs_gw        max_err=", MAX_ERR
ELSE
  WRITE(*,'(A)') "PASS  test1_analytic_vs_gw"
END IF

MAX_ERR = MAXVAL(ABS(TEST2A - TEST2B))
IF (MAX_ERR > 1.0E-12_dp) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  test2_analytic_vs_vsh       max_err=", MAX_ERR
ELSE
  WRITE(*,'(A)') "PASS  test2_analytic_vs_vsh"
END IF

MAX_ERR = MAXVAL(ABS(TEST2A - TEST2C))
IF (MAX_ERR > 1.0E-12_dp) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  test2_analytic_vs_gw        max_err=", MAX_ERR
ELSE
  WRITE(*,'(A)') "PASS  test2_analytic_vs_gw"
END IF

! ── Batch consistency tests ──────────────────────────────────────────────────
OPEN(UNIT=21, FILE="./validation/batch_alm_cons.dat")
CALL BATCH_ALM_CONS(8, 20, 21, STATUS)
CLOSE(21)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_alm_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_alm_cons"
END IF

OPEN(UNIT=22, FILE="./validation/batch_dalm_cons.dat")
CALL BATCH_DALM_CONS(8, 20, 22, STATUS)
CLOSE(22)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_dalm_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_dalm_cons"
END IF

OPEN(UNIT=23, FILE="./validation/batch_ssh_cons.dat")
CALL BATCH_SSH_CONS(5, 15, 23, STATUS)
CLOSE(23)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_ssh_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_ssh_cons"
END IF

OPEN(UNIT=24, FILE="./validation/ssh_ortho.dat")
CALL SSH_ORTHO(6, 10000, 24, STATUS)
CLOSE(24)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  ssh_ortho"
ELSE
  WRITE(*,'(A)') "PASS  ssh_ortho"
END IF

OPEN(UNIT=31, FILE="./validation/batch_grad_ssh_cons.dat")
CALL BATCH_GRAD_SSH_CONS(5, 15, 31, STATUS)
CLOSE(31)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_grad_ssh_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_grad_ssh_cons"
END IF

OPEN(UNIT=32, FILE="./validation/batch_l_ssh_cons.dat")
CALL BATCH_L_SSH_CONS(5, 15, 32, STATUS)
CLOSE(32)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_l_ssh_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_l_ssh_cons"
END IF

OPEN(UNIT=33, FILE="./validation/batch_pvsh_rad_cons.dat")
CALL BATCH_PVSH_RAD_CONS(5, 15, 33, STATUS)
CLOSE(33)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_pvsh_rad_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_pvsh_rad_cons"
END IF

OPEN(UNIT=34, FILE="./validation/batch_pvsh_pol_cons.dat")
CALL BATCH_PVSH_POL_CONS(5, 15, 34, STATUS)
CLOSE(34)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_pvsh_pol_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_pvsh_pol_cons"
END IF

OPEN(UNIT=35, FILE="./validation/batch_pvsh_tor_cons.dat")
CALL BATCH_PVSH_TOR_CONS(5, 15, 35, STATUS)
CLOSE(35)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_pvsh_tor_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_pvsh_tor_cons"
END IF

OPEN(UNIT=36, FILE="./validation/batch_vsh_tor_cons.dat")
CALL BATCH_VSH_TOR_CONS(5, 15, 36, STATUS)
CLOSE(36)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_vsh_tor_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_vsh_tor_cons"
END IF

OPEN(UNIT=37, FILE="./validation/batch_vsh_pol_up_cons.dat")
CALL BATCH_VSH_POL_UP_CONS(5, 15, 37, STATUS)
CLOSE(37)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_vsh_pol_up_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_vsh_pol_up_cons"
END IF

OPEN(UNIT=38, FILE="./validation/batch_vsh_pol_dn_cons.dat")
CALL BATCH_VSH_POL_DN_CONS(5, 15, 38, STATUS)
CLOSE(38)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  batch_vsh_pol_dn_cons"
ELSE
  WRITE(*,'(A)') "PASS  batch_vsh_pol_dn_cons"
END IF

OPEN(UNIT=39, FILE="./validation/pvsh_pol_tor_ortho.dat")
CALL PVSH_POL_TOR_ORTHO(5, 15, 39, STATUS)
CLOSE(39)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  pvsh_pol_tor_ortho"
ELSE
  WRITE(*,'(A)') "PASS  pvsh_pol_tor_ortho"
END IF

OPEN(UNIT=40, FILE="./validation/vsh_pol_inversion.dat")
CALL VSH_POL_INVERSION(5, 15, 40, STATUS)
CLOSE(40)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  vsh_pol_inversion"
ELSE
  WRITE(*,'(A)') "PASS  vsh_pol_inversion"
END IF

! ── Wigner-D / rotation validation ──────────────────────────────────────────────
OPEN(UNIT=41, FILE="./validation/wigner_d_small_identity.dat")
CALL WIGNER_D_SMALL_IDENTITY(6, 41, STATUS)
CLOSE(41)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  wigner_d_small_identity"
ELSE
  WRITE(*,'(A)') "PASS  wigner_d_small_identity"
END IF

OPEN(UNIT=42, FILE="./validation/wigner_d_small_symmetry.dat")
CALL WIGNER_D_SMALL_SYMMETRY(6, 42, STATUS)
CLOSE(42)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  wigner_d_small_symmetry"
ELSE
  WRITE(*,'(A)') "PASS  wigner_d_small_symmetry"
END IF

OPEN(UNIT=43, FILE="./validation/wigner_d_identity.dat")
CALL WIGNER_D_IDENTITY(6, 43, STATUS)
CLOSE(43)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  wigner_d_identity"
ELSE
  WRITE(*,'(A)') "PASS  wigner_d_identity"
END IF

OPEN(UNIT=44, FILE="./validation/wigner_d_unitarity.dat")
CALL WIGNER_D_UNITARITY(6, 44, STATUS)
CLOSE(44)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  wigner_d_unitarity"
ELSE
  WRITE(*,'(A)') "PASS  wigner_d_unitarity"
END IF

OPEN(UNIT=45, FILE="./validation/wigner_d_composition.dat")
CALL WIGNER_D_COMPOSITION(6, 45, STATUS)
CLOSE(45)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  wigner_d_composition"
ELSE
  WRITE(*,'(A)') "PASS  wigner_d_composition"
END IF

OPEN(UNIT=46, FILE="./validation/rotate_ssh_inversion.dat")
CALL ROTATE_SSH_INVERSION(6, 46, STATUS)
CLOSE(46)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  rotate_ssh_inversion"
ELSE
  WRITE(*,'(A)') "PASS  rotate_ssh_inversion"
END IF

OPEN(UNIT=47, FILE="./validation/rotate_pvsh_inversion.dat")
CALL ROTATE_PVSH_INVERSION(6, 47, STATUS)
CLOSE(47)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  rotate_pvsh_inversion"
ELSE
  WRITE(*,'(A)') "PASS  rotate_pvsh_inversion"
END IF

OPEN(UNIT=48, FILE="./validation/rotate_vsh_std_inversion.dat")
CALL ROTATE_VSH_STD_INVERSION(6, 48, STATUS)
CLOSE(48)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  rotate_vsh_std_inversion"
ELSE
  WRITE(*,'(A)') "PASS  rotate_vsh_std_inversion"
END IF

OPEN(UNIT=49, FILE="./validation/rotate_pvsh_spectrum_gt.dat")
CALL ROTATE_PVSH_SPECTRUM_GT(3, 200, 32, 49, STATUS)
CLOSE(49)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  rotate_pvsh_spectrum_gt"
ELSE
  WRITE(*,'(A)') "PASS  rotate_pvsh_spectrum_gt"
END IF

OPEN(UNIT=50, FILE="./validation/rotate_vsh_std_spectrum_gt.dat")
CALL ROTATE_VSH_STD_SPECTRUM_GT(3, 200, 32, 50, STATUS)
CLOSE(50)
IF (STATUS /= 0) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  rotate_vsh_std_spectrum_gt"
ELSE
  WRITE(*,'(A)') "PASS  rotate_vsh_std_spectrum_gt"
END IF

OPEN(UNIT=51, FILE="./validation/wigner_d_spotcheck.dat")
CALL WIGNER_D_SPOTCHECK(51)
CLOSE(51)
WRITE(*,'(A)') "WROTE wigner_d_spotcheck (cross-check via py/sympy_reference.py)"

! ── Summary ──────────────────────────────────────────────────────────────────
IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " test(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - all tests"
END IF

END PROGRAM MAIN
