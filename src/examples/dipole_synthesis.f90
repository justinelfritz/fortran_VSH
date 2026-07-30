!> Worked example (2a): synthesize the classic exterior vacuum magnetic
!> dipole (\( \ell=1,\,m=0 \)) from FORTVSH's poloidal basis functions and
!> confirm (i) it reproduces the textbook closed form
!> \( B_r=2\mu\cos\theta/r^3,\,B_\theta=\mu\sin\theta/r^3 \), and (ii) it
!> is numerically pure -- its overlap with every other \( \ell \) (at
!> \( m=0 \)) is zero to quadrature precision. No induction equation or
!> time evolution is involved: this is a static field-synthesis and
!> orthogonal-projection check only.
!>
!> Poloidal decomposition used (derived directly from FORTVSH's own
!> normalization, [[PVSH_RAD]]'s radial component is \( Y_\ell^m \) and
!> [[PVSH_POL]]'s \( \hat\theta \) component is
!> \( \partial_\theta Y_\ell^m/\sqrt{\ell(\ell+1)} \)):
!> $$ B_r = \frac{\ell(\ell+1)}{r^2}\,S(r)\,\mathrm{[[PVSH_RAD]]}(\ell,m,\theta,\phi)_r, $$
!> $$ B_\theta = \frac{\sqrt{\ell(\ell+1)}}{r}\,\frac{dS}{dr}\,
!>    \mathrm{[[PVSH_POL]]}(\ell,m,\theta,\phi)_\theta. $$
!> For \( \ell=1 \) this reproduces the exterior vacuum dipole with the
!> radial stream function \( S(r)=C/r,\,C=\mu/\sqrt{3/(4\pi)} \).
PROGRAM DIPOLE_SYNTHESIS
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: SSH, PVSH_RAD, PVSH_POL
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: NTHETA = 181     !- field-profile grid
INTEGER(KIND=i4), PARAMETER :: NQUAD  = 20000   !- purity-check quadrature
INTEGER(KIND=i4), PARAMETER :: LTESTMIN = 0, LTESTMAX = 3
REAL(KIND=dp),    PARAMETER :: MU = 1.0_dp      !- dipole moment (normalized units)
REAL(KIND=dp),    PARAMETER :: R0 = 1.0_dp      !- field-profile evaluation radius
REAL(KIND=dp),    PARAMETER :: SYNTH_TOL = 1.0E-12_dp
REAL(KIND=dp),    PARAMETER :: PURITY_TOL = 1.0E-4_dp

REAL(KIND=dp) :: THETA, PHI, U, DU
REAL(KIND=dp) :: BR_ANALYTIC, BTH_ANALYTIC, BR_LIB, BTH_LIB
REAL(KIND=dp) :: C_STREAM, DSDR, MAX_ERR
COMPLEX(KIND=dp), DIMENSION(3) :: RADVEC, POLVEC
COMPLEX(KIND=dp) :: CRAD(LTESTMIN:LTESTMAX), CPOL(LTESTMIN:LTESTMAX)
INTEGER(KIND=i4) :: I, L, IQ, N_FAIL

PHI = 0.0_dp   !- m=0: field independent of phi
N_FAIL = 0

! ── Field profile: library synthesis vs. closed-form dipole ─────────────────
C_STREAM = MU / DSQRT(3.0_dp/(4.0_dp*pi))
DSDR = -C_STREAM / R0**2   !- d/dr (C/r) at r=R0

OPEN(UNIT=61, FILE="./examples/dipole_synthesis/field_profile.dat")
WRITE(61,'(A)') '# theta  Br_lib  Btheta_lib  Br_analytic  Btheta_analytic'

MAX_ERR = 0.0_dp
DO I = 0, NTHETA-1
  THETA = I*pi/(NTHETA-1)
  BR_ANALYTIC  = 2.0_dp*MU*DCOS(THETA)/R0**3
  BTH_ANALYTIC =         MU*DSIN(THETA)/R0**3

  RADVEC = PVSH_RAD(1, 0, THETA, PHI)
  POLVEC = PVSH_POL(1, 0, THETA, PHI)
  BR_LIB  = (2.0_dp*C_STREAM/R0**2) * DBLE(RADVEC(1))
  BTH_LIB = (DSQRT(2.0_dp)/R0 * DSDR) * DBLE(POLVEC(2))

  MAX_ERR = MAX(MAX_ERR, ABS(BR_LIB-BR_ANALYTIC), ABS(BTH_LIB-BTH_ANALYTIC))
  WRITE(61,*) THETA, BR_LIB, BTH_LIB, BR_ANALYTIC, BTH_ANALYTIC
ENDDO
CLOSE(61)

IF (MAX_ERR > SYNTH_TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  dipole_synthesis_accuracy  max_err=", MAX_ERR
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  dipole_synthesis_accuracy  max_err=", MAX_ERR
END IF

! ── Purity check: project the analytic dipole onto l=0..3 at m=0 ────────────
! Integrate over u=cos(theta) in [-1,1] via midpoint rule (matches the
! quadrature already used and validated by SSH_ORTHO in src/tests.f90),
! which absorbs the sin(theta) solid-angle weight; multiply by 2*pi for
! the trivial phi integral (m=0: no phi dependence).
DU = 2.0_dp / NQUAD
CRAD = DCMPLX(0.0_dp, 0.0_dp)
CPOL = DCMPLX(0.0_dp, 0.0_dp)

DO IQ = 1, NQUAD
  U = -1.0_dp + (IQ - 0.5_dp) * DU
  THETA = DACOS(U)
  BR_ANALYTIC  = 2.0_dp*MU*DCOS(THETA)/R0**3
  BTH_ANALYTIC =         MU*DSIN(THETA)/R0**3
  DO L = LTESTMIN, LTESTMAX
    CRAD(L) = CRAD(L) + BR_ANALYTIC  * CONJG(SSH(L, 0, THETA, PHI))
    POLVEC  = PVSH_POL(L, 0, THETA, PHI)
    CPOL(L) = CPOL(L) + BTH_ANALYTIC * CONJG(POLVEC(2))
  ENDDO
ENDDO
CRAD = CRAD * DU * 2.0_dp*pi
CPOL = CPOL * DU * 2.0_dp*pi

OPEN(UNIT=62, FILE="./examples/dipole_synthesis/mode_overlap.dat")
WRITE(62,'(A)') '# L  |overlap_rad|  |overlap_pol|'
DO L = LTESTMIN, LTESTMAX
  WRITE(62,*) L, ABS(CRAD(L)), ABS(CPOL(L))
  IF (L /= 1) THEN
    IF (ABS(CRAD(L)) > PURITY_TOL .OR. ABS(CPOL(L)) > PURITY_TOL) THEN
      N_FAIL = N_FAIL + 1
      WRITE(*,'(A,I0,A,ES10.3,A,ES10.3)') &
        "FAIL  dipole_purity l=", L, "  |rad|=", ABS(CRAD(L)), &
        "  |pol|=", ABS(CPOL(L))
    ELSE
      WRITE(*,'(A,I0)') "PASS  dipole_purity l=", L
    END IF
  END IF
ENDDO
CLOSE(62)

! ── Summary ──────────────────────────────────────────────────────────────────
IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - dipole_synthesis"
END IF

END PROGRAM DIPOLE_SYNTHESIS
