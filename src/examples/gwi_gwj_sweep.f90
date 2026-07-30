!     Worked example (2e): Geppert-Wiebicke (GWI/GWJ) coefficients swept
!     across a range of mode-number combinations, each checked against
!     an independent formula (not just the closed-form GWI/GWJ code
!     itself evaluated a second way).
!
!     GWI check: GWI_DP's own formula,
!       sqrt((2J1+1)(2J2+1)/(2L+1)/4pi) * CG(J1,0,J2,0,L,0) * CG(J1,M1,J2,M2,L,M)
!     is algebraically the standard closed-form "Gaunt coefficient" --
!     the textbook expression for the triple-scalar-harmonic integral
!         GWI(J1,M1,J2,M2,L,M) = int Y_J1^M1 Y_J2^M2 conj(Y_L^M) dOmega
!     This was verified numerically against direct 2D quadrature (SHGLQ
!     in theta x uniform trapezoidal in phi, the same technique
!     introduced in vsh_decomposition.f90) before writing this sweep,
!     confirming the identification to ~1e-16 across several mode
!     combinations, including a selection-rule-forbidden case that
!     correctly evaluates to exactly zero on both sides.
!
!     GWJ check: GWJ's defining relationship is the coefficient in the
!     expansion of a poloidal-toroidal VSH dot product in scalar
!     harmonics (the same pattern src/tests.f90's TEST2_GW/TEST2_VSH
!     already validate for one hardcoded case):
!       DOT(PVSH_POL(J2,M2),PVSH_TOR(J1,M1))
!         = -i/sqrt(J1(J1+1)*J2(J2+1)) * sum_L GWJ(J1,M1,J2,M2,L,M1+M2)*SSH(L,M1+M2)
!     summed over L=|J1-J2|..J1+J2 (terms outside the valid triangle/
!     parity range correctly evaluate to zero after the GWJ_DP guard
!     fix). This is a pointwise identity (no integration needed) and was
!     independently verified against several new mode combinations
!     (beyond the one hardcoded in TEST2_GW) before writing this sweep.

PROGRAM GWI_GWJ_SWEEP
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi, j
USE VSH,     ONLY: GWI, GWJ, SSH, PVSH_POL, PVSH_TOR, DOT, SHGLQ
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: JMAX = 3
INTEGER(KIND=i4), PARAMETER :: LMAX_QUAD = 20
INTEGER(KIND=i4), PARAMETER :: NPHI = 16
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-8_dp

INTEGER(KIND=i4), PARAMETER :: N_POINTS = 2
REAL(KIND=dp), PARAMETER :: THETA_PTS(N_POINTS) = (/0.7d0, 1.9d0/)
REAL(KIND=dp), PARAMETER :: PHI_PTS(N_POINTS)   = (/1.3d0, 0.4d0/)

INTEGER(KIND=i4) :: N_FAIL

N_FAIL = 0
CALL RUN_GWI_SWEEP(N_FAIL)
CALL RUN_GWJ_SWEEP(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - gwi_gwj_sweep"
END IF

CONTAINS

!     GWI vs. the brute-force Gaunt-coefficient integral, swept over
!     J1,J2=0..JMAX, all valid M1,M2,L,M.
  SUBROUTINE RUN_GWI_SWEEP(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ZERO(LMAX_QUAD+1), WTHETA(LMAX_QUAD+1)
  REAL(KIND=dp) :: DPHI, THETA, PHI, MAX_ERR, ABSDIFF
  COMPLEX(KIND=dp) :: CLOSED, BRUTE
  INTEGER(KIND=i4) :: J1, J2, M1, M2, L, M, IQ, IP

  CALL SHGLQ(ZERO, WTHETA, LMAX_QUAD)
  DPHI = 2.d0*pi/NPHI
  MAX_ERR = 0.0_dp

  OPEN(UNIT=81, FILE="./examples/gwi_gwj_sweep/gwi_check.dat")
  WRITE(81,'(A)') '# J1  M1  J2  M2  L  M  '// &
                  'Re(closed)  Im(closed)  Re(brute)  Im(brute)  absdiff'

  DO J1 = 0, JMAX
    DO J2 = 0, JMAX
      DO M1 = -J1, J1
        DO M2 = -J2, J2
          DO L = ABS(J1-J2), J1+J2
            M = M1+M2
            IF (ABS(M) > L) CYCLE

            CLOSED = GWI(J1,M1,J2,M2,L,M)

            BRUTE = DCMPLX(0.d0, 0.d0)
            DO IQ = 1, LMAX_QUAD+1
              THETA = DACOS(ZERO(IQ))
              DO IP = 0, NPHI-1
                PHI = IP*DPHI
                BRUTE = BRUTE + WTHETA(IQ)*DPHI * SSH(J1,M1,THETA,PHI) * &
                        SSH(J2,M2,THETA,PHI) * CONJG(SSH(L,M,THETA,PHI))
              ENDDO
            ENDDO

            ABSDIFF = ABS(CLOSED-BRUTE)
            MAX_ERR = MAX(MAX_ERR, ABSDIFF)
            WRITE(81,*) J1, M1, J2, M2, L, M, DBLE(CLOSED), AIMAG(CLOSED), &
                        DBLE(BRUTE), AIMAG(BRUTE), ABSDIFF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CLOSE(81)

  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  gwi_gaunt_sweep  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  gwi_gaunt_sweep  max_err=", MAX_ERR
  END IF
  END SUBROUTINE RUN_GWI_SWEEP

!     GWJ vs. the DOT(PVSH_POL,PVSH_TOR) pointwise identity, swept over
!     J1,J2=0..JMAX, all valid M1,M2, at N_POINTS fixed (theta,phi).
  SUBROUTINE RUN_GWJ_SWEEP(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  COMPLEX(KIND=dp) :: LHS, RHS, PREFACTOR
  REAL(KIND=dp) :: MAX_ERR, ABSDIFF, THETA, PHI
  INTEGER(KIND=i4) :: J1, J2, M1, M2, L, M, IP

  MAX_ERR = 0.0_dp

  OPEN(UNIT=82, FILE="./examples/gwi_gwj_sweep/gwj_check.dat")
  WRITE(82,'(A)') '# J1  M1  J2  M2  theta  phi  '// &
                  'Re(lhs)  Im(lhs)  Re(rhs)  Im(rhs)  absdiff'

  DO J1 = 1, JMAX
    DO J2 = 1, JMAX
      PREFACTOR = -j/DSQRT(DBLE(J1)*(J1+1.d0)*DBLE(J2)*(J2+1.d0))
      DO M1 = -J1, J1
        DO M2 = -J2, J2
          M = M1+M2
          DO IP = 1, N_POINTS
            THETA = THETA_PTS(IP)
            PHI   = PHI_PTS(IP)

            LHS = DOT(PVSH_POL(J2,M2,THETA,PHI), PVSH_TOR(J1,M1,THETA,PHI))

            RHS = DCMPLX(0.d0, 0.d0)
            DO L = ABS(J1-J2), J1+J2
              IF (ABS(M) <= L) THEN
                RHS = RHS + GWJ(J1,M1,J2,M2,L,M)*SSH(L,M,THETA,PHI)
              END IF
            ENDDO
            RHS = PREFACTOR * RHS

            ABSDIFF = ABS(LHS-RHS)
            MAX_ERR = MAX(MAX_ERR, ABSDIFF)
            WRITE(82,*) J1, M1, J2, M2, THETA, PHI, DBLE(LHS), AIMAG(LHS), &
                        DBLE(RHS), AIMAG(RHS), ABSDIFF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CLOSE(82)

  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  gwj_identity_sweep  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  gwj_identity_sweep  max_err=", MAX_ERR
  END IF
  END SUBROUTINE RUN_GWJ_SWEEP

END PROGRAM GWI_GWJ_SWEEP
