!     Worked example (2d): arbitrary-field VSH spectral decomposition,
!     reconstruction, and convergence.
!
!     Demonstrates the manuscript Introduction's own claim that VSH
!     "constitute a complete orthonormal basis" enabling "spectral
!     expansion of MHD vector fields with arbitrary angular
!     configuration" -- using a genuine 2D (theta,phi) numerical
!     quadrature built from SHGLQ (Gauss-Legendre in x=cos(theta), exact
!     for the polynomial degrees involved here) times a uniform
!     trapezoidal sum in phi (spectrally exact for finite Fourier
!     content). SHGLQ is defined in the library but was not called
!     anywhere before this example.
!
!     Part A: a field built from equal unit coefficients on exactly the
!     eight modes (l,m) = (1,-1),(1,0),(1,1),(2,-2),(2,-1),(2,0),(2,1),
!     (2,2), applied to all three basis families (radial, poloidal,
!     toroidal). Decompose and confirm the recovered coefficients equal
!     1 at those eight modes and ~0 (quadrature floor) everywhere else --
!     the first genuinely non-axisymmetric (m/=0) purity check in this
!     repository (the dipole examples are all m=0).
!
!     Part B: a smooth field with effectively unbounded spectral content
!     (a Gaussian angular bump with slight phi modulation) is
!     reconstructed at increasing Lmax, and the reconstruction error
!     (checked against the true closed-form field on an independent
!     evaluation grid, not the quadrature grid) is shown to shrink as
!     more modes are retained -- the classic spectral truncation
!     convergence curve.

PROGRAM VSH_DECOMPOSITION
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: SHGLQ, PVSH_RAD, PVSH_POL, PVSH_TOR, &
                    PVSH_RAD_ALL, PVSH_POL_ALL, PVSH_TOR_ALL, &
                    YLM_INDEX, DOT
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_KNOWN = 8
INTEGER(KIND=i4), PARAMETER :: L_KNOWN(N_KNOWN) = (/1,1,1,2,2,2,2,2/)
INTEGER(KIND=i4), PARAMETER :: M_KNOWN(N_KNOWN) = (/-1,0,1,-2,-1,0,1,2/)
INTEGER(KIND=i4), PARAMETER :: LMAX_A = 4
REAL(KIND=dp),    PARAMETER :: TOL_A = 1.0E-8_dp

INTEGER(KIND=i4), PARAMETER :: N_LMAX_B = 6
INTEGER(KIND=i4), PARAMETER :: LMAX_B_LIST(N_LMAX_B) = (/2,4,8,16,24,32/)
REAL(KIND=dp),    PARAMETER :: SIGMA_B = 0.4_dp

INTEGER(KIND=i4) :: N_FAIL

N_FAIL = 0

CALL RUN_PART_A(N_FAIL)
CALL RUN_PART_B(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - vsh_decomposition"
END IF

CONTAINS

!     Target field for Part A: equal unit weight on all three basis
!     families at each of the eight known modes; zero elsewhere. A
!     genuine finite sum over eight independent modes, not a shortcut.
  FUNCTION FIELD_KNOWN(THETA, PHI) RESULT(F)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: THETA, PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: F
  INTEGER(KIND=i4) :: K
  F = DCMPLX(0.d0, 0.d0)
  DO K = 1, N_KNOWN
    F = F + PVSH_RAD(L_KNOWN(K), M_KNOWN(K), THETA, PHI) &
          + PVSH_POL(L_KNOWN(K), M_KNOWN(K), THETA, PHI) &
          + PVSH_TOR(L_KNOWN(K), M_KNOWN(K), THETA, PHI)
  ENDDO
  END FUNCTION FIELD_KNOWN

!     Target field for Part B: a smooth angular Gaussian bump with slight
!     phi modulation (carries m=0 and m=+/-1 content, on top of
!     effectively unbounded l content from the Gaussian profile in
!     theta). Purely radial, for simplicity.
!
!     The phi modulation is weighted by sin(theta) so the field stays
!     single-valued at the poles: at theta=0 (or pi), every phi maps to
!     the same physical point, so any legitimate smooth field must be
!     phi-independent there. An earlier version of this field used a
!     bare cos(phi) term with no such weighting -- it was multi-valued
!     at the pole (1.5 approaching along phi=0, 0.5 along phi=pi), and
!     no finite spherical harmonic series can converge to a discontinuous
!     function; the reconstruction error plateaued at exactly 0.5 (the
!     size of the discontinuity) instead of shrinking with Lmax. Caught
!     by comparing this driver's convergence curve against an isolated,
!     single-Lmax reconstruction check, which did not show the plateau.
  FUNCTION FIELD_GAUSSIAN(THETA, PHI) RESULT(F)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: THETA, PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: F
  REAL(KIND=dp) :: BUMP
  BUMP = DEXP(-THETA**2/(2.d0*SIGMA_B**2)) * &
         (1.d0 + 0.5d0*DSIN(THETA)*DCOS(PHI))
  F(1) = DCMPLX(BUMP, 0.d0)
  F(2) = DCMPLX(0.d0, 0.d0)
  F(3) = DCMPLX(0.d0, 0.d0)
  END FUNCTION FIELD_GAUSSIAN

!     Decompose FIELD_FN into (a_rad,a_pol,a_tor) for 0<=l<=LMAX,
!     -l<=m<=l, via 2D Gauss-Legendre(theta) x uniform-trapezoidal(phi)
!     quadrature. LMAX_QUAD sets the theta quadrature order (LMAX_QUAD+1
!     nodes); NPHI sets the number of phi points.
  SUBROUTINE DECOMPOSE(FIELD_FN, LMAX, LMAX_QUAD, NPHI, A_RAD, A_POL, A_TOR)
  IMPLICIT NONE
  INTERFACE
    FUNCTION FIELD_FN(THETA, PHI) RESULT(F)
      IMPORT :: dp
      REAL(KIND=dp), INTENT(IN) :: THETA, PHI
      COMPLEX(KIND=dp), DIMENSION(3) :: F
    END FUNCTION FIELD_FN
  END INTERFACE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX, LMAX_QUAD, NPHI
  COMPLEX(KIND=dp), INTENT(OUT) :: A_RAD((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(OUT) :: A_POL((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(OUT) :: A_TOR((LMAX+1)**2)
  REAL(KIND=dp), ALLOCATABLE :: ZERO(:), WTHETA(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: RADALL(:,:), POLALL(:,:), TORALL(:,:)
  COMPLEX(KIND=dp), DIMENSION(3) :: FVEC
  REAL(KIND=dp) :: THETA, PHI, DPHI, WEIGHT
  INTEGER(KIND=i4) :: NYLM, IQ, IP, IDX

  NYLM = (LMAX+1)**2
  ALLOCATE(ZERO(LMAX_QUAD+1), WTHETA(LMAX_QUAD+1))
  ALLOCATE(RADALL(3,NYLM), POLALL(3,NYLM), TORALL(3,NYLM))
  CALL SHGLQ(ZERO, WTHETA, LMAX_QUAD)

  DPHI = 2.d0*pi/NPHI
  A_RAD = DCMPLX(0.d0, 0.d0)
  A_POL = DCMPLX(0.d0, 0.d0)
  A_TOR = DCMPLX(0.d0, 0.d0)

  DO IQ = 1, LMAX_QUAD+1
    THETA = DACOS(ZERO(IQ))
    DO IP = 0, NPHI-1
      PHI = IP*DPHI
      FVEC = FIELD_FN(THETA, PHI)
      CALL PVSH_RAD_ALL(RADALL, LMAX, THETA, PHI)
      CALL PVSH_POL_ALL(POLALL, LMAX, THETA, PHI)
      CALL PVSH_TOR_ALL(TORALL, LMAX, THETA, PHI)
      WEIGHT = WTHETA(IQ) * DPHI
      DO IDX = 1, NYLM
        A_RAD(IDX) = A_RAD(IDX) + WEIGHT * DOT(FVEC, CONJG(RADALL(:,IDX)))
        A_POL(IDX) = A_POL(IDX) + WEIGHT * DOT(FVEC, CONJG(POLALL(:,IDX)))
        A_TOR(IDX) = A_TOR(IDX) + WEIGHT * DOT(FVEC, CONJG(TORALL(:,IDX)))
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE(ZERO, WTHETA, RADALL, POLALL, TORALL)
  END SUBROUTINE DECOMPOSE

!     Reconstruct the field at (THETA,PHI) from coefficients up to LMAX.
  FUNCTION RECONSTRUCT(LMAX, A_RAD, A_POL, A_TOR, THETA, PHI) RESULT(F)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  COMPLEX(KIND=dp), INTENT(IN) :: A_RAD((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(IN) :: A_POL((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(IN) :: A_TOR((LMAX+1)**2)
  REAL(KIND=dp), INTENT(IN) :: THETA, PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: F
  COMPLEX(KIND=dp), ALLOCATABLE :: RADALL(:,:), POLALL(:,:), TORALL(:,:)
  INTEGER(KIND=i4) :: NYLM, IDX
  NYLM = (LMAX+1)**2
  ALLOCATE(RADALL(3,NYLM), POLALL(3,NYLM), TORALL(3,NYLM))
  CALL PVSH_RAD_ALL(RADALL, LMAX, THETA, PHI)
  CALL PVSH_POL_ALL(POLALL, LMAX, THETA, PHI)
  CALL PVSH_TOR_ALL(TORALL, LMAX, THETA, PHI)
  F = DCMPLX(0.d0, 0.d0)
  DO IDX = 1, NYLM
    F = F + A_RAD(IDX)*RADALL(:,IDX) + A_POL(IDX)*POLALL(:,IDX) &
          + A_TOR(IDX)*TORALL(:,IDX)
  ENDDO
  DEALLOCATE(RADALL, POLALL, TORALL)
  END FUNCTION RECONSTRUCT

!     Part A: decompose FIELD_KNOWN, confirm exact recovery of the eight
!     known unit coefficients and near-zero everywhere else.
  SUBROUTINE RUN_PART_A(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: LMAX_QUAD_A = 30
  INTEGER(KIND=i4), PARAMETER :: NPHI_A = 32
  COMPLEX(KIND=dp) :: A_RAD((LMAX_A+1)**2), A_POL((LMAX_A+1)**2), &
                       A_TOR((LMAX_A+1)**2)
  REAL(KIND=dp) :: EXPECTED, MAX_ERR
  INTEGER(KIND=i4) :: L, M, IDX, K

  CALL DECOMPOSE(FIELD_KNOWN, LMAX_A, LMAX_QUAD_A, NPHI_A, A_RAD, A_POL, A_TOR)

  OPEN(UNIT=71, FILE="./examples/vsh_decomposition/known_modes.dat")
  WRITE(71,'(A)') '# L  M  Re(a_rad)  Im(a_rad)  Re(a_pol)  Im(a_pol)  '// &
                  'Re(a_tor)  Im(a_tor)  expected'
  MAX_ERR = 0.0_dp
  DO L = 0, LMAX_A
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      EXPECTED = 0.0_dp
      DO K = 1, N_KNOWN
        IF (L_KNOWN(K) == L .AND. M_KNOWN(K) == M) EXPECTED = 1.0_dp
      ENDDO
      MAX_ERR = MAX(MAX_ERR, ABS(A_RAD(IDX)-EXPECTED), &
                    ABS(A_POL(IDX)-EXPECTED), ABS(A_TOR(IDX)-EXPECTED))
      WRITE(71,*) L, M, DBLE(A_RAD(IDX)), AIMAG(A_RAD(IDX)), &
                  DBLE(A_POL(IDX)), AIMAG(A_POL(IDX)), &
                  DBLE(A_TOR(IDX)), AIMAG(A_TOR(IDX)), EXPECTED
    ENDDO
  ENDDO
  CLOSE(71)

  IF (MAX_ERR > TOL_A) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  known_mode_recovery  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  known_mode_recovery  max_err=", MAX_ERR
  END IF
  END SUBROUTINE RUN_PART_A

!     Part B: decompose+reconstruct FIELD_GAUSSIAN at increasing Lmax,
!     tracking max reconstruction error against the true closed-form
!     field on an independent evaluation grid. Self-check: error at the
!     largest Lmax must be much smaller than at the smallest (genuine
!     convergence, not a flat/noisy curve -- which would indicate a
!     quadrature-order-vs-Lmax scaling bug).
  SUBROUTINE RUN_PART_B(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: NEVAL = 40
  REAL(KIND=dp),    PARAMETER :: CONVERGENCE_RATIO_TOL = 0.01_dp
  INTEGER(KIND=i4) :: IB, LMAX_TEST, LMAX_QUAD, NPHI, ITH, IPH
  COMPLEX(KIND=dp), ALLOCATABLE :: A_RAD(:), A_POL(:), A_TOR(:)
  COMPLEX(KIND=dp), DIMENSION(3) :: F_TRUE, F_RECON
  REAL(KIND=dp) :: THETA, PHI, MAX_ERR, FIRST_ERR, LAST_ERR

  OPEN(UNIT=72, FILE="./examples/vsh_decomposition/convergence.dat")
  WRITE(72,'(A)') '# Lmax  max_reconstruction_error'

  DO IB = 1, N_LMAX_B
    LMAX_TEST = LMAX_B_LIST(IB)
    LMAX_QUAD = 3*LMAX_TEST + 10
    NPHI = 4*LMAX_TEST + 16

    ALLOCATE(A_RAD((LMAX_TEST+1)**2), A_POL((LMAX_TEST+1)**2), &
             A_TOR((LMAX_TEST+1)**2))
    CALL DECOMPOSE(FIELD_GAUSSIAN, LMAX_TEST, LMAX_QUAD, NPHI, &
                   A_RAD, A_POL, A_TOR)

    MAX_ERR = 0.0_dp
    DO ITH = 0, NEVAL-1
      THETA = ITH*pi/(NEVAL-1)
      DO IPH = 0, NEVAL-1
        PHI = IPH*2.d0*pi/NEVAL
        F_TRUE  = FIELD_GAUSSIAN(THETA, PHI)
        F_RECON = RECONSTRUCT(LMAX_TEST, A_RAD, A_POL, A_TOR, THETA, PHI)
        MAX_ERR = MAX(MAX_ERR, ABS(F_RECON(1)-F_TRUE(1)), &
                      ABS(F_RECON(2)-F_TRUE(2)), ABS(F_RECON(3)-F_TRUE(3)))
      ENDDO
    ENDDO

    IF (IB == 1) FIRST_ERR = MAX_ERR
    IF (IB == N_LMAX_B) LAST_ERR = MAX_ERR

    WRITE(72,*) LMAX_TEST, MAX_ERR
    WRITE(*,'(A,I0,A,ES10.3)') "  Lmax=", LMAX_TEST, &
      "  max_reconstruction_error=", MAX_ERR

    DEALLOCATE(A_RAD, A_POL, A_TOR)
  ENDDO
  CLOSE(72)

  IF (LAST_ERR > CONVERGENCE_RATIO_TOL*FIRST_ERR) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3,A,ES10.3)') "FAIL  spectral_convergence  first=", &
      FIRST_ERR, "  last=", LAST_ERR
  ELSE
    WRITE(*,'(A,ES10.3,A,ES10.3)') "PASS  spectral_convergence  first=", &
      FIRST_ERR, "  last=", LAST_ERR
  END IF
  END SUBROUTINE RUN_PART_B

END PROGRAM VSH_DECOMPOSITION
