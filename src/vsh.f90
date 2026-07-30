!     Function definitions for scalar spherical harmonics (SSH),
!     vector spherical harmonics (VSH), and related
!     angular functions

MODULE VSH
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi, j
IMPLICIT NONE
PRIVATE
PUBLIC :: &
  ! Legendre polynomials and associated Legendre functions
  LEGENDRE, DDX_LEGENDRE, &
  ASSOC_LEGENDRE, DDX_ASSOC_LEGENDRE, &
  ASSOC_LEGENDRE_ALL, DDX_ASSOC_LEGENDRE_ALL, &
  ASSOC_LEGENDRE_NORM_ALL, DDX_ASSOC_LEGENDRE_NORM_ALL, &
  ! Index helpers and quadrature
  PLM_INDEX, YLM_INDEX, LOG_FACT, SHGLQ, &
  ! Scalar spherical harmonics
  SSH, SSH_ALL, &
  ! SSH gradient and L operator
  GRAD_SSH, GRAD_SSH_ALL, L_SSH, L_SSH_ALL, &
  ! Polar vector spherical harmonics
  PVSH_RAD, PVSH_RAD_ALL, &
  PVSH_TOR, PVSH_TOR_ALL, &
  PVSH_POL, PVSH_POL_ALL, &
  ! Standard vector spherical harmonics
  VSH_TOR, VSH_TOR_ALL, &
  VSH_POL_DN, VSH_POL_DN_ALL, &
  VSH_POL_UP, VSH_POL_UP_ALL, &
  ! Coupling coefficients and integrals
  CGCOEFF, SYMBOL3J, DOT, GWI, GWJ

  ! ── Precision-generic interfaces ──────────────────────────────────────────
  ! Each INTERFACE exposes a generic name that currently dispatches to the
  ! dp-precision specific.  When an sp variant is added, insert its name as a
  ! second MODULE PROCEDURE line.  The PUBLIC list and call sites are
  ! unaffected.

  ! -- Scalar Legendre functions --
  INTERFACE LEGENDRE
    MODULE PROCEDURE LEGENDRE_DP
  END INTERFACE LEGENDRE

  INTERFACE DDX_LEGENDRE
    MODULE PROCEDURE DDX_LEGENDRE_DP
  END INTERFACE DDX_LEGENDRE

  INTERFACE ASSOC_LEGENDRE
    MODULE PROCEDURE ASSOC_LEGENDRE_DP
  END INTERFACE ASSOC_LEGENDRE

  INTERFACE DDX_ASSOC_LEGENDRE
    MODULE PROCEDURE DDX_ASSOC_LEGENDRE_DP
  END INTERFACE DDX_ASSOC_LEGENDRE

  ! -- Batch Legendre --
  INTERFACE ASSOC_LEGENDRE_ALL
    MODULE PROCEDURE ASSOC_LEGENDRE_ALL_DP
  END INTERFACE ASSOC_LEGENDRE_ALL

  INTERFACE DDX_ASSOC_LEGENDRE_ALL
    MODULE PROCEDURE DDX_ASSOC_LEGENDRE_ALL_DP
  END INTERFACE DDX_ASSOC_LEGENDRE_ALL

  INTERFACE ASSOC_LEGENDRE_NORM_ALL
    MODULE PROCEDURE ASSOC_LEGENDRE_NORM_ALL_DP
  END INTERFACE ASSOC_LEGENDRE_NORM_ALL

  INTERFACE DDX_ASSOC_LEGENDRE_NORM_ALL
    MODULE PROCEDURE DDX_ASSOC_LEGENDRE_NORM_ALL_DP
  END INTERFACE DDX_ASSOC_LEGENDRE_NORM_ALL

  ! -- Quadrature and log-factorial --
  INTERFACE SHGLQ
    MODULE PROCEDURE SHGLQ_DP
  END INTERFACE SHGLQ

  INTERFACE LOG_FACT
    MODULE PROCEDURE LOG_FACT_DP
  END INTERFACE LOG_FACT

  ! -- Scalar harmonics (single-point) --
  INTERFACE SSH
    MODULE PROCEDURE SSH_DP
  END INTERFACE SSH

  INTERFACE GRAD_SSH
    MODULE PROCEDURE GRAD_SSH_DP
  END INTERFACE GRAD_SSH

  INTERFACE L_SSH
    MODULE PROCEDURE L_SSH_DP
  END INTERFACE L_SSH

  ! -- Scalar harmonics (batch) --
  INTERFACE SSH_ALL
    MODULE PROCEDURE SSH_ALL_DP
  END INTERFACE SSH_ALL

  INTERFACE GRAD_SSH_ALL
    MODULE PROCEDURE GRAD_SSH_ALL_DP
  END INTERFACE GRAD_SSH_ALL

  INTERFACE L_SSH_ALL
    MODULE PROCEDURE L_SSH_ALL_DP
  END INTERFACE L_SSH_ALL

  ! -- Poloidal/toroidal VSH (single-point) --
  INTERFACE PVSH_RAD
    MODULE PROCEDURE PVSH_RAD_DP
  END INTERFACE PVSH_RAD

  INTERFACE PVSH_TOR
    MODULE PROCEDURE PVSH_TOR_DP
  END INTERFACE PVSH_TOR

  INTERFACE PVSH_POL
    MODULE PROCEDURE PVSH_POL_DP
  END INTERFACE PVSH_POL

  ! -- Poloidal/toroidal VSH (batch) --
  INTERFACE PVSH_RAD_ALL
    MODULE PROCEDURE PVSH_RAD_ALL_DP
  END INTERFACE PVSH_RAD_ALL

  INTERFACE PVSH_TOR_ALL
    MODULE PROCEDURE PVSH_TOR_ALL_DP
  END INTERFACE PVSH_TOR_ALL

  INTERFACE PVSH_POL_ALL
    MODULE PROCEDURE PVSH_POL_ALL_DP
  END INTERFACE PVSH_POL_ALL

  ! -- VSH (single-point) --
  INTERFACE VSH_TOR
    MODULE PROCEDURE VSH_TOR_DP
  END INTERFACE VSH_TOR

  INTERFACE VSH_POL_DN
    MODULE PROCEDURE VSH_POL_DN_DP
  END INTERFACE VSH_POL_DN

  INTERFACE VSH_POL_UP
    MODULE PROCEDURE VSH_POL_UP_DP
  END INTERFACE VSH_POL_UP

  ! -- VSH (batch) --
  INTERFACE VSH_TOR_ALL
    MODULE PROCEDURE VSH_TOR_ALL_DP
  END INTERFACE VSH_TOR_ALL

  INTERFACE VSH_POL_DN_ALL
    MODULE PROCEDURE VSH_POL_DN_ALL_DP
  END INTERFACE VSH_POL_DN_ALL

  INTERFACE VSH_POL_UP_ALL
    MODULE PROCEDURE VSH_POL_UP_ALL_DP
  END INTERFACE VSH_POL_UP_ALL

  ! -- Coupling coefficients and dot product --
  INTERFACE CGCOEFF
    MODULE PROCEDURE CGCOEFF_DP
  END INTERFACE CGCOEFF

  INTERFACE SYMBOL3J
    MODULE PROCEDURE SYMBOL3J_DP
  END INTERFACE SYMBOL3J

  INTERFACE DOT
    MODULE PROCEDURE DOT_DP
  END INTERFACE DOT

  INTERFACE GWI
    MODULE PROCEDURE GWI_DP
  END INTERFACE GWI

  INTERFACE GWJ
    MODULE PROCEDURE GWJ_DP
  END INTERFACE GWJ

  ! -- Private generics for internal helpers --
  INTERFACE ASSOC_LEGENDRE_AND_DERIV
    MODULE PROCEDURE ASSOC_LEGENDRE_AND_DERIV_DP
  END INTERFACE ASSOC_LEGENDRE_AND_DERIV

  INTERFACE VSH_CORE
    MODULE PROCEDURE VSH_CORE_DP
  END INTERFACE VSH_CORE

  INTERFACE FACTORIAL
    MODULE PROCEDURE FACTORIAL_DP
  END INTERFACE FACTORIAL

CONTAINS
!     Compute Factorial k! for k>=0
  FUNCTION FACTORIAL_DP(K) RESULT(FACTORIAL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp) :: FACTORIAL
  INTEGER(KIND=i4) :: I
  FACTORIAL = 1.0d0
  DO I = 1,K
      FACTORIAL = I * FACTORIAL
  END DO
  RETURN
  END FUNCTION FACTORIAL_DP


!     Compute Legendre polynomial P_l(x) order l at -1<=x<=+1
  FUNCTION LEGENDRE_DP(L,X) RESULT(LEGENDRE)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: X
  INTEGER(KIND=I4), INTENT(IN) :: L
  REAL(KIND=dp) :: P0,PIM1,PIM2
  REAL(KIND=dp) :: FI,LEGENDRE
  INTEGER(KIND=I4) :: I
!     CHECK FOR VALID VALUES OF N AND X HERE BEFORE PROCEEDING C 
!     PROCEEDING WITH CALCULATIONS       
  IF (L.EQ.0) THEN         
      LEGENDRE=1.d0       
  ELSEIF (L.EQ.1) THEN         
      LEGENDRE=X       
  ELSE 
!     USE RECURRENCE RELATIONS       
      PIM1=1.d0         
      P0=X         
      DO I=2,L    
          FI=I        
          PIM2=PIM1         
          PIM1=P0         
          P0=((2.d0*I-1.d0)*X*PIM1-(I-1.d0)*PIM2)/FI    
      ENDDO
      LEGENDRE=P0
  ENDIF       
  RETURN
  END FUNCTION LEGENDRE_DP

!     Compute d/dx of Legendre polynomial P_l(x) order l at -1<=x<=+1
  FUNCTION DDX_LEGENDRE_DP(L,X) RESULT(DDX_LEGENDRE)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: X
  INTEGER(KIND=i4), INTENT(IN) :: L
  REAL(KIND=dp) :: DDX_LEGENDRE
  IF (L.EQ.0) THEN
      DDX_LEGENDRE=0.d0
  ELSEIF (L.EQ.1) THEN
      DDX_LEGENDRE=1.d0
  ELSEIF (L.GT.1 .AND. ABS(X).LT.1.d0) THEN
      DDX_LEGENDRE=L*(X*LEGENDRE(L,X)- &
      LEGENDRE(L-1,X))/(X**2-1.d0)
  ELSEIF (L.GT.1 .AND. X.EQ.1.d0) THEN
      DDX_LEGENDRE=L*(L+1)/2.d0
  ELSEIF (L.GT.1 .AND. X.EQ.-1.d0) THEN
      DDX_LEGENDRE=(L*(L+1)/2.d0)*(-1)**L
  ENDIF
  RETURN
  END FUNCTION DDX_LEGENDRE_DP

  FUNCTION ASSOC_LEGENDRE_DP(L, K, X) RESULT(ASSOC_LEGENDRE)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L, K
  REAL(KIND=dp), INTENT(IN) :: X
  INTEGER(KIND=i4) :: ABSK, I
  REAL(KIND=dp) :: ASSOC_LEGENDRE
  REAL(KIND=dp) :: Pmm, Pmm1, Plm, somx2

  ABSK = ABS(K)

  IF (ABSK > L .OR. ABS(X) > 1.0_dp) THEN
    ASSOC_LEGENDRE = 0.0_dp
    RETURN
  END IF

  Pmm = 1.0_dp
  IF (ABSK > 0) THEN
    somx2 = sqrt((1.0_dp - X) * (1.0_dp + X))
    DO I = 1, ABSK
      Pmm = -Pmm * (2*I - 1) * somx2
    END DO
  END IF

  IF (L == ABSK) THEN
    Plm = Pmm
  ELSE
    Pmm1 = X * (2 * ABSK + 1) * Pmm
    IF (L == ABSK + 1) THEN
      Plm = Pmm1
    ELSE
      DO I = ABSK + 2, L
        Plm = (X*(2*I-1)*Pmm1-(I+ABSK-1)*Pmm) / (I-ABSK)
        Pmm = Pmm1
        Pmm1 = Plm
      END DO
    END IF
  END IF

  IF (K < 0) THEN
    Plm = ((-1.0_dp)**ABSK)* &
      EXP(LOG_FACT(L-ABSK)-LOG_FACT(L+ABSK))*Plm
  END IF

  ASSOC_LEGENDRE = Plm
  END FUNCTION ASSOC_LEGENDRE_DP

!     Compute P_l^k(x) and d/dx P_l^k(x) in a single upward recurrence.
!     Eliminates the two separate ASSOC_LEGENDRE calls in DDX_ASSOC_LEGENDRE
!     and the redundant SSH + DDX_ASSOC_LEGENDRE pair in GRAD_SSH / L_SSH.
!     Returns PLK=0, DPLK=0 for |k|>l or |x|>1; DPLK=0 at the poles (|x|=1).
  SUBROUTINE ASSOC_LEGENDRE_AND_DERIV_DP(L, K, X, PLK, DPLK)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: L, K
  REAL(KIND=dp),    INTENT(IN)  :: X
  REAL(KIND=dp),    INTENT(OUT) :: PLK, DPLK
  INTEGER(KIND=i4) :: ABSK, I
  REAL(KIND=dp)    :: Pmm, Pmm1, Plm, Plm1, somx2, DOM, SIGN_K
  ABSK = ABS(K)
  IF (ABSK .GT. L .OR. ABS(X) .GT. 1.0_dp) THEN
    PLK  = 0.0_dp
    DPLK = 0.0_dp
    RETURN
  END IF
  DOM  = 1.0_dp - X**2
  Pmm  = 1.0_dp
  Plm1 = 0.0_dp
  IF (ABSK .GT. 0) THEN
    somx2 = SQRT((1.0_dp - X) * (1.0_dp + X))
    DO I = 1, ABSK
      Pmm = -Pmm * (2*I - 1) * somx2
    END DO
  END IF
  IF (L .EQ. ABSK) THEN
    Plm = Pmm
  ELSE
    Pmm1 = X * (2 * ABSK + 1) * Pmm
    IF (L .EQ. ABSK + 1) THEN
      Plm  = Pmm1
      Plm1 = Pmm
    ELSE
      DO I = ABSK + 2, L
        Plm  = (X*(2*I-1)*Pmm1-(I+ABSK-1)*Pmm) / (I-ABSK)
        Pmm  = Pmm1
        Pmm1 = Plm
      END DO
      Plm1 = Pmm
    END IF
  END IF
  IF (K .LT. 0) THEN
    SIGN_K = (-1.0_dp)**ABSK
    Plm  = SIGN_K*EXP(LOG_FACT(L-ABSK)-LOG_FACT(L+ABSK))*Plm
    IF (L .GT. ABSK) &
      Plm1 = SIGN_K* &
             EXP(LOG_FACT(L-1-ABSK)-LOG_FACT(L-1+ABSK))*Plm1
  END IF
  PLK = Plm
  IF (L .EQ. 0 .OR. ABS(X) .GE. 1.0_dp) THEN
    DPLK = 0.0_dp
  ELSE
    DPLK = ((L + K) * Plm1 - L * X * Plm) / DOM
  END IF
  END SUBROUTINE ASSOC_LEGENDRE_AND_DERIV_DP

  FUNCTION DDX_ASSOC_LEGENDRE_DP(L, K, X) RESULT(RES)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L, K
  REAL(KIND=dp), INTENT(IN)    :: X
  REAL(KIND=dp)                :: RES, PLK_UNUSED
  CALL ASSOC_LEGENDRE_AND_DERIV(L, K, X, PLK_UNUSED, RES)
  END FUNCTION DDX_ASSOC_LEGENDRE_DP

!     Compute scalar spherical harmonic Y_l^k(theta,phi) at 0<=theta<=pi, 0<=phi<=2pi
!     Ylk = sqrt{(2l+1)/4pi}sqrt{(l-k)!/(l+k)!}P_l^k*exp(i k phi)
  FUNCTION SSH_DP(L,K,THETA,PHI) RESULT(SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: NORM
  COMPLEX(KIND=dp) :: SSH
  NORM = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
      EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
  SSH = NORM*ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)
  RETURN
  END FUNCTION SSH_DP

!     Compute angular gradient of scalar spherical harmonic Ylk at 0<=theta<=pi, 0<=phi<=2pi
  FUNCTION GRAD_SSH_DP(L,K,THETA,PHI) RESULT(GRAD_SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK
  COMPLEX(KIND=dp) :: EPHIM, YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: GRAD_SSH
  SINTH = DSIN(THETA)
  NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
      EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
  CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
  EPHIM = EXP(j*K*PHI)
  YLM   = NORM * PLK * EPHIM
  GRAD_SSH(1) = CMPLX(0.d0, 0.d0)
  GRAD_SSH(2) = -SINTH * NORM * DPLK * EPHIM
  IF (SINTH .NE. 0.d0) THEN
    GRAD_SSH(3) = j*K*YLM/SINTH
  ELSE
    GRAD_SSH(3) = DCMPLX(0.d0, 0.d0)
  END IF
  RETURN
  END FUNCTION GRAD_SSH_DP


!     Compute rhat X angular gradient of scalar spherical harmonic Ylk at 0<=theta<=pi, 0<=phi<=2pi
  FUNCTION L_SSH_DP(L,K,THETA,PHI) RESULT(L_SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK
  COMPLEX(KIND=dp) :: EPHIM, YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: L_SSH
  SINTH = DSIN(THETA)
  NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
      EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
  CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
  EPHIM = EXP(j*K*PHI)
  YLM   = NORM * PLK * EPHIM
  L_SSH(1) = DCMPLX(0.d0, 0.d0)
  IF (SINTH .NE. 0.d0) THEN
    L_SSH(2) = -j*K*YLM/SINTH
  ELSE
    L_SSH(2) = DCMPLX(0.d0, 0.d0)
  END IF
  L_SSH(3) = -SINTH * NORM * DPLK * EPHIM
  RETURN
  END FUNCTION L_SSH_DP


!     Compute longitudinal/radial component of polar vector spherical harmonic Y_(L,K)^(-1) = rhat.Ylm
  FUNCTION PVSH_RAD_DP(L,K,THETA,PHI) RESULT(PVSH_RAD)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  COMPLEX(KIND=dp) :: YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_RAD
  YLM = SSH(L,K,THETA,PHI)
  PVSH_RAD(1) = YLM
  PVSH_RAD(2) = DCMPLX(0.d0, 0.d0)
  PVSH_RAD(3) = DCMPLX(0.d0, 0.d0)
  RETURN
  END FUNCTION PVSH_RAD_DP


!     Compute transverse/toroidal component of polar vector spherical harmonic Y_(L,K)^(0) ~ L(SSH)
  FUNCTION PVSH_TOR_DP(L,K,THETA,PHI) RESULT(PVSH_TOR)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK,SCALE
  COMPLEX(KIND=dp) :: EPHIM,YLM,GTH,GPH
  COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_TOR
  PVSH_TOR(1) = DCMPLX(0.d0, 0.d0)
  IF (L .GT. 0) THEN
    SINTH = DSIN(THETA)
    NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
            EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
    CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
    EPHIM = EXP(j*K*PHI)
    YLM   = NORM * PLK * EPHIM
    GTH   = -SINTH * NORM * DPLK * EPHIM
    IF (SINTH .NE. 0.d0) THEN
      GPH = j*K*YLM/SINTH
    ELSE
      GPH = DCMPLX(0.d0, 0.d0)
    END IF
    SCALE = 1.d0/DSQRT(L*(L+1.d0))
    PVSH_TOR(2) =  j*GPH*SCALE
    PVSH_TOR(3) = -j*GTH*SCALE
  ELSE
    PVSH_TOR(2) = DCMPLX(0.d0, 0.d0)
    PVSH_TOR(3) = DCMPLX(0.d0, 0.d0)
  END IF
  RETURN
  END FUNCTION PVSH_TOR_DP


!     Compute transverse/polar component of polar vector spherical harmonic Y_(L,K)^(+1) ~ GRADSSH
  FUNCTION PVSH_POL_DP(L,K,THETA,PHI) RESULT(PVSH_POL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK,SCALE
  COMPLEX(KIND=dp) :: EPHIM,YLM,GTH,GPH
  COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_POL
  PVSH_POL(1) = DCMPLX(0.d0, 0.d0)
  IF (L .GT. 0) THEN
    SINTH = DSIN(THETA)
    NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
            EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
    CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
    EPHIM = EXP(j*K*PHI)
    YLM   = NORM * PLK * EPHIM
    GTH   = -SINTH * NORM * DPLK * EPHIM
    IF (SINTH .NE. 0.d0) THEN
      GPH = j*K*YLM/SINTH
    ELSE
      GPH = DCMPLX(0.d0, 0.d0)
    END IF
    SCALE = 1.d0/DSQRT(L*(L+1.d0))
    PVSH_POL(2) = GTH*SCALE
    PVSH_POL(3) = GPH*SCALE
  ELSE
    PVSH_POL(2) = DCMPLX(0.d0, 0.d0)
    PVSH_POL(3) = DCMPLX(0.d0, 0.d0)
  END IF
  RETURN
  END FUNCTION PVSH_POL_DP


!     Compute transverse/toroidal component of standard vector spherical harmonic Y_(L,K)^L
  FUNCTION VSH_TOR_DP(L,K,THETA,PHI) RESULT(VSH_TOR)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  COMPLEX(KIND=dp) :: YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: VSH_TOR
  VSH_TOR = PVSH_TOR(L,K,THETA,PHI)
  RETURN
  END FUNCTION VSH_TOR_DP


!     Compute transverse component of standard vector spherical harmonic Y_(L,K)^L-1
  FUNCTION VSH_POL_DN_DP(L,K,THETA,PHI) RESULT(VSH_POL_DN)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK,SC1,SC23
  COMPLEX(KIND=dp) :: EPHIM,YLM,GTH,GPH
  COMPLEX(KIND=dp), DIMENSION(3) :: VSH_POL_DN
  VSH_POL_DN(1) = DCMPLX(0.d0, 0.d0)
  VSH_POL_DN(2) = DCMPLX(0.d0, 0.d0)
  VSH_POL_DN(3) = DCMPLX(0.d0, 0.d0)
  IF (L .GT. 0) THEN
    SINTH = DSIN(THETA)
    NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
            EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
    CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
    EPHIM = EXP(j*K*PHI)
    YLM   = NORM * PLK * EPHIM
    GTH   = -SINTH * NORM * DPLK * EPHIM
    IF (SINTH .NE. 0.d0) THEN
      GPH = j*K*YLM/SINTH
    ELSE
      GPH = DCMPLX(0.d0, 0.d0)
    END IF
    SC1  = DSQRT(DBLE(L)/(2*L+1.d0))
    SC23 = 1.d0/DSQRT((2*L+1.d0)*DBLE(L))
    VSH_POL_DN(1) = SC1  * YLM
    VSH_POL_DN(2) = SC23 * GTH
    VSH_POL_DN(3) = SC23 * GPH
  END IF
  RETURN
  END FUNCTION VSH_POL_DN_DP


!     Compute transverse component of polar vector spherical harmonic Y_(L,K)^L+1
  FUNCTION VSH_POL_UP_DP(L,K,THETA,PHI) RESULT(VSH_POL_UP)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L,K
  REAL(KIND=dp), INTENT(IN) :: THETA,PHI
  REAL(KIND=dp) :: SINTH,NORM,PLK,DPLK,SC1,SC23
  COMPLEX(KIND=dp) :: EPHIM,YLM,GTH,GPH
  COMPLEX(KIND=dp), DIMENSION(3) :: VSH_POL_UP
  SINTH = DSIN(THETA)
  NORM  = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
          EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
  CALL ASSOC_LEGENDRE_AND_DERIV(L, K, DCOS(THETA), PLK, DPLK)
  EPHIM = EXP(j*K*PHI)
  YLM   = NORM * PLK * EPHIM
  SC1   = -DSQRT((L+1.d0)/(2*L+1.d0))
  VSH_POL_UP(1) = SC1 * YLM
  IF (L .GT. 0) THEN
    GTH = -SINTH * NORM * DPLK * EPHIM
    IF (SINTH .NE. 0.d0) THEN
      GPH = j*K*YLM/SINTH
    ELSE
      GPH = DCMPLX(0.d0, 0.d0)
    END IF
    SC23 = 1.d0/DSQRT((2*L+1.d0)*(L+1.d0))
    VSH_POL_UP(2) = SC23 * GTH
    VSH_POL_UP(3) = SC23 * GPH
  ELSE
    VSH_POL_UP(2) = DCMPLX(0.d0, 0.d0)
    VSH_POL_UP(3) = DCMPLX(0.d0, 0.d0)
  END IF
  RETURN
  END FUNCTION VSH_POL_UP_DP

!     Compute Clebsch-Gordan coefficient
!     Adapted from David Simpson (NASA GSFC)
  FUNCTION CGCOEFF_DP(J1,M1,J2,M2,J3,M3) RESULT(CGCOEFF)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1,M1,J2,M2,J3,M3
  INTEGER(KIND=i4) :: K,KMIN,KMAX
  REAL(KIND=dp) :: SUMK,TERM,CGCOEFF

  IF ((J3 .LT. ABS(J1-J2)) .OR. &
      (J3 .GT. (J1+J2))    .OR. &
      (ABS(M1) .GT. J1)    .OR. &
      (ABS(M2) .GT. J2)    .OR. &
      (ABS(M3) .GT. J3)  .OR. &
      ((M1+M2) .NE. M3))  THEN
      CGCOEFF = 0.0d0
  ELSE
      CGCOEFF = DSQRT(2.0_dp*J3+1.0_dp)* &
       EXP(-0.5_dp*LOG_FACT(J1+J2+J3+1))
      CGCOEFF = CGCOEFF * EXP(0.5_dp*(LOG_FACT(J1+J2-J3)+ &
                LOG_FACT(J2+J3-J1)+LOG_FACT(J3+J1-J2)))

      CGCOEFF = CGCOEFF * EXP(0.5_dp*(LOG_FACT(J1+M1)+ &
                LOG_FACT(J1-M1)+LOG_FACT(J2+M2)+LOG_FACT(J2-M2)+ &
                LOG_FACT(J3+M3)+LOG_FACT(J3-M3)))

      SUMK = 0.0_dp
      KMIN = MAX(0, J1-J3+M2, J2-J3-M1)
      KMAX = MIN(J1+J2-J3, J1-M1, J2+M2)
      IF (KMIN .LE. KMAX) THEN
        TERM = LOG_FACT(J1+J2-J3-KMIN)+LOG_FACT(J3-J1-M2+KMIN)+ &
               LOG_FACT(J3-J2+M1+KMIN)+LOG_FACT(J1-M1-KMIN)+ &
               LOG_FACT(J2+M2-KMIN)+LOG_FACT(KMIN)
        DO K = KMIN, KMAX
          IF (MOD(K,2) == 1) THEN
            SUMK = SUMK - EXP(-TERM)
          ELSE
            SUMK = SUMK + EXP(-TERM)
          END IF
          IF (K .LT. KMAX) THEN
            TERM = TERM &
                 - LOG(DBLE(J1+J2-J3-K)) &
                 + LOG(DBLE(J3-J1-M2+K+1)) &
                 + LOG(DBLE(J3-J2+M1+K+1)) &
                 - LOG(DBLE(J1-M1-K)) &
                 - LOG(DBLE(J2+M2-K)) &
                 + LOG(DBLE(K+1))
          END IF
        END DO
      END IF
      CGCOEFF = CGCOEFF*SUMK
  ENDIF
  RETURN
  END FUNCTION CGCOEFF_DP

!     Compute log(factorial()) for stability  
  FUNCTION LOG_FACT_DP(N) RESULT(LOG_FACT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: N
  REAL(KIND=dp) :: LOG_FACT
  IF (N <= 1) THEN
    LOG_FACT = 0.0_dp
  ELSE
    ! ln(n!) = ln(Gamma(n+1))
    LOG_FACT = LOG_GAMMA(REAL(N + 1, KIND=dp))
  END IF
  END FUNCTION LOG_FACT_DP

!     Compute Gauss-Legendre quadrature nodes and weights for exact
!     integration of degree <= 2*LMAX+1 using LMAX+1 points.
!     ZERO(i): GL nodes as cos(theta), ordered from -1 to +1.
!     W(i)   : corresponding weights for integral int_{-1}^{1} f(x) dx.
!     Nodes are zeros of P_{LMAX+1}(x); weights via w=2/((1-x^2)*P'(x)^2).
!     Newton iteration converges to machine precision in <10 steps.
  SUBROUTINE SHGLQ_DP(ZERO, W, LMAX)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(OUT) :: ZERO(LMAX+1)
  REAL(KIND=dp),    INTENT(OUT) :: W(LMAX+1)
  INTEGER(KIND=i4) :: N, I, J, K, NHALF
  REAL(KIND=dp)    :: X, P0, P1, P2, DPN, DX, TOL
  TOL = 1.0e-15_dp

  N     = LMAX + 1
  NHALF = (N + 1) / 2

  DO I = 1, NHALF
    X = DCOS(pi * (I - 0.25_dp) / (N + 0.5_dp))
    DO J = 1, 100
      P0 = 1.0_dp
      P1 = X
      DO K = 2, N
        P2  = ((2*K - 1) * X * P1 - (K-1) * P0) / K
        P0  = P1
        P1  = P2
      END DO
      DPN = N * (P0 - X * P1) / (1.0_dp - X**2)
      DX  = P1 / DPN
      X   = X - DX
      IF (ABS(DX) .LE. TOL) EXIT
    END DO
    ZERO(I)         = -X
    ZERO(N + 1 - I) =  X
    W(I)            = 2.0_dp / ((1.0_dp - X**2) * DPN**2)
    W(N + 1 - I)    = W(I)
  END DO
  END SUBROUTINE SHGLQ_DP

!     Compute Wigner 3j symbol
  FUNCTION SYMBOL3J_DP(J1,M1,J2,M2,J3,M3) RESULT(SYMBOL3J)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1,M1,J2,M2,J3,M3
  INTEGER(KIND=i4) :: K
  REAL(KIND=dp) :: TERM,CG,SYMBOL3J
  TERM = (-1)**(J3+M3+NINT(2.d0*J1))/DSQRT(2.d0*J3+1.d0)
  CG = CGCOEFF(J1,-M1,J2,-M2,J3,M3)
  SYMBOL3J = TERM*CG
  RETURN
  END FUNCTION SYMBOL3J_DP

!     Compute Wigner 6j symbol
!-----Placeholder for future extension

!     Compute Wigner 9j symbol
!-----Placeholder for future extension

!     Compute dot product of two VSH evaluated at theta,phi coordinates
  FUNCTION DOT_DP(VSH1, VSH2) RESULT(DOT)
  IMPLICIT NONE
  COMPLEX(KIND=dp), DIMENSION(3), INTENT(IN) :: VSH1, VSH2
  COMPLEX(KIND=dp) :: DOT
  DOT = VSH1(1)*VSH2(1) + VSH1(2)*VSH2(2) + VSH1(3)*VSH2(3)
  RETURN
  END FUNCTION DOT_DP


  FUNCTION GWI_DP(J1,M1,J2,M2,L,M) RESULT(GWI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1,M1,J2,M2,L,M
  COMPLEX(KIND=dp) :: GWI
  GWI=DSQRT((2.d0*J1+1.d0)*(2.d0*J2+1.d0)/(2.d0*L+1.d0)/4.d0/pi)* &
    CGCOEFF(J1,0,J2,0,L,0)*CGCOEFF(J1,M1,J2,M2,L,M)
  RETURN
  END FUNCTION GWI_DP


  FUNCTION GWJ_DP(J1,M1,J2,M2,L,M) RESULT(GWJ)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1,M1,J2,M2,L,M
  REAL(KIND=dp) :: TRIANGLE_ARG
  COMPLEX(KIND=dp) :: GWJ
  ! Guard the sqrt() argument below the same way CGCOEFF guards its own
  ! triangle inequality: this product is negative outside the valid
  ! (J1+1,J2,L) triangle range, which would otherwise silently produce
  ! NaN (or trap under -ffpe-trap builds) instead of the correct zero.
  TRIANGLE_ARG = (J1+J2+L+2.d0)*(J2+L-J1)* &
    (J1+J2-L+1.d0)*(J1-J2+L+1.d0)
  IF (TRIANGLE_ARG < 0.d0) THEN
    GWJ = DCMPLX(0.d0, 0.d0)
    RETURN
  END IF
  GWJ=-DSQRT((2.d0*J1+1.d0)*(2.d0*J2+1.d0)/(2.d0*L+1.d0)/4.d0/pi)* &
    CGCOEFF(J1+1,0,J2,0,L,0)*CGCOEFF(J1,M1,J2,M2,L,M)* &
    DCMPLX(0.d0, 1.d0)*DSQRT(TRIANGLE_ARG)/2.d0
  RETURN
  END FUNCTION GWJ_DP


!     Return 1D index for (L,M) in ASSOC_LEGENDRE_ALL output array
!     SHTOOLS PlmIndex convention: l*(l+1)/2 + m + 1, requires 0<=m<=l
  FUNCTION PLM_INDEX(L, M)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L, M
  INTEGER(KIND=i4) :: PLM_INDEX
  PLM_INDEX = L*(L+1)/2 + M + 1
  END FUNCTION PLM_INDEX

!     Return 1D index for (L,M) in SSH_ALL output array, -l<=m<=l
  FUNCTION YLM_INDEX(L, M)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L, M
  INTEGER(KIND=i4) :: YLM_INDEX
  YLM_INDEX = L**2 + L + M + 1
  END FUNCTION YLM_INDEX

!     Compute all P_l^m(x) for 0<=l<=lmax, 0<=m<=l using Bonnet recurrence
!     P indexed by PLM_INDEX(l,m) = l*(l+1)/2 + m + 1
!     Condon-Shortley convention; output size (lmax+1)*(lmax+2)/2
  SUBROUTINE ASSOC_LEGENDRE_ALL_DP(P, LMAX, X)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  REAL(KIND=dp), INTENT(IN) :: X
  REAL(KIND=dp), INTENT(OUT) :: P((LMAX+1)*(LMAX+2)/2)
  INTEGER(KIND=i4) :: L, M
  REAL(KIND=dp) :: SOMX2

  SOMX2 = SQRT((1.0_dp - X) * (1.0_dp + X))
  P(PLM_INDEX(0, 0)) = 1.0_dp

  DO M = 0, LMAX
    IF (M .GT. 0) &
      P(PLM_INDEX(M, M)) = -(2*M - 1) * SOMX2 * &
        P(PLM_INDEX(M-1, M-1))
      IF (M .LT. LMAX) &
        P(PLM_INDEX(M+1, M)) = X * (2*M + 1) * &
            P(PLM_INDEX(M, M))
      DO L = M+2, LMAX
        P(PLM_INDEX(L, M)) = &
            (X*(2*L-1)*P(PLM_INDEX(L-1, M)) - &
            (L+M-1)*P(PLM_INDEX(L-2, M))) / (L-M)
      END DO
  END DO
  END SUBROUTINE ASSOC_LEGENDRE_ALL_DP

!     Compute all d/dx P_l^m(x) for 0<=l<=lmax, 0<=m<=l
!     Requires precomputed P from ASSOC_LEGENDRE_ALL
!     DP_OUT indexed by PLM_INDEX(l,m); returns 0 at poles (|x|>=1)
  SUBROUTINE DDX_ASSOC_LEGENDRE_ALL_DP(DP_OUT, P, LMAX, X)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  REAL(KIND=dp), INTENT(IN) :: X
  REAL(KIND=dp), INTENT(IN) :: P((LMAX+1)*(LMAX+2)/2)
  REAL(KIND=dp), INTENT(OUT) :: DP_OUT((LMAX+1)*(LMAX+2)/2)
  INTEGER(KIND=i4) :: L, M
  REAL(KIND=dp) :: DOM

  DOM = 1.0_dp - X**2

  DO M = 0, LMAX
    DO L = M, LMAX
      IF (ABS(X) .GE. 1.0_dp .OR. L .EQ. 0) THEN
        DP_OUT(PLM_INDEX(L, M)) = 0.0_dp
      ELSE IF (L .EQ. M) THEN
        DP_OUT(PLM_INDEX(L, M)) = &
            -L * X * P(PLM_INDEX(L, M)) / DOM
      ELSE
        DP_OUT(PLM_INDEX(L, M)) = &
            ((L+M)*P(PLM_INDEX(L-1, M)) - &
            L*X*P(PLM_INDEX(L, M))) / DOM
      END IF
    END DO
  END DO
  END SUBROUTINE DDX_ASSOC_LEGENDRE_ALL_DP

!     Compute all 4pi-normalized P_l^m(cos theta) for 0<=l<=lmax, 0<=m<=l
!     using the Holmes & Featherstone (2002) modified forward-column recurrence.
!     Normalization N_lm = sqrt((2l+1)/(4pi)*(l-m)!/(l+m)!) is folded into
!     the recurrence coefficients, keeping all intermediates O(1/sqrt(4pi))
!     and extending numerical stability to l~2700 vs ~1400 for Bonnet.
!     PNORM indexed by PLM_INDEX(l,m); output size (lmax+1)*(lmax+2)/2.
  SUBROUTINE ASSOC_LEGENDRE_NORM_ALL_DP(PNORM, LMAX, X)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: X
  REAL(KIND=dp),    INTENT(OUT) :: PNORM((LMAX+1)*(LMAX+2)/2)
  INTEGER(KIND=i4) :: L, M
  REAL(KIND=dp)    :: SINTH

  SINTH = SQRT((1.0_dp - X) * (1.0_dp + X))

  PNORM(PLM_INDEX(0, 0)) = 1.0_dp / DSQRT(4.0_dp*pi)

  DO M = 0, LMAX
    IF (M .GT. 0) &
      PNORM(PLM_INDEX(M, M)) = &
        -DSQRT((2*M+1.d0)/(2*M)) * SINTH * &
        PNORM(PLM_INDEX(M-1, M-1))
    IF (M .LT. LMAX) &
      PNORM(PLM_INDEX(M+1, M)) = &
        DSQRT(2*M+3.d0) * X * PNORM(PLM_INDEX(M, M))
    DO L = M+2, LMAX
      PNORM(PLM_INDEX(L, M)) = &
        DSQRT((4.d0*L**2-1.d0)/(L**2-M**2)) * X * &
        PNORM(PLM_INDEX(L-1, M)) - &
        DSQRT((2*L+1.d0)*(L-M-1.d0)*(L+M-1.d0)/ &
              ((2*L-3.d0)*(L**2-M**2))) * &
        PNORM(PLM_INDEX(L-2, M))
    END DO
  END DO
  END SUBROUTINE ASSOC_LEGENDRE_NORM_ALL_DP


!     Compute d/dx of normalized P_l^m for 0<=l<=lmax, 0<=m<=l.
!     Requires precomputed PNORM from ASSOC_LEGENDRE_NORM_ALL.
!     Formula: (sqrt((2l+1)(l^2-m^2)/(2l-1))*PNORM(l-1,m) - l*x*PNORM(l,m))
!              / (1-x^2); returns 0 at poles (|x|>=1).
  SUBROUTINE DDX_ASSOC_LEGENDRE_NORM_ALL_DP(DPNORM, PNORM, LMAX, X)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: X
  REAL(KIND=dp),    INTENT(IN)  :: PNORM((LMAX+1)*(LMAX+2)/2)
  REAL(KIND=dp),    INTENT(OUT) :: DPNORM((LMAX+1)*(LMAX+2)/2)
  INTEGER(KIND=i4) :: L, M
  REAL(KIND=dp)    :: DOM, COEFF

  DOM = 1.0_dp - X**2

  DO M = 0, LMAX
    DO L = M, LMAX
      IF (ABS(X) .GE. 1.0_dp .OR. L .EQ. 0) THEN
        DPNORM(PLM_INDEX(L, M)) = 0.0_dp
      ELSE IF (L .EQ. M) THEN
        DPNORM(PLM_INDEX(L, M)) = &
          -L * X * PNORM(PLM_INDEX(L, M)) / DOM
      ELSE
        COEFF = DSQRT((2*L+1.d0)*(L**2-M**2)/(2*L-1.d0))
        DPNORM(PLM_INDEX(L, M)) = &
          (COEFF * PNORM(PLM_INDEX(L-1, M)) - &
           L * X * PNORM(PLM_INDEX(L, M))) / DOM
      END IF
    END DO
  END DO
  END SUBROUTINE DDX_ASSOC_LEGENDRE_NORM_ALL_DP


!     Compute all Y_l^m(theta,phi) for 0<=l<=lmax, -l<=m<=l
!     YLM indexed by YLM_INDEX(l,m) = l**2 + l + m + 1
!     Y_l^{-m} = (-1)^m * conj(Y_l^m) via conjugate symmetry
  SUBROUTINE SSH_ALL_DP(YLM, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  REAL(KIND=dp), INTENT(IN) :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: YLM((LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, PSIZE
  REAL(KIND=dp), ALLOCATABLE :: PNORM(:)
  REAL(KIND=dp) :: SIGN_M

  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(PNORM(PSIZE))
  CALL ASSOC_LEGENDRE_NORM_ALL(PNORM, LMAX, DCOS(THETA))

  DO L = 0, LMAX
    DO M = 0, L
        YLM(YLM_INDEX(L, M)) = &
            PNORM(PLM_INDEX(L, M)) * EXP(j*M*PHI)
    END DO
    DO M = 1, L
        SIGN_M = 1.0_dp - 2.0_dp * MOD(M, 2)
        YLM(YLM_INDEX(L, -M)) = SIGN_M * &
            CONJG(YLM(YLM_INDEX(L, M)))
    END DO
  END DO

  DEALLOCATE(PNORM)
  END SUBROUTINE SSH_ALL_DP

!     Helper: precompute YLM_OUT, GSSH_TH, GSSH_PH for all (l,m),
!     0<=l<=LMAX, -l<=m<=l, indexed by YLM_INDEX.
!     GSSH_TH/GPH are the theta/phi components of GRAD_SSH.
  SUBROUTINE VSH_CORE_DP(YLM_OUT, GSSH_TH, GSSH_PH, &
                      LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: YLM_OUT((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(OUT) :: GSSH_TH((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(OUT) :: GSSH_PH((LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, PSIZE
  REAL(KIND=dp), ALLOCATABLE :: PNORM(:), DPNORM(:)
  REAL(KIND=dp)    :: SINTH, SIGN_M
  COMPLEX(KIND=dp) :: EPHIM, YP, GTH_P, GPH_P
  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(PNORM(PSIZE), DPNORM(PSIZE))
  CALL ASSOC_LEGENDRE_NORM_ALL(PNORM, LMAX, DCOS(THETA))
  CALL DDX_ASSOC_LEGENDRE_NORM_ALL(DPNORM, PNORM, LMAX, &
                                  DCOS(THETA))
  SINTH = DSIN(THETA)
  DO L = 0, LMAX
    DO M = 0, L
      EPHIM = EXP(j*M*PHI)
      YP    = PNORM(PLM_INDEX(L,M)) * EPHIM
      GTH_P = -SINTH * DPNORM(PLM_INDEX(L,M)) * EPHIM
      IF (SINTH .NE. 0.d0) THEN
        GPH_P = j*M*YP/SINTH
      ELSE
        GPH_P = DCMPLX(0.d0, 0.d0)
      END IF
      YLM_OUT(YLM_INDEX(L, M)) = YP
      GSSH_TH(YLM_INDEX(L, M)) = GTH_P
      GSSH_PH(YLM_INDEX(L, M)) = GPH_P
      IF (M .GT. 0) THEN
        SIGN_M = 1.0_dp - 2.0_dp * MOD(M, 2)
        YLM_OUT(YLM_INDEX(L,-M)) = SIGN_M * CONJG(YP)
        GSSH_TH(YLM_INDEX(L,-M)) = SIGN_M * CONJG(GTH_P)
        GSSH_PH(YLM_INDEX(L,-M)) = SIGN_M * CONJG(GPH_P)
      END IF
    END DO
  END DO
  DEALLOCATE(PNORM, DPNORM)
  END SUBROUTINE VSH_CORE_DP


!     Compute all GRAD_SSH for 0<=l<=LMAX, -l<=m<=l
!     OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE GRAD_SSH_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: NYLM
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT(1,:) = DCMPLX(0.d0, 0.d0)
  OUT(2,:) = GTH
  OUT(3,:) = GPH
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE GRAD_SSH_ALL_DP


!     Compute all L_SSH (rhat x grad SSH) for 0<=l<=LMAX, -l<=m<=l
!     OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE L_SSH_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: NYLM
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT(1,:) = DCMPLX(0.d0, 0.d0)
  OUT(2,:) = -GPH
  OUT(3,:) =  GTH
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE L_SSH_ALL_DP


!     Compute all PVSH_RAD (radial polar VSH) for 0<=l<=LMAX, -l<=m<=l
!     OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE PVSH_RAD_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: NYLM
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT(1,:) = YLM
  OUT(2,:) = DCMPLX(0.d0, 0.d0)
  OUT(3,:) = DCMPLX(0.d0, 0.d0)
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE PVSH_RAD_ALL_DP


!     Compute all PVSH_POL (poloidal polar VSH) for 0<=l<=LMAX, -l<=m<=l
!     Zero for l=0. OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE PVSH_POL_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, NYLM, IDX
  REAL(KIND=dp)    :: SCALE
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT = DCMPLX(0.d0, 0.d0)
  DO L = 1, LMAX
    SCALE = 1.d0/DSQRT(L*(L+1.d0))
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      OUT(2,IDX) = GTH(IDX) * SCALE
      OUT(3,IDX) = GPH(IDX) * SCALE
    END DO
  END DO
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE PVSH_POL_ALL_DP


!     Compute all PVSH_TOR (toroidal polar VSH) for 0<=l<=LMAX, -l<=m<=l
!     Zero for l=0. OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE PVSH_TOR_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, NYLM, IDX
  REAL(KIND=dp)    :: SCALE
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT = DCMPLX(0.d0, 0.d0)
  DO L = 1, LMAX
    SCALE = 1.d0/DSQRT(L*(L+1.d0))
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      OUT(2,IDX) =  j*GPH(IDX)*SCALE
      OUT(3,IDX) = -j*GTH(IDX)*SCALE
    END DO
  END DO
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE PVSH_TOR_ALL_DP


!     Compute all VSH_TOR for 0<=l<=LMAX, -l<=m<=l (identical to PVSH_TOR)
!     OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE VSH_TOR_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  CALL PVSH_TOR_ALL(OUT, LMAX, THETA, PHI)
  END SUBROUTINE VSH_TOR_ALL_DP


!     Compute all VSH_POL_UP (J=L+1 polar VSH) for 0<=l<=LMAX, -l<=m<=l
!     sqrt(L/(2L+1))*PVSH_POL - sqrt((L+1)/(2L+1))*PVSH_RAD
!     OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE VSH_POL_UP_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, NYLM, IDX
  REAL(KIND=dp)    :: SC_TH, SC_R
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT = DCMPLX(0.d0, 0.d0)
  OUT(1, YLM_INDEX(0,0)) = -YLM(YLM_INDEX(0,0))
  DO L = 1, LMAX
    SC_TH = 1.d0/DSQRT((2*L+1.d0)*(L+1.d0))
    SC_R  = -DSQRT((L+1.d0)/(2*L+1.d0))
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      OUT(1,IDX) = SC_R  * YLM(IDX)
      OUT(2,IDX) = SC_TH * GTH(IDX)
      OUT(3,IDX) = SC_TH * GPH(IDX)
    END DO
  END DO
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE VSH_POL_UP_ALL_DP


!     Compute all VSH_POL_DN (J=L-1 polar VSH) for 0<=l<=LMAX, -l<=m<=l
!     sqrt((L+1)/(2L+1))*PVSH_POL + sqrt(L/(2L+1))*PVSH_RAD
!     Zero for l=0. OUT(c, YLM_INDEX(l,m)), c=1(r), 2(theta), 3(phi)
  SUBROUTINE VSH_POL_DN_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA, PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  INTEGER(KIND=i4) :: L, M, NYLM, IDX
  REAL(KIND=dp)    :: SC_TH, SC_R
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:), GTH(:), GPH(:)
  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM), GTH(NYLM), GPH(NYLM))
  CALL VSH_CORE(YLM, GTH, GPH, LMAX, THETA, PHI)
  OUT = DCMPLX(0.d0, 0.d0)
  DO L = 1, LMAX
    SC_TH = 1.d0/DSQRT((2*L+1.d0)*DBLE(L))
    SC_R  = DSQRT(DBLE(L)/(2*L+1.d0))
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      OUT(1,IDX) = SC_R  * YLM(IDX)
      OUT(2,IDX) = SC_TH * GTH(IDX)
      OUT(3,IDX) = SC_TH * GPH(IDX)
    END DO
  END DO
  DEALLOCATE(YLM, GTH, GPH)
  END SUBROUTINE VSH_POL_DN_ALL_DP


END MODULE VSH
