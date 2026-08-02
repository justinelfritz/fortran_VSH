!> Scalar and vector spherical harmonics (VSH), Legendre functions, and
!> the angular-momentum coupling coefficients (Clebsch-Gordan, Wigner
!> 3-j, and Geppert-Wiebicke) used to build them.
!>
!> Two vector bases are provided over the same underlying scalar
!> spherical harmonic \( Y_\ell^m \): the "polar" basis
!> (`PVSH_RAD`/`PVSH_POL`/`PVSH_TOR`, indexed by \( (\ell,m) \) directly)
!> familiar from the poloidal/toroidal decomposition of divergence-free
!> vector fields in MHD, and the "standard" total-angular-momentum-coupled
!> basis (`VSH_TOR`/`VSH_POL_DN`/`VSH_POL_UP`, indexed by the coupled
!> total angular momentum \( J=\ell,\ell{-}1,\ell{+}1 \)) familiar from the
!> quantum-mechanical vector-spherical-harmonic literature (e.g. Barrera
!> et al., 1985). Every routine has a single-point form (one mode at one
!> \( (\theta,\phi) \)) and a batch `_ALL` form (every
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \) mode at one point in a
!> single call, reusing shared recurrences internally for efficiency).
!>
!> All angles are in radians, with \( \theta \) the colatitude
!> (\( 0\le\theta\le\pi \), from the north pole) and \( \phi \) the
!> longitude (\( 0\le\phi<2\pi \)). Every spherical harmonic here uses the
!> fully-normalized convention (\( \int |Y_\ell^m|^2\,d\Omega=1 \)) with
!> the Condon-Shortley phase, matching common physics/quantum-mechanics
!> conventions (e.g. Edmonds, 1960).
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

!> Factorial \( k! \) for \( k \ge 0 \), computed by direct multiplication.
!>
!> @param K Non-negative integer argument.
!> Returns: \( k! \) as a double-precision real (avoids integer overflow
!>   for moderate `K`; see [[LOG_FACT]] for large-`K` stability).
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


!> Legendre polynomial \( P_\ell(x) \), degree \( \ell \ge 0 \), via the
!> standard three-term recurrence
!> $$ (\ell{+}1)\,P_{\ell+1}(x) = (2\ell{+}1)\,x\,P_\ell(x) - \ell\,P_{\ell-1}(x). $$
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param X Argument, \( -1 \le x \le 1 \) (typically \( x=\cos\theta \)).
!> Returns: \( P_\ell(x) \).
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

!> Derivative \( dP_\ell/dx \) of the Legendre polynomial, via
!> \( (x^2{-}1)\,P_\ell'(x) = \ell\,[x P_\ell(x) - P_{\ell-1}(x)] \) away
!> from the poles, and the closed-form pole values
!> \( P_\ell'(\pm1) = (\pm1)^{\ell+1}\,\ell(\ell{+}1)/2 \) at \( x=\pm1 \).
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param X Argument, \( -1 \le x \le 1 \).
!> Returns: \( dP_\ell/dx \) evaluated at `X`.
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

!> Associated Legendre function \( P_\ell^k(x) \), Condon-Shortley phase
!> convention, computed by the classic Bonnet-type upward recurrence
!> starting from \( P_m^m(x) = (-1)^m(2m{-}1)!!\,(1{-}x^2)^{m/2} \). For
!> negative order, \( P_\ell^{-k}(x) = (-1)^k\,\frac{(\ell-k)!}{(\ell+k)!}
!> P_\ell^k(x) \). Unnormalized (unlike [[ASSOC_LEGENDRE_NORM_ALL]]); loses
!> numerical stability well before the intermediate \( P_m^m \) values
!> literally overflow, since the later upward-in-\(\ell\) recurrence
!> involves a subtraction of two comparably huge terms (catastrophic
!> cancellation).
!>
!> @warning **Numerical stability**: verified accurate (relative error
!>   `< 3e-14`) for \( \ell\le150 \) at any \( m,\theta \) -- the worst
!>   case, \( \ell=151 \), occurs near the equator (\( \theta\approx90 \)
!>   degrees) at \( m=\ell \); safety improves toward either pole. This
!>   supersedes an earlier "\( \ell\sim1400 \)" claim, verified too
!>   optimistic by direct comparison against an arbitrary-precision
!>   reference. Full methodology and per-degree table:
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2 (regenerable via the
!>   opt-in `VSH_BUILD_STABILITY` CMake target's `vsh_unnorm_stability`
!>   executable).
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param X Argument, \( -1 \le x \le 1 \) (typically \( x=\cos\theta \)).
!> Returns: \( P_\ell^k(x) \).
  FUNCTION ASSOC_LEGENDRE_DP(L, K, X) RESULT(ASSOC_LEGENDRE)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
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

!> Private helper: \( P_\ell^k(x) \) and \( dP_\ell^k/dx \) together, from a
!> single upward recurrence pass. Avoids the redundant separate
!> [[ASSOC_LEGENDRE]] calls that [[DDX_ASSOC_LEGENDRE]] would otherwise
!> need, and the repeated [[SSH]] + [[DDX_ASSOC_LEGENDRE]] evaluation that
!> [[GRAD_SSH]] and [[L_SSH]] both build on.
!>
!> @warning Inherits [[ASSOC_LEGENDRE]]'s stability limit (safe for
!>   \( \ell\le150 \) at any \( m,\theta \)) -- see `STABILITY_FINDINGS.md`
!>   at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param X Argument, \( -1 \le x \le 1 \).
!> @param PLK Output \( P_\ell^k(x) \); 0 for \( |k|>\ell \) or \( |x|>1 \).
!> @param DPLK Output \( dP_\ell^k/dx \); 0 at the poles (\( |x|=1 \)).
  SUBROUTINE ASSOC_LEGENDRE_AND_DERIV_DP(L, K, X, PLK, DPLK)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: L
  INTEGER(KIND=i4), INTENT(IN)  :: K
  REAL(KIND=dp),    INTENT(IN)  :: X
  REAL(KIND=dp),    INTENT(OUT) :: PLK
  REAL(KIND=dp),    INTENT(OUT) :: DPLK
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

!> Derivative \( dP_\ell^k/dx \) of the associated Legendre function,
!> via [[ASSOC_LEGENDRE_AND_DERIV]] (thin wrapper that discards \(
!> P_\ell^k(x) \) itself).
!>
!> @warning Inherits [[ASSOC_LEGENDRE]]'s stability limit (safe for
!>   \( \ell\le150 \) at any \( m,\theta \)) -- see `STABILITY_FINDINGS.md`
!>   at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param X Argument, \( -1 \le x \le 1 \).
!> Returns: \( dP_\ell^k/dx \) evaluated at `X`; 0 at the poles (\( |x|=1
!>   \)) or for \( |k|>\ell \).
  FUNCTION DDX_ASSOC_LEGENDRE_DP(L, K, X) RESULT(RES)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN)    :: X
  REAL(KIND=dp)                :: RES, PLK_UNUSED
  CALL ASSOC_LEGENDRE_AND_DERIV(L, K, X, PLK_UNUSED, RES)
  END FUNCTION DDX_ASSOC_LEGENDRE_DP

!> Scalar spherical harmonic \( Y_\ell^k(\theta,\phi) \), fully normalized
!> with the Condon-Shortley phase:
!> $$ Y_\ell^k(\theta,\phi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-k)!}{(\ell+k)!}}\,
!>    P_\ell^k(\cos\theta)\,e^{ik\phi}. $$
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE]] -- inherits its
!>   stability limit (safe for \( \ell\le150 \) at any \( m,\theta \)),
!>   *not* [[SSH_ALL]]'s \( \ell\le2000 \) (that routine uses the
!>   normalized recurrence instead). See `STABILITY_FINDINGS.md` at the
!>   repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( Y_\ell^k(\theta,\phi) \), complex.
  FUNCTION SSH_DP(L,K,THETA,PHI) RESULT(SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
  REAL(KIND=dp) :: NORM
  COMPLEX(KIND=dp) :: SSH
  NORM = DSQRT((2.d0*L+1.d0)/(4.d0*pi))* &
      EXP(0.5_dp*(LOG_FACT(L-K)-LOG_FACT(L+K)))
  SSH = NORM*ASSOC_LEGENDRE(L,K,DCOS(THETA))*EXP(j*K*PHI)
  RETURN
  END FUNCTION SSH_DP

!> Angular (surface) gradient of the scalar spherical harmonic,
!> \( \nabla_\perp Y_\ell^k = \hat\theta\,\partial_\theta Y_\ell^k +
!> \hat\phi\,\frac{1}{\sin\theta}\partial_\phi Y_\ell^k \). Returned as a
!> 3-vector \( (\hat r,\hat\theta,\hat\phi) \) with a zero radial
!> component, since the gradient of a purely angular function is purely
!> tangential.
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[GRAD_SSH_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( (0,\,\partial_\theta Y_\ell^k,\,\frac{1}{\sin\theta}
!>   \partial_\phi Y_\ell^k) \), complex 3-vector.
  FUNCTION GRAD_SSH_DP(L,K,THETA,PHI) RESULT(GRAD_SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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


!> The angular-momentum-operator field \( \hat r \times \nabla_\perp
!> Y_\ell^k \) (proportional to the quantum-mechanical orbital
!> angular-momentum operator \( \mathbf{L}\,Y_\ell^k \) acting on the
!> scalar harmonic). Simply [[GRAD_SSH]] rotated 90 degrees within the
!> tangent plane: \( \hat r\times(\hat\theta\,G_\theta+\hat\phi\,G_\phi)
!> = \hat\theta\,(-G_\phi) + \hat\phi\,G_\theta \).
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[L_SSH_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( (0,\,-\frac{1}{\sin\theta}\partial_\phi Y_\ell^k,\,
!>   \partial_\theta Y_\ell^k) \), complex 3-vector.
  FUNCTION L_SSH_DP(L,K,THETA,PHI) RESULT(L_SSH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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


!> Radial member of the polar vector spherical harmonic basis,
!> \( \mathbf{Y}_{\ell m}^{(-1)} = \hat r\,Y_\ell^m(\theta,\phi) \). Purely
!> radial by construction (zero horizontal components) -- combined with a
!> radial stream function \( f(r) \), this is the \( B_r \)-generating
!> piece of a poloidal magnetic (or any divergence-free) field,
!> \( B_r = \frac{\ell(\ell+1)}{r^2}\,f(r)\,Y_\ell^m \).
!>
!> @warning Built on [[SSH]] -- inherits the unnormalized
!>   [[ASSOC_LEGENDRE]]'s stability limit (safe for \( \ell\le150 \) at
!>   any \( m,\theta \)), *not* [[PVSH_RAD_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( (Y_\ell^m,\,0,\,0) \), complex 3-vector.
  FUNCTION PVSH_RAD_DP(L,K,THETA,PHI) RESULT(PVSH_RAD)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
  COMPLEX(KIND=dp) :: YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: PVSH_RAD
  YLM = SSH(L,K,THETA,PHI)
  PVSH_RAD(1) = YLM
  PVSH_RAD(2) = DCMPLX(0.d0, 0.d0)
  PVSH_RAD(3) = DCMPLX(0.d0, 0.d0)
  RETURN
  END FUNCTION PVSH_RAD_DP


!> Toroidal member of the polar vector spherical harmonic basis,
!> \( \mathbf{Y}_{\ell m}^{(0)} = \frac{1}{\sqrt{\ell(\ell+1)}}\,
!> \hat r\times\nabla_\perp Y_\ell^m \) -- purely horizontal (zero radial
!> component), identical in form to [[L_SSH]] but normalized by
!> \( 1/\sqrt{\ell(\ell+1)} \). This is the horizontal-field piece of a
!> toroidal magnetic field, \( \mathbf{B}_{tor} = T(r)\,
!> \mathbf{Y}_{\ell m}^{(0)} \) for a toroidal stream function
!> \( T(r) \). Zero for \( \ell=0 \) (no toroidal monopole).
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[PVSH_TOR_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( (0,\,\frac{i}{\sqrt{\ell(\ell+1)}}\frac{1}{\sin\theta}
!>   \partial_\phi Y_\ell^m,\,-\frac{i}{\sqrt{\ell(\ell+1)}}\partial_\theta
!>   Y_\ell^m) \), complex 3-vector.
  FUNCTION PVSH_TOR_DP(L,K,THETA,PHI) RESULT(PVSH_TOR)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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


!> Poloidal (horizontal) member of the polar vector spherical harmonic
!> basis, \( \mathbf{Y}_{\ell m}^{(+1)} = \frac{1}{\sqrt{\ell(\ell+1)}}\,
!> \nabla_\perp Y_\ell^m \) -- purely horizontal (zero radial component),
!> i.e. [[GRAD_SSH]] normalized by \( 1/\sqrt{\ell(\ell+1)} \). Combined
!> with a radial stream function \( f(r) \), this is the horizontal-field
!> piece of a poloidal field,
!> \( \mathbf{B}_{\theta,\phi} = \frac{1}{r}\frac{d(rf)}{dr}\,
!> \sqrt{\ell(\ell+1)}\,\mathbf{Y}_{\ell m}^{(+1)} \). Zero for
!> \( \ell=0 \) (a constant has no horizontal gradient).
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[PVSH_POL_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( (0,\,\frac{1}{\sqrt{\ell(\ell+1)}}\partial_\theta Y_\ell^m,\,
!>   \frac{1}{\sqrt{\ell(\ell+1)}}\frac{1}{\sin\theta}\partial_\phi
!>   Y_\ell^m) \), complex 3-vector.
  FUNCTION PVSH_POL_DP(L,K,THETA,PHI) RESULT(PVSH_POL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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


!> Toroidal member of the "standard" (total-angular-momentum-coupled)
!> vector spherical harmonic basis, \( \mathbf{Y}_{\ell m}^{J=\ell} \).
!> Identical to [[PVSH_TOR]] -- the toroidal member is common to both the
!> polar and standard bases, since \( \hat r\times\nabla_\perp Y_\ell^m \)
!> already carries total angular momentum \( J=\ell \) with no
!> \( J=\ell\pm1 \) admixture.
!>
!> @warning Identical to [[PVSH_TOR]] -- inherits the unnormalized
!>   [[ASSOC_LEGENDRE_AND_DERIV]]'s stability limit (safe for
!>   \( \ell\le150 \) at any \( m,\theta \)), *not* [[VSH_TOR_ALL]]'s
!>   \( \ell\le2000 \). See `STABILITY_FINDINGS.md` at the repo root,
!>   Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( \mathbf{Y}_{\ell m}^{J=\ell}(\theta,\phi) \), complex
!>   3-vector (identical to [[PVSH_TOR]]).
  FUNCTION VSH_TOR_DP(L,K,THETA,PHI) RESULT(VSH_TOR)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
  COMPLEX(KIND=dp) :: YLM
  COMPLEX(KIND=dp), DIMENSION(3) :: VSH_TOR
  VSH_TOR = PVSH_TOR(L,K,THETA,PHI)
  RETURN
  END FUNCTION VSH_TOR_DP


!> The \( J=\ell{-}1 \) member of the standard (total-angular-momentum
!> coupled) vector spherical harmonic basis, built from the polar basis by
!> $$ \mathbf{Y}_{\ell m}^{\ell-1} = \sqrt{\tfrac{\ell}{2\ell+1}}\,
!>    \mathbf{Y}_{\ell m}^{(-1)} + \sqrt{\tfrac{\ell+1}{2\ell+1}}\,
!>    \mathbf{Y}_{\ell m}^{(+1)}, $$
!> i.e. [[PVSH_RAD]] and [[PVSH_POL]] rotated together by an
!> \( \ell \)-dependent angle (see [[VSH_POL_UP]] for the companion
!> \( J=\ell{+}1 \) rotation). Zero for \( \ell=0 \) (no \( J=-1 \)).
!>
!> @warning At \( \ell=1 \) this does *not* reproduce the classic 2:1
!>   magnetic-dipole radial:horizontal ratio -- `VSH_POL_DN(1,0,...)` is
!>   proportional to \( (\cos\theta,\,-\sin\theta,\,0) \), a 1:-1 ratio.
!>   [[VSH_POL_UP]] is the one that reproduces the dipole pattern.
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[VSH_POL_DN_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( \mathbf{Y}_{\ell m}^{\ell-1}(\theta,\phi) \), complex
!>   3-vector.
  FUNCTION VSH_POL_DN_DP(L,K,THETA,PHI) RESULT(VSH_POL_DN)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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


!> The \( J=\ell{+}1 \) member of the standard (total-angular-momentum
!> coupled) vector spherical harmonic basis, built from the polar basis by
!> $$ \mathbf{Y}_{\ell m}^{\ell+1} = \sqrt{\tfrac{\ell+1}{2\ell+1}}\,
!>    \mathbf{Y}_{\ell m}^{(-1)} - \sqrt{\tfrac{\ell}{2\ell+1}}\,
!>    \mathbf{Y}_{\ell m}^{(+1)}, $$
!> the companion rotation to [[VSH_POL_DN]]'s \( J=\ell{-}1 \). Unlike
!> [[VSH_POL_DN]] this is well-defined (nonzero radial part) at
!> \( \ell=0 \), since \( J=1 \) still exists when \( \ell=0 \).
!>
!> @note At \( \ell=1 \), `VSH_POL_UP(1,0,...)` is proportional to
!>   \( (2\cos\theta,\,\sin\theta,\,0) \) -- the classic 2:1 magnetic-dipole
!>   radial:horizontal ratio.
!>
!> @warning Built on the unnormalized [[ASSOC_LEGENDRE_AND_DERIV]] --
!>   inherits its stability limit (safe for \( \ell\le150 \) at any
!>   \( m,\theta \)), *not* [[VSH_POL_UP_ALL]]'s \( \ell\le2000 \). See
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2.
!>
!> @param L Degree, \( \ell \ge 0 \).
!> @param K Order, \( -\ell \le k \le \ell \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: \( \mathbf{Y}_{\ell m}^{\ell+1}(\theta,\phi) \), complex
!>   3-vector.
  FUNCTION VSH_POL_UP_DP(L,K,THETA,PHI) RESULT(VSH_POL_UP)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: K
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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

!> Clebsch-Gordan coefficient \( \langle j_1 m_1 j_2 m_2 | j_3 m_3\rangle \)
!> for coupling two angular momenta \( j_1,j_2 \) to a total \( j_3 \),
!> via the closed-form Racah sum (evaluated in log-factorial form for
!> numerical stability at large arguments -- see [[LOG_FACT]]):
!> $$ \langle j_1 m_1 j_2 m_2|j_3 m_3\rangle = \delta_{m_3,m_1+m_2}
!>    \sqrt{2j_3+1}\,\Delta(j_1j_2j_3)\times $$
!> $$ \sqrt{(j_1{+}m_1)!(j_1{-}m_1)!(j_2{+}m_2)!(j_2{-}m_2)!(j_3{+}m_3)!(j_3{-}m_3)!}
!>    \sum_k \frac{(-1)^k}{k!\,(j_1{+}j_2{-}j_3{-}k)!(j_3{-}j_1{-}m_2{+}k)!
!>    (j_3{-}j_2{+}m_1{+}k)!(j_1{-}m_1{-}k)!(j_2{+}m_2{-}k)!}, $$
!> where \( \Delta(j_1j_2j_3) = \sqrt{\frac{(j_1+j_2-j_3)!(j_2+j_3-j_1)!
!> (j_3+j_1-j_2)!}{(j_1+j_2+j_3+1)!}} \). Adapted from David Simpson (NASA
!> GSFC). Automatically returns 0 for any invalid combination -- violated
!> triangle inequality \( |j_1-j_2|\le j_3\le j_1+j_2 \), \( |m_i|>j_i \>
!> \), or \( m_1+m_2\ne m_3 \) -- so it is always safe to call with an
!> arbitrary sextuple of integers, e.g. when sweeping a parameter range.
!>
!> @param J1 First angular momentum, \( j_1\ge0 \).
!> @param M1 First projection, \( -j_1\le m_1\le j_1 \).
!> @param J2 Second angular momentum, \( j_2\ge0 \).
!> @param M2 Second projection, \( -j_2\le m_2\le j_2 \).
!> @param J3 Coupled angular momentum, \( j_3\ge0 \).
!> @param M3 Coupled projection, \( -j_3\le m_3\le j_3 \).
!> Returns: \( \langle j_1 m_1 j_2 m_2|j_3 m_3\rangle \), real (0 if
!>   selection rules are violated).
  FUNCTION CGCOEFF_DP(J1,M1,J2,M2,J3,M3) RESULT(CGCOEFF)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1
  INTEGER(KIND=i4), INTENT(IN) :: M1
  INTEGER(KIND=i4), INTENT(IN) :: J2
  INTEGER(KIND=i4), INTENT(IN) :: M2
  INTEGER(KIND=i4), INTENT(IN) :: J3
  INTEGER(KIND=i4), INTENT(IN) :: M3
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

!> \( \ln(n!) \), via \( \ln(n!) = \ln\Gamma(n{+}1) \) (the intrinsic
!> `LOG_GAMMA`). Used throughout [[CGCOEFF]], [[SSH]], and
!> [[ASSOC_LEGENDRE]] to evaluate factorial ratios like
!> \( \sqrt{(\ell-k)!/(\ell+k)!} \) as
!> \( \exp(\tfrac12[\texttt{LOG\_FACT}(\ell{-}k)-\texttt{LOG\_FACT}(\ell{+}k)]) \)
!> without overflowing the individual factorials for large \( \ell \).
!>
!> @param N Non-negative integer argument (\( n\le1 \) returns 0).
!> Returns: \( \ln(n!) \).
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

!> Gauss-Legendre quadrature nodes and weights, exact for integrating any
!> polynomial of degree \( \le 2\ell_{max}{+}1 \) using only
!> \( \ell_{max}{+}1 \) points -- the natural quadrature for spherical
!> harmonic transforms, since \( P_\ell^m(\cos\theta) \) is a polynomial in
!> \( x=\cos\theta \) times \( (1{-}x^2)^{|m|/2} \), and the transform
!> integral separates exactly into this \( x \)-quadrature times a
!> uniform sum over \( \phi \) (itself exact for finite Fourier content).
!> Nodes are the zeros of \( P_{\ell_{max}+1}(x) \), found by Newton
!> iteration (converges to machine precision in under 10 steps); weights
!> are \( w_i = 2/[(1-x_i^2)\,P_{\ell_{max}+1}'(x_i)^2] \).
!>
!> @param ZERO Output nodes \( x_i=\cos\theta_i \), size `LMAX+1`,
!>   ascending from -1 to +1.
!> @param W Output weights for \( \int_{-1}^{1}f(x)\,dx \approx
!>   \sum_i w_i f(x_i) \), size `LMAX+1`.
!> @param LMAX Requested exactness is degree \( 2\,\texttt{LMAX}{+}1 \);
!>   uses `LMAX+1` quadrature points.
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

!> Wigner 3-j symbol \( \begin{pmatrix}j_1&j_2&j_3\\m_1&m_2&m_3\end{pmatrix} \),
!> via the standard conversion from [[CGCOEFF]]:
!> $$ \begin{pmatrix}j_1&j_2&j_3\\m_1&m_2&m_3\end{pmatrix} =
!>    \frac{(-1)^{j_1-j_2-m_3}}{\sqrt{2j_3+1}}\,
!>    \langle j_1,-m_1,j_2,-m_2|j_3,-m_3\rangle. $$
!> Inherits [[CGCOEFF]]'s selection-rule guards -- safe to call with any
!> integer sextuple, returning 0 for invalid combinations.
!>
!> @param J1 First angular momentum, \( j_1\ge0 \).
!> @param M1 First projection, \( -j_1\le m_1\le j_1 \).
!> @param J2 Second angular momentum, \( j_2\ge0 \).
!> @param M2 Second projection, \( -j_2\le m_2\le j_2 \).
!> @param J3 Third angular momentum, \( j_3\ge0 \).
!> @param M3 Third projection, \( -j_3\le m_3\le j_3 \).
!> Returns: The Wigner 3-j symbol, real (0 if selection rules are
!>   violated).
  FUNCTION SYMBOL3J_DP(J1,M1,J2,M2,J3,M3) RESULT(SYMBOL3J)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1
  INTEGER(KIND=i4), INTENT(IN) :: M1
  INTEGER(KIND=i4), INTENT(IN) :: J2
  INTEGER(KIND=i4), INTENT(IN) :: M2
  INTEGER(KIND=i4), INTENT(IN) :: J3
  INTEGER(KIND=i4), INTENT(IN) :: M3
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

!> Dot product of two 3-vectors (typically two VSH values evaluated at
!> the same \( (\theta,\phi) \)): \( \mathbf{V}_1\cdot\mathbf{V}_2 =
!> \sum_{c=1}^3 V_{1,c}V_{2,c} \). Note this is a plain component-wise
!> product, *not* a Hermitian inner product -- callers wanting
!> \( \mathbf{V}_1\cdot\mathbf{V}_2^{*} \) must conjugate `VSH2`
!> themselves first (as done throughout `src/tests.f90` and the
!> `examples/` drivers).
!>
!> @param VSH1 First complex 3-vector \( (\hat r,\hat\theta,\hat\phi) \).
!> @param VSH2 Second complex 3-vector \( (\hat r,\hat\theta,\hat\phi) \).
!> Returns: \( \mathbf{VSH1}\cdot\mathbf{VSH2} \), complex scalar.
  FUNCTION DOT_DP(VSH1, VSH2) RESULT(DOT)
  IMPLICIT NONE
  COMPLEX(KIND=dp), DIMENSION(3), INTENT(IN) :: VSH1
  COMPLEX(KIND=dp), DIMENSION(3), INTENT(IN) :: VSH2
  COMPLEX(KIND=dp) :: DOT
  DOT = VSH1(1)*VSH2(1) + VSH1(2)*VSH2(2) + VSH1(3)*VSH2(3)
  RETURN
  END FUNCTION DOT_DP


!> Geppert-Wiebicke axisymmetric coupling coefficient \(
!> I_{j_1 m_1 j_2 m_2}^{\ell m} \) (Geppert & Wiebicke, 1991), a compact
!> shorthand for the generation/amplification amplitude of a two-mode
!> interaction in the Hall-MHD induction equation:
!> $$ I_{j_1 m_1 j_2 m_2}^{\ell m} = \sqrt{\frac{(2j_1+1)(2j_2+1)}
!>    {4\pi(2\ell+1)}}\,\langle j_1 0\,j_2 0|\ell 0\rangle\,
!>    \langle j_1 m_1\,j_2 m_2|\ell m\rangle. $$
!> This is algebraically the standard Gaunt coefficient -- the closed
!> form for \( \int Y_{j_1}^{m_1}Y_{j_2}^{m_2}(Y_\ell^m)^{*}\,d\Omega \) --
!> written in terms of two [[CGCOEFF]] evaluations. Built entirely from
!> [[CGCOEFF]], which guards all selection rules, so `GWI` is safe to call
!> with any integer sextuple and returns 0 for invalid combinations
!> (unlike [[GWJ]], see its own selection-rule note).
!>
!> @param J1 First coupled degree, \( j_1\ge0 \).
!> @param M1 First coupled order, \( -j_1\le m_1\le j_1 \).
!> @param J2 Second coupled degree, \( j_2\ge0 \).
!> @param M2 Second coupled order, \( -j_2\le m_2\le j_2 \).
!> @param L Resulting degree, \( \ell\ge0 \).
!> @param M Resulting order, \( -\ell\le m\le\ell \).
!> Returns: \( I_{j_1 m_1 j_2 m_2}^{\ell m} \), complex (real-valued in
!>   practice, but returned as `COMPLEX` for a uniform interface with
!>   [[GWJ]]).
  FUNCTION GWI_DP(J1,M1,J2,M2,L,M) RESULT(GWI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1
  INTEGER(KIND=i4), INTENT(IN) :: M1
  INTEGER(KIND=i4), INTENT(IN) :: J2
  INTEGER(KIND=i4), INTENT(IN) :: M2
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: M
  COMPLEX(KIND=dp) :: GWI
  GWI=DSQRT((2.d0*J1+1.d0)*(2.d0*J2+1.d0)/(2.d0*L+1.d0)/4.d0/pi)* &
    CGCOEFF(J1,0,J2,0,L,0)*CGCOEFF(J1,M1,J2,M2,L,M)
  RETURN
  END FUNCTION GWI_DP


!> Geppert-Wiebicke non-axisymmetric coupling coefficient \(
!> J_{j_1 m_1 j_2 m_2}^{\ell m} \) (Geppert & Wiebicke, 1991), the
!> companion to [[GWI]] encoding the amplitude of a two-mode Hall-MHD
!> interaction that generates or amplifies a non-axisymmetric component:
!> $$ J_{j_1 m_1 j_2 m_2}^{\ell m} = -\frac{i}{2}
!>    \sqrt{\frac{(2j_1+1)(2j_2+1)}{4\pi(2\ell+1)}}\,
!>    \sqrt{(j_1{+}j_2{+}\ell{+}2)(j_2{+}\ell{-}j_1)(j_1{+}j_2{-}\ell{+}1)
!>    (j_1{-}j_2{+}\ell{+}1)}\; \times $$
!> $$ \langle j_1{+}1,0\,j_2 0|\ell 0\rangle\,
!>    \langle j_1 m_1\,j_2 m_2|\ell m\rangle. $$
!> Cross-validated against an independent VSH-dot-product expansion in
!> `TEST2_GW` (`src/tests.f90`).
!>
!> @warning Unlike [[GWI]], the extra square-root factor above is *not*
!>   internally selection-rule-guarded by [[CGCOEFF]] alone -- it can go
!>   negative for out-of-range `(J1,J2,L)` triples. This routine guards it
!>   explicitly (returning 0, matching the mathematically correct value)
!>   so it is still safe to call with any integer sextuple.
!>
!> @param J1 First coupled degree, \( j_1\ge0 \).
!> @param M1 First coupled order, \( -j_1\le m_1\le j_1 \).
!> @param J2 Second coupled degree, \( j_2\ge0 \).
!> @param M2 Second coupled order, \( -j_2\le m_2\le j_2 \).
!> @param L Resulting degree, \( \ell\ge0 \).
!> @param M Resulting order, \( -\ell\le m\le\ell \).
!> Returns: \( J_{j_1 m_1 j_2 m_2}^{\ell m} \), complex.
  FUNCTION GWJ_DP(J1,M1,J2,M2,L,M) RESULT(GWJ)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: J1
  INTEGER(KIND=i4), INTENT(IN) :: M1
  INTEGER(KIND=i4), INTENT(IN) :: J2
  INTEGER(KIND=i4), INTENT(IN) :: M2
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: M
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


!> 1D packed-triangular index for \( (\ell,m) \) into the output arrays of
!> [[ASSOC_LEGENDRE_ALL]], [[DDX_ASSOC_LEGENDRE_ALL]],
!> [[ASSOC_LEGENDRE_NORM_ALL]], and [[DDX_ASSOC_LEGENDRE_NORM_ALL]] (all
!> of which only need \( 0\le m\le\ell \), by conjugate/parity symmetry).
!> Matches the SHTOOLS `PlmIndex` convention.
!>
!> @param L Degree, \( \ell\ge0 \).
!> @param M Order, \( 0\le m\le\ell \).
!> Returns: Packed index \( \ell(\ell+1)/2+m+1 \), \( 1\le
!>   \texttt{PLM\_INDEX}\le(\ell_{max}{+}1)(\ell_{max}{+}2)/2 \).
  FUNCTION PLM_INDEX(L, M)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: M
  INTEGER(KIND=i4) :: PLM_INDEX
  PLM_INDEX = L*(L+1)/2 + M + 1
  END FUNCTION PLM_INDEX

!> 1D index for \( (\ell,m) \) into the output arrays of every batch
!> `_ALL` routine that needs the full \( -\ell\le m\le\ell \) range
!> ([[SSH_ALL]], [[GRAD_SSH_ALL]], [[L_SSH_ALL]], and every batch VSH
!> routine).
!>
!> @param L Degree, \( \ell\ge0 \).
!> @param M Order, \( -\ell\le m\le\ell \).
!> Returns: Index \( \ell^2+\ell+m+1 \), \( 1\le\texttt{YLM\_INDEX}\le
!>   (\ell_{max}{+}1)^2 \).
  FUNCTION YLM_INDEX(L, M)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: L
  INTEGER(KIND=i4), INTENT(IN) :: M
  INTEGER(KIND=i4) :: YLM_INDEX
  YLM_INDEX = L**2 + L + M + 1
  END FUNCTION YLM_INDEX

!> All unnormalized associated Legendre functions \( P_\ell^m(x) \) for
!> \( 0\le\ell\le\ell_{max},\,0\le m\le\ell \), via the Bonnet recurrence,
!> in one pass (amortizing shared recurrence terms across every \( \ell
!> \) at fixed \( m \) -- far cheaper than \( (\ell_{max}{+}1)^2/2 \)
!> independent [[ASSOC_LEGENDRE]] calls). Condon-Shortley phase; same
!> recurrence and same stability limit as [[ASSOC_LEGENDRE]] (see
!> [[ASSOC_LEGENDRE_NORM_ALL]] for the far-more-stable normalized
!> alternative, itself safe to \( \ell\le2000 \)).
!>
!> @warning **Numerical stability**: verified accurate (relative error
!>   `< 3e-14`) for \( \ell\le150 \) at any \( m,\theta \) -- the worst
!>   case, \( \ell=151 \), occurs near the equator (\( \theta\approx90 \)
!>   degrees) at \( m=\ell \); safety improves toward either pole. This
!>   supersedes an earlier "\( \ell\sim1400 \)" claim, verified too
!>   optimistic by direct comparison against an arbitrary-precision
!>   reference. Full methodology and per-degree table:
!>   `STABILITY_FINDINGS.md` at the repo root, Part 2 (regenerable via the
!>   opt-in `VSH_BUILD_STABILITY` CMake target's `vsh_unnorm_stability`
!>   executable).
!>
!> @param P Output, size `(LMAX+1)*(LMAX+2)/2`, indexed by
!>   [[PLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param X Argument, \( -1\le x\le1 \) (typically \( x=\cos\theta \)).
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

!> All derivatives \( dP_\ell^m/dx \) for \( 0\le\ell\le\ell_{max},\,
!> 0\le m\le\ell \), reusing a precomputed [[ASSOC_LEGENDRE_ALL]] table
!> rather than recomputing \( P_\ell^m \) from scratch.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_ALL]]'s stability limit (safe for
!>   \( \ell\le150 \) at any \( m,\theta \)) -- see `STABILITY_FINDINGS.md`
!>   at the repo root, Part 2.
!>
!> @param DP_OUT Output derivatives, size `(LMAX+1)*(LMAX+2)/2`, indexed
!>   by [[PLM_INDEX]](l,m); 0 at the poles (\( |x|\ge1 \)).
!> @param P Precomputed \( P_\ell^m(x) \) table from
!>   [[ASSOC_LEGENDRE_ALL]], same size/indexing.
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \) (must match `P`'s).
!> @param X Argument, \( -1\le x\le1 \) (must match the `X` used to build
!>   `P`).
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

!> All \( 4\pi \)-normalized associated Legendre functions
!> \( N_\ell^m P_\ell^m(\cos\theta) \), \( N_\ell^m=\sqrt{\frac{2\ell+1}
!> {4\pi}\frac{(\ell-m)!}{(\ell+m)!}} \) (i.e. the real-valued radial part
!> of [[SSH]], for \( 0\le\ell\le\ell_{max},\,0\le m\le\ell \)), using the
!> Holmes & Featherstone (2002) modified forward-column recurrence. The
!> normalization is folded directly into the recurrence coefficients,
!> keeping every intermediate value \( O(1/\sqrt{4\pi}) \) rather than
!> letting \( P_m^m \) itself grow/shrink combinatorially like the
!> unnormalized [[ASSOC_LEGENDRE_ALL]]/Bonnet recurrence does. This is
!> the recurrence [[SSH_ALL]] and every batch VSH routine build on.
!>
!> @warning **Numerical stability**: verified safe for \( \ell\le2000 \)
!>   at any \( m,\theta \). Beyond that, safety depends on \( \theta \):
!>   this recurrence has a distinct failure mode near the classical
!>   turning point of the associated Legendre ODE (\( \sin\theta\approx
!>   m/\ell \)) -- for \( \theta \) within about 30-56 degrees of either
!>   pole, values can blow up by 100+ orders of magnitude once \( \ell \)
!>   exceeds a \(\theta\)-dependent onset (as low as \( \ell=2100 \) in
!>   the worst case, \( \theta\approx30 \) degrees). Near the poles
!>   (\( \theta<30 \) degrees) or equator (\( 57\le\theta\le123 \)
!>   degrees) it remains accurate at least to \( \ell=5000 \) (untested
!>   beyond that). This supersedes an earlier, overly optimistic
!>   "stable to \( \ell\sim2700 \)" claim that didn't account for the
!>   \(\theta\)-dependence. Full methodology, per-degree onset table, and
!>   a `(\ell,\theta)` heatmap: `STABILITY_FINDINGS.md` at the repo root
!>   (regenerable via the opt-in `VSH_BUILD_STABILITY` CMake target).
!>
!> @param PNORM Output, size `(LMAX+1)*(LMAX+2)/2`, indexed by
!>   [[PLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param X Argument, \( -1\le x\le1 \) (typically \( x=\cos\theta \)).
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


!> All derivatives \( d(N_\ell^m P_\ell^m)/dx \) of the normalized
!> associated Legendre functions, for \( 0\le\ell\le\ell_{max},\,0\le
!> m\le\ell \), reusing a precomputed [[ASSOC_LEGENDRE_NORM_ALL]] table:
!> $$ \frac{d(N_\ell^mP_\ell^m)}{dx} = \frac{\sqrt{\frac{(2\ell+1)
!>    (\ell^2-m^2)}{2\ell-1}}\,N_{\ell-1}^mP_{\ell-1}^m - \ell\,x\,
!>    N_\ell^mP_\ell^m}{1-x^2}. $$
!> This is the recurrence [[VSH_CORE]] (and therefore every batch VSH
!> routine) uses for the \( \hat\theta \) component.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param DPNORM Output derivatives, size `(LMAX+1)*(LMAX+2)/2`, indexed
!>   by [[PLM_INDEX]](l,m); 0 at the poles (\( |x|\ge1 \)).
!> @param PNORM Precomputed normalized table from
!>   [[ASSOC_LEGENDRE_NORM_ALL]], same size/indexing.
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \) (must match `PNORM`'s).
!> @param X Argument, \( -1\le x\le1 \) (must match the `X` used to build
!>   `PNORM`).
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


!> All scalar spherical harmonics \( Y_\ell^m(\theta,\phi) \) for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), in one call. Computes
!> [[ASSOC_LEGENDRE_NORM_ALL]] once, then fills \( m\ge0 \) directly and
!> \( m<0 \) via the conjugate symmetry \( Y_\ell^{-m} =
!> (-1)^m\,(Y_\ell^m)^{*} \) rather than a second recurrence pass.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param YLM Output, size `(LMAX+1)**2`, indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE SSH_ALL_DP(YLM, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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

!> Private helper shared by every batch VSH routine
!> ([[GRAD_SSH_ALL]], [[L_SSH_ALL]], [[PVSH_RAD_ALL]], [[PVSH_POL_ALL]],
!> [[PVSH_TOR_ALL]], [[VSH_TOR_ALL]], [[VSH_POL_UP_ALL]],
!> [[VSH_POL_DN_ALL]]): computes \( Y_\ell^m \) and its two angular
!> gradient components ([[GRAD_SSH]]'s \( \hat\theta,\hat\phi \) parts)
!> for every \( (\ell,m) \) in one pass, so each caller only has to apply
!> its own linear combination/normalization on top.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param YLM_OUT Output \( Y_\ell^m \), size `(LMAX+1)**2`, indexed by
!>   [[YLM_INDEX]](l,m).
!> @param GSSH_TH Output \( \hat\theta \) component of \( \nabla_\perp
!>   Y_\ell^m \), same size/indexing.
!> @param GSSH_PH Output \( \hat\phi \) component of \( \nabla_\perp
!>   Y_\ell^m \), same size/indexing.
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE VSH_CORE_DP(YLM_OUT, GSSH_TH, GSSH_PH, &
                      LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All angular gradients \( \nabla_\perp Y_\ell^m \) ([[GRAD_SSH]]) for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one
!> [[VSH_CORE]] pass.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\) (always 0), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode
!>   axis indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE GRAD_SSH_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All angular-momentum-operator fields \( \hat r\times\nabla_\perp
!> Y_\ell^m \) ([[L_SSH]]) for \( 0\le\ell\le\ell_{max},\,-\ell\le
!> m\le\ell \), from one [[VSH_CORE]] pass.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\) (always 0), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode
!>   axis indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE L_SSH_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All radial polar VSH members [[PVSH_RAD]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one [[VSH_CORE]]
!> pass.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\)=\(Y_\ell^m\), `2`=\(\hat\theta\) (always 0),
!>   `3`=\(\hat\phi\) (always 0); mode axis indexed by
!>   [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE PVSH_RAD_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All poloidal (horizontal) polar VSH members [[PVSH_POL]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one [[VSH_CORE]]
!> pass. Zero for \( \ell=0 \).
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\) (always 0), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode
!>   axis indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE PVSH_POL_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All toroidal polar VSH members [[PVSH_TOR]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one [[VSH_CORE]]
!> pass. Zero for \( \ell=0 \).
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\) (always 0), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode
!>   axis indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE PVSH_TOR_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All toroidal standard-VSH members [[VSH_TOR]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \). Identical to
!> [[PVSH_TOR_ALL]] (thin wrapper) -- the toroidal member is common to
!> both bases.
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\) (always 0), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode
!>   axis indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE VSH_TOR_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
  COMPLEX(KIND=dp), INTENT(OUT) :: OUT(3,(LMAX+1)**2)
  CALL PVSH_TOR_ALL(OUT, LMAX, THETA, PHI)
  END SUBROUTINE VSH_TOR_ALL_DP


!> All \( J=\ell{+}1 \) standard-VSH members [[VSH_POL_UP]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one [[VSH_CORE]]
!> pass: \( -\sqrt{(\ell{+}1)/(2\ell{+}1)}\,\mathbf{Y}_{\ell m}^{(-1)} +
!> \sqrt{\ell/(2\ell{+}1)}\,\mathbf{Y}_{\ell m}^{(+1)} \). Well-defined
!> (nonzero) at \( \ell=0 \) -- see [[VSH_POL_UP]] for the dipole-ratio
!> note at \( \ell=1 \).
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode axis
!>   indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE VSH_POL_UP_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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


!> All \( J=\ell{-}1 \) standard-VSH members [[VSH_POL_DN]] for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), from one [[VSH_CORE]]
!> pass: \( \sqrt{\ell/(2\ell{+}1)}\,\mathbf{Y}_{\ell m}^{(-1)} +
!> \sqrt{(\ell{+}1)/(2\ell{+}1)}\,\mathbf{Y}_{\ell m}^{(+1)} \). Zero for
!> \( \ell=0 \) -- see [[VSH_POL_DN]] for the dipole-ratio caveat at
!> \( \ell=1 \).
!>
!> @warning Inherits [[ASSOC_LEGENDRE_NORM_ALL]]'s \(\theta\)-dependent
!>   stability limit (safe for \( \ell\le2000 \) unconditionally) --
!>   see `STABILITY_FINDINGS.md` at the repo root.
!>
!> @param OUT Output, shape `(3,(LMAX+1)**2)`; component axis is
!>   `1`=\(\hat r\), `2`=\(\hat\theta\), `3`=\(\hat\phi\); mode axis
!>   indexed by [[YLM_INDEX]](l,m).
!> @param LMAX Maximum degree, \( \ell_{max}\ge0 \).
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE VSH_POL_DN_ALL_DP(OUT, LMAX, THETA, PHI)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  REAL(KIND=dp),    INTENT(IN)  :: THETA
  REAL(KIND=dp),    INTENT(IN)  :: PHI
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
