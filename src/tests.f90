!> The numerical validation suite backing FORTVSH's accuracy claims (see
!> the "Numerical Validation" appendix of the accompanying manuscript).
!> Three families of checks:
!>
!> 1. **Cross-validation** ([[TEST1_ANALYTIC]]/[[TEST1_VSH]]/[[TEST1_GW]]
!>    and their `TEST2_*` counterparts): the same closed-form quantity
!>    computed three independent ways -- a hand-derived analytic formula,
!>    a direct [[DOT]] product of VSH values, and the
!>    [[GWI]]/[[GWJ]]-coefficient expansion -- must agree to machine
!>    precision.
!> 2. **Batch consistency** (`BATCH_*_CONS`): every `_ALL` batch routine
!>    must agree pointwise with its single-mode counterpart, looped over
!>    the same \( (\ell,m) \) range.
!> 3. **Mathematical identities** ([[SSH_ORTHO]],
!>    [[PVSH_POL_TOR_ORTHO]], [[VSH_POL_INVERSION]]): orthonormality,
!>    orthogonality, and the [[VSH_POL_UP]]/[[VSH_POL_DN]] rotation back
!>    to [[PVSH_POL]]/[[PVSH_RAD]], each checked against its known exact
!>    value.
!>
!> Every `_CONS`/`_ORTHO`/`_INVERSION` routine returns a `STATUS` (0 =
!> pass, 1 = fail against its stated tolerance) and writes its full
!> per-point data to `OUTUNIT`, driven from `src/main.f90` -- the same
!> data that feeds the manuscript's convergence figures via
!> `py/compute_validation.py`.
MODULE TESTS
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi, j
USE VSH,     ONLY: &
  ASSOC_LEGENDRE, DDX_ASSOC_LEGENDRE, &
  ASSOC_LEGENDRE_ALL, DDX_ASSOC_LEGENDRE_ALL, &
  PLM_INDEX, YLM_INDEX, LOG_FACT, &
  SSH, SSH_ALL, &
  GRAD_SSH, GRAD_SSH_ALL, L_SSH, L_SSH_ALL, &
  PVSH_RAD, PVSH_RAD_ALL, &
  PVSH_TOR, PVSH_TOR_ALL, &
  PVSH_POL, PVSH_POL_ALL, &
  VSH_TOR, VSH_TOR_ALL, &
  VSH_POL_DN, VSH_POL_DN_ALL, &
  VSH_POL_UP, VSH_POL_UP_ALL, &
  DOT, GWI, GWJ
IMPLICIT NONE
PRIVATE
PUBLIC :: &
  ! Analytic and numerical test functions
  TEST1_ANALYTIC, TEST1_VSH, TEST1_GW, &
  TEST2_ANALYTIC, TEST2_VSH, TEST2_GW, &
  ! Batch consistency checks
  BATCH_ALM_CONS, BATCH_DALM_CONS, BATCH_SSH_CONS, &
  BATCH_GRAD_SSH_CONS, BATCH_L_SSH_CONS, &
  BATCH_PVSH_RAD_CONS, BATCH_PVSH_POL_CONS, BATCH_PVSH_TOR_CONS, &
  BATCH_VSH_TOR_CONS, BATCH_VSH_POL_UP_CONS, BATCH_VSH_POL_DN_CONS, &
  ! Orthogonality and inversion validation
  SSH_ORTHO, PVSH_POL_TOR_ORTHO, VSH_POL_INVERSION

CONTAINS

!> Hand-derived closed form for Test 1's target quantity (a specific
!> product of VSH inner products), evaluated over a colatitude grid at
!> \( \phi=0 \) for comparison against [[TEST1_VSH]] and [[TEST1_GW]].
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, the analytic Test 1 values.
  FUNCTION TEST1_ANALYTIC(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp) :: TH
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_ANALYTIC

  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEST1_ANALYTIC(I)=3.d0*DSQRT(5.d0)/64.d0/pi**2* &
      (DSIN(TH)*DCOS(TH))**2
  ENDDO
  RETURN
  END FUNCTION TEST1_ANALYTIC

!> Test 1's target quantity computed directly from [[DOT]] products of
!> [[VSH_TOR]]/[[VSH_POL_DN]]/[[VSH_POL_UP]] values -- must agree with
!> [[TEST1_ANALYTIC]] to machine precision.
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, matching [[TEST1_ANALYTIC]]'s
!>   grid.
  FUNCTION TEST1_VSH(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp) :: TH,PH
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_VSH

  PH = 0.123d0   !- arbitrary choice
  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEST1_VSH(I) = DOT(VSH_TOR(2,0,TH,PH),VSH_TOR(1,1,TH,PH))* &
         DOT(VSH_POL_DN(1,-1,TH,PH),CONJG(VSH_POL_UP(0,0,TH,PH)))
  ENDDO
  RETURN
  END FUNCTION TEST1_VSH

!> Test 1's target quantity computed via the [[GWI]]-coefficient expansion
!> in scalar spherical harmonics, rather than direct VSH [[DOT]] products
!> -- must agree with [[TEST1_ANALYTIC]] and [[TEST1_VSH]] to machine
!> precision. This is the check that specifically exercises [[GWI]] (see
!> also [[TEST2_GW]] for the [[GWJ]] counterpart).
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, matching [[TEST1_ANALYTIC]]'s
!>   grid.
  FUNCTION TEST1_GW(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I,I1,I2
  REAL(KIND=dp) :: TH,PH
  COMPLEX(KIND=dp) :: TEMP1, TEMP2
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST1_GW

  PH = 0.123d0   !- arbitrary choice
  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEMP1 = DCMPLX(0.d0,0.d0)
    TEMP2 = DCMPLX(0.d0,0.d0)
    DO I1=1,3,2
      TEMP1 = TEMP1-((8.d0-I1*(I1+1.d0))/2.d0/DSQRT(12.d0))* &
        GWI(2,0,1,1,I1,1)*SSH(I1,1,TH,PH)
    ENDDO
    I2 = 1
    TEMP2 = -(1.d0/DSQRT(3.d0))*GWI(1,-1,1,1,0,0)* &
        CONJG(SSH(1,1,TH,PH))

    TEST1_GW(I) = TEMP1*TEMP2
  ENDDO
  RETURN
  END FUNCTION TEST1_GW

!> Hand-derived closed form for Test 2's target quantity (a different
!> product of VSH inner products than Test 1, with nonzero \( m \)
!> content), evaluated over a colatitude grid at a fixed \( \phi \) for
!> comparison against [[TEST2_VSH]] and [[TEST2_GW]].
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, the analytic Test 2 values.
  FUNCTION TEST2_ANALYTIC(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp) :: TH,PH
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_ANALYTIC
  PH = 1.006d0   !- arbitrary choice
  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEST2_ANALYTIC(I)=15.d0*DSQRT(15.d0*7.d0)/128.d0/pi**2* &
      (DSIN(TH)*DCOS(TH))**3*DSIN(TH)/DSQRT(2.d0)* &
      exp(4.d0*j*PH)
  ENDDO
  RETURN
  END FUNCTION TEST2_ANALYTIC

!> Test 2's target quantity computed directly from [[DOT]] products of
!> [[PVSH_POL]]/[[PVSH_TOR]]/[[PVSH_RAD]] values -- must agree with
!> [[TEST2_ANALYTIC]] to machine precision.
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, matching [[TEST2_ANALYTIC]]'s
!>   grid.
  FUNCTION TEST2_VSH(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp) :: TH,PH
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_VSH

  PH = 1.006d0   !- arbitrary choice
  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEST2_VSH(I) = DOT(PVSH_POL(2,0,TH,PH),PVSH_TOR(3,2,TH,PH))* &
         DOT(PVSH_RAD(1,0,TH,PH),CONJG(PVSH_RAD(2,-2,TH,PH)))
  ENDDO
  RETURN
  END FUNCTION TEST2_VSH

!> Test 2's target quantity computed via the [[GWJ]]-coefficient expansion
!> (combined with a [[GWI]] term) in scalar spherical harmonics -- must
!> agree with [[TEST2_ANALYTIC]] and [[TEST2_VSH]] to machine precision.
!> This is the check that specifically exercises [[GWJ]].
!>
!> @param NTH Number of colatitude samples over \( [0,\pi] \).
!> Returns: Complex array of size `NTH`, matching [[TEST2_ANALYTIC]]'s
!>   grid.
  FUNCTION TEST2_GW(NTH)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4) :: I,I1,I2
  REAL(KIND=dp) :: TH,PH
  COMPLEX(KIND=dp) :: TEMP1, TEMP2
  COMPLEX(KIND=dp), DIMENSION(NTH) :: TEST2_GW

  PH = 1.006d0   !- arbitrary choice
  DO I=1,NTH
    TH = (I-1)*pi/(NTH-1)
    TEMP1 = DCMPLX(0.d0,0.d0)
    TEMP2 = DCMPLX(0.d0,0.d0)
    DO I1=2,4,2
      TEMP1 = TEMP1+GWJ(3,2,2,0,I1,2)*SSH(I1,2,TH,PH)
    ENDDO
    TEMP1=-TEMP1*j/DSQRT(6.d0*12.d0)
    DO I2=1,3,2
      TEMP2 = TEMP2+GWI(1,0,I2,-2,2,-2)* &
        CONJG(SSH(I2,-2,TH,PH))
    ENDDO
    TEST2_GW(I) = TEMP1*TEMP2
  ENDDO
  RETURN
  END FUNCTION TEST2_GW


!> Batch-consistency check: [[ASSOC_LEGENDRE_ALL]] vs. looped
!> [[ASSOC_LEGENDRE]] calls, for \( 0\le\ell\le\ell_{max},\,0\le m\le\ell
!> \), swept over `NX` points \( x\in[-0.9,0.9] \). Writes every
!> `(l,m,x)` triple's batch value, single-mode value, and absolute
!> difference to `OUTUNIT`; fails (`STATUS=1`) if any relative error
!> exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NX Number of \( x \) points swept over \( [-0.9,0.9] \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   X P_batch P_single abs_diff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some relative error
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_ALM_CONS(LMAX, NX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, IX, PSIZE
  REAL(KIND=dp) :: X, DX, P_SINGLE, SCALE, REL_ERR
  REAL(KIND=dp), ALLOCATABLE :: P_BATCH(:)

  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(P_BATCH(PSIZE))
  STATUS = 0
  DX = 1.8_dp / (NX - 1)

  WRITE(OUTUNIT, '(A)') '# L  M  X  P_batch  P_single  abs_diff'
  DO IX = 0, NX-1
    X = -0.9_dp + IX * DX
    CALL ASSOC_LEGENDRE_ALL(P_BATCH, LMAX, X)
    DO L = 0, LMAX
      DO M = 0, L
        P_SINGLE = ASSOC_LEGENDRE(L, M, X)
        WRITE(OUTUNIT, *) L, M, X, &
            P_BATCH(PLM_INDEX(L, M)), P_SINGLE, &
            ABS(P_BATCH(PLM_INDEX(L, M)) - P_SINGLE)
        SCALE = MAX(ABS(P_SINGLE), 1.0_dp)
        REL_ERR = ABS(P_BATCH(PLM_INDEX(L, M)) - P_SINGLE) / SCALE
        IF (REL_ERR > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO

  DEALLOCATE(P_BATCH)
  END SUBROUTINE BATCH_ALM_CONS


!> Batch-consistency check: [[DDX_ASSOC_LEGENDRE_ALL]] vs. looped
!> [[DDX_ASSOC_LEGENDRE]] calls, for \( 0\le\ell\le\ell_{max},\,0\le
!> m\le\ell \), swept over `NX` points \( x\in[-0.9,0.9] \). Same
!> tolerance/output convention as [[BATCH_ALM_CONS]].
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NX Number of \( x \) points swept over \( [-0.9,0.9] \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   X dP_batch dP_single abs_diff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some relative error
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_DALM_CONS(LMAX, NX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NX, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, IX, PSIZE
  REAL(KIND=dp) :: X, DX, DP_SINGLE, SCALE, REL_ERR
  REAL(KIND=dp), ALLOCATABLE :: P_BATCH(:), DP_BATCH(:)

  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(P_BATCH(PSIZE), DP_BATCH(PSIZE))
  STATUS = 0
  DX = 1.8_dp / (NX - 1)

  WRITE(OUTUNIT, '(A)') '# L  M  X  dP_batch  dP_single  abs_diff'
  DO IX = 0, NX-1
    X = -0.9_dp + IX * DX
    CALL ASSOC_LEGENDRE_ALL(P_BATCH, LMAX, X)
    CALL DDX_ASSOC_LEGENDRE_ALL(DP_BATCH, P_BATCH, LMAX, X)
    DO L = 0, LMAX
      DO M = 0, L
        DP_SINGLE = DDX_ASSOC_LEGENDRE(L, M, X)
        WRITE(OUTUNIT, *) L, M, X, &
            DP_BATCH(PLM_INDEX(L, M)), DP_SINGLE, &
            ABS(DP_BATCH(PLM_INDEX(L, M)) - DP_SINGLE)
        SCALE = MAX(ABS(DP_SINGLE), 1.0_dp)
        REL_ERR = ABS(DP_BATCH(PLM_INDEX(L, M)) - DP_SINGLE) / SCALE
        IF (REL_ERR > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO

  DEALLOCATE(P_BATCH, DP_BATCH)
  END SUBROUTINE BATCH_DALM_CONS


!> Batch-consistency check: [[SSH_ALL]] vs. looped [[SSH]] calls, for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over `NTH`
!> colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`) if any
!> \( |Y_{batch}-Y_{single}| \) exceeds \( 10^{-13} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi re(batch) im(batch) re(single) im(single)`).
!> @param STATUS Output: 0 = pass, 1 = fail (some
!>   \( |Y_{batch}-Y_{single}|>10^{-13} \)).
  SUBROUTINE BATCH_SSH_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp) :: Y_SINGLE
  COMPLEX(KIND=dp), ALLOCATABLE :: YLM(:)

  NYLM = (LMAX+1)**2
  ALLOCATE(YLM(NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0

  WRITE(OUTUNIT, '(A)') &
      '# L  M  theta  phi  re(batch)  im(batch)'// &
      '  re(single)  im(single)'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL SSH_ALL(YLM, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        Y_SINGLE = SSH(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            DREAL(YLM(YLM_INDEX(L, M))), &
            DIMAG(YLM(YLM_INDEX(L, M))), &
            DREAL(Y_SINGLE), DIMAG(Y_SINGLE)
        IF (ABS(YLM(YLM_INDEX(L,M)) - Y_SINGLE) > 1.0E-13_dp) STATUS = 1
      END DO
    END DO
  END DO

  DEALLOCATE(YLM)
  END SUBROUTINE BATCH_SSH_CONS


!> Verifies scalar spherical harmonic orthonormality,
!> \( \int Y_{\ell_1}^m(Y_{\ell_2}^m)^{*}\,d\Omega = \delta_{\ell_1\ell_2}
!> \), for every \( 0\le m\le\ell_1\le\ell_2\le\ell_{max} \). The
!> \( \phi \) integral is done analytically (\( 2\pi \) for matching
!> \( m \), by construction 0 otherwise, so only same-\( m \) pairs are
!> checked); the \( \theta \) integral (in \( u=\cos\theta \)) uses a
!> plain `NQUAD`-point midpoint rule via [[ASSOC_LEGENDRE_ALL]], not
!> [[SHGLQ]] -- sufficient here since only a scalar pass/fail diagonal
!> check is needed, not a high-order transform.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NQUAD Number of midpoint-rule quadrature points over
!>   \( u\in[-1,1] \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L1 M
!>   L2 inner_product expected(0_or_1)`, upper triangle \( \ell_1\le
!>   \ell_2 \) only).
!> @param STATUS Output: 0 = pass, 1 = fail (any diagonal off from 1, or
!>   off-diagonal off from 0, by more than \( 10^{-3} \)).
  SUBROUTINE SSH_ORTHO(LMAX, NQUAD, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NQUAD, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L1, L2, M, IQ, PSIZE, EXPECTED
  REAL(KIND=dp) :: U, DU, N_L1M, N_L2M, FULL_INT
  REAL(KIND=dp), ALLOCATABLE :: P(:), ORTHO(:,:,:)

  PSIZE = (LMAX+1)*(LMAX+2)/2
  ALLOCATE(P(PSIZE))
  ALLOCATE(ORTHO(0:LMAX, 0:LMAX, 0:LMAX))
  STATUS = 0
  ORTHO = 0.0_dp
  DU = 2.0_dp / NQUAD

  !     Accumulate int_{-1}^{1} P_l1^m(u) * P_l2^m(u) du via midpoint rule
  DO IQ = 1, NQUAD
    U = -1.0_dp + (IQ - 0.5_dp) * DU
    CALL ASSOC_LEGENDRE_ALL(P, LMAX, U)
    DO M = 0, LMAX
      DO L1 = M, LMAX
        DO L2 = L1, LMAX
          ORTHO(L1, L2, M) = ORTHO(L1, L2, M) + &
              P(PLM_INDEX(L1, M)) * P(PLM_INDEX(L2, M))
        END DO
      END DO
    END DO
  END DO
  ORTHO = ORTHO * DU

  WRITE(OUTUNIT, '(A)') &
      '# L1  M  L2  inner_product  expected(0_or_1)'
  DO M = 0, LMAX
    DO L1 = M, LMAX
      DO L2 = L1, LMAX
        N_L1M = DSQRT((2.d0*L1+1.d0)/(4.d0*pi)) * &
            EXP(0.5_dp*(LOG_FACT(L1-M) - LOG_FACT(L1+M)))
        N_L2M = DSQRT((2.d0*L2+1.d0)/(4.d0*pi)) * &
            EXP(0.5_dp*(LOG_FACT(L2-M) - LOG_FACT(L2+M)))
        FULL_INT = ORTHO(L1, L2, M) * 2.d0*pi * N_L1M * N_L2M
        IF (L1 .EQ. L2) THEN
          EXPECTED = 1
        ELSE
          EXPECTED = 0
        END IF
        WRITE(OUTUNIT, *) L1, M, L2, FULL_INT, EXPECTED
        IF (L1 .EQ. L2) THEN
          IF (ABS(FULL_INT - 1.0_dp) > 1.0E-3_dp) STATUS = 1
        ELSE
          IF (ABS(FULL_INT) > 1.0E-3_dp) STATUS = 1
        END IF
      END DO
    END DO
  END DO

  DEALLOCATE(P, ORTHO)
  END SUBROUTINE SSH_ORTHO


!> Batch-consistency check: [[GRAD_SSH_ALL]] vs. looped [[GRAD_SSH]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`)
!> if any component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_GRAD_SSH_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL GRAD_SSH_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = GRAD_SSH(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_GRAD_SSH_CONS


!> Batch-consistency check: [[L_SSH_ALL]] vs. looped [[L_SSH]] calls, for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over `NTH`
!> colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`) if any
!> component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_L_SSH_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL L_SSH_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = L_SSH(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_L_SSH_CONS


!> Batch-consistency check: [[PVSH_RAD_ALL]] vs. looped [[PVSH_RAD]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`)
!> if any component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_PVSH_RAD_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL PVSH_RAD_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = PVSH_RAD(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_PVSH_RAD_CONS


!> Batch-consistency check: [[PVSH_POL_ALL]] vs. looped [[PVSH_POL]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`)
!> if any component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_PVSH_POL_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL PVSH_POL_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = PVSH_POL(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_PVSH_POL_CONS


!> Batch-consistency check: [[PVSH_TOR_ALL]] vs. looped [[PVSH_TOR]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \) (both routines agree
!> at \( \ell=0 \) too -- [[PVSH_TOR]] explicitly branches to return zero
!> there rather than evaluating a 0/0 form). Fails (`STATUS=1`) if any
!> component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_PVSH_TOR_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL PVSH_TOR_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = PVSH_TOR(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_PVSH_TOR_CONS


!> Batch-consistency check: [[VSH_TOR_ALL]] vs. looped [[VSH_TOR]] calls,
!> for \( 1\le\ell\le\ell_{max},\,-\ell\le m\le\ell \) (starts at
!> \( \ell=1 \); [[VSH_TOR]] is identical to [[PVSH_TOR]], already
!> covered at \( \ell=0 \) by [[BATCH_PVSH_TOR_CONS]]), swept over `NTH`
!> colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`) if any
!> component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge1 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_VSH_TOR_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL VSH_TOR_ALL(OUT, LMAX, THETA, PHI)
    DO L = 1, LMAX
      DO M = -L, L
        VS = VSH_TOR(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_VSH_TOR_CONS


!> Batch-consistency check: [[VSH_POL_UP_ALL]] vs. looped [[VSH_POL_UP]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`)
!> if any component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_VSH_POL_UP_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL VSH_POL_UP_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = VSH_POL_UP(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_VSH_POL_UP_CONS


!> Batch-consistency check: [[VSH_POL_DN_ALL]] vs. looped [[VSH_POL_DN]]
!> calls, for \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), swept over
!> `NTH` colatitude points at fixed \( \phi=\pi/4 \). Fails (`STATUS=1`)
!> if any component's absolute difference exceeds \( 10^{-12} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi absdiff_r absdiff_th absdiff_ph`).
!> @param STATUS Output: 0 = pass, 1 = fail (some component difference
!>   \( >10^{-12} \)).
  SUBROUTINE BATCH_VSH_POL_DN_CONS(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp), DIMENSION(3) :: VS
  COMPLEX(KIND=dp), ALLOCATABLE :: OUT(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(OUT(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  absdiff_r  absdiff_th  absdiff_ph'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL VSH_POL_DN_ALL(OUT, LMAX, THETA, PHI)
    DO L = 0, LMAX
      DO M = -L, L
        VS = VSH_POL_DN(L, M, THETA, PHI)
        WRITE(OUTUNIT, *) L, M, THETA, PHI, &
            ABS(OUT(1,YLM_INDEX(L,M))-VS(1)), &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)), &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3))
        IF (ABS(OUT(1,YLM_INDEX(L,M))-VS(1)) > 1.0E-12_dp .OR. &
            ABS(OUT(2,YLM_INDEX(L,M))-VS(2)) > 1.0E-12_dp .OR. &
            ABS(OUT(3,YLM_INDEX(L,M))-VS(3)) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(OUT)
  END SUBROUTINE BATCH_VSH_POL_DN_CONS


!> Verifies pointwise orthogonality of the poloidal and toroidal polar
!> VSH members, \( \mathbf{Y}_{\ell m}^{(+1)}\cdot\mathbf{Y}_{\ell
!> m}^{(0)} = 0 \) for every \( \ell\ge1,\,m,\,\theta,\,\phi \) -- an
!> exact analytic identity (poloidal is a pure gradient, toroidal a pure
!> \( \hat r\times \)gradient, and \( \nabla_\perp f\cdot\hat r\times
!> \nabla_\perp f=0 \) always), so any nonzero [[DOT]] value here is
!> purely numerical error.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge1 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi |DOT(PVSH_POL,PVSH_TOR)|`).
!> @param STATUS Output: 0 = pass, 1 = fail (some \( |{\rm DOT}|>10^{-12}
!>   \)).
  SUBROUTINE PVSH_POL_TOR_ORTHO(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, ITH, NYLM, IDX
  REAL(KIND=dp) :: THETA, PHI, DTH
  COMPLEX(KIND=dp) :: DOTVAL
  COMPLEX(KIND=dp), ALLOCATABLE :: POL(:,:), TOR(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(POL(3, NYLM), TOR(3, NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  |DOT(PVSH_POL,PVSH_TOR)|'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL PVSH_POL_ALL(POL, LMAX, THETA, PHI)
    CALL PVSH_TOR_ALL(TOR, LMAX, THETA, PHI)
    DO L = 1, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L, M)
        DOTVAL = DOT(POL(:,IDX), TOR(:,IDX))
        WRITE(OUTUNIT, *) L, M, THETA, PHI, ABS(DOTVAL)
        IF (ABS(DOTVAL) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(POL, TOR)
  END SUBROUTINE PVSH_POL_TOR_ORTHO


!> Verifies the orthogonal-rotation identity linking the standard and
!> polar VSH bases can be inverted exactly:
!> $$ \sqrt{\tfrac{\ell}{2\ell+1}}\,\mathbf{Y}_{\ell m}^{\ell+1} +
!>    \sqrt{\tfrac{\ell+1}{2\ell+1}}\,\mathbf{Y}_{\ell m}^{\ell-1} =
!>    \mathbf{Y}_{\ell m}^{(+1)}, $$
!> $$ -\sqrt{\tfrac{\ell+1}{2\ell+1}}\,\mathbf{Y}_{\ell m}^{\ell+1} +
!>    \sqrt{\tfrac{\ell}{2\ell+1}}\,\mathbf{Y}_{\ell m}^{\ell-1} =
!>    \mathbf{Y}_{\ell m}^{(-1)}, $$
!> i.e. rotating [[VSH_POL_UP]]/[[VSH_POL_DN]] back by the same angle used
!> to build them from [[PVSH_POL]]/[[PVSH_RAD]] recovers the originals
!> exactly.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTH Number of colatitude points swept over \( (0,\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   theta phi max_absdiff_pol max_absdiff_rad`, each the max over the 3
!>   vector components).
!> @param STATUS Output: 0 = pass, 1 = fail (either max difference
!>   \( >10^{-12} \)).
  SUBROUTINE VSH_POL_INVERSION(LMAX, NTH, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX, NTH, OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, C, ITH, NYLM, IDX
  REAL(KIND=dp) :: THETA, PHI, DTH, SC_POL, SC_RAD
  REAL(KIND=dp) :: DIFF_POL, DIFF_RAD
  COMPLEX(KIND=dp), ALLOCATABLE :: UP(:,:), DN(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: POL(:,:), RAD(:,:)
  NYLM = (LMAX+1)**2
  ALLOCATE(UP(3,NYLM), DN(3,NYLM), POL(3,NYLM), RAD(3,NYLM))
  STATUS = 0
  DTH = pi / (NTH + 1)
  PHI = pi / 4.d0
  WRITE(OUTUNIT,'(A)') &
      '# L  M  theta  phi  max_absdiff_pol  max_absdiff_rad'
  DO ITH = 1, NTH
    THETA = ITH * DTH
    CALL VSH_POL_UP_ALL(UP,  LMAX, THETA, PHI)
    CALL VSH_POL_DN_ALL(DN,  LMAX, THETA, PHI)
    CALL PVSH_POL_ALL(POL, LMAX, THETA, PHI)
    CALL PVSH_RAD_ALL(RAD, LMAX, THETA, PHI)
    DO L = 0, LMAX
      SC_POL = DSQRT(DBLE(L)/(2*L+1.d0))
      SC_RAD = DSQRT((L+1.d0)/(2*L+1.d0))
      DO M = -L, L
        IDX = YLM_INDEX(L, M)
        DIFF_POL = 0.d0
        DIFF_RAD = 0.d0
        DO C = 1, 3
          DIFF_POL = MAX(DIFF_POL, ABS( &
              SC_POL*UP(C,IDX) + SC_RAD*DN(C,IDX) - POL(C,IDX)))
          DIFF_RAD = MAX(DIFF_RAD, ABS( &
             -SC_RAD*UP(C,IDX) + SC_POL*DN(C,IDX) - RAD(C,IDX)))
        END DO
        WRITE(OUTUNIT, *) L, M, THETA, PHI, DIFF_POL, DIFF_RAD
        IF (DIFF_POL > 1.0E-12_dp .OR. DIFF_RAD > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(UP, DN, POL, RAD)
  END SUBROUTINE VSH_POL_INVERSION


END MODULE TESTS
