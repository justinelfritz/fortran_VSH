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
!>
!> 4. **Wigner-D rotation infrastructure** ([[WIGNER_D_SMALL_IDENTITY]]
!>    through [[ROTATE_VSH_STD_SPECTRUM_GT]]): algebraic identities of
!>    the small-d/D-matrices themselves (identity angle, symmetries,
!>    unitarity, group composition), coefficient-array round-trip
!>    inversion for all three rotation entry points
!>    ([[ROTATE_SSH_ALL]]/[[ROTATE_PVSH_ALL]]/[[ROTATE_VSH_STD_ALL]]),
!>    and -- the strongest check in this file -- an end-to-end
!>    from-scratch ground truth
!>    ([[ROTATE_PVSH_SPECTRUM_GT]]/[[ROTATE_VSH_STD_SPECTRUM_GT]])
!>    that physically rotates a synthesized VSH field in Cartesian
!>    space (independent 3x3 rotation matrix, independent
!>    spherical/Cartesian reprojection, no Wigner-D machinery in the
!>    comparison path) and re-decomposes it by direct quadrature,
!>    confirming that coefficient-space rotation via the library's own
!>    engine agrees with physically rotating the field it represents.
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
  DOT, GWI, GWJ, &
  WIGNER_D_SMALL, WIGNER_D, &
  ROTATE_SSH_ALL, ROTATE_PVSH_ALL, ROTATE_VSH_STD_ALL
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
  SSH_ORTHO, PVSH_POL_TOR_ORTHO, VSH_POL_INVERSION, &
  ! Wigner-D / rotation validation
  WIGNER_D_SMALL_IDENTITY, WIGNER_D_SMALL_SYMMETRY, &
  WIGNER_D_IDENTITY, WIGNER_D_UNITARITY, WIGNER_D_COMPOSITION, &
  ROTATE_SSH_INVERSION, ROTATE_PVSH_INVERSION, ROTATE_VSH_STD_INVERSION, &
  ROTATE_PVSH_SPECTRUM_GT, ROTATE_VSH_STD_SPECTRUM_GT, &
  WIGNER_D_SPOTCHECK

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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NQUAD
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTH
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
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


! ── Rotation test scaffolding (module-internal only) ──────────────────────────
! Small helpers shared by the Wigner-D/rotation checks below. None of these
! reuse any rotation logic from [[VSH]] -- the whole point of
! [[ROTATE_PVSH_SPECTRUM_GT]]/[[ROTATE_VSH_STD_SPECTRUM_GT]] is an
! independent, from-scratch ground truth, so the Cartesian rotation matrix
! and the spherical/Cartesian reprojection here are built directly from
! their own definitions, not delegated to [[WIGNER_D]]/[[WIGNER_D_SMALL]].

!> Builds the Cartesian rotation matrix for the \( z\text{-}y\text{-}z \)
!> Euler convention used throughout [[VSH]]'s Wigner-D machinery,
!> \( R=R_z(\alpha)R_y(\beta)R_z(\gamma) \). Independent of
!> [[WIGNER_D]]/[[WIGNER_D_SMALL]] -- exists purely so
!> [[ROTATE_PVSH_SPECTRUM_GT]]/[[ROTATE_VSH_STD_SPECTRUM_GT]] can
!> physically rotate a Cartesian vector field without going anywhere near
!> the coefficient-space machinery being tested.
!>
!> @param ALPHA First Euler angle in radians.
!> @param BETA Second Euler angle in radians.
!> @param GAMMA Third Euler angle in radians.
!> Returns: The \( 3\times3 \) real rotation matrix \( R \).
  FUNCTION EULER_ROTATION_MATRIX(ALPHA, BETA, GAMMA) RESULT(R)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: ALPHA
  REAL(KIND=dp), INTENT(IN) :: BETA
  REAL(KIND=dp), INTENT(IN) :: GAMMA
  REAL(KIND=dp) :: CA, SA, CB, SB, CG, SG
  REAL(KIND=dp) :: RZ1(3,3), RY(3,3), RZ2(3,3), R(3,3)
  CA = DCOS(ALPHA); SA = DSIN(ALPHA)
  CB = DCOS(BETA);  SB = DSIN(BETA)
  CG = DCOS(GAMMA); SG = DSIN(GAMMA)
  RZ1 = RESHAPE([CA,SA,0.0_dp, -SA,CA,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
  RY  = RESHAPE([CB,0.0_dp,-SB, 0.0_dp,1.0_dp,0.0_dp, SB,0.0_dp,CB], [3,3])
  RZ2 = RESHAPE([CG,SG,0.0_dp, -SG,CG,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
  R = MATMUL(RZ1, MATMUL(RY, RZ2))
  END FUNCTION EULER_ROTATION_MATRIX

!> Local spherical unit vectors \( (\hat r,\hat\theta,\hat\phi) \) in
!> Cartesian components at a point \( (\theta,\phi) \) -- the same
!> convention as every VSH basis function's \( (r,\theta,\phi) \)
!> component ordering, used here only to move between that local basis
!> and global Cartesian coordinates for the ground-truth rotation checks.
!>
!> @param THETA Colatitude in radians.
!> @param PHI Longitude in radians.
!> @param RHAT Output radial unit vector, Cartesian.
!> @param THHAT Output colatitude unit vector, Cartesian.
!> @param PHHAT Output longitude unit vector, Cartesian.
  SUBROUTINE SPH_UNIT_VECTORS(THETA, PHI, RHAT, THHAT, PHHAT)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN)  :: THETA
  REAL(KIND=dp), INTENT(IN)  :: PHI
  REAL(KIND=dp), INTENT(OUT) :: RHAT(3)
  REAL(KIND=dp), INTENT(OUT) :: THHAT(3)
  REAL(KIND=dp), INTENT(OUT) :: PHHAT(3)
  REAL(KIND=dp) :: ST, CT, SP, CP
  ST = DSIN(THETA); CT = DCOS(THETA)
  SP = DSIN(PHI);   CP = DCOS(PHI)
  RHAT  = [ST*CP, ST*SP, CT]
  THHAT = [CT*CP, CT*SP, -ST]
  PHHAT = [-SP, CP, 0.0_dp]
  END SUBROUTINE SPH_UNIT_VECTORS

!> Inverse of [[SPH_UNIT_VECTORS]]'s point map: recovers
!> \( (\theta,\phi) \) from a (not-necessarily-unit) Cartesian vector,
!> with \( \phi \) wrapped to \( [0,2\pi) \) to match [[VSH]]'s own
!> convention.
!>
!> @param V Input Cartesian vector (normalized internally).
!> @param THETA Output colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Output longitude in radians, \( 0\le\phi<2\pi \).
  SUBROUTINE CART_TO_SPH(V, THETA, PHI)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN)  :: V(3)
  REAL(KIND=dp), INTENT(OUT) :: THETA
  REAL(KIND=dp), INTENT(OUT) :: PHI
  REAL(KIND=dp) :: VN(3), NORM
  NORM = DSQRT(V(1)**2 + V(2)**2 + V(3)**2)
  VN = V / NORM
  THETA = DACOS(MAX(-1.0_dp, MIN(1.0_dp, VN(3))))
  PHI = DATAN2(VN(2), VN(1))
  IF (PHI < 0.0_dp) PHI = PHI + 2.0_dp*pi
  END SUBROUTINE CART_TO_SPH

!> Applies a real \( 3\times3 \) matrix to a complex 3-vector by explicit
!> summation (deliberately not the `MATMUL` intrinsic, to avoid relying
!> on mixed real/complex argument promotion across compilers -- this
!> file's minimum-toolchain target is gfortran 9, see `README.md`).
!>
!> @param RMAT Real \( 3\times3 \) matrix.
!> @param V Complex 3-vector.
!> Returns: \( \texttt{RMAT}\cdot V \), complex 3-vector.
  FUNCTION ROTATE_CVEC3(RMAT, V) RESULT(RV)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: RMAT(3,3)
  COMPLEX(KIND=dp), INTENT(IN) :: V(3)
  COMPLEX(KIND=dp) :: RV(3)
  INTEGER(KIND=i4) :: I
  DO I = 1, 3
    RV(I) = RMAT(I,1)*V(1) + RMAT(I,2)*V(2) + RMAT(I,3)*V(3)
  END DO
  END FUNCTION ROTATE_CVEC3

!> \( (-1)^N \), computed via `MOD(ABS(N),2)` so negative `N` (common
!> throughout the \( m\)-index alternating-sign identities below) is
!> handled portably rather than relying on integer exponentiation to a
!> negative power.
!>
!> @param N Integer exponent (any sign).
!> Returns: \( +1 \) if `N` is even, \( -1 \) if `N` is odd.
  FUNCTION ALT_SIGN(N) RESULT(S)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: N
  REAL(KIND=dp) :: S
  IF (MOD(ABS(N), 2) .EQ. 0) THEN
    S = 1.0_dp
  ELSE
    S = -1.0_dp
  END IF
  END FUNCTION ALT_SIGN


! ── Wigner-D / rotation validation ─────────────────────────────────────────────

!> Verifies \( d^\ell_{m'm}(0)=\delta_{m'm} \) -- rotating by zero angle
!> about \( y \) is the identity -- for every \( 0\le\ell\le\ell_{max},\,
!> -\ell\le m',m\le\ell \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L MP
!>   M value expected(0_or_1)`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `|value-expected| >
!>   1e-12`).
  SUBROUTINE WIGNER_D_SMALL_IDENTITY(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, MP, M, EXPECTED
  REAL(KIND=dp) :: VAL
  STATUS = 0
  WRITE(OUTUNIT,'(A)') '# L  MP  M  value  expected'
  DO L = 0, LMAX
    DO MP = -L, L
      DO M = -L, L
        VAL = WIGNER_D_SMALL(L, MP, M, 0.0_dp)
        IF (MP .EQ. M) THEN
          EXPECTED = 1
        ELSE
          EXPECTED = 0
        END IF
        WRITE(OUTUNIT, *) L, MP, M, VAL, EXPECTED
        IF (ABS(VAL - EXPECTED) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_SMALL_IDENTITY


!> Verifies the small-d matrix's exact symmetries,
!> \( d^\ell_{m'm}(\beta)=(-1)^{m'-m}d^\ell_{m,m'}(\beta)=
!> d^\ell_{-m,-m'}(\beta) \), swept over five representative
!> \( \beta \) (including \( \beta=\pi \), where the additional special
!> value \( d^\ell_{m'm}(\pi)=(-1)^{\ell-m}\delta_{m',-m} \) is also
!> checked).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L MP
!>   M beta absdiff_flip absdiff_negate absdiff_pi_special`).
!> @param STATUS Output: 0 = pass, 1 = fail (some difference \( >10^{-12}
!>   \)).
  SUBROUTINE WIGNER_D_SMALL_SYMMETRY(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, MP, M, IB
  REAL(KIND=dp) :: BETA, D1, D2, D3
  REAL(KIND=dp) :: DIFF_FLIP, DIFF_NEG, DIFF_PI
  REAL(KIND=dp), DIMENSION(5) :: BETAS
  BETAS = [0.3_dp, 1.0_dp, 2.0_dp, pi/2.0_dp, pi]
  STATUS = 0
  WRITE(OUTUNIT,'(A)') &
      '# L  MP  M  beta  absdiff_flip  absdiff_negate  absdiff_pi_special'
  DO IB = 1, 5
    BETA = BETAS(IB)
    DO L = 0, LMAX
      DO MP = -L, L
        DO M = -L, L
          D1 = WIGNER_D_SMALL(L, MP, M, BETA)
          D2 = WIGNER_D_SMALL(L, M, MP, BETA) * ALT_SIGN(MP-M)
          D3 = WIGNER_D_SMALL(L, -M, -MP, BETA)
          DIFF_FLIP = ABS(D1 - D2)
          DIFF_NEG  = ABS(D1 - D3)
          DIFF_PI = 0.0_dp
          IF (ABS(BETA-pi) < 1.0E-12_dp) THEN
            IF (MP .EQ. -M) THEN
              DIFF_PI = ABS(D1 - ALT_SIGN(L-M))
            ELSE
              DIFF_PI = ABS(D1)
            END IF
          END IF
          WRITE(OUTUNIT, *) L, MP, M, BETA, DIFF_FLIP, DIFF_NEG, DIFF_PI
          IF (DIFF_FLIP > 1.0E-12_dp .OR. DIFF_NEG > 1.0E-12_dp .OR. &
              DIFF_PI > 1.0E-12_dp) STATUS = 1
        END DO
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_SMALL_SYMMETRY


!> Verifies \( D^\ell_{m'm}(0,0,0)=\delta_{m'm} \) -- the zero rotation is
!> the identity element of the full complex D-matrix -- for every
!> \( 0\le\ell\le\ell_{max} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L MP
!>   M re(value) im(value) expected(0_or_1)`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `|value-expected| >
!>   1e-12`).
  SUBROUTINE WIGNER_D_IDENTITY(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, MP, M, EXPECTED
  COMPLEX(KIND=dp) :: VAL
  STATUS = 0
  WRITE(OUTUNIT,'(A)') '# L  MP  M  re(value)  im(value)  expected'
  DO L = 0, LMAX
    DO MP = -L, L
      DO M = -L, L
        VAL = WIGNER_D(L, MP, M, 0.0_dp, 0.0_dp, 0.0_dp)
        IF (MP .EQ. M) THEN
          EXPECTED = 1
        ELSE
          EXPECTED = 0
        END IF
        WRITE(OUTUNIT, *) L, MP, M, DREAL(VAL), DIMAG(VAL), EXPECTED
        IF (ABS(VAL - EXPECTED) > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_IDENTITY


!> Verifies \( D^\ell(\alpha,\beta,\gamma) \) is unitary,
!> \( \sum_{m'}\overline{D^\ell_{m'm_1}}\,D^\ell_{m'm_2}=\delta_{m_1m_2}
!> \), for four representative \( (\alpha,\beta,\gamma) \) triples
!> (including \( \beta=\pi \)) and every \( 0\le\ell\le\ell_{max} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M1
!>   M2 alpha beta gamma absdiff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-12`).
  SUBROUTINE WIGNER_D_UNITARITY(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M1, M2, MP, IT, EXPECTED
  REAL(KIND=dp) :: ALPHA, BETA, GAMMA, DIFF
  REAL(KIND=dp), DIMENSION(4) :: ALPHAS, BETAS, GAMMAS
  COMPLEX(KIND=dp) :: ACC
  ALPHAS = [0.4_dp, 1.9_dp, pi/5.0_dp, 0.0_dp]
  BETAS  = [0.7_dp, 2.6_dp, pi/3.0_dp, pi]
  GAMMAS = [1.1_dp, 0.2_dp, pi/7.0_dp, 0.5_dp]
  STATUS = 0
  WRITE(OUTUNIT,'(A)') '# L  M1  M2  alpha  beta  gamma  absdiff'
  DO IT = 1, 4
    ALPHA = ALPHAS(IT); BETA = BETAS(IT); GAMMA = GAMMAS(IT)
    DO L = 0, LMAX
      DO M1 = -L, L
        DO M2 = -L, L
          ACC = DCMPLX(0.d0, 0.d0)
          DO MP = -L, L
            ACC = ACC + CONJG(WIGNER_D(L,MP,M1,ALPHA,BETA,GAMMA)) * &
                        WIGNER_D(L,MP,M2,ALPHA,BETA,GAMMA)
          END DO
          IF (M1 .EQ. M2) THEN
            EXPECTED = 1
          ELSE
            EXPECTED = 0
          END IF
          DIFF = ABS(ACC - EXPECTED)
          WRITE(OUTUNIT, *) L, M1, M2, ALPHA, BETA, GAMMA, DIFF
          IF (DIFF > 1.0E-12_dp) STATUS = 1
        END DO
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_UNITARITY


!> Verifies \( D^\ell \) satisfies the rotation-group composition law,
!> \( D^\ell(\text{combined})_{m'm}=\sum_k D^\ell(\text{first})_{m'k}\,
!> D^\ell(\text{second})_{km} \), for two special cases chosen so the
!> "combined" angle is known in closed form without needing a general
!> Euler-angle composition formula as a test dependency: pure-\(\alpha\)
!> rotations (diagonal; angles add directly) and pure-\(\beta\) rotations
!> (small-d; angles add directly, \( R_y(\beta_1)R_y(\beta_2)=
!> R_y(\beta_1+\beta_2) \)).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L MP
!>   M absdiff_alpha_composition absdiff_beta_composition`).
!> @param STATUS Output: 0 = pass, 1 = fail (some difference
!>   \( >10^{-11} \)).
  SUBROUTINE WIGNER_D_COMPOSITION(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  REAL(KIND=dp), PARAMETER :: ALPHA1 = 0.4_dp, ALPHA2 = 0.7_dp
  REAL(KIND=dp), PARAMETER :: BETA1 = 0.5_dp, BETA2 = 0.9_dp
  INTEGER(KIND=i4) :: L, MP, M, K
  COMPLEX(KIND=dp) :: ACC_A, DIRECT_A
  REAL(KIND=dp) :: ACC_B, DIRECT_B, DIFF_A, DIFF_B
  STATUS = 0
  WRITE(OUTUNIT,'(A)') &
      '# L  MP  M  absdiff_alpha_composition  absdiff_beta_composition'
  DO L = 0, LMAX
    DO MP = -L, L
      DO M = -L, L
        ACC_A = DCMPLX(0.d0, 0.d0)
        ACC_B = 0.0_dp
        DO K = -L, L
          ACC_A = ACC_A + WIGNER_D(L,MP,K,ALPHA1,0.0_dp,0.0_dp) * &
                          WIGNER_D(L,K,M,ALPHA2,0.0_dp,0.0_dp)
          ACC_B = ACC_B + WIGNER_D_SMALL(L,MP,K,BETA1) * &
                          WIGNER_D_SMALL(L,K,M,BETA2)
        END DO
        DIRECT_A = WIGNER_D(L,MP,M,ALPHA1+ALPHA2,0.0_dp,0.0_dp)
        DIRECT_B = WIGNER_D_SMALL(L,MP,M,BETA1+BETA2)
        DIFF_A = ABS(ACC_A - DIRECT_A)
        DIFF_B = ABS(ACC_B - DIRECT_B)
        WRITE(OUTUNIT, *) L, MP, M, DIFF_A, DIFF_B
        IF (DIFF_A > 1.0E-11_dp .OR. DIFF_B > 1.0E-11_dp) STATUS = 1
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_COMPOSITION


!> Verifies that rotating a full SSH coefficient array by
!> \( (\alpha,\beta,\gamma) \) and then by the inverse Euler triple
!> \( (-\gamma,-\beta,-\alpha) \) exactly recovers the input, for three
!> representative angle triples (a generic case, a near-\(\beta=0\)
!> case, and an exact \( \beta=\pi \) case). Direct structural analogue
!> of [[VSH_POL_INVERSION]], extended to a full forward/backward
!> round-trip through [[ROTATE_SSH_ALL]] rather than a single
!> fixed-angle relation.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns:
!>   `alpha beta gamma L M absdiff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-12`).
  SUBROUTINE ROTATE_SSH_INVERSION(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, IT, NYLM, IDX
  REAL(KIND=dp), DIMENSION(3) :: ALPHAS, BETAS, GAMMAS
  REAL(KIND=dp) :: ALPHA, BETA, GAMMA, DIFF
  COMPLEX(KIND=dp), ALLOCATABLE :: ALM_IN(:), ALM_MID(:), ALM_OUT(:)
  ALPHAS = [pi/5.0_dp, 0.1_dp, 0.0_dp]
  BETAS  = [pi/3.0_dp, 0.05_dp, pi]
  GAMMAS = [pi/7.0_dp, -1.3_dp, 0.3_dp]
  NYLM = (LMAX+1)**2
  ALLOCATE(ALM_IN(NYLM), ALM_MID(NYLM), ALM_OUT(NYLM))
  STATUS = 0
  WRITE(OUTUNIT,'(A)') '# alpha  beta  gamma  L  M  absdiff'
  DO IT = 1, 3
    ALPHA = ALPHAS(IT); BETA = BETAS(IT); GAMMA = GAMMAS(IT)
    DO L = 0, LMAX
      DO M = -L, L
        ALM_IN(YLM_INDEX(L,M)) = &
            DCMPLX(0.1d0*L+0.05d0*M, 0.2d0*L-0.03d0*M)
      END DO
    END DO
    CALL ROTATE_SSH_ALL(ALM_MID, ALM_IN, LMAX, ALPHA, BETA, GAMMA)
    CALL ROTATE_SSH_ALL(ALM_OUT, ALM_MID, LMAX, -GAMMA, -BETA, -ALPHA)
    DO L = 0, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L,M)
        DIFF = ABS(ALM_OUT(IDX) - ALM_IN(IDX))
        WRITE(OUTUNIT, *) ALPHA, BETA, GAMMA, L, M, DIFF
        IF (DIFF > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(ALM_IN, ALM_MID, ALM_OUT)
  END SUBROUTINE ROTATE_SSH_INVERSION


!> Same forward/backward round-trip as [[ROTATE_SSH_INVERSION]], applied
!> to [[ROTATE_PVSH_ALL]] with all three polar-basis families
!> ([[PVSH_RAD]]/[[PVSH_POL]]/[[PVSH_TOR]]) seeded simultaneously with
!> distinct coefficient formulas -- catches any accidental cross-family
!> mixing that a single-family test could miss, in addition to the basic
!> inversion identity.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns:
!>   `alpha beta gamma L M absdiff_rad absdiff_pol absdiff_tor`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-12`).
  SUBROUTINE ROTATE_PVSH_INVERSION(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, IT, NYLM, IDX
  REAL(KIND=dp), DIMENSION(3) :: ALPHAS, BETAS, GAMMAS
  REAL(KIND=dp) :: ALPHA, BETA, GAMMA
  REAL(KIND=dp) :: DIFF_RAD, DIFF_POL, DIFF_TOR
  COMPLEX(KIND=dp), ALLOCATABLE :: RAD_IN(:), POL_IN(:), TOR_IN(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: RAD_MID(:), POL_MID(:), TOR_MID(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: RAD_OUT(:), POL_OUT(:), TOR_OUT(:)
  ALPHAS = [pi/5.0_dp, 0.1_dp, 0.0_dp]
  BETAS  = [pi/3.0_dp, 0.05_dp, pi]
  GAMMAS = [pi/7.0_dp, -1.3_dp, 0.3_dp]
  NYLM = (LMAX+1)**2
  ALLOCATE(RAD_IN(NYLM), POL_IN(NYLM), TOR_IN(NYLM))
  ALLOCATE(RAD_MID(NYLM), POL_MID(NYLM), TOR_MID(NYLM))
  ALLOCATE(RAD_OUT(NYLM), POL_OUT(NYLM), TOR_OUT(NYLM))
  STATUS = 0
  WRITE(OUTUNIT,'(A)') &
      '# alpha  beta  gamma  L  M  absdiff_rad  absdiff_pol  absdiff_tor'
  DO IT = 1, 3
    ALPHA = ALPHAS(IT); BETA = BETAS(IT); GAMMA = GAMMAS(IT)
    DO L = 0, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L,M)
        RAD_IN(IDX) = DCMPLX(0.1d0*L+0.05d0*M,  0.2d0*L-0.03d0*M)
        POL_IN(IDX) = DCMPLX(0.07d0*L-0.02d0*M, 0.15d0*L+0.04d0*M)
        TOR_IN(IDX) = DCMPLX(0.12d0*L+0.06d0*M,-0.09d0*L+0.01d0*M)
      END DO
    END DO
    CALL ROTATE_PVSH_ALL(RAD_MID, POL_MID, TOR_MID, &
                          RAD_IN, POL_IN, TOR_IN, LMAX, ALPHA, BETA, GAMMA)
    CALL ROTATE_PVSH_ALL(RAD_OUT, POL_OUT, TOR_OUT, &
                          RAD_MID, POL_MID, TOR_MID, &
                          LMAX, -GAMMA, -BETA, -ALPHA)
    DO L = 0, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L,M)
        DIFF_RAD = ABS(RAD_OUT(IDX) - RAD_IN(IDX))
        DIFF_POL = ABS(POL_OUT(IDX) - POL_IN(IDX))
        DIFF_TOR = ABS(TOR_OUT(IDX) - TOR_IN(IDX))
        WRITE(OUTUNIT, *) ALPHA, BETA, GAMMA, L, M, &
            DIFF_RAD, DIFF_POL, DIFF_TOR
        IF (DIFF_RAD > 1.0E-12_dp .OR. DIFF_POL > 1.0E-12_dp .OR. &
            DIFF_TOR > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(RAD_IN, POL_IN, TOR_IN, RAD_MID, POL_MID, TOR_MID)
  DEALLOCATE(RAD_OUT, POL_OUT, TOR_OUT)
  END SUBROUTINE ROTATE_PVSH_INVERSION


!> Same forward/backward round-trip as [[ROTATE_SSH_INVERSION]], applied
!> to [[ROTATE_VSH_STD_ALL]] with all three standard-basis families
!> ([[VSH_TOR]]/[[VSH_POL_UP]]/[[VSH_POL_DN]]) seeded simultaneously with
!> distinct coefficient formulas.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param OUTUNIT Fortran unit number to write results to (columns:
!>   `alpha beta gamma L M absdiff_tor absdiff_up absdiff_dn`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-12`).
  SUBROUTINE ROTATE_VSH_STD_INVERSION(LMAX, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  INTEGER(KIND=i4) :: L, M, IT, NYLM, IDX
  REAL(KIND=dp), DIMENSION(3) :: ALPHAS, BETAS, GAMMAS
  REAL(KIND=dp) :: ALPHA, BETA, GAMMA
  REAL(KIND=dp) :: DIFF_TOR, DIFF_UP, DIFF_DN
  COMPLEX(KIND=dp), ALLOCATABLE :: TOR_IN(:), UP_IN(:), DN_IN(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: TOR_MID(:), UP_MID(:), DN_MID(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: TOR_OUT(:), UP_OUT(:), DN_OUT(:)
  ALPHAS = [pi/5.0_dp, 0.1_dp, 0.0_dp]
  BETAS  = [pi/3.0_dp, 0.05_dp, pi]
  GAMMAS = [pi/7.0_dp, -1.3_dp, 0.3_dp]
  NYLM = (LMAX+1)**2
  ALLOCATE(TOR_IN(NYLM), UP_IN(NYLM), DN_IN(NYLM))
  ALLOCATE(TOR_MID(NYLM), UP_MID(NYLM), DN_MID(NYLM))
  ALLOCATE(TOR_OUT(NYLM), UP_OUT(NYLM), DN_OUT(NYLM))
  STATUS = 0
  WRITE(OUTUNIT,'(A)') &
      '# alpha  beta  gamma  L  M  absdiff_tor  absdiff_up  absdiff_dn'
  DO IT = 1, 3
    ALPHA = ALPHAS(IT); BETA = BETAS(IT); GAMMA = GAMMAS(IT)
    DO L = 0, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L,M)
        TOR_IN(IDX) = DCMPLX(0.1d0*L+0.05d0*M,  0.2d0*L-0.03d0*M)
        UP_IN(IDX)  = DCMPLX(0.07d0*L-0.02d0*M, 0.15d0*L+0.04d0*M)
        DN_IN(IDX)  = DCMPLX(0.12d0*L+0.06d0*M,-0.09d0*L+0.01d0*M)
      END DO
    END DO
    CALL ROTATE_VSH_STD_ALL(TOR_MID, UP_MID, DN_MID, &
                             TOR_IN, UP_IN, DN_IN, LMAX, ALPHA, BETA, GAMMA)
    CALL ROTATE_VSH_STD_ALL(TOR_OUT, UP_OUT, DN_OUT, &
                             TOR_MID, UP_MID, DN_MID, &
                             LMAX, -GAMMA, -BETA, -ALPHA)
    DO L = 0, LMAX
      DO M = -L, L
        IDX = YLM_INDEX(L,M)
        DIFF_TOR = ABS(TOR_OUT(IDX) - TOR_IN(IDX))
        DIFF_UP  = ABS(UP_OUT(IDX)  - UP_IN(IDX))
        DIFF_DN  = ABS(DN_OUT(IDX)  - DN_IN(IDX))
        WRITE(OUTUNIT, *) ALPHA, BETA, GAMMA, L, M, &
            DIFF_TOR, DIFF_UP, DIFF_DN
        IF (DIFF_TOR > 1.0E-12_dp .OR. DIFF_UP > 1.0E-12_dp .OR. &
            DIFF_DN > 1.0E-12_dp) STATUS = 1
      END DO
    END DO
  END DO
  DEALLOCATE(TOR_IN, UP_IN, DN_IN, TOR_MID, UP_MID, DN_MID)
  DEALLOCATE(TOR_OUT, UP_OUT, DN_OUT)
  END SUBROUTINE ROTATE_VSH_STD_INVERSION


!> End-to-end ground-truth check for [[ROTATE_PVSH_ALL]]: seeds a
!> [[PVSH_TOR]]-only coefficient array, rotates it through the public
!> engine, and independently verifies the result by physically rotating
!> the Cartesian vector field the *input* coefficients represent --
!> \( \mathbf F'(\mathbf n)=R^{-1}\mathbf F(R\mathbf n) \) at each
!> quadrature point (the physical-rotation convention this identity was
!> empirically pinned to match [[WIGNER_D]]'s own sign convention; see
!> the derivation in the accompanying rotation-infrastructure design
!> notes) -- and re-decomposing it by direct numerical quadrature
!> (midpoint rule in \( \cos\theta \), uniform in \( \phi \), no
!> Wigner-D or [[ROTATE_SSH_ALL]] machinery anywhere in this path). Also
!> confirms exactly zero leakage into the untouched
!> [[PVSH_RAD]]/[[PVSH_POL]] families, the practical form of the
!> "no \( L/\lambda \) mixing under rotation" claim [[ROTATE_PVSH_ALL]]
!> is built on.
!>
!> @warning The `1e-3` tolerance reflects the quadrature's own \(
!>   O(1/\texttt{NTHETA}^2) \) truncation error (same midpoint-rule
!>   character as [[SSH_ORTHO]]), not the accuracy of the rotation
!>   engine itself -- at `NTHETA=200,NPHI=32` the observed error is
!>   \( \sim10^{-4} \), and the engine's own error (checked separately,
!>   not in this file, against the same ground truth at tighter
!>   quadrature) is \( \sim10^{-11} \).
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \). Kept small
!>   (\( \lesssim5 \)) -- runtime scales as
!>   \( O(\texttt{NTHETA}\cdot\texttt{NPHI}\cdot\ell_{max}^2) \).
!> @param NTHETA Number of midpoint-rule quadrature points over
!>   \( \cos\theta\in[-1,1] \).
!> @param NPHI Number of uniform quadrature points over
!>   \( \phi\in[0,2\pi) \) (trapezoid rule, spectrally accurate for the
!>   periodic integrand here, so a modest value suffices).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   re(engine) im(engine) re(ground_truth) im(ground_truth) absdiff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-3`, or
!>   any cross-family leakage `> 1e-12`).
  SUBROUTINE ROTATE_PVSH_SPECTRUM_GT(LMAX, NTHETA, NPHI, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTHETA
  INTEGER(KIND=i4), INTENT(IN) :: NPHI
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  REAL(KIND=dp), PARAMETER :: ALPHA = 3.14159265358979324_dp/5.0_dp
  REAL(KIND=dp), PARAMETER :: BETA  = 3.14159265358979324_dp/3.0_dp
  REAL(KIND=dp), PARAMETER :: GAMMA = 3.14159265358979324_dp/7.0_dp
  INTEGER(KIND=i4) :: L, M, ITH, IPH, NYLM, IDX
  REAL(KIND=dp) :: THETA, PHI, THETA1, PHI1, U, DU, DPHI, MAXDIFF, DIFF
  REAL(KIND=dp) :: RMAT(3,3), RMAT_T(3,3)
  REAL(KIND=dp) :: RHAT0(3), THHAT0(3), PHHAT0(3)
  REAL(KIND=dp) :: RHAT1(3), THHAT1(3), PHHAT1(3), FORWARD(3)
  COMPLEX(KIND=dp) :: F0_CART(3), FPRIME_CART(3), F_TH, F_PH
  COMPLEX(KIND=dp), ALLOCATABLE :: C_RAD_IN(:), C_POL_IN(:), C_TOR_IN(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: C_RAD_OUT(:), C_POL_OUT(:), C_TOR_OUT(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: GT_TOR(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: BASIS0(:,:), BASIS1(:,:)

  NYLM = (LMAX+1)**2
  ALLOCATE(C_RAD_IN(NYLM), C_POL_IN(NYLM), C_TOR_IN(NYLM))
  ALLOCATE(C_RAD_OUT(NYLM), C_POL_OUT(NYLM), C_TOR_OUT(NYLM))
  ALLOCATE(GT_TOR(NYLM), BASIS0(3,NYLM), BASIS1(3,NYLM))
  STATUS = 0

  C_RAD_IN = DCMPLX(0.d0, 0.d0)
  C_POL_IN = DCMPLX(0.d0, 0.d0)
  DO L = 0, LMAX
    DO M = -L, L
      C_TOR_IN(YLM_INDEX(L,M)) = DCMPLX(0.1d0*L+0.05d0*M, 0.2d0*L-0.03d0*M)
    END DO
  END DO

  CALL ROTATE_PVSH_ALL(C_RAD_OUT, C_POL_OUT, C_TOR_OUT, &
                        C_RAD_IN, C_POL_IN, C_TOR_IN, &
                        LMAX, ALPHA, BETA, GAMMA)

  RMAT = EULER_ROTATION_MATRIX(ALPHA, BETA, GAMMA)
  RMAT_T = TRANSPOSE(RMAT)
  GT_TOR = DCMPLX(0.d0, 0.d0)
  DU = 2.0_dp / NTHETA
  DPHI = 2.0_dp*pi / NPHI

  DO ITH = 1, NTHETA
    U = -1.0_dp + (ITH - 0.5_dp) * DU
    THETA = DACOS(U)
    DO IPH = 1, NPHI
      PHI = (IPH - 1) * DPHI
      CALL SPH_UNIT_VECTORS(THETA, PHI, RHAT0, THHAT0, PHHAT0)
      FORWARD = MATMUL(RMAT, RHAT0)
      CALL CART_TO_SPH(FORWARD, THETA1, PHI1)
      CALL SPH_UNIT_VECTORS(THETA1, PHI1, RHAT1, THHAT1, PHHAT1)
      CALL PVSH_TOR_ALL(BASIS1, LMAX, THETA1, PHI1)
      F0_CART = DCMPLX(0.d0, 0.d0)
      DO IDX = 1, NYLM
        F0_CART = F0_CART + C_TOR_IN(IDX) * &
            (BASIS1(2,IDX)*THHAT1 + BASIS1(3,IDX)*PHHAT1)
      END DO
      FPRIME_CART = ROTATE_CVEC3(RMAT_T, F0_CART)
      F_TH = FPRIME_CART(1)*THHAT0(1) + FPRIME_CART(2)*THHAT0(2) + &
             FPRIME_CART(3)*THHAT0(3)
      F_PH = FPRIME_CART(1)*PHHAT0(1) + FPRIME_CART(2)*PHHAT0(2) + &
             FPRIME_CART(3)*PHHAT0(3)
      CALL PVSH_TOR_ALL(BASIS0, LMAX, THETA, PHI)
      DO IDX = 1, NYLM
        GT_TOR(IDX) = GT_TOR(IDX) + &
            (F_TH*CONJG(BASIS0(2,IDX)) + F_PH*CONJG(BASIS0(3,IDX))) * &
            DU * DPHI
      END DO
    END DO
  END DO

  MAXDIFF = 0.0_dp
  WRITE(OUTUNIT,'(A)') &
      '# L  M  re(engine)  im(engine)  re(ground_truth)'// &
      '  im(ground_truth)  absdiff'
  DO L = 0, LMAX
    DO M = -L, L
      IDX = YLM_INDEX(L,M)
      DIFF = ABS(C_TOR_OUT(IDX) - GT_TOR(IDX))
      MAXDIFF = MAX(MAXDIFF, DIFF)
      WRITE(OUTUNIT, *) L, M, DREAL(C_TOR_OUT(IDX)), DIMAG(C_TOR_OUT(IDX)), &
          DREAL(GT_TOR(IDX)), DIMAG(GT_TOR(IDX)), DIFF
    END DO
  END DO
  IF (MAXDIFF > 1.0E-3_dp) STATUS = 1

  DO IDX = 1, NYLM
    IF (ABS(C_RAD_OUT(IDX)) > 1.0E-12_dp .OR. &
        ABS(C_POL_OUT(IDX)) > 1.0E-12_dp) STATUS = 1
  END DO

  DEALLOCATE(C_RAD_IN, C_POL_IN, C_TOR_IN, C_RAD_OUT, C_POL_OUT, C_TOR_OUT)
  DEALLOCATE(GT_TOR, BASIS0, BASIS1)
  END SUBROUTINE ROTATE_PVSH_SPECTRUM_GT


!> Companion to [[ROTATE_PVSH_SPECTRUM_GT]] for [[ROTATE_VSH_STD_ALL]]:
!> same from-scratch ground-truth methodology, but seeds a
!> [[VSH_POL_UP]]-only coefficient array instead of a toroidal one, so
!> the check additionally exercises a *nonzero radial component*
!> (the \( J=\ell+1 \) fixed linear combination of [[PVSH_RAD]] and
!> [[PVSH_POL]]) rather than a purely tangential field. Also confirms
!> exactly zero leakage into the untouched [[VSH_TOR]]/[[VSH_POL_DN]]
!> families.
!>
!> @warning Same quadrature-truncation caveat as
!>   [[ROTATE_PVSH_SPECTRUM_GT]]'s `@warning` applies here.
!>
!> @param LMAX Maximum degree tested, \( \ell_{max}\ge0 \).
!> @param NTHETA Number of midpoint-rule quadrature points over
!>   \( \cos\theta\in[-1,1] \).
!> @param NPHI Number of uniform quadrature points over
!>   \( \phi\in[0,2\pi) \).
!> @param OUTUNIT Fortran unit number to write results to (columns: `L M
!>   re(engine) im(engine) re(ground_truth) im(ground_truth) absdiff`).
!> @param STATUS Output: 0 = pass, 1 = fail (some `absdiff > 1e-3`, or
!>   any cross-family leakage `> 1e-12`).
  SUBROUTINE ROTATE_VSH_STD_SPECTRUM_GT(LMAX, NTHETA, NPHI, OUTUNIT, STATUS)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  INTEGER(KIND=i4), INTENT(IN) :: NTHETA
  INTEGER(KIND=i4), INTENT(IN) :: NPHI
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), INTENT(OUT) :: STATUS
  REAL(KIND=dp), PARAMETER :: ALPHA = 3.14159265358979324_dp/5.0_dp
  REAL(KIND=dp), PARAMETER :: BETA  = 3.14159265358979324_dp/3.0_dp
  REAL(KIND=dp), PARAMETER :: GAMMA = 3.14159265358979324_dp/7.0_dp
  INTEGER(KIND=i4) :: L, M, ITH, IPH, NYLM, IDX
  REAL(KIND=dp) :: THETA, PHI, THETA1, PHI1, U, DU, DPHI, MAXDIFF, DIFF
  REAL(KIND=dp) :: RMAT(3,3), RMAT_T(3,3)
  REAL(KIND=dp) :: RHAT0(3), THHAT0(3), PHHAT0(3)
  REAL(KIND=dp) :: RHAT1(3), THHAT1(3), PHHAT1(3), FORWARD(3)
  COMPLEX(KIND=dp) :: F0_CART(3), FPRIME_CART(3), F_R, F_TH, F_PH
  COMPLEX(KIND=dp), ALLOCATABLE :: C_TOR_IN(:), C_UP_IN(:), C_DN_IN(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: C_TOR_OUT(:), C_UP_OUT(:), C_DN_OUT(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: GT_UP(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: BASIS0(:,:), BASIS1(:,:)

  NYLM = (LMAX+1)**2
  ALLOCATE(C_TOR_IN(NYLM), C_UP_IN(NYLM), C_DN_IN(NYLM))
  ALLOCATE(C_TOR_OUT(NYLM), C_UP_OUT(NYLM), C_DN_OUT(NYLM))
  ALLOCATE(GT_UP(NYLM), BASIS0(3,NYLM), BASIS1(3,NYLM))
  STATUS = 0

  C_TOR_IN = DCMPLX(0.d0, 0.d0)
  C_DN_IN  = DCMPLX(0.d0, 0.d0)
  DO L = 0, LMAX
    DO M = -L, L
      C_UP_IN(YLM_INDEX(L,M)) = &
          DCMPLX(0.07d0*L+0.04d0*M, 0.15d0*L-0.02d0*M)
    END DO
  END DO

  CALL ROTATE_VSH_STD_ALL(C_TOR_OUT, C_UP_OUT, C_DN_OUT, &
                           C_TOR_IN, C_UP_IN, C_DN_IN, &
                           LMAX, ALPHA, BETA, GAMMA)

  RMAT = EULER_ROTATION_MATRIX(ALPHA, BETA, GAMMA)
  RMAT_T = TRANSPOSE(RMAT)
  GT_UP = DCMPLX(0.d0, 0.d0)
  DU = 2.0_dp / NTHETA
  DPHI = 2.0_dp*pi / NPHI

  DO ITH = 1, NTHETA
    U = -1.0_dp + (ITH - 0.5_dp) * DU
    THETA = DACOS(U)
    DO IPH = 1, NPHI
      PHI = (IPH - 1) * DPHI
      CALL SPH_UNIT_VECTORS(THETA, PHI, RHAT0, THHAT0, PHHAT0)
      FORWARD = MATMUL(RMAT, RHAT0)
      CALL CART_TO_SPH(FORWARD, THETA1, PHI1)
      CALL SPH_UNIT_VECTORS(THETA1, PHI1, RHAT1, THHAT1, PHHAT1)
      CALL VSH_POL_UP_ALL(BASIS1, LMAX, THETA1, PHI1)
      F0_CART = DCMPLX(0.d0, 0.d0)
      DO IDX = 1, NYLM
        F0_CART = F0_CART + C_UP_IN(IDX) * &
            (BASIS1(1,IDX)*RHAT1 + BASIS1(2,IDX)*THHAT1 + &
             BASIS1(3,IDX)*PHHAT1)
      END DO
      FPRIME_CART = ROTATE_CVEC3(RMAT_T, F0_CART)
      F_R  = FPRIME_CART(1)*RHAT0(1)  + FPRIME_CART(2)*RHAT0(2)  + &
             FPRIME_CART(3)*RHAT0(3)
      F_TH = FPRIME_CART(1)*THHAT0(1) + FPRIME_CART(2)*THHAT0(2) + &
             FPRIME_CART(3)*THHAT0(3)
      F_PH = FPRIME_CART(1)*PHHAT0(1) + FPRIME_CART(2)*PHHAT0(2) + &
             FPRIME_CART(3)*PHHAT0(3)
      CALL VSH_POL_UP_ALL(BASIS0, LMAX, THETA, PHI)
      DO IDX = 1, NYLM
        GT_UP(IDX) = GT_UP(IDX) + &
            (F_R*CONJG(BASIS0(1,IDX)) + F_TH*CONJG(BASIS0(2,IDX)) + &
             F_PH*CONJG(BASIS0(3,IDX))) * DU * DPHI
      END DO
    END DO
  END DO

  MAXDIFF = 0.0_dp
  WRITE(OUTUNIT,'(A)') &
      '# L  M  re(engine)  im(engine)  re(ground_truth)'// &
      '  im(ground_truth)  absdiff'
  DO L = 0, LMAX
    DO M = -L, L
      IDX = YLM_INDEX(L,M)
      DIFF = ABS(C_UP_OUT(IDX) - GT_UP(IDX))
      MAXDIFF = MAX(MAXDIFF, DIFF)
      WRITE(OUTUNIT, *) L, M, DREAL(C_UP_OUT(IDX)), DIMAG(C_UP_OUT(IDX)), &
          DREAL(GT_UP(IDX)), DIMAG(GT_UP(IDX)), DIFF
    END DO
  END DO
  IF (MAXDIFF > 1.0E-3_dp) STATUS = 1

  DO IDX = 1, NYLM
    IF (ABS(C_TOR_OUT(IDX)) > 1.0E-12_dp .OR. &
        ABS(C_DN_OUT(IDX)) > 1.0E-12_dp) STATUS = 1
  END DO

  DEALLOCATE(C_TOR_IN, C_UP_IN, C_DN_IN, C_TOR_OUT, C_UP_OUT, C_DN_OUT)
  DEALLOCATE(GT_UP, BASIS0, BASIS1)
  END SUBROUTINE ROTATE_VSH_STD_SPECTRUM_GT


!> Writes a representative table of \( D^\ell_{m'm}(\alpha,\beta,\gamma) \)
!> values (\( 0\le\ell\le4 \), every \( m',m \), at three angle triples)
!> to `OUTUNIT` for cross-checking against `py/sympy_reference.py` --
!> `sympy.physics.quantum.spin.Rotation.D`, a genuinely independent
!> library/algorithm, not just a higher-precision re-evaluation of the
!> same closed form. Data-generation only, no pass/fail `STATUS` -- the
!> comparison itself is manual/investigative (run
!> `py/compute_validation.py` after the suite), matching
!> `py/mpmath_reference.py`'s own self-check convention rather than a
!> ctest gate.
!>
!> @param OUTUNIT Fortran unit number to write results to (columns: `L
!>   MP M alpha beta gamma re(D) im(D)`).
  SUBROUTINE WIGNER_D_SPOTCHECK(OUTUNIT)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: OUTUNIT
  INTEGER(KIND=i4), PARAMETER :: LMAX = 4
  INTEGER(KIND=i4) :: L, MP, M, IT
  REAL(KIND=dp) :: ALPHA, BETA, GAMMA
  REAL(KIND=dp), DIMENSION(3) :: ALPHAS, BETAS, GAMMAS
  COMPLEX(KIND=dp) :: D
  ALPHAS = [0.4_dp, pi/5.0_dp, 1.9_dp]
  BETAS  = [0.9_dp, pi/3.0_dp, 2.6_dp]
  GAMMAS = [1.3_dp, pi/7.0_dp, 0.2_dp]
  WRITE(OUTUNIT,'(A)') '# L  MP  M  alpha  beta  gamma  re(D)  im(D)'
  DO IT = 1, 3
    ALPHA = ALPHAS(IT); BETA = BETAS(IT); GAMMA = GAMMAS(IT)
    DO L = 0, LMAX
      DO MP = -L, L
        DO M = -L, L
          D = WIGNER_D(L, MP, M, ALPHA, BETA, GAMMA)
          WRITE(OUTUNIT, *) L, MP, M, ALPHA, BETA, GAMMA, &
              DREAL(D), DIMAG(D)
        END DO
      END DO
    END DO
  END DO
  END SUBROUTINE WIGNER_D_SPOTCHECK


END MODULE TESTS
