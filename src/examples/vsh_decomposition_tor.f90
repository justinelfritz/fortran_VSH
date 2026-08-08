!> Worked example (2f): pure-toroidal VSH spectral decomposition,
!> reconstruction, and convergence -- companion to `vsh_decomposition.f90`
!> (example 2d), which decomposes onto the radial family. This example
!> decomposes onto the *toroidal* family only ([[PVSH_TOR]], coefficients
!> `a_tor`) instead.
!>
!> [[PVSH_RAD]](l,m) = \( (Y_\ell^m,0,0) \), so any scalar pattern assigned
!> directly to the radial slot is automatically pure-radial (`a_pol`,
!> `a_tor` vanish identically) -- no extra care needed, which is why
!> `vsh_decomposition.f90`'s `FIELD_GAUSSIAN` could just assign a bump
!> shape straight to the radial component. [[PVSH_TOR]] does not have this
!> property: it has *both* \( \theta \) and \( \phi \) components in
!> general (\( \mathrm{[[PVSH_TOR]]}=-i\,\mathrm{[[L_SSH]]}/
!> \sqrt{\ell(\ell+1)} \), nonzero in both components for \( m\ne0 \)), so
!> naively assigning a pattern to \( B_\phi \) alone with \( B_r=B_\theta=0
!> \) would leak into `a_pol` too. The fix used here: for *any* smooth
!> scalar generating function \( T(\theta,\phi) \) (not just a single
!> \( Y_\ell^m \)), \( \mathbf{B}=\hat r\times\nabla_\perp T \) is
!> *exactly* pure-toroidal -- \( a_{rad}=a_{pol}=0 \) identically -- since
!> \( \hat r\times\nabla_\perp Y_\ell^m \) is proportional to
!> \( \mathrm{[[PVSH_TOR]]}(\ell,m) \) for every \( (\ell,m) \), and by
!> linearity that holds for any \( T \) expressible in the \( Y_\ell^m \)
!> basis. Concretely \( B_\theta=-(1/\sin\theta)\,\partial_\phi T \),
!> \( B_\phi=\partial_\theta T \) -- the same relationship [[L_SSH]]
!> already implements per-mode, applied here to an arbitrary function.
!>
!> **Part A**: the same eight-mode exact-recovery test as example 2d
!> (\( (\ell,m)=(1,{-}1),(1,0),(1,1),(2,{-}2),(2,{-}1),(2,0),(2,1),(2,2) \)),
!> but built entirely from [[PVSH_TOR]] calls (unit coefficient each) --
!> `a_tor` should recover 1 at those eight modes and ~0 elsewhere, *and*
!> `a_rad`/`a_pol` should be ~0 at every mode (the direct check that this
!> field really is purely toroidal, not just that the known coefficients
!> come back right).
!>
!> **Part B**: a "quadrupole" axisymmetric toroidal background
!> (\( A_{quad}\,\mathrm{[[PVSH_TOR]]}(2,0,\theta,\phi) \), exact by
!> construction) superimposed with a broadband toroidal field built from
!> an equator-centered Gaussian bump. The background already vanishes at
!> both poles with no help needed: its \( \phi \)-component is
!> \( \propto\sin\theta\cos\theta \), zero at \( \theta=0,\pi \) for any
!> axisymmetric toroidal mode. The bump does need help. A generating
!> function of the form \( \exp(\cdots)\,(1+c\sin\theta\cos\phi) \) (the
!> shape used in an earlier version of this example, and in example 2d's
!> `FIELD_GAUSSIAN`) keeps \( T \) itself single-valued at the poles, but
!> not its *gradient*: \( \sin\theta\sim\theta \) near \( \theta=0 \), so
!> \( B_\theta=-(1/\sin\theta)\partial_\phi T \) comes out finite but
!> \( \phi \)-dependent right at the pole (physically fine for genuine
!> \( m=1 \) content, but not what "vanish at the poles" asks for here).
!> Forcing the *derivatives* to vanish, not just \( T \), needs the
!> \( \phi \)-dependent term to go to zero *quadratically* in \( \theta \)
!> near either pole -- hence \( \sin^2\theta \) rather than \( \sin\theta \)
!> below. The generating function used here,
!> $$ T(\theta,\phi) = \exp\!\left(-\frac{(\theta-\pi/2)^2}{2\sigma^2}
!>    \right)\sin^2\theta\,(1+0.5\cos\phi), $$
!> centers the Gaussian envelope at the equator (\( \theta=\pi/2 \)) as
!> requested, and its closed-form partial derivatives (verified against
!> both a finite-difference check and direct evaluation at the poles)
!> give
!> $$ B_\theta = 0.5\,e^{-(\theta-\pi/2)^2/2\sigma^2}\sin\theta\sin\phi,
!> $$
!> $$ B_\phi = e^{-(\theta-\pi/2)^2/2\sigma^2}\left[-\frac{\theta-\pi/2}
!>    {\sigma^2}\sin^2\theta + \sin2\theta\right](1+0.5\cos\phi), $$
!> both *exactly* zero at \( \theta=0,\pi \) -- every term carries an
!> explicit \( \sin\theta \) or \( \sin^2\theta \) or \( \sin2\theta \)
!> factor, each individually zero there, so the vanishing holds
!> identically, not just numerically. Decomposed/reconstructed at
!> increasing \( \ell_{max} \), same convergence check as example 2d, plus
!> a purity check that `a_rad`/`a_pol` stay ~0 throughout the sweep.
PROGRAM VSH_DECOMPOSITION_TOR
USE KINDS,   ONLY: dp, i4
USE GLOBALS, ONLY: pi
USE VSH,     ONLY: SHGLQ, PVSH_TOR, &
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
REAL(KIND=dp),    PARAMETER :: A_QUAD = 1.0_dp

INTEGER(KIND=i4) :: N_FAIL

N_FAIL = 0

CALL RUN_PART_A(N_FAIL)
CALL RUN_PART_B(N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - vsh_decomposition_tor"
END IF

CONTAINS

!> Target field for Part A: equal unit weight on [[PVSH_TOR]] alone at
!> each of the eight known modes; zero elsewhere. Purely toroidal by
!> construction (each term is a direct [[PVSH_TOR]] call).
!>
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: The target vector field \( \mathbf{F}(\theta,\phi) \), complex
!>   3-vector.
  FUNCTION FIELD_KNOWN_TOR(THETA, PHI) RESULT(F)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: F
  INTEGER(KIND=i4) :: K
  F = DCMPLX(0.d0, 0.d0)
  DO K = 1, N_KNOWN
    F = F + PVSH_TOR(L_KNOWN(K), M_KNOWN(K), THETA, PHI)
  ENDDO
  END FUNCTION FIELD_KNOWN_TOR

!> Target field for Part B: an exact axisymmetric \( (\ell,m)=(2,0) \)
!> "quadrupole" toroidal background plus a broadband toroidal field built
!> from \( \mathbf{B}=\hat r\times\nabla_\perp T \) applied to an
!> equator-centered Gaussian bump generating function (see module
!> docstring for the closed-form derivatives and why \( \sin^2\theta \),
!> not \( \sin\theta \), is needed to force \( B_\theta \)/\( B_\phi \) to
!> vanish exactly at both poles) -- purely toroidal by construction in
!> both pieces.
!>
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: The target vector field \( \mathbf{F}(\theta,\phi) \), complex
!>   3-vector.
  FUNCTION FIELD_QUAD_BUMP(THETA, PHI) RESULT(F)
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: F
  REAL(KIND=dp) :: ENV, SIN2TH, B_THETA, B_PHI
  COMPLEX(KIND=dp), DIMENSION(3) :: QUAD

  QUAD = A_QUAD * PVSH_TOR(2, 0, THETA, PHI)

  ENV = DEXP(-(THETA-0.5d0*pi)**2/(2.d0*SIGMA_B**2))
  SIN2TH = DSIN(THETA)**2
  B_THETA = 0.5d0 * ENV * DSIN(THETA) * DSIN(PHI)
  B_PHI   = ENV * ( -((THETA-0.5d0*pi)/SIGMA_B**2)*SIN2TH + &
                     DSIN(2.d0*THETA) ) * (1.d0 + 0.5d0*DCOS(PHI))

  F(1) = QUAD(1)
  F(2) = QUAD(2) + DCMPLX(B_THETA, 0.d0)
  F(3) = QUAD(3) + DCMPLX(B_PHI, 0.d0)
  END FUNCTION FIELD_QUAD_BUMP

!> Decomposes `FIELD_FN` into \( (a_{rad},a_{pol},a_{tor}) \) for
!> \( 0\le\ell\le\ell_{max},\,-\ell\le m\le\ell \), via 2D
!> Gauss-Legendre(\( \theta \)) times uniform-trapezoidal(\( \phi \))
!> quadrature. Identical to `vsh_decomposition.f90`'s `DECOMPOSE` --
!> duplicated rather than shared, matching this repository's convention of
!> fully independent worked-example programs (see `CMakeLists.txt`).
!>
!> @param FIELD_FN Target field function, `(theta,phi) -> COMPLEX(3)`.
!> @param LMAX Maximum degree to decompose, \( \ell_{max}\ge0 \).
!> @param LMAX_QUAD Theta-quadrature order ([[SHGLQ]] uses `LMAX_QUAD+1`
!>   nodes); must scale with `LMAX` for an accurate result.
!> @param NPHI Number of uniform phi-quadrature points.
!> @param A_RAD Output radial coefficients, size `(LMAX+1)**2`, indexed by
!>   [[YLM_INDEX]](l,m).
!> @param A_POL Output poloidal coefficients, same size/indexing.
!> @param A_TOR Output toroidal coefficients, same size/indexing.
  SUBROUTINE DECOMPOSE(FIELD_FN, LMAX, LMAX_QUAD, NPHI, A_RAD, A_POL, A_TOR)
  IMPLICIT NONE
  INTERFACE
    FUNCTION FIELD_FN(THETA, PHI) RESULT(F)
      IMPORT :: dp
      REAL(KIND=dp), INTENT(IN) :: THETA, PHI
      COMPLEX(KIND=dp), DIMENSION(3) :: F
    END FUNCTION FIELD_FN
  END INTERFACE
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX
  INTEGER(KIND=i4), INTENT(IN)  :: LMAX_QUAD
  INTEGER(KIND=i4), INTENT(IN)  :: NPHI
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

!> Reconstructs the field at \( (\theta,\phi) \) from coefficients up to
!> `LMAX`. Identical to `vsh_decomposition.f90`'s `RECONSTRUCT`.
!>
!> @param LMAX Maximum degree retained, \( \ell_{max}\ge0 \) (must match
!>   the coefficient arrays' size).
!> @param A_RAD Radial coefficients from [[DECOMPOSE]], size
!>   `(LMAX+1)**2`.
!> @param A_POL Poloidal coefficients from [[DECOMPOSE]], same size.
!> @param A_TOR Toroidal coefficients from [[DECOMPOSE]], same size.
!> @param THETA Colatitude in radians, \( 0\le\theta\le\pi \).
!> @param PHI Longitude in radians, \( 0\le\phi<2\pi \).
!> Returns: Reconstructed field \( \mathbf{F}(\theta,\phi) \), complex
!>   3-vector.
  FUNCTION RECONSTRUCT(LMAX, A_RAD, A_POL, A_TOR, THETA, PHI) RESULT(F)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(IN) :: LMAX
  COMPLEX(KIND=dp), INTENT(IN) :: A_RAD((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(IN) :: A_POL((LMAX+1)**2)
  COMPLEX(KIND=dp), INTENT(IN) :: A_TOR((LMAX+1)**2)
  REAL(KIND=dp), INTENT(IN) :: THETA
  REAL(KIND=dp), INTENT(IN) :: PHI
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

!> Part A: [[DECOMPOSE]] [[FIELD_KNOWN_TOR]], confirm exact recovery of
!> the eight known unit `a_tor` coefficients and near-zero everywhere
!> else, *and* confirm `a_rad`/`a_pol` are near-zero at every mode (the
!> direct purity check). Writes
!> `examples/vsh_decomposition_tor/known_modes.dat`.
!>
!> @param N_FAIL Running failure counter, incremented if either the max
!>   `a_tor` coefficient error or the max `a_rad`/`a_pol` purity leakage
!>   exceeds `TOL_A`.
  SUBROUTINE RUN_PART_A(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: LMAX_QUAD_A = 30
  INTEGER(KIND=i4), PARAMETER :: NPHI_A = 32
  COMPLEX(KIND=dp) :: A_RAD((LMAX_A+1)**2), A_POL((LMAX_A+1)**2), &
                       A_TOR((LMAX_A+1)**2)
  REAL(KIND=dp) :: EXPECTED, MAX_ERR, MAX_PURITY
  INTEGER(KIND=i4) :: L, M, IDX, K

  CALL DECOMPOSE(FIELD_KNOWN_TOR, LMAX_A, LMAX_QUAD_A, NPHI_A, &
                 A_RAD, A_POL, A_TOR)

  OPEN(UNIT=73, FILE="./examples/vsh_decomposition_tor/known_modes.dat")
  WRITE(73,'(A)') '# L  M  Re(a_rad)  Im(a_rad)  Re(a_pol)  Im(a_pol)  '// &
                  'Re(a_tor)  Im(a_tor)  expected'
  MAX_ERR = 0.0_dp
  MAX_PURITY = 0.0_dp
  DO L = 0, LMAX_A
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      EXPECTED = 0.0_dp
      DO K = 1, N_KNOWN
        IF (L_KNOWN(K) == L .AND. M_KNOWN(K) == M) EXPECTED = 1.0_dp
      ENDDO
      MAX_ERR = MAX(MAX_ERR, ABS(A_TOR(IDX)-EXPECTED))
      MAX_PURITY = MAX(MAX_PURITY, ABS(A_RAD(IDX)), ABS(A_POL(IDX)))
      WRITE(73,*) L, M, DBLE(A_RAD(IDX)), AIMAG(A_RAD(IDX)), &
                  DBLE(A_POL(IDX)), AIMAG(A_POL(IDX)), &
                  DBLE(A_TOR(IDX)), AIMAG(A_TOR(IDX)), EXPECTED
    ENDDO
  ENDDO
  CLOSE(73)

  IF (MAX_ERR > TOL_A) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  known_mode_recovery  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  known_mode_recovery  max_err=", MAX_ERR
  END IF

  IF (MAX_PURITY > TOL_A) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  toroidal_purity  max_leakage=", MAX_PURITY
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  toroidal_purity  max_leakage=", MAX_PURITY
  END IF
  END SUBROUTINE RUN_PART_A

!> Part B: [[DECOMPOSE]]+[[RECONSTRUCT]] [[FIELD_QUAD_BUMP]] at increasing
!> \( \ell_{max} \) (quadrature order scaled with each: `LMAX_QUAD =
!> 3*LMAX_TEST+30`, `NPHI = 4*LMAX_TEST+32` -- a higher floor than example
!> 2d needs, since the equator-centered bump's sharper angular structure
!> otherwise leaves a small quadrature-truncation leak in `a_rad`/`a_pol`
!> at the smallest `LMAX_TEST`, shrinking with quadrature order rather
!> than indicating an impure construction), tracking max reconstruction
!> error against the true closed-form field on an independent evaluation
!> grid, plus max `a_rad`/`a_pol` leakage across the sweep. Self-check:
!> error at the largest \( \ell_{max} \) must be much smaller than at the
!> smallest (genuine convergence). Writes
!> `examples/vsh_decomposition_tor/convergence.dat`.
!>
!> @note The evaluation grid still avoids landing exactly on the poles
!>   (`THETA = (ITH+0.5)*pi/NEVAL`, not `ITH*pi/(NEVAL-1)` as example 2d
!>   uses), a holdover from an earlier version of this example where it
!>   was load-bearing: [[PVSH_TOR_ALL]]/[[PVSH_POL_ALL]] guard against
!>   dividing by \( \sin\theta=0 \) by returning exactly zero for every
!>   \( m\ne0 \) basis function at \( \theta=0,\pi \), and that earlier
!>   bump had a nonzero, perfectly well-defined true limit right at the
!>   pole (its \( 1/\sin\theta \) cancelled analytically before being
!>   coded up), so comparing exactly at the pole produced a hard plateau
!>   at reconstruction error \( \approx0.5 \) that never shrank with
!>   \( \ell_{max} \) -- not an unphysical field, but a coordinate-singular
!>   comparison the library's own safety guard could never match. The
!>   current [[FIELD_QUAD_BUMP]] no longer has that problem -- its true
!>   value *is* exactly zero at both poles (see module docstring), matching
!>   what [[PVSH_TOR_ALL]]/[[PVSH_POL_ALL]]'s guard already returns there --
!>   but the off-pole grid is kept anyway as ordinary numerical hygiene
!>   near a coordinate singularity.
!>
!> @param N_FAIL Running failure counter, incremented if the error at the
!>   largest \( \ell_{max} \) is not far smaller than at the smallest, or
!>   if `a_rad`/`a_pol` leakage exceeds `TOL_A` anywhere in the sweep.
  SUBROUTINE RUN_PART_B(N_FAIL)
  IMPLICIT NONE
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  INTEGER(KIND=i4), PARAMETER :: NEVAL = 40
  REAL(KIND=dp),    PARAMETER :: CONVERGENCE_RATIO_TOL = 0.01_dp
  INTEGER(KIND=i4) :: IB, LMAX_TEST, LMAX_QUAD, NPHI, ITH, IPH
  COMPLEX(KIND=dp), ALLOCATABLE :: A_RAD(:), A_POL(:), A_TOR(:)
  COMPLEX(KIND=dp), DIMENSION(3) :: F_TRUE, F_RECON
  REAL(KIND=dp) :: THETA, PHI, MAX_ERR, FIRST_ERR, LAST_ERR, MAX_PURITY

  MAX_PURITY = 0.0_dp

  OPEN(UNIT=74, FILE="./examples/vsh_decomposition_tor/convergence.dat")
  WRITE(74,'(A)') '# Lmax  max_reconstruction_error'

  DO IB = 1, N_LMAX_B
    LMAX_TEST = LMAX_B_LIST(IB)
    LMAX_QUAD = 3*LMAX_TEST + 30
    NPHI = 4*LMAX_TEST + 32

    ALLOCATE(A_RAD((LMAX_TEST+1)**2), A_POL((LMAX_TEST+1)**2), &
             A_TOR((LMAX_TEST+1)**2))
    CALL DECOMPOSE(FIELD_QUAD_BUMP, LMAX_TEST, LMAX_QUAD, NPHI, &
                   A_RAD, A_POL, A_TOR)

    MAX_PURITY = MAX(MAX_PURITY, MAXVAL(ABS(A_RAD)), MAXVAL(ABS(A_POL)))

    MAX_ERR = 0.0_dp
    DO ITH = 0, NEVAL-1
      THETA = (ITH+0.5d0)*pi/NEVAL
      DO IPH = 0, NEVAL-1
        PHI = IPH*2.d0*pi/NEVAL
        F_TRUE  = FIELD_QUAD_BUMP(THETA, PHI)
        F_RECON = RECONSTRUCT(LMAX_TEST, A_RAD, A_POL, A_TOR, THETA, PHI)
        MAX_ERR = MAX(MAX_ERR, ABS(F_RECON(1)-F_TRUE(1)), &
                      ABS(F_RECON(2)-F_TRUE(2)), ABS(F_RECON(3)-F_TRUE(3)))
      ENDDO
    ENDDO

    IF (IB == 1) FIRST_ERR = MAX_ERR
    IF (IB == N_LMAX_B) LAST_ERR = MAX_ERR

    WRITE(74,*) LMAX_TEST, MAX_ERR
    WRITE(*,'(A,I0,A,ES10.3)') "  Lmax=", LMAX_TEST, &
      "  max_reconstruction_error=", MAX_ERR

    DEALLOCATE(A_RAD, A_POL, A_TOR)
  ENDDO
  CLOSE(74)

  IF (LAST_ERR > CONVERGENCE_RATIO_TOL*FIRST_ERR) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3,A,ES10.3)') "FAIL  spectral_convergence  first=", &
      FIRST_ERR, "  last=", LAST_ERR
  ELSE
    WRITE(*,'(A,ES10.3,A,ES10.3)') "PASS  spectral_convergence  first=", &
      FIRST_ERR, "  last=", LAST_ERR
  END IF

  IF (MAX_PURITY > TOL_A) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  toroidal_purity_sweep  max_leakage=", MAX_PURITY
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  toroidal_purity_sweep  max_leakage=", MAX_PURITY
  END IF
  END SUBROUTINE RUN_PART_B

END PROGRAM VSH_DECOMPOSITION_TOR
