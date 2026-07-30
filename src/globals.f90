!> Physical/mathematical constants shared across the library: \( \pi \)
!> and the imaginary unit \( i \), both at [[KINDS:dp]] precision so every
!> expression built from them (spherical harmonics, VSH, coupling
!> coefficients) stays in double precision without repeated literal
!> casts.
MODULE GLOBALS
  USE KINDS, ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pi, j
  !> \( \pi \), double precision.
  REAL(KIND=dp), PARAMETER :: pi=3.14159265358979324d0
  !> Imaginary unit \( i=\sqrt{-1} \), double-precision complex.
  COMPLEX(KIND=dp), PARAMETER :: j = (0,1)
END MODULE GLOBALS
