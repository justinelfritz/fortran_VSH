Fortran module `vsh` containing functions for Legendre polynomials and derivatives,
scalar spherical harmonics and angular derivatives, vector spherical harmonics and 
angular derivatives, and related functions. 

Independent modules for global declarations and datatype definitions are used
to enable greater portability and convenience. These are `globals` and `kinds`, 
respectively.

Please review the LaTeX whitepaper (in `/tex/`) for important information regarding 
normalization choices and other standards that propagate through the numerical 
routines. 

---

To build a test program:
* Ensure working directory contains makefile, `src/` containing all `*.f` and corresponding
`*.mod` files, an empty `obj/` directory, and an empty `out/` directory.
* Update makefile with your desired program name
* Update `main.f` to perform the desired test computations
* Run `make`

WIP - 17 July 2024