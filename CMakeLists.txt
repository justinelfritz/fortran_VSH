cmake_minimum_required(VERSION 3.14)
project(FORTVSH
  VERSION   1.0.0
  LANGUAGES Fortran C
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
add_compile_options(-O3)

option(VSH_BUILD_SHARED    "Build shared library (libvsh.so) in addition to static" OFF)
option(VSH_INSTALL_TESTS   "Install the vsh_test executable"                         OFF)
option(VSH_BUILD_BENCHMARK "Build the vsh_benchmark executable (batch vs. naive-loop timing)" OFF)
option(VSH_BUILD_EXAMPLES  "Build the poloidal-dipole worked-example executables"    OFF)
option(VSH_BUILD_CONVERGENCE "Build the vsh_convergence executable (batch-consistency accuracy sweep up to Lmax=300)" OFF)
option(VSH_BUILD_STABILITY "Build the vsh_stability executable (ASSOC_LEGENDRE_NORM_ALL stability scan over l,m,theta)" OFF)
option(VSH_BUILD_ROTATION_STABILITY "Build the vsh_rotation_stability executable (Wigner-D small-d stability scan over l,beta)" OFF)

# ── Generate version Fortran module ──────────────────────────────────────────
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/vsh_version.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/src/vsh_version.f90"
  @ONLY
)

# ── Library object compilation (compile once for both static + shared) ────────
add_library(vsh_objects OBJECT
  src/kinds.f90
  src/globals.f90
  src/vsh.f90
  "${CMAKE_CURRENT_BINARY_DIR}/src/vsh_version.f90"
)
set_target_properties(vsh_objects PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY  "${CMAKE_CURRENT_BINARY_DIR}/mod"
)

# ── Static library (always built) ─────────────────────────────────────────────
add_library(vsh STATIC $<TARGET_OBJECTS:vsh_objects>)
add_library(FORTVSH::vsh ALIAS vsh)
# Propagate the module directory to consumers (don't link vsh_objects as a
# PUBLIC dependency — OBJECT libraries cannot be exported)
target_include_directories(vsh INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/vsh>
)
set_target_properties(vsh PROPERTIES OUTPUT_NAME vsh)

# ── Shared library (opt-in) ───────────────────────────────────────────────────
if(VSH_BUILD_SHARED)
  add_library(vsh_shared SHARED $<TARGET_OBJECTS:vsh_objects>)
  add_library(FORTVSH::vsh_shared ALIAS vsh_shared)
  target_include_directories(vsh_shared INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/vsh>
  )
  set_target_properties(vsh_shared PROPERTIES
    OUTPUT_NAME vsh
    VERSION     ${PROJECT_VERSION}
    SOVERSION   ${PROJECT_VERSION_MAJOR}
  )
endif()

# ── Test executable ───────────────────────────────────────────────────────────
add_executable(vsh_test src/tests.f90 src/main.f90)
target_link_libraries(vsh_test PRIVATE vsh)
set_target_properties(vsh_test PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME create_validation_dir
  COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/validation"
)
set_tests_properties(create_validation_dir PROPERTIES
  FIXTURES_SETUP vsh_validation_dirs
)

add_test(NAME vsh_validation
  COMMAND vsh_test
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
set_tests_properties(vsh_validation PROPERTIES
  FIXTURES_REQUIRED vsh_validation_dirs
)

option(VSH_REGRESSION_TEST
  "Compare output .dat files against committed reference files in validation/reference/"
  OFF
)
if(VSH_REGRESSION_TEST)
  find_package(Python3 REQUIRED)
  add_test(NAME vsh_regression
    COMMAND ${Python3_EXECUTABLE}
      "${CMAKE_SOURCE_DIR}/cmake/check_regression.py"
      "${CMAKE_SOURCE_DIR}/validation/reference"
      "${CMAKE_SOURCE_DIR}/validation"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  )
  set_tests_properties(vsh_regression PROPERTIES
    DEPENDS vsh_validation
    FIXTURES_REQUIRED vsh_validation_dirs
  )
endif()

# ── Benchmark executable (opt-in; not part of ctest) ──────────────────────────
# Timing is machine/compiler-dependent, so it is a manually-run build target
# rather than a pass/fail test. Run the resulting executable from the repo
# root so it can create ./benchmark/*.dat, mirroring how vsh_test writes to
# ./validation/.
if(VSH_BUILD_BENCHMARK)
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/benchmark")
  add_executable(vsh_benchmark src/benchmark.f90 src/bench_main.f90)
  target_link_libraries(vsh_benchmark PRIVATE vsh)
  set_target_properties(vsh_benchmark PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif()

# ── Convergence executable (opt-in; not part of ctest) ────────────────────────
# Calls the same TESTS batch-consistency routines main.f90 uses for the
# correctness gate, but at Lmax=300 instead of the small fixed values that
# gate matter -- error legitimately grows with Lmax, so this is a figure
# data source, not a pass/fail check. Compiles src/tests.f90 a second time
# (same pattern as vsh_test) rather than restructuring TESTS into its own
# library target. Writes to ./convergence/, kept separate from
# ./validation/ so it can never interact with VSH_REGRESSION_TEST.
if(VSH_BUILD_CONVERGENCE)
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/convergence")
  add_executable(vsh_convergence src/tests.f90 src/convergence_main.f90)
  target_link_libraries(vsh_convergence PRIVATE vsh)
  set_target_properties(vsh_convergence PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj_convergence"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif()

# ── Stability-scan executable (opt-in; not part of ctest) ─────────────────────
# Maps where ASSOC_LEGENDRE_NORM_ALL's Holmes-Featherstone recurrence stays
# accurate vs. blows up across (l, m, theta) -- see STABILITY_FINDINGS.md.
# Standalone (doesn't need TESTS), writes to ./stability/.
if(VSH_BUILD_STABILITY)
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/stability")
  add_executable(vsh_stability src/stability_main.f90)
  target_link_libraries(vsh_stability PRIVATE vsh)
  set_target_properties(vsh_stability PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj_stability"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  # Companion scan for the unnormalized Bonnet recurrence (ASSOC_LEGENDRE/
  # _ALL and everything single-mode built on it) -- separate PROGRAM unit,
  # same opt-in flag/output directory.
  add_executable(vsh_unnorm_stability src/unnorm_stability_main.f90)
  target_link_libraries(vsh_unnorm_stability PRIVATE vsh)
  set_target_properties(vsh_unnorm_stability PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj_unnorm_stability"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif()

# ── Rotation (Wigner-D) stability-scan executable (opt-in; not part of
# ctest) ────────────────────────────────────────────────────────────────────
# Separate from VSH_BUILD_STABILITY: that option's two executables are about
# the Legendre recurrences; WIGNER_D_SMALL is a different closed-form family
# (Wigner's explicit finite sum, same technique as CGCOEFF). See
# src/rotation_stability_main.f90 for the exact |d|<=1 bound this scan uses.
if(VSH_BUILD_ROTATION_STABILITY)
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/rotation_stability")
  add_executable(vsh_rotation_stability src/rotation_stability_main.f90)
  target_link_libraries(vsh_rotation_stability PRIVATE vsh)
  set_target_properties(vsh_rotation_stability PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj_rotation_stability"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif()

# ── Worked-example executables (opt-in; not part of ctest) ────────────────────
# Three fully independent poloidal-dipole demonstrations (no shared module
# between them, matching the manuscript's "kept separated" worked-example
# structure). Physics demonstrations, not the correctness gate, so — like
# the benchmark target — they are a manually-run build option rather than
# a ctest, though each prints its own PASS/FAIL self-check. Run from the
# repo root so each can create its own ./examples/<name>/ subdirectory.
if(VSH_BUILD_EXAMPLES)
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/dipole_synthesis")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/dipole_uniform_sphere")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/dipole_flux_expulsion")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/vsh_decomposition")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/gwi_gwj_sweep")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/vsh_decomposition_tor")

  foreach(EXAMPLE_NAME dipole_synthesis dipole_uniform_sphere dipole_flux_expulsion
                       vsh_decomposition gwi_gwj_sweep vsh_decomposition_tor)
    add_executable(${EXAMPLE_NAME} src/examples/${EXAMPLE_NAME}.f90)
    target_link_libraries(${EXAMPLE_NAME} PRIVATE vsh)
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/obj"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
  endforeach()
endif()

# ── FortranCInterface: name-mangling header for future C/Python bindings ──────
include(FortranCInterface)
FortranCInterface_HEADER(
  "${CMAKE_CURRENT_BINARY_DIR}/FortranCInterface.h"
  MACRO_NAMESPACE "VSH_"
)

# ── pkg-config ────────────────────────────────────────────────────────────────
configure_file(fortvsh.pc.in "${CMAKE_CURRENT_BINARY_DIR}/fortvsh.pc" @ONLY)

# ── CMake package config (enables find_package(FORTVSH)) ─────────────────────
configure_package_config_file(
  "cmake/fortvsh-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/fortvsh-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FORTVSH"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/fortvsh-config-version.cmake"
  VERSION     ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ── Install rules ─────────────────────────────────────────────────────────────
install(TARGETS vsh
  EXPORT  fortvsh-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
if(VSH_BUILD_SHARED)
  install(TARGETS vsh_shared
    EXPORT  fortvsh-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/vsh"
  FILES_MATCHING PATTERN "*.mod"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/fortvsh.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)
install(EXPORT fortvsh-targets
  FILE      fortvsh-targets.cmake
  NAMESPACE FORTVSH::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FORTVSH"
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/fortvsh-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/fortvsh-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FORTVSH"
)
if(VSH_INSTALL_TESTS)
  install(TARGETS vsh_test RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
