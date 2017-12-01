# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.0)
project(linalg C CXX Fortran)

# Define version information
set(linalg_VERSION_MAJOR 1)
set(linalg_VERSION_MINOR 4)
set(linalg_VERSION_PATCH 0)

# Define output directories, if undefined
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(STATUS "Library output directories undefined.  Using default directories.")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    # Debug Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug/mod)
    endif ()
  elseif (CMAKE_BUILD_TYPE MATCHES Release)
    # Release Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release/mod)
    endif ()
  else ()
    # Default Condition
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug/mod)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug/mod)
    endif ()
  endif ()
endif()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -Wno-c-binding-type")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    # ifort (untested)
    set(CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
endif()

# Define the output directories
set(ExternalProjectCMakeArgs
  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  -CMAKE_Fortran_MODULE_DIRECTORY=${CMAKE_Fortran_MODULE_DIRECTORY}
  -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

# Locate the files
add_subdirectory(src)

# ------------------------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------------------------
option(BUILD_EXAMPLES "Build examples?" OFF)
if (BUILD_EXAMPLES)
  # Inform the user we're building the examples
  message(STATUS "Building LINALG examples.")

  # LU Example
  add_executable(lu_example examples/linalg_lu_example.f90)
  target_link_libraries(lu_example linalg)

  # QR Example
  add_executable(qr_example examples/linalg_qr_example.f90)
  target_link_libraries(qr_example linalg)

  # Overdetermined Example
  add_executable(od_example examples/linalg_od_example.f90)
  target_link_libraries(od_example linalg)

  # Eigen Example
  add_executable(eigen_example examples/linalg_eigen_example.f90)
  target_link_libraries(eigen_example linalg)

# --------------------
  include_directories(${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_SOURCE_DIR}/src/external/ferror/include)

  # LU Example
  add_executable(lu_example_c examples/linalg_lu_c_example.c)
  target_link_libraries(lu_example_c linalg)

  # QR Example
  add_executable(qr_example_c examples/linalg_qr_c_example.c)
  target_link_libraries(qr_example_c linalg)

  # Overdetermined Example
  add_executable(od_example_c examples/linalg_od_c_example.c)
  target_link_libraries(od_example_c linalg)

  # Eigen Example
  add_executable(eigen_example_c examples/linalg_eigen_c_example.c)
  target_link_libraries(eigen_example_c linalg m)
endif()

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
set(linalg_test_sources tests/linalg_test.f90 tests/test_core.f90 
    tests/test_qr.f90 tests/test_svd.f90 tests/test_cholesky.f90 
    tests/test_eigen.f90 tests/test_mtx_inverse.f90 tests/test_misc.f90 
    tests/test_lu.f90 tests/test_sort.f90)
add_executable(linalg_test ${linalg_test_sources})
target_link_libraries(linalg_test linalg)
add_custom_command(OUTPUT linalg_tests DEPENDS ${linalg_test_sources}
  COMMAND linalg_test)
add_custom_target(run_linalg_tests ALL DEPENDS linalg_tests)

# C Test
include_directories(${PROJECT_SOURCE_DIR}/tests)
set(c_linalg_test_sources tests/linalg_c_test.c tests/c_test_core.c
  tests/c_test_misc.c tests/c_test_lu.c tests/c_test_qr.c tests/c_test_svd.c
  tests/c_test_inv.c tests/c_test_eigen.c tests/c_test_sort.c)
add_executable(linalg_c_test ${c_linalg_test_sources})
target_link_libraries(linalg_c_test linalg m)
add_custom_command(OUTPUT linalg_c_tests DEPENDS ${c_linalg_test_sources}
  COMMAND linalg_c_test)
add_custom_target(run_linalg_c_tests ALL DEPENDS linalg_c_tests)
