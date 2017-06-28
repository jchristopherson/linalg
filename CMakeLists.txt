# Master CMAKE Build Script
cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
project(linalg)

# Define version information
set(linalg_VERSION_MAJOR 1)
set(linalg_VERSION_MINOR 0)
set(linalg_VERSION_PATCH 0)
set(linalg_VERSION_TWEAK 1)

# Specify output directories
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

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  # Options: https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Specify source files
if (BUILD_SHARED_LIBS)
    file(GLOB_RECURSE sources src/*.f90 src/*f src/QRUPDATE/*f ferror/src/*f90)
else()
    # Include everything into the static library
    file(GLOB_RECURSE sources src/*.f90 src/*f src/QRUPDATE/*f 
        BLAS/src/*f LAPACK/src/*f ferror/src/*f90)
endif()

# Straighten out dependencies
# if (BUILD_SHARED_LIBS)
#   # See if the user wishes to use the reference BLAS
#   option(use_ref_blas "Use reference BLAS?" ON)
#   if (use_ref_blas)
#     # Add the BLAS source files to the list
#     file(GLOB_RECURSE blas_sources BLAS/src*f)
#     list(APPEND sources ${blas_sources})
#   else()
#     # TO DO: Ask user where to find BLAS library
#   endif()

#   # See if the user wishes to use the supplied LAPACK
#   option(use_lapack "Use supplied LAPACK?" ON)
#   if (use_lapack)
#     # Add the LAPACK source files to the list
#     file(GLOB_RECURSE lapack_sources LAPACK/src/*f)
#     list(APPEND sources ${lapack_sources})
#   else()
#     # TO DO: Ask user where to find LAPACK library
#   endif()
# endif()
if (BUILD_SHARED_LIBS)
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
  # add_external(lapack "")
  # add_dependencies(linalg lapack)
  # add_dependencies(linalg_test lapack)
endif ()

# Define the testing files
set(testing_sources ${sources})
list(APPEND testing_sources tests/test.f90 tests/test_core.f90 tests/test_qr.f90 
    tests/test_svd.f90 tests/test_cholesky.f90 tests/test_eigen.f90 
    tests/test_mtx_inverse.f90 tests/test_misc.f90 tests/test_lu.f90)

# If BUILD_SHARED_LIBS is defined, then a shared library (DLL) will be built;
# else, a static library will be built.
add_library(linalg ${sources})
target_link_libraries(linalg lapack)

# Build a test executable
add_executable(linalg_test ${testing_sources})
target_link_libraries(linalg_test lapack)

