# Master CMAKE Build Script
cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
project(linalg)

# Include the LAPACK Directory - If a shared build
if (BUILD_SHARED_LIBS)
    add_subdirectory(${CMAKE_SOURCE_DIR}/LAPACK)
endif()

# Specify output directories
if (CMAKE_BUILD_TYPE MATCHES Debug)
  # Debug Build
  if (BUILD_SHARED_LIBS)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin/mod)
  else()
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib/mod)
  endif ()
elseif (CMAKE_BUILD_TYPE MATCHES Release)
  # Release Build
  if (BUILD_SHARED_LIBS)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/bin)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/bin)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Release/bin/mod)
  else()
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Release/lib)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Release/lib/mod)
  endif ()
else ()
  # Default Condition
  if (BUILD_SHARED_LIBS)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/bin/mod)
  else()
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/Debug/lib/mod)
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

# Specify source files for the various builds
if (BUILD_SHARED_LIBS)
    file(GLOB_RECURSE sources src/*.f90 src/*f src/QRUPDATE/*f ferror/src/*f90)
else()
    # Include everything into the static library
    file(GLOB_RECURSE sources src/*.f90 src/*f src/QRUPDATE/*f 
        LAPACK/BLAS/src/*f LAPACK/src/*f ferror/src/*f90)
endif()

# Remove testing items from the library build
set(testing_sources ${sources})
list(REMOVE_ITEM sources src/test.f90 src/test_core.f90)

# If BUILD_SHARED_LIBS is defined, then a shared library (DLL) will be built;
# else, a static library will be built.
add_library(linalg ${sources})

# Build a test executable, only for debug builds
if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_executable(linalg_test ${testing_sources})
endif()

if (BUILD_SHARED_LIBS)
    # Establish Dependicies - Only for a shared library build
    target_link_libraries(linalg lapack)
endif()
