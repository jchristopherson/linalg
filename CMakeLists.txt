# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
  linalg
  LANGUAGES Fortran C
  VERSION 1.7.4
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Build the C API?
option(BUILD_C_API "Build C API?" OFF)

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(BLAS)
find_package(LAPACK)
add_subdirectory(dependencies)

# Source
add_subdirectory(src)
add_fortran_library(
  ${PROJECT_NAME}
  ${PROJECT_INCLUDE_DIR}
  ${CMAKE_INSTALL_INCLUDEDIR}
  ${PROJECT_VERSION}
  ${PROJECT_VERSION_MAJOR}
  ${LINALG_SOURCES}
)
target_link_libraries(
  ${PROJECT_NAME} 
  ${BLAS_LIBRARIES} 
  ${LAPACK_LIBRARIES}
)
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})

# Installation
add_subdirectory(install)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(tests_c)
endif()

# Examples
option(BUILD_LINALG_EXAMPLES "Build LINALG examples")
message(STATUS "Build LINALG examples: ${BUILD_LINALG_EXAMPLES}")
if (BUILD_LINALG_EXAMPLES)
  add_subdirectory(examples)
endif()