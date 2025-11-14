cmake_minimum_required(VERSION 3.24)
project(
  linalg
  LANGUAGES Fortran C
  VERSION 1.8.6
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)

# Build the C API?
option(BUILD_LINALG_C_API "Build LINALG C API?" OFF)

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)

# Installation Instructions
install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (${BUILD_TESTING})
  message(STATUS "Building LINALG tests.")
  enable_testing()
  add_subdirectory(tests)
endif()

# Examples
option(BUILD_LINALG_EXAMPLES "Build LINALG examples.")
message(STATUS "Build LINALG examples: ${BUILD_LINALG_EXAMPLES}")
if (${BUILD_LINALG_EXAMPLES})
  add_subdirectory(examples)
endif()