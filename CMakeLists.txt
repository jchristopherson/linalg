# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.17)
project(
  linalg
  LANGUAGES Fortran C
  VERSION 1.6.1
)

# Utilize the GNU installation structure
include(GNUInstallDirs)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(BLAS)
find_package(LAPACK)
find_package(ferror 1.4.0 QUIET)
add_subdirectory(dependencies)

# print_all_variables()

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
target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})

# Installation

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

# Examples