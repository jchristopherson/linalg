# Define the QRUPDATE source files
set(qrupdate_sources
    external/qrupdate/dch1dn.f
    external/qrupdate/dch1up.f
    external/qrupdate/dlup1up.f
    external/qrupdate/dqhqr.f
    external/qrupdate/dqr1up.f
    external/qrupdate/dqrot.f
    external/qrupdate/dqrqh.f
    external/qrupdate/dqrtv1.f
)

# Define the LINALG source files
set(linalg_sources
    linalg_core.f90
    linalg_constants.f90
    linalg_eigen.f90
    linalg_factor.f90
    linalg_solve.f90
    linalg_sorting.f90
    linalg_basic.f90
    linalg_immutable.f90
)
list(APPEND linalg_sources ${qrupdate_sources})

# Find LAPACK
option(USE_EXISTING_LAPACK "Use an existing LAPACK build?" OFF)
if (USE_EXISTING_LAPACK)
  find_package(lapack)
else()
  message(STATUS "Using the default LAPACK implementation.")
  add_subdirectory(external/lapack)
endif()

# Find FERROR
option(USE_EXISTING_FERROR "Use an existing FERROR build?" OFF)
if (USE_EXISTING_FERROR)
    find_package(ferror 1.3.0)
else()
    message(STATUS "Using the default FERROR implementation.")
    add_subdirectory(external/ferror)
    set(ferror_INCLUDE_DIRS src/external/ferror/include)
endif()
include_directories(${ferror_INCLUDE_DIRS})

# Build the library
add_library(linalg ${linalg_sources})
target_link_libraries(linalg blas lapack ferror)

# https://stackoverflow.com/questions/30987311/exporting-an-imported-library

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET linalg PROPERTY VERSION ${linalg_VERSION})
set_property(TARGET linalg PROPERTY SOVERSION ${linalg_VERSION_PATCH})
set_property(TARGET linalg PROPERTY INTERFACE_linalg_MAJOR_VERSION ${linalg_VERSION_MINOR})
set_property(TARGET linalg APPEND PROPERTY COMPATIBLE_INTERFACE_STRING linalg_VERSION_MAJOR)

# Locate the "include" directories
set(linalg_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS linalg ${BLAS_LIBRARY}
    EXPORT linalgTargets
    RUNTIME DESTINATION linalg/bin
    LIBRARY DESTINATION linalg/lib
    ARCHIVE DESTINATION linalg/lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION linalg/include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}/linalg
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION linalg/doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/linalgConfigVersion.cmake"
    VERSION ${linalg_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT linalgTargets
    FILE "${CMAKE_BINARY_DIR}/linalgTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/linalgConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/linalgConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation linalg/lib/cmake/linalg)
install(
    EXPORT linalgTargets
    FILE linalgTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/linalgConfig.cmake"
        "${CMAKE_BINARY_DIR}/linalgConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)
