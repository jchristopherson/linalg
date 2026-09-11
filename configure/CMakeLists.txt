# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# By default, static library
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Export all symbols on Windows when building libraries
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Utilize the GNU installation structure
include(GNUInstallDirs)

# Locate the local include directory
if(NOT DEFINED CMAKE_INSTALL_MODULEDIR)
  set(
    CMAKE_INSTALL_MODULEDIR
    "${CMAKE_INSTALL_INCLUDEDIR}"
    CACHE
    STRING
    "Directory in prefix to install generated module files"
  )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

# Locate the dependencies.  This must happen before the package config file is
# generated as that file must know whether a dependency was built here.
include(FetchContent)

find_package(BLAS)
find_package(LAPACK)
find_package(arpackng)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "BLAS/LAPACK could not be found.  A reference version will be employed.")
  FetchContent_Declare(
      lapack
      GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
  )
  FetchContent_MakeAvailable(lapack)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
  set(BUILD_LAPACK_BLAS TRUE)

  # ARPACK-NG performs its own REQUIRED search unless these targets already
  # exist, which would either fail or link it against a different BLAS/LAPACK
  if (NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED GLOBAL)
    set_target_properties(BLAS::BLAS PROPERTIES
      INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
  endif()
  if (NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED GLOBAL)
    set_target_properties(LAPACK::LAPACK PROPERTIES
      INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
  endif()
else()
  set(BUILD_LAPACK_BLAS FALSE)
endif()

if (NOT arpackng_FOUND)
  message(STATUS "ARPACK could not be found.  A reference version will be employed.")
  FetchContent_Declare(
      arpack
      GIT_REPOSITORY "https://github.com/opencollab/arpack-ng.git"
  )
  FetchContent_MakeAvailable(arpack)
endif()
# Both the system package and arpack-ng itself export this target name
set(ARPACK_LIBRARIES arpack)

set(BLAS_LIBRARIES "${BLAS_LIBRARIES}" PARENT_SCOPE)
set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" PARENT_SCOPE)
set(ARPACK_LIBRARIES "${ARPACK_LIBRARIES}" PARENT_SCOPE)
set(BUILD_LAPACK_BLAS "${BUILD_LAPACK_BLAS}" PARENT_SCOPE)

# Export a pkg-config file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

# Export CMake package file
include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

if(BUILD_SHARED_LIBS OR PROJECT_VERSION_MAJOR EQUAL 0)
  # Due to the uncertain ABI compatibility of Fortran shared libraries
  # limit compatibility for dynamic linking to same minor version.
  set(COMPATIBILITY SameMinorVersion)
else()
  # Require API compatibility via semantic versioning for static linking.
  set(COMPATIBILITY SameMajorVersion)
endif()

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY ${COMPATIBILITY}
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
