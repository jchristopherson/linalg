@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

# The exported targets reference these dependencies by their unqualified target
# names, so the packages defining them must be loaded first
if(NOT TARGET arpack)
  find_dependency(arpackng)
endif()

if(@BUILD_LAPACK_BLAS@ AND NOT TARGET LAPACK::lapack)
  # The reference version was built and installed alongside LINALG.  Its own
  # config file loads lapack-targets.cmake but never blas-targets.cmake, which
  # leaves LAPACK::lapack referencing an undefined BLAS::blas, so load both.
  file(GLOB _linalg_lapack_dirs
    "${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_LIBDIR@/cmake/lapack-*")
  foreach(_linalg_dir IN LISTS _linalg_lapack_dirs)
    foreach(_linalg_file blas-targets.cmake lapack-targets.cmake)
      if(EXISTS "${_linalg_dir}/${_linalg_file}")
        include("${_linalg_dir}/${_linalg_file}")
      endif()
    endforeach()
  endforeach()
  unset(_linalg_file)
  unset(_linalg_dir)
  unset(_linalg_lapack_dirs)
endif()

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()