include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get FORTRAN_TEST_HELPER
add_subdirectory(fortran_test_helper)

# Fortran API tests source code
set(linalg_test_sources
    linalg_test.f90
    test_core.f90
    test_qr.f90
    test_svd.f90
    test_cholesky.f90
    test_eigen.f90
    test_mtx_inverse.f90
    test_misc.f90
    test_lu.f90
    test_sort.f90
    test_lq.f90
)

# Build the Fortran API tests
add_executable(linalg_test ${linalg_test_sources})
link_library(linalg_test linalg ${PROJECT_INCLUDE_DIR})
link_library(linalg_test ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME linalg_test
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:linalg_test>
)
