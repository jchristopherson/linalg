include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

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
add_test(
    NAME linalg_test
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:linalg_test>
)

# C API tests
if (${BUILD_C_API})
    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${ferror_INCLUDE_DIR}
    )

    set(linalg_c_test_sources
        c_linalg_test.c
        c_test_core.c
        c_linalg_test_misc.c
        c_linalg_test_factor.c
        c_linalg_test_eigen.c
    )

    # Build the C API tests
    add_executable(linalg_c_test ${linalg_c_test_sources})
    link_library(linalg_c_test linalg ${PROJECT_INCLUDE_DIR})
    add_test(
        NAME linalg_c_test
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:linalg_c_test>
    )
endif()