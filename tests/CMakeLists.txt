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
    test_sparse.f90
    test_rz.f90
)

# Dependencies
find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    set(_linalg_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(fortran_test_helper)
    set(BUILD_TESTING ${_linalg_build_testing})
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

# Build the tests
add_executable(linalg_test ${linalg_test_sources})
target_link_libraries(
    linalg_test
    linalg 
    ${FTH_LIBRARIES}
)
add_test(
    NAME linalg_test
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:linalg_test>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET linalg_test
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_RUNTIME_DLLS:linalg_test> 
        $<TARGET_FILE_DIR:linalg_test>
        COMMAND_EXPAND_LISTS
    )
endif()
