# Fortran API tests source code
set(linalg_test_sources
    linalg_test.f90
    test_core.f90
    test_qr.f90
    test_svd.f90
    test_least_squares.f90
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

# Pure function APIs report invalid input with ERROR STOP in isolated processes.
foreach(validation_case IN ITEMS
    csr_direct_bad_rhs
    csr_direct_bad_storage
    csr_lu_bad_rhs
    csr_lu_bad_pointer
    csr_gmres_bad_restart
    csr_gmres_bad_storage
    operator_gmres_bad_restart
    csr_operator_bad_storage
    ilu_operator_bad_pointer
)
    add_test(
        NAME ${validation_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:linalg_test> ${validation_case}
    )
    set_tests_properties(${validation_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(test_case IN ITEMS
    basic_utility_edges
    pseudoinverse_edge_cases
    complex_rank_deficient_pseudoinverse
    zero_matrix_pseudoinverse
    least_squares_shapes
    least_squares_rank_deficient
    rz_square
    sort_index_edges
)
    add_test(
        NAME ${test_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:linalg_test> ${test_case}
    )
endforeach()

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
