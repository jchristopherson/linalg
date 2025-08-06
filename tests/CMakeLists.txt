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
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

# Build the Fortran API tests
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

if (${BUILD_LINALG_C_API})
    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/C
        ${CMAKE_SOURCE_DIR}/include
    )

    set(linalg_c_test_sources
        C/c_linalg_test.c
        C/c_test_core.c
        C/c_linalg_test_misc.c
        C/c_linalg_test_factor.c
        C/c_linalg_test_eigen.c
    )

    add_executable(linalg_c_test ${linalg_c_test_sources})
    target_link_libraries(linalg_c_test linalg)

    add_test(
        NAME linalg_c_test
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:linalg_c_test>
    )
endif()

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
