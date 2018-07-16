# Fortran Tests
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
    test_immutable.f90
)
add_executable(linalg_test ${linalg_test_sources})
target_link_libraries(linalg_test linalg)

# Copy the FERROR shared object to the runtime directory on Windows
if (WIN32)
    add_custom_command(
        TARGET linalg_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:linalg_test>
    )
endif()

# Run the tests
add_custom_command(OUTPUT linalg_tests DEPENDS ${linalg_test_sources}
    COMMAND linalg_test)
add_custom_target(run_linalg_tests ALL DEPENDS linalg_tests)
