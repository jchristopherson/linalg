# LU Example
add_executable(lu_example linalg_lu_example.f90)
target_link_libraries(lu_example linalg)

# QR Example
add_executable(qr_example linalg_qr_example.f90)
target_link_libraries(qr_example linalg)

# Overdetermined Example
add_executable(od_example linalg_od_example.f90)
target_link_libraries(od_example linalg)

# Eigen Example
add_executable(eigen_example linalg_eigen_example.f90)
target_link_libraries(eigen_example linalg)

# LU Example - Full Method
add_executable(lu_full_example linalg_lu_full_example.f90)
target_link_libraries(lu_full_example linalg)

# QR Example - Full Method, No Pivoting
add_executable(qr_full_example linalg_qr_full_example.f90)
target_link_libraries(qr_full_example linalg)

# QR Example - QR Multiplication Example, No Pivoting
add_executable(qr_mult_example linalg_qr_mult_example.f90)
target_link_libraries(qr_mult_example linalg)

# QR Update Example
add_executable(qr_update_example linalg_qr_update_example.f90)
target_link_libraries(qr_update_example linalg)

# Cholesky Example
add_executable(cholesky_example linalg_cholesky_example.f90)
target_link_libraries(cholesky_example linalg)

# Cholesky Update Example
add_executable(cholesky_update_example linalg_cholesky_update_example.f90)
target_link_libraries(cholesky_update_example linalg)

# Cholesky Downdate Example
add_executable(cholesky_downdate_example linalg_cholesky_downdate_example.f90)
target_link_libraries(cholesky_downdate_example linalg)

# Matrix Inverse Example
add_executable(inverse_example linalg_inverse_example.f90)
target_link_libraries(inverse_example linalg)

# Moore-Penrose Example
add_executable(pinverse_example linalg_pinverse_example.f90)
target_link_libraries(pinverse_example linalg)

# SVD Example
add_executable(svd_example linalg_svd_example.f90)
target_link_libraries(svd_example linalg)

# --------------------
# Immutable LU Example
add_executable(lu_example_immutable linalg_lu_example_immutable.f90)
target_link_libraries(lu_example_immutable linalg)

# --------------------
# Copy FERROR, BLAS, and LAPACK to the runtime directory
if (WIN32)
    add_custom_command(
        TARGET svd_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:svd_example>
    )
    add_custom_command(
        TARGET svd_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${blas_LibLocation} $<TARGET_FILE_DIR:svd_example>
    )
    add_custom_command(
        TARGET svd_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${lapack_LibLocation} $<TARGET_FILE_DIR:svd_example>
    )
endif()
