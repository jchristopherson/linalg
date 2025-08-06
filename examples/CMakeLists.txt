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

# Full LQ Example
add_executable(lq_full_example linalg_lq_full_example.f90)
target_link_libraries(lq_full_example linalg)

# LQ Multiplication Example
add_executable(lq_mult_example linalg_lq_mult_example.f90)
target_link_libraries(lq_mult_example linalg)

# LQ Example
add_executable(lq_example linalg_lq_example.f90)
target_link_libraries(lq_example linalg)

# Sparse Direct Solve Example
add_executable(sparse_direct_example linalg_sparse_direct_example.f90)
target_link_libraries(sparse_direct_example linalg)

# --------------------
# C API Eigenvalue Example
if (${BUILD_LINALG_C_API})
    include_directories(${PROJECT_SOURCE_DIR}/include)
    add_executable(c_eigen_example c_linalg_eigen_example.c)
    target_link_libraries(c_eigen_example linalg)

    add_executable(c_lu_example c_linalg_lu_example.c)
    target_link_libraries(c_lu_example linalg)
endif()

# ------------------------------------------------------------------------------
if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET lu_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:lu_example> 
        $<TARGET_FILE_DIR:lu_example>
        COMMAND_EXPAND_LISTS
    )
endif()
