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

# --------------------
include_directories(${PROJECT_SOURCE_DIR}/include 
${PROJECT_SOURCE_DIR}/src/external/ferror/include)

# LU Example
add_executable(lu_example_c linalg_lu_c_example.c)
target_link_libraries(lu_example_c linalg)

# QR Example
add_executable(qr_example_c linalg_qr_c_example.c)
target_link_libraries(qr_example_c linalg)

# Overdetermined Example
add_executable(od_example_c linalg_od_c_example.c)
target_link_libraries(od_example_c linalg)

# Eigen Example
add_executable(eigen_example_c linalg_eigen_c_example.c)
target_link_libraries(eigen_example_c linalg m)