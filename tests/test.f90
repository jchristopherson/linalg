! test.f90

! A testing application.
program main
    ! Imported Modules
    use test_qr
    use test_svd_ops
    use test_mtx_inverse
    use test_cholesky
    use test_eigen
    use test_misc
    use test_lu

    ! Introduce the testing application
    print '(A)', "Hello from the LINALG test application."

    ! Misc. Item Tests
    call test_diagonal_mtx_mult()
    call test_rank1_update()
    call test_rank()

    ! LU Factorization Tests
    call test_lu_factor()
    call test_lu_solve()

    ! QR Factorization Tests
    call test_qr_factor()
    call test_qr_factor_od()
    call test_qr_factor_ud()
    call test_qr_mult()
    call test_qr_mult_od()
    call test_qr_mult_ud()
    call test_qr_mult_right()
    call test_qr_mult_right_od()
    call test_qr_mult_right_ud()
    call test_qr_mult_vector()
    call test_qr_solve_no_pivot()
    call test_qr_solve_pivot()
    call test_qr_solve_pivot_ud()
    call test_qr_update_1()

    ! SVD Tests
    call test_svd()
    call test_svd_od()
    call test_svd_ud()

    ! Matrix Inverse Tests
    call test_pinv()
    call test_pinv_od()
    call test_inv()

    ! Cholesky Factorization Tests
    call test_cholesky_factor()
    call test_cholesky_rank1_update()
    call test_cholesky_rank1_downdate()

    ! Eigenvalue/Eigenvector Tests
    call test_eigen_symm()
    call test_eigen_asymm()
    call test_eigen_gen()

end program
