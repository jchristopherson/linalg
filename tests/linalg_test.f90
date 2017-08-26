! linalg_test.f90

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
    use test_sort

    ! Local Variables
    logical :: rst, overall

    ! Initialization
    overall = .true.

    ! Misc. Item Tests
    rst = test_diagonal_mtx_mult()
    if (.not.rst) overall = .false.

    rst = test_rank1_update()
    if (.not.rst) overall = .false.
    
    rst = test_rank()
    if (.not.rst) overall = .false.
    
    rst = test_tri_mtx_mult_1()
    if (.not.rst) overall = .false.
    
    rst = test_tri_mtx_mult_2()
    if (.not.rst) overall = .false.
    

    ! LU Factorization Tests
    rst = test_lu_factor()
    if (.not.rst) overall = .false.
    
    rst = test_lu_solve()
    if (.not.rst) overall = .false.
    

    ! QR Factorization Tests
    rst = test_qr_factor()
    if (.not.rst) overall = .false.
    
    rst = test_qr_factor_od()
    if (.not.rst) overall = .false.
    
    rst = test_qr_factor_ud()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_od()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_ud()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_right()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_right_od()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_right_ud()
    if (.not.rst) overall = .false.
    
    rst = test_qr_mult_vector()
    if (.not.rst) overall = .false.
    
    rst = test_qr_solve_no_pivot()
    if (.not.rst) overall = .false.
    
    rst = test_qr_solve_pivot()
    if (.not.rst) overall = .false.

    rst = test_qr_solve_pivot_od()
    if (.not.rst) overall = .false.
    
    rst = test_qr_solve_pivot_ud()
    if (.not.rst) overall = .false.
    
    rst = test_qr_update_1()
    if (.not.rst) overall = .false.
    
    ! SVD Tests
    rst = test_svd()
    if (.not.rst) overall = .false.
    
    rst = test_svd_od()
    if (.not.rst) overall = .false.
    
    rst = test_svd_ud()
    if (.not.rst) overall = .false.
    
    ! Matrix Inverse Tests
    rst = test_pinv()
    if (.not.rst) overall = .false.
    
    rst = test_pinv_od()
    if (.not.rst) overall = .false.
    
    rst = test_inv()
    if (.not.rst) overall = .false.
    
    ! Cholesky Factorization Tests
    rst =  test_cholesky_factor()
    if (.not.rst) overall = .false.
    
    rst = test_cholesky_rank1_update()
    if (.not.rst) overall = .false.
    
    rst = test_cholesky_rank1_downdate()
    if (.not.rst) overall = .false.
    
    ! Eigenvalue/Eigenvector Tests
    rst = test_eigen_symm()
    if (.not.rst) overall = .false.
    
    rst = test_eigen_asymm()
    if (.not.rst) overall = .false.
    
    rst = test_eigen_gen()
    if (.not.rst) overall = .false.

    ! Sorting Tests
    rst = test_dbl_ascend_sort()
    if (.not.rst) overall = .false.

    rst = test_dbl_descend_sort()
    if (.not.rst) overall = .false.
    
    ! End
    if (overall) then
        print '(A)', "LINALG TEST STATUS: PASS"
    else
        print '(A)', "LINALG TEST STATUS: FAILED"
    end if
end program
