! test.f90

! A testing application.
program main
    ! Imported Modules
    use test_qr
    use test_svd_ops

    ! Introduce the testing application
    print '(A)', "Hello from the LINALG test application."

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
end program
