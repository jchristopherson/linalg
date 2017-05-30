! test.f90

! A testing application.
program main
    ! Imported Modules
    use test_qr

    ! Introduce the testing application
    print '(A)', "Hello from the LINALG test application."

    ! Tests
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
end program
