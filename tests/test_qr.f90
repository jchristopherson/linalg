! test_qr.f90

! Tests for QR factorization/solution operations.
module test_qr
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! QR FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_qr_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r1, r2
        real(real64), dimension(m, m) :: q1, q2
        real(real64), dimension(n, n) :: p2
        real(real64), dimension(n) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Factorization Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Factorization Test 2, Part C"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r1, r2
        complex(real64), dimension(m, m) :: q1, q2
        complex(real64), dimension(n, n) :: p2
        complex(real64), dimension(n) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Factorization Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Factorization Test 2, Part C"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_factor_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, r1, r2
        real(real64), dimension(m, m) :: q1, q2
        real(real64), dimension(n, n) :: p2
        real(real64), dimension(n) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Test 2, Part C"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_factor_od_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        integer(int32) :: i, j
        complex(real64), dimension(m, n) :: a, r1, r2
        complex(real64), dimension(m, m) :: q1, q2
        complex(real64), dimension(n, n) :: p2
        complex(real64), dimension(n) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Overdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Overdetermined QR Test 2, Part C"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_factor_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r1, r2
        real(real64), dimension(m, m) :: q1, q2
        real(real64), dimension(n, n) :: p2
        real(real64), dimension(m) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Test 2, Part C"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_factor_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r1, r2
        complex(real64), dimension(m, m) :: q1, q2
        complex(real64), dimension(n, n) :: p2
        complex(real64), dimension(m) :: tau1, tau2
        integer(int32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.assert(a, matmul(q1, r1), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.assert(matmul(a, p2), matmul(q2, r2), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Test 2, Part C"
        end if
    end function

! ******************************************************************************
! QR MULTIPLICATION TEST
! ------------------------------------------------------------------------------
    ! LEFT SIDE
    function test_qr_mult() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r, c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r, c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Multiplication Test 1"
        end if

        ! Compute C = Q**H * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**H * C
        ans = matmul(conjg(transpose(q)), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! OVERDETERMINED - LEFT
    function test_qr_mult_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, r, c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_od_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r, c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**H * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**H * C
        ans = matmul(conjg(transpose(q)), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Overdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! UNDERDETERMINED - LEFT
    function test_qr_mult_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r, c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r, c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = Q * C
        call mult_qr(.true., .false., r, tau, c1)

        ! Compute ANS = Q * C
        ans = matmul(q, c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**H * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**H * C
        ans = matmul(conjg(transpose(q)), c2)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! RIGHT SIDE
    function test_qr_mult_right() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r
        real(real64), dimension(n, m) :: c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_right_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r
        complex(real64), dimension(n, m) :: c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**H
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**H
        ans = matmul(c2, conjg(transpose(q)))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! OVERDETERMINED - RIGHT SIDE
    function test_qr_mult_right_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, r
        real(real64), dimension(n, m) :: c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right Overdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_right_od_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r
        complex(real64), dimension(n, m) :: c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**H
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**H
        ans = matmul(c2, conjg(transpose(q)))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right Overdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! UNDERDETERMINED - RIGHT SIDE
    function test_qr_mult_right_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r
        real(real64), dimension(n, m) :: c1, c2, ans
        real(real64), dimension(m, m) :: q
        real(real64), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Right Underdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_right_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r
        complex(real64), dimension(n, m) :: c1, c2, ans
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Generate the QR factorization of A
        call qr_factor(a, tau)
        r = a
        call form_qr(a, tau, q)

        ! Compute C = C * Q
        call mult_qr(.false., .false., r, tau, c1)

        ! Compute ANS = C * Q
        ans = matmul(c2, q)

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**H
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**H
        ans = matmul(c2, conjg(transpose(q)))

        ! Test
        if (.not.assert(c1, ans, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Right Underdetermined QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    ! VECTOR
    function test_qr_mult_vector() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a, r
        real(real64), dimension(m, m) :: q
        real(real64), dimension(m) :: c1, c2, ans
        real(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the QR factorization
        call qr_factor(a, tau)
        r = a
        call form_qr(r, tau, q)

        ! Compute C1 = Q * C1
        call mult_qr(.false., a, tau, c1)

        ! Compute ANS = Q * C2
        ans = matmul(q, c2)

        ! Compare
        if (.not.assert(ans, c1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Vector QR Multiplication Test 1"
        end if

        ! Compute C1 = Q**T * C1
        c1 = c2
        call mult_qr(.true., a, tau, c1)

        ! Compute ANS = Q * C2
        ans = matmul(transpose(q), c2)

        ! Compare
        if (.not.assert(ans, c1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Vector QR Multiplication Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_mult_vector_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60

        ! Local Variables
        complex(real64), dimension(m, n) :: a, r
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(m) :: c1, c2, ans
        complex(real64), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(c1)
        c2 = c1

        ! Compute the QR factorization
        call qr_factor(a, tau)
        r = a
        call form_qr(r, tau, q)

        ! Compute C1 = Q * C1
        call mult_qr(.false., a, tau, c1)

        ! Compute ANS = Q * C2
        ans = matmul(q, c2)

        ! Compare
        if (.not.assert(ans, c1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Vector QR Multiplication Test 1"
        end if

        ! Compute C1 = Q**H * C1
        c1 = c2
        call mult_qr(.true., a, tau, c1)

        ! Compute ANS = Q * C2
        ans = matmul(conjg(transpose(q)), c2)

        ! Compare
        if (.not.assert(ans, c1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Vector QR Multiplication Test 2"
        end if
    end function

! ******************************************************************************
! QR SOLUTION TEST
! ------------------------------------------------------------------------------
    function test_qr_solve_no_pivot() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(m, nrhs) :: b, b1, ans1
        real(real64), dimension(n, nrhs) :: x1
        real(real64), dimension(n) :: tau
        real(real64), dimension(m) :: b2a, b2, ans2
        real(real64), dimension(n) :: x2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2a)
        a1 = a
        b1 = b
        b2 = b2a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau)

        ! Solve the system of equations
        call solve_qr(a1, tau, b1)

        ! Get X1 from B1
        x1 = b1(1:n,:)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 1, No Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 2, No Pivoting"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_no_pivot_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 60
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(m, nrhs) :: b, b1, ans1
        complex(real64), dimension(n, nrhs) :: x1
        complex(real64), dimension(n) :: tau
        complex(real64), dimension(m) :: b2a, b2, ans2
        complex(real64), dimension(n) :: x2
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2a)
        a1 = a
        b1 = b
        b2 = b2a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau)

        ! Solve the system of equations
        call solve_qr(a1, tau, b1)

        ! Get X1 from B1
        x1 = b1(1:n,:)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Solution Test 1, No Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Solution Test 2, No Pivoting"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1
        real(real64), dimension(m, nrhs) :: b, b1, ans1
        real(real64), dimension(n, nrhs) :: x1
        real(real64), dimension(n) :: tau
        real(real64), dimension(m) :: b2a, b2, ans2
        real(real64), dimension(n) :: x2
        integer(int32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2a)
        a1 = a
        b1 = b
        b2 = b2a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau, pvt)

        ! Solve the system of equations
        call solve_qr(a1, tau, pvt, b1)

        ! Get X1 from B1
        x1 = b1(1:n,:)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, pvt, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 2, With Pivoting"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1
        complex(real64), dimension(m, nrhs) :: b, b1, ans1
        complex(real64), dimension(n, nrhs) :: x1
        complex(real64), dimension(n) :: tau
        complex(real64), dimension(m) :: b2a, b2, ans2
        complex(real64), dimension(n) :: x2
        integer(int32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2a)
        a1 = a
        b1 = b
        b2 = b2a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau, pvt)

        ! Solve the system of equations
        call solve_qr(a1, tau, pvt, b1)

        ! Get X1 from B1
        x1 = b1(1:n,:)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, pvt, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2a, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued QR Solution Test 2, With Pivoting"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot_od() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 200
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: a1, a2
        real(real64), dimension(m, nrhs) :: b1, b2
        real(real64), dimension(n) :: tau
        real(real64), allocatable, dimension(:) :: work
        real(real64) :: temp(1), rcond
        integer(int32), dimension(n) :: pvt
        integer(int32) :: lwork, info, rnk
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a1)
        call create_random_array(b1)
        a2 = a1
        b2 = b1

         ! Compute the solution via DGELSY
         rcond = epsilon(rcond)
         call dgelsy(m, n, nrhs, a1, m, b1, m, pvt, rcond, rnk, temp, -1, info)
         lwork = int(temp(1))
         allocate(work(lwork))
         call dgelsy(m, n, nrhs, a1, m, b1, m, pvt, rcond, rnk, work, lwork, &
            info)
        
         ! Compute the solution via QR factorization
         call qr_factor(a2, tau, pvt)
         call solve_qr(a2, tau, pvt, b2)

         ! Test
         if (.not.assert(b1(1:n,:), b2(1:n,:), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Solution Test, With Pivoting"
         end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot_od_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 200
        integer(int32), parameter :: n = 100
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(m, n) :: a1, a2
        complex(real64), dimension(m, nrhs) :: b1, b2
        complex(real64), dimension(n) :: tau
        complex(real64), allocatable, dimension(:) :: work
        real(real64), dimension(2 * n) :: rwork
        complex(real64) :: temp(1)
        real(real64) :: rcond
        integer(int32), dimension(n) :: pvt
        integer(int32) :: lwork, info, rnk
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a1)
        call create_random_array(b1)
        a2 = a1
        b2 = b1

         ! Compute the solution via DGELSY
         rcond = epsilon(rcond)
         call zgelsy(m, n, nrhs, a1, m, b1, m, pvt, rcond, rnk, temp, -1, &
            rwork, info)
         lwork = int(temp(1))
         allocate(work(lwork))
         call zgelsy(m, n, nrhs, a1, m, b1, m, pvt, rcond, rnk, work, lwork, &
            rwork, info)
        
         ! Compute the solution via QR factorization
         call qr_factor(a2, tau, pvt)
         call solve_qr(a2, tau, pvt, b2)

         ! Test
         if (.not.assert(b1(1:n,:), b2(1:n,:), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Overdetermined QR Solution Test, With Pivoting"
         end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot_ud() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 50
        integer(int32), parameter :: n = 60
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, a2
        real(real64), dimension(m, nrhs) :: b, ans1
        real(real64), dimension(n, nrhs) :: x1
        real(real64), dimension(m) :: tau
        real(real64), dimension(m) :: b2, ans2
        real(real64), dimension(n) :: x2
        integer(int32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2)
        a1 = a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau, pvt)
        a2 = a1

        ! Solve the system of equations
        x1(1:m,:) = b
        call solve_qr(a1, tau, pvt, x1)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        x2(1:m) = b2
        call solve_qr(a2, tau, pvt, x2)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Solution Test 2, With Pivoting"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_solve_pivot_ud_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 5
        integer(int32), parameter :: n = 6
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1, a2
        complex(real64), dimension(m, nrhs) :: b, ans1
        complex(real64), dimension(n, nrhs) :: x1
        complex(real64), dimension(m) :: tau
        complex(real64), dimension(m) :: b2, ans2
        complex(real64), dimension(n) :: x2
        integer(int32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call create_random_array(a)
        call create_random_array(b)
        call create_random_array(b2)
        a1 = a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau, pvt)
        a2 = a1

        ! Solve the system of equations
        x1(1:m,:) = b
        call solve_qr(a1, tau, pvt, x1)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.assert(ans1, b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        x2(1:m) = b2
        call solve_qr(a2, tau, pvt, x2)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.assert(ans2, b2, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Underdetermined QR Solution Test 2, With Pivoting"
        end if
    end function

! ******************************************************************************
! QR UPDATE TEST
! ------------------------------------------------------------------------------
    function test_qr_update_1() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50

        ! Local Variables
        real(real64), dimension(m, n) :: a, a1, r
        real(real64), dimension(m, m) :: q
        real(real64), dimension(m) :: u
        real(real64), dimension(n) :: v, tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(u)
        call create_random_array(v)

        ! Compute the QR factorization of A
        r = a
        call qr_factor(r, tau)

        ! Form Q and R
        call form_qr(r, tau, q)

        ! Compute A1 = A + u * v**T
        a1 = a
        call rank1_update(1.0d0, u, v, a1)

        ! Use the QR update to update the original R and Q
        call qr_rank1_update(q, r, u, v)

        ! Test
        if (.not.assert(a1, matmul(q, r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Rank 1 QR Update"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_qr_update_1_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 60
        integer(int32), parameter :: n = 50
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a, a1, r
        complex(real64), dimension(m, m) :: q
        complex(real64), dimension(m) :: u
        complex(real64), dimension(n) :: v, tau
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(u)
        call create_random_array(v)

        ! Compute the QR factorization of A
        r = a
        call qr_factor(r, tau)

        ! Form Q and R
        call form_qr(r, tau, q)

        ! Compute A1 = A + u * v**H
        a1 = a
        call rank1_update(one, u, v, a1)

        ! Use the QR update to update the original R and Q
        call qr_rank1_update(q, r, u, v)

        ! Test
        if (.not.assert(a1, matmul(q, r), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Rank 1 QR Update"
        end if
    end function

! ------------------------------------------------------------------------------
end module
