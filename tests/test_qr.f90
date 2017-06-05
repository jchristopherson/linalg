! test_qr.f90

! Tests for QR factorization/solution operations.
module test_qr
    use linalg_constants
    use test_core
    use linalg_core, only : rank1_update
    use linalg_factor, only : qr_factor, mult_qr, form_qr, qr_rank1_update
    use linalg_solve, only : solve_qr
    implicit none
contains
! ******************************************************************************
! QR FACTORIZATION TEST
! ------------------------------------------------------------------------------
    subroutine test_qr_factor()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r1, r2
        real(dp), dimension(m, m) :: q1, q2
        real(dp), dimension(n, n) :: p2
        real(dp), dimension(n) :: tau1, tau2
        integer(i32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.is_mtx_equal(a, matmul(q1, r1), tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Factorization Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.is_mtx_equal(matmul(a, p2), matmul(q2, r2), tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Factorization Test 2, Part C"
        end if
        if (rst) print '(A)', "Test Passed: QR Factorization"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_qr_factor_od()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 50
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r1, r2
        real(dp), dimension(m, m) :: q1, q2
        real(dp), dimension(n, n) :: p2
        real(dp), dimension(n) :: tau1, tau2
        integer(i32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.is_mtx_equal(a, matmul(q1, r1), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.is_mtx_equal(matmul(a, p2), matmul(q2, r2), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Test 2, Part C"
        end if
        if (rst) print '(A)', "Test Passed: Overdetermined QR Factorization"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_qr_factor_ud()
        ! Parameters
        integer(i32), parameter :: m = 50
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r1, r2
        real(dp), dimension(m, m) :: q1, q2
        real(dp), dimension(n, n) :: p2
        real(dp), dimension(m) :: tau1, tau2
        integer(i32), dimension(n) :: pvt2
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        r1 = a
        r2 = a

        ! Compute the QR factorization of A
        call qr_factor(r1, tau1)

        ! Extract Q and R, and then check that Q * R = A
        call form_qr(r1, tau1, q1)
        if (.not.is_mtx_equal(a, matmul(q1, r1), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Test 1, Part C"
        end if

        ! Compute the QR factorization of A with pivoting
        pvt2 = 0
        call qr_factor(r2, tau2, pvt2)

        ! Extract Q, R, and P, and then check that Q * R = A * P
        call form_qr(r2, tau2, pvt2, q2, p2)
        if (.not.is_mtx_equal(matmul(a, p2), matmul(q2, r2), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Test 2, Part C"
        end if
        if (rst) print '(A)', "Test Passed: Underdetermined QR Factorization"
    end subroutine

! ******************************************************************************
! QR MULTIPLICATION TEST
! ------------------------------------------------------------------------------
    ! LEFT SIDE
    subroutine test_qr_mult()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r, c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! OVERDETERMINED - LEFT
    subroutine test_qr_mult_od()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 50
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r, c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Overdetermined QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! UNDERDETERMINED - LEFT
    subroutine test_qr_mult_ud()
        ! Parameters
        integer(i32), parameter :: m = 50
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r, c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = Q**T * C
        c1 = c2
        call mult_qr(.true., .true., r, tau, c1)

        ! Compute ANS = Q**T * C
        ans = matmul(transpose(q), c2)

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Underdetermined QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! RIGHT SIDE
    subroutine test_qr_mult_right()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r
        real(dp), dimension(n, m) :: c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Right QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! OVERDETERMINED - RIGHT SIDE
    subroutine test_qr_mult_right_od()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 50
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r
        real(dp), dimension(n, m) :: c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right Overdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right Overdetermined QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Right Overdetermined QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! UNDERDETERMINED - RIGHT SIDE
    subroutine test_qr_mult_right_ud()
        ! Parameters
        integer(i32), parameter :: m = 50
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r
        real(dp), dimension(n, m) :: c1, c2, ans
        real(dp), dimension(m, m) :: q
        real(dp), dimension(m) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right Underdetermined QR Multiplication Test 1"
        end if

        ! Compute C = C * Q**T
        c1 = c2
        call mult_qr(.false., .true., r, tau, c1)

        ! Compute ANS = C * Q**T
        ans = matmul(c2, transpose(q))

        ! Test
        if (.not.is_mtx_equal(c1, ans, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Right Underdetermined QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Right Underdetermined QR Multiplication"
    end subroutine

! ------------------------------------------------------------------------------
    ! VECTOR
    subroutine test_qr_mult_vector()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, r
        real(dp), dimension(m, m) :: q
        real(dp), dimension(m) :: c1, c2, ans
        real(dp), dimension(n) :: tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(c1)
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
        if (.not.is_mtx_equal(ans, c1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Vector QR Multiplication Test 1"
        end if

        ! Compute C1 = Q**T * C1
        c1 = c2
        call mult_qr(.true., a, tau, c1)

        ! Compute ANS = Q * C2
        ans = matmul(transpose(q), c2)

        ! Compare
        if (.not.is_mtx_equal(ans, c1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Vector QR Multiplication Test 2"
        end if
        if (rst) print '(A)', "Test Passed: Vector QR Multiplication"
    end subroutine

! ******************************************************************************
! QR SOLUTION TEST
! ------------------------------------------------------------------------------
    subroutine test_qr_solve_no_pivot()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        integer(i32), parameter :: nrhs = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1
        real(dp), dimension(m, nrhs) :: b, b1, ans1
        real(dp), dimension(n, nrhs) :: x1
        real(dp), dimension(n) :: tau
        real(dp), dimension(m) :: b2a, b2, ans2
        real(dp), dimension(n) :: x2
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)
        call random_number(b2a)
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
        if (.not.is_mtx_equal(ans1, b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 1, No Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.is_mtx_equal(ans2, b2a, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 2, No Pivoting"
        end if
        if (rst) print '(A)', "Test Passed: QR Solution, No Pivoting"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_qr_solve_pivot()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 60
        integer(i32), parameter :: nrhs = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1
        real(dp), dimension(m, nrhs) :: b, b1, ans1
        real(dp), dimension(n, nrhs) :: x1
        real(dp), dimension(n) :: tau
        real(dp), dimension(m) :: b2a, b2, ans2
        real(dp), dimension(n) :: x2
        integer(i32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call random_number(a)
        call random_number(b)
        call random_number(b2a)
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
        if (.not.is_mtx_equal(ans1, b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        call solve_qr(a1, tau, pvt, b2)

        ! Get X2 from B2
        x2 = b2(1:n)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.is_mtx_equal(ans2, b2a, tol)) then
            rst = .false.
            print '(A)', "Test Failed: QR Solution Test 2, With Pivoting"
        end if
        if (rst) print '(A)', "Test Passed: QR Solution, With Pivoting"
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_qr_solve_pivot_ud()
        ! Parameters
        integer(i32), parameter :: m = 5
        integer(i32), parameter :: n = 6
        integer(i32), parameter :: nrhs = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1, a2
        real(dp), dimension(m, nrhs) :: b, ans1
        real(dp), dimension(n, nrhs) :: x1
        real(dp), dimension(m) :: tau
        real(dp), dimension(m) :: b2, ans2
        real(dp), dimension(n) :: x2
        integer(i32), dimension(n) :: pvt
        logical :: rst

        ! Initialization
        rst = .true.
        pvt = 0
        call random_number(a)
        call random_number(b)
        call random_number(b2)
        a1 = a

        ! Compute the QR factorization of A
        call qr_factor(a1, tau, pvt)
        a2 = a1

        ! Solve the system of equations
        x1(1:m,:) = b
        call solve_qr(a1, tau, pvt, x1)

        ! Test
        ans1 = matmul(a, x1)
        if (.not.is_mtx_equal(ans1, b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Solution Test 1, With Pivoting"
        end if

        ! Solve the system of equations - vector
        x2(1:m) = b2
        call solve_qr(a2, tau, pvt, x2)

        ! Test
        ans2 = matmul(a, x2)
        if (.not.is_mtx_equal(ans2, b2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Underdetermined QR Solution Test 2, With Pivoting"
        end if
        if (rst) print '(A)', "Test Passed: Underdetermined QR Solution, With Pivoting"
    end subroutine

! ******************************************************************************
! QR UPDATE TEST
! ------------------------------------------------------------------------------
    subroutine test_qr_update_1()
        ! Parameters
        integer(i32), parameter :: m = 60
        integer(i32), parameter :: n = 50
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(m, n) :: a, a1, r
        real(dp), dimension(m, m) :: q
        real(dp), dimension(m) :: u
        real(dp), dimension(n) :: v, tau
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(u)
        call random_number(v)

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
        if (.not.is_mtx_equal(a1, matmul(q, r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Rank 1 QR Update"
        end if
        if (rst) print '(A)', "Test Passed: Rank 1 QR Update"
    end subroutine

end module
