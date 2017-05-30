! test_qr.f90

! Tests for QR factorization/solution operations.
module test_qr
    use linalg_constants
    use qr
    use test_core
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


end module
