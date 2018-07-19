! test_immutable.f90

module test_immutable
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg_immutable
    implicit none
contains
! ******************************************************************************
! LU FACTORIZATION
! ------------------------------------------------------------------------------
    function test_im_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(n, n) :: a
        type(lu_results) :: lu_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a) 

        ! Compute the factorization
        lu_rst = mat_lu(a)

        ! Ensuure P * A = L * U
        if (.not.is_mtx_equal(matmul(lu_rst%p, a), &
                matmul(lu_rst%l, lu_rst%u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_lu_solve() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(n, n), b(n, nrhs), x(n, nrhs), y(n, nrhs)
        type(lu_results) :: lu_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)

        ! Compute the factorization
        lu_rst = mat_lu(a)

        ! Solve L * Y = P * B, for Y
        y = mat_solve_lower_tri(lu_rst%l, matmul(lu_rst%p, b))

        ! Solve U * X = Y, for X
        x = mat_solve_upper_tri(lu_rst%u, y)

        ! Ensure A * X = B
        if (.not.is_mtx_equal(matmul(a, x), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable LU Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_qr_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n)
        type(qr_results) :: qr_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Compute the factorization
        qr_rst = mat_qr(a)

        ! Test to see if A = Q * R
        if (.not.is_mtx_equal(a, matmul(qr_rst%q, qr_rst%r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable QR Factorization Test"
        end if
    end function 

! ------------------------------------------------------------------------------
    function test_im_qr_factor_pvt() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n)
        type(qr_results) :: qr_rst
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Compute the factorization
        qr_rst = mat_qr(a, .true.)

        ! Test to see if A * P = Q * R
        if (.not.is_mtx_equal(matmul(a, qr_rst%p), matmul(qr_rst%q, qr_rst%r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable QR Factorization Test with Pivoting"
        end if
    end function 

! ------------------------------------------------------------------------------
    function test_im_cholesky() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(n, n), r(n, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Ensure A is positive definite
        a = matmul(a, transpose(a))

        ! Compute the factorization
        r = mat_cholesky(a)

        ! Ensure A = R**T * R
        if (.not.is_mtx_equal(a, matmul(transpose(r), r), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable Cholesky Factorization"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_svd() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n)
        type(svd_results) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Compute the factorization
        x = mat_svd(a)

        ! Ensure U * S * V**T = A
        if (.not.is_mtx_equal(matmul(x%u, matmul(x%s, x%vt)), a, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable SVD"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_inverse() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(n, n), i(n, n), ainv(n, n)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a) 
        i = identity(n)

        ! Compute the inverse
        ainv = mat_inverse(a)

        ! Ensure ainv * a == I
        if (.not.is_mtx_equal(matmul(ainv, a), i, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable Inverse"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_pinverse() result(rst)
        ! Parameters
        integer(int32), parameter :: m = 100
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(m, n), ainv(n, m)
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Compute the inverse
        ainv = mat_pinverse(a)

        ! Ensure A * A+ * A = A
        if (.not.is_mtx_equal(matmul(a, matmul(ainv, a)), a, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable Pseudo-Inverse"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_im_eigen() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 100
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64) :: a(n, n)
        type(eigen_results) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)

        ! Compute the eigenvalues and eigenvectors
        x = mat_eigen(a)

        ! Ensure A * v = lambda * v
        if (.not.is_mtx_equal(matmul(a, x%vectors), mat_mult_diag(x%vectors, x%values), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Immutable Eigen Analysis"
        end if
    end function
end module
