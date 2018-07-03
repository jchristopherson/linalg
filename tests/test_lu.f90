! test_lu.f90

! Tests the LU factorization routines
module test_lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg_core
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1, l, u, p
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a  ! Forces us to keep a copy of the original matrix

        ! Compute the factorization
        call lu_factor(a1, ipvt)

        ! Extract L, U, and P to determine if P * A = L * U
        l = a1
        call form_lu(l, ipvt, u, p)
        if (.not.is_mtx_equal(matmul(p, a), matmul(l, u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1
        real(real64), dimension(n, nrhs) :: b, x
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)
        a1 = a
        x = b

        ! Factor A
        call lu_factor(a1, ipvt)

        ! Solve for X
        call solve_lu(a1, ipvt, x)

        ! Test by determining if A * X = B
        if (.not.is_mtx_equal(matmul(a, x), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        complex(real64), dimension(n, n) :: a, a1, l, u
        real(real64), dimension(n, n) :: p
        real(real64) :: temp1, temp2
        integer(int32) :: i, j
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                call random_number(temp1)
                call random_number(temp2)
                a(i,j) = cmplx(temp1, temp2, real64)
                a1(i,j) = a(i,j) ! Store a copy of the original matrix
            end do
        end do

        ! Compute the factorization
        call lu_factor(a1, ipvt)

        ! Extract L, U, and P to determine if P * A = L * U
        l = a1
        call form_lu(l, ipvt, u, p)
        if (.not.is_mtx_equal(matmul(p, a), matmul(l, u), tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        complex(real64), dimension(n, n) :: a, a1
        complex(real64), dimension(n, nrhs) :: b, x
        real(real64) :: temp1, temp2
        integer(int32) :: i, j
        integer(int32), dimension(n) :: ipvt
        logical :: rst

        ! Initialization
        rst = .true.
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                call random_number(temp1)
                call random_number(temp2)
                a(i,j) = cmplx(temp1, temp2, real64)
                a1(i,j) = a(i,j)
            end do
        end do

        do j = 1, size(b, 2)
            do i = 1, size(b, 1)
                call random_number(temp1)
                call random_number(temp2)
                b(i,j) = cmplx(temp1, temp2, real64)
                x(i,j) = b(i,j)
            end do
        end do

        ! Factor A
        call lu_factor(a1, ipvt)

        ! Solve for X
        call solve_lu(a1, ipvt, x)

        ! Test by determining if A * X = B
        if (.not.is_mtx_equal(matmul(a, x), b, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization & Solution Test"
        end if
    end function

end module
