! test_lu.f90

! Tests the LU factorization routines
module test_lu
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
contains
! ******************************************************************************
! LU FACTORIZATION TEST
! ------------------------------------------------------------------------------
    function test_lu_factor() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75

        ! Local Variables
        real(real64), dimension(n, n) :: a
        real(real64), allocatable, dimension(:,:) :: l, u, p
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the factorization
        call lu_factor(a, l = l, u = u, p = p)

        ! Determine if P * A = L * U
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(n, n) :: a
        real(real64), dimension(n, nrhs) :: b, x
        integer(int32), allocatable :: ipvt(:)
        real(real64), allocatable :: lu(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Factor A
        call lu_factor(a, ipvt = ipvt, lu = lu)

        ! Solve for X
        x = solve_lu(lu, ipvt, b)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factor_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75

        ! Local Variables
        complex(real64), dimension(n, n) :: a
        complex(real64), allocatable, dimension(:,:) :: l, u
        real(real64), allocatable, dimension(:,:) :: p
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the factorization
        call lu_factor(a, l = l, u = u, p = p)

        ! Determine if P * A = L * U
        if (.not.assert(matmul(p, a), matmul(l, u), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_solve_cmplx() result(rst)
        ! Parameters
        integer(int32), parameter :: n = 75
        integer(int32), parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(n, n) :: a
        complex(real64), dimension(n, nrhs) :: b, x
        integer(int32) :: i, j
        integer(int32), allocatable :: ipvt(:)
        complex(real64), allocatable :: lu(:,:)
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Factor A
        call lu_factor(a, ipvt = ipvt, lu = lu)

        ! Solve for X
        x = solve_lu(lu, ipvt, b)

        ! Test by determining if A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued LU Factorization & Solution Test"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_lu_factor_pure() result(rst)
        use linear_algebra

        ! Parameters
        integer(int32), parameter :: n = 75

        ! Local Variables
        real(real64) :: a(n,n)
        type(lu_factors) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the factorization
        x = lu_factor(a)

        ! Tests
        if (.not.assert(matmul(x%P, a), matmul(x%L, x%U), tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_lu_factor_pure"
        end if
    end function
end module
