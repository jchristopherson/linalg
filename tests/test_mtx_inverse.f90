! test_mtx_inverse.f90
! SPDX-License-Identifier: MIT

! Tests matrix inversion routines
module test_mtx_inverse
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use test_core
    use linalg
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! PSEUDO-INVERSE TESTS
! ------------------------------------------------------------------------------
    function test_pinv() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 60
        integer, parameter :: nrhs = 20

        ! Local Variables
        real(real64), dimension(m, n) :: a
        real(real64), dimension(n, m) :: ainv
        real(real64), dimension(m, nrhs) :: b
        real(real64), dimension(n, nrhs) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Compute the inverse
        ainv = mtx_pinverse(a)

        ! Compute X = inv(A) * B
        x = matmul(ainv, b)

        ! Test: A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_od() result(rst)
        ! Parameters
        integer, parameter :: m = 80
        integer, parameter :: n = 60

        ! Local Variables
        real(real64), dimension(m, n) :: a
        real(real64), dimension(n, m) :: ainv
        real(real64), dimension(n, n) :: identity, check
        logical :: rst
        integer :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)

        identity = 0.0d0
        do i = 1, size(identity, 1)
            identity(i,i) = 1.0d0
        end do

        ! Compute the inverse
        ainv = mtx_pinverse(a)

        ! Compute A+ * A - should = I
        check = matmul(ainv, a)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Overdetermined Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_inv() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        real(real64), dimension(n, n) :: a, a1, ainv
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)

        ! Compute the inverse of A
        ainv = mtx_inverse(a)

        ! Compute the inverse using the already tested pseudo-inverse
        a1 = mtx_pinverse(a)

        ! Test
        if (.not.assert(ainv, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Matrix Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 60
        integer, parameter :: nrhs = 20

        ! Local Variables
        complex(real64), dimension(m, n) :: a
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, nrhs) :: b
        complex(real64), dimension(n, nrhs) :: x
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a)
        call create_random_array(b)

        ! Compute the inverse
        ainv = mtx_pinverse(a)

        ! Compute X = inv(A) * B
        x = matmul(ainv, b)

        ! Test: A * X = B
        if (.not.assert(matmul(a, x), b, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_od_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 80
        integer, parameter :: n = 60

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(n, n) :: identity, check
        logical :: rst
        integer :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)

        identity = zero
        do i = 1, size(identity, 1)
            identity(i,i) = one
        end do
        
        ! Compute the inverse
        ainv = mtx_pinverse(a)

        ! Compute A+ * A - should = I
        check = matmul(ainv, a)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pinv_ud_cmplx() result(rst)
        ! Parameters
        integer, parameter :: m = 60
        integer, parameter :: n = 80

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        complex(real64), dimension(m, n) :: a
        complex(real64), dimension(n, m) :: ainv
        complex(real64), dimension(m, m) :: identity, check
        logical :: rst
        integer :: i

        ! Initialization
        rst = .true.
        call create_random_array(a)

        identity = zero
        do i = 1, size(identity, 1)
            identity(i,i) = one
        end do
        
        ! Compute the inverse
        ainv = mtx_pinverse(a)

        ! Compute A * A+ - should = I
        check = matmul(a, ainv)
        if (.not.assert(check, identity, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex Pseudo-Inverse Test 3"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_inv_cmplx() result(rst)
        ! Parameters
        integer, parameter :: n = 100

        ! Local Variables
        complex(real64), dimension(n, n) :: a, ainv, a1
        logical :: rst

        ! Initialization
        rst = .true.
        call create_random_array(a, mtype = SYMMETRIC_MATRIX)

        ! Compute the inverse of A
        ainv = mtx_inverse(a)

        ! Compute the inverse using the already tested pseudo-inverse
        a1 = mtx_pinverse(a)

        ! Test
        if (.not.assert(ainv, a1, tol = REAL64_TOL)) then
            rst = .false.
            print '(A)', "Test Failed: Complex-Valued Matrix Inverse Test 1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_pseudoinverse_edge_cases() result(rst)
        real(real64) :: a(3,2), scaled(2,2), expected(2,2)
        real(real64) :: projector_left(3,3), projector_right(2,2)
        real(real64), allocatable :: pinv(:,:), ainv(:,:)
        complex(real64), allocatable :: ainvc(:,:)
        logical :: rst

        rst = .true.
        a = reshape([1.0_real64, 2.0_real64, 3.0_real64, &
            2.0_real64, 4.0_real64, 6.0_real64], [3,2])
        pinv = mtx_pinverse(a)
        if (.not.assert(matmul(matmul(a, pinv), a), a, tol = REAL64_TOL) .or. &
            .not.assert(matmul(matmul(pinv, a), pinv), pinv, tol = REAL64_TOL)) then
            rst = .false.
        end if
        projector_left = matmul(a, pinv)
        if (.not.assert(projector_left, transpose(projector_left), &
            tol = REAL64_TOL)) rst = .false.
        projector_right = matmul(pinv, a)
        if (.not.assert(projector_right, transpose(projector_right), &
            tol = REAL64_TOL)) rst = .false.

        ! A supplied cutoff must remove small singular directions exactly.
        scaled = 0.0_real64
        scaled(1,1) = 1.0_real64
        scaled(2,2) = 1.0e-4_real64
        expected = 0.0_real64
        expected(1,1) = 1.0_real64
        pinv = mtx_pinverse(scaled, tol = 1.0e-3_real64)
        if (.not.assert(pinv, expected, tol = REAL64_TOL)) rst = .false.

        ! The ordinary inverse documents NaN output for singular matrices.
        scaled = 0.0_real64
        scaled(1,1) = 1.0_real64
        ainv = mtx_inverse(scaled)
        if (.not.all(ieee_is_nan(ainv))) rst = .false.
        ainvc = mtx_inverse(cmplx(scaled, 0.0_real64, kind = real64))
        if (.not.all(ieee_is_nan(real(ainvc, real64))) .or. &
            .not.all(ieee_is_nan(aimag(ainvc)))) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: Pseudoinverse Edge Cases"
    end function

! ------------------------------------------------------------------------------
    function test_complex_rank_deficient_pseudoinverse() result(rst)
        complex(real64) :: a(3,2), projector_left(3,3), projector_right(2,2)
        complex(real64), allocatable :: pinv(:,:)
        logical :: rst

        rst = .true.
        a(:,1) = [(1.0_real64, 1.0_real64), (2.0_real64, -1.0_real64), &
            (0.0_real64, 0.0_real64)]
        a(:,2) = 2.0_real64 * a(:,1)
        pinv = mtx_pinverse(a)
        if (.not.assert(matmul(matmul(a, pinv), a), a, tol = REAL64_TOL) .or. &
            .not.assert(matmul(matmul(pinv, a), pinv), pinv, tol = REAL64_TOL)) then
            rst = .false.
        end if
        projector_left = matmul(a, pinv)
        if (.not.assert(projector_left, transpose(conjg(projector_left)), &
            tol = REAL64_TOL)) rst = .false.
        projector_right = matmul(pinv, a)
        if (.not.assert(projector_right, transpose(conjg(projector_right)), &
            tol = REAL64_TOL)) rst = .false.
        if (.not.rst) print '(A)', "Test Failed: Complex Rank-Deficient Pseudoinverse"
    end function

! ------------------------------------------------------------------------------
    function test_zero_matrix_pseudoinverse() result(rst)
        real(real64) :: zero(2,2)
        complex(real64) :: czero(2,2)
        real(real64), allocatable :: pinv(:,:)
        complex(real64), allocatable :: pinvc(:,:)
        logical :: rst

        rst = .true.
        zero = 0.0_real64
        czero = (0.0_real64, 0.0_real64)
        pinv = mtx_pinverse(zero)
        if (any(pinv /= 0.0_real64)) rst = .false.
        pinvc = mtx_pinverse(czero)
        if (any(pinvc /= (0.0_real64, 0.0_real64))) rst = .false.
        if (.not.rst) print '(A)', "Test Failed: Zero Matrix Pseudoinverse"
    end function

! ------------------------------------------------------------------------------
end module
