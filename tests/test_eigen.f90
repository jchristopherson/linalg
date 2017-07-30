! test_eigen.f90

! Tests the eigenvalue/eigenvector routines
module test_eigen
    use linalg_constants
    use linalg_core, only : diag_mtx_mult
    use test_core
    use linalg_eigen
    implicit none
contains
! ******************************************************************************
! EIGEN ANALYSIS TESTS
! ------------------------------------------------------------------------------
    function test_eigen_symm() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, vecs, x, y
        real(dp), dimension(n) :: vals
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a = matmul(a, transpose(a))

        ! Compute the eigenvalues and eigenvectors of A
        vecs = a
        call eigen(.true., vecs, vals)

        ! Compute vecs * vals, where vals is a diagonal matrix
        call diag_mtx_mult(.false., .false., 1.0d0, vals, vecs, 0.0d0, x)

        ! Test
        y = matmul(a, vecs)
        if (.not.is_mtx_equal(x, y, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Symmetric Eigen Values"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eigen_asymm() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1
        complex(dp), dimension(n, n) :: vecs, vmtx, x, y
        complex(dp), dimension(n) :: vals, vals1
        integer(i32) :: i
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        a1 = a
        vmtx = cmplx(0.0d0, 0.0d0, dp)

        ! Compute the eigenvalues and eigenvectors of A
        call eigen(a1, vals, vecs)

        ! Compute vecs * vals, where vals is a diagonal matrix
        do i = 1, n
            vmtx(i,i) = vals(i)
        end do
        x = matmul(vecs, vmtx)

        ! Test 1
        y = matmul(a, vecs)
        if (.not.is_mtx_equal(x, y, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Asymmetric Eigen Values Test 1"
        end if

        ! Compute just the eigenvalues
        call eigen(a, vals1)
        if (.not.is_mtx_equal(vals, vals1, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Asymmetric Eigen Values Test 2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_eigen_gen() result(rst)
        ! Parameters
        integer(i32), parameter :: n = 100
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        real(dp), dimension(n, n) :: a, a1, b, b1
        complex(dp), dimension(n) :: vals, vals2
        complex(dp), dimension(n, n) :: vecs, x, y
        logical :: rst

        ! Initialization
        rst = .true.
        call random_number(a)
        call random_number(b)
        a1 = a
        b1 = b

        ! Test 1
        call eigen(a1, b1, vals, vecs = vecs)

        ! Compute X = A * VECS
        x = matmul(a, vecs)

        ! And Y = B * (VECS * LAMBDA)
        call diag_mtx_mult(.false., .false., 1.0d0, vals, vecs, 0.0d0, y)
        y = matmul(b, y)

        ! Check
        if (.not.is_mtx_equal(x, y, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Generalized Eigen Values Test 1"
        end if

        ! Test 2 - Eigenvalues Only
        call eigen(a, b, vals2)
        if (.not.is_mtx_equal(vals, vals2, tol)) then
            rst = .false.
            print '(A)', "Test Failed: Generalized Eigen Values Test 2"
        end if
    end function


end module
