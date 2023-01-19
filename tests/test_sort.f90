! test_sort.f90

! A module containing routines to support testing of the sorting routines.
module test_sort
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use test_core
    use linalg
    use fortran_test_helper
    implicit none

contains
! ******************************************************************************
! DOUBLE PRECISION TESTS
! ------------------------------------------------------------------------------
function test_dbl_ascend_sort() result(rst)
    ! Parameters
    integer(int32), parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer(int32) :: i, ind(n)
    real(real64) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .true.
    call create_random_array(x1)
    do i = 1, n
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Use the LAPACK-based sorting routine as a benchmark
    call sort(x1, ascend)

    ! Now, track the indices
    call sort(x2, ind, ascend)

    ! Compare the two arrays - we really don't have a good means of checking
    ! ind at this point
    if (.not.assert(x1, x2, tol = REAL64_TOL)) then
        rst = .false.
        print '(A)', "Test Failed: Ascending sort of a double-precision array."
    end if

    ! Ensure the array is in ascending order
    do i = 2, n
        if (x1(i) < x1(i-1)) then
            rst = .false.
            print '(A)', "Test Failed: Sorted array is not in ascending order."
            exit
        end if
    end do

    ! End
end function

! ------------------------------------------------------------------------------
function test_ascend_sort_cmplx() result(rst)
    ! Parameters
    integer(int32), parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer(int32) :: i, ind(n)
    complex(real64) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .true.
    call create_random_array(x1)
    do i = 1, n
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Use the LAPACK-based sorting routine as a benchmark
    call sort(x1, ascend)

    ! Now, track the indices
    call sort(x2, ind, ascend)

    ! Compare the two arrays - we really don't have a good means of checking
    ! ind at this point
    if (.not.assert(x1, x2, tol = REAL64_TOL)) then
        rst = .false.
        print '(A)', "Test Failed: Complex-Valued Ascending sort of a double-precision array."
    end if

    ! End
end function

! ------------------------------------------------------------------------------
function test_dbl_descend_sort() result(rst)
    ! Parameters
    integer(int32), parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer(int32) :: i, ind(n)
    real(real64) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .false.
    call create_random_array(x1)
    do i = 1, n
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Use the LAPACK-based sorting routine as a benchmark
    call sort(x1, ascend)

    ! Now, track the indices
    call sort(x2, ind, ascend)

    ! Compare the two arrays - we really don't have a good means of checking
    ! ind at this point
    if (.not.assert(x1, x2, tol = REAL64_TOL)) then
        rst = .false.
        print '(A)', "Test Failed: Descending sort of a double-precision array."
    end if

    ! Ensure the array is in descending order
    do i = 2, n
        if (x1(i) > x1(i-1)) then
            rst = .false.
            print '(A)', "Test Failed: Sorted array is not in descending order."
            exit
        end if
    end do

    ! End
end function

! ------------------------------------------------------------------------------
function test_descend_sort_cmplx() result(rst)
    ! Parameters
    integer(int32), parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer(int32) :: i, ind(n)
    complex(real64) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .false.
    call create_random_array(x1)
    do i = 1, n
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Use the LAPACK-based sorting routine as a benchmark
    call sort(x1, ascend)

    ! Now, track the indices
    call sort(x2, ind, ascend)

    ! Compare the two arrays - we really don't have a good means of checking
    ! ind at this point
    if (.not.assert(x1, x2, tol = REAL64_TOL)) then
        rst = .false.
        print '(A)', "Test Failed: Complex-Valued descending sort of a double-precision array."
    end if

    ! End
end function

! ------------------------------------------------------------------------------
end module
