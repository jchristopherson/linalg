! test_sort.f90
! SPDX-License-Identifier: MIT

! A module containing routines to support testing of the sorting routines.
module test_sort
    use, intrinsic :: iso_fortran_env, only : real64
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
    integer, parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer :: i, ind(n)
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
    integer, parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer :: i, ind(n)
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
    integer, parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer :: i, ind(n)
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
    integer, parameter :: n = 200

    ! Local Variables
    logical :: rst, ascend
    integer :: i, ind(n)
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
function test_int32_ascend_sort() result(rst)
    ! Parameters
    integer, parameter :: n = 200
    integer, parameter :: first = -10000
    integer, parameter :: last = 10000

    ! Local Variables
    logical :: rst, ascend
    integer :: i, ind(n), x1(n), x2(n), ri
    real(real64) :: rv

    ! Initialization
    rst = .true.
    ascend = .true.
    do i = 1, n
        call random_number(rv)
        ri = first + floor((last + 1 - first) * rv)
        x1(i) = ri + 1
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Sort into ascending order
    call sort(x1, ascend)

    ! Ensure it's increasing
    do i = 2, n
        if (x1(i) < x1(i-1)) then
            rst = .false.
            print "(A)", "Test Failed: Ascending sort of 32-bit integer array."
            exit
        end if
    end do

    ! Sort the other array
    call sort(x2, ind, ascend)

    ! Compare
    if (.not.assert(x1, x2)) then
        rst = .false.
        print "(A)", "Test Failed: Sorted integer arrays do not match."
    end if
end function

! ------------------------------------------------------------------------------
function test_sort_index_edges() result(rst)
    real(real64) :: original(6), values(6)
    complex(real64) :: coriginal(6), cvalues(6)
    integer :: indices(6), i, direction
    logical :: seen(6), ascend, rst

    rst = .true.
    original = [3.0_real64, 1.0_real64, 2.0_real64, &
        1.0_real64, -1.0_real64, 3.0_real64]
    coriginal = [(3.0_real64, 1.0_real64), (1.0_real64, -2.0_real64), &
        (2.0_real64, 0.5_real64), (1.0_real64, 4.0_real64), &
        (-1.0_real64, 0.0_real64), (3.0_real64, -1.0_real64)]

    do direction = 1, 2
        ascend = direction == 1
        values = original
        indices = [(i, i = 1, size(indices))]
        call sort(values, indices, ascend)
        seen = .false.
        do i = 1, size(values)
            if (indices(i) < 1 .or. indices(i) > size(original)) then
                rst = .false.
                cycle
            end if
            if (seen(indices(i)) .or. values(i) /= original(indices(i))) rst = .false.
            seen(indices(i)) = .true.
            if (i > 1) then
                if ((ascend .and. values(i) < values(i-1)) .or. &
                    (.not.ascend .and. values(i) > values(i-1))) rst = .false.
            end if
        end do
        if (.not.all(seen)) rst = .false.

        cvalues = coriginal
        indices = [(i, i = 1, size(indices))]
        call sort(cvalues, indices, ascend)
        seen = .false.
        do i = 1, size(cvalues)
            if (indices(i) < 1 .or. indices(i) > size(coriginal)) then
                rst = .false.
                cycle
            end if
            if (seen(indices(i)) .or. cvalues(i) /= coriginal(indices(i))) rst = .false.
            seen(indices(i)) = .true.
            if (i > 1) then
                if ((ascend .and. real(cvalues(i)) < real(cvalues(i-1))) .or. &
                    (.not.ascend .and. real(cvalues(i)) > real(cvalues(i-1)))) rst = .false.
            end if
        end do
        if (.not.all(seen)) rst = .false.
    end do

    if (.not.rst) print '(A)', "Test Failed: Sorting Indices with Duplicate Values"
end function

! ------------------------------------------------------------------------------
end module
