! test_sort.f90

! A module containing routines to support testing of the sorting routines.
module test_sort
    use test_core
    use linalg_constants
    use linalg_sorting
    implicit none

contains
! ******************************************************************************
! DOUBLE PRECISION TESTS
! ------------------------------------------------------------------------------
function test_dbl_ascend_sort() result(rst)
    ! Parameters
    integer(i32), parameter :: n = 200
    real(dp), parameter :: tol = 1.0d-8

    ! Local Variables
    logical :: rst, ascend
    integer(i32) :: i, ind(n)
    real(dp) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .true.
    call random_number(x1)
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
    if (.not.is_mtx_equal(x1, x2, tol)) then
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
function test_dbl_descend_sort() result(rst)
    ! Parameters
    integer(i32), parameter :: n = 200
    real(dp), parameter :: tol = 1.0d-8

    ! Local Variables
    logical :: rst, ascend
    integer(i32) :: i, ind(n)
    real(dp) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    ascend = .false.
    call random_number(x1)
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
    if (.not.is_mtx_equal(x1, x2, tol)) then
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

end module
