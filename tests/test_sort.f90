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
    logical :: rst
    integer(i32) :: i, ind(n)
    real(dp) :: x1(n), x2(n)

    ! Initialization
    rst = .true.
    call random_number(x1)
    do i = 1, n
        x2(i) = x1(i)
        ind(i) = i
    end do

    ! Use the LAPACK-based sorting routine as a benchmark
    call sort(x1, .true.)

    ! Now, track the indices
    call sort(x2, ind, .true.)

    ! Compare the two arrays - we really don't have a good means of checking
    ! ind at this point
    if (.not.is_mtx_equal(x1, x2, tol)) then
        rst = .false.
        print '(A)', "Test Failed: Ascending sort of a double-precision array."
    end if

    ! End
end function

end module
