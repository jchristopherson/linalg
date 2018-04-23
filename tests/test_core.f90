! test_core.f90

! A module containing routines to support basic testing operations.
module test_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    implicit none
    private
    public :: is_mtx_equal

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface is_mtx_equal
        module procedure :: is_mtx_equal_double
        module procedure :: is_vec_equal_double
        module procedure :: is_vec_equal_i32
        module procedure :: is_mtx_equal_complex
        module procedure :: is_vec_equal_complex
    end interface

contains
! ******************************************************************************
! MATRIX COMPARISON TESTS
! ------------------------------------------------------------------------------
function is_mtx_equal_double(x, y, tol) result(check)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: x, y
    real(real64), intent(in) :: tol
    logical :: check

    ! Local Variables
    integer(int32) :: i, j, m, n

    ! Initialization
    m = size(x, 1)
    n = size(x, 2)
    check = .true.

    ! Process
    if (size(y, 1) /= m .or. size(y, 2) /= n) then
        check = .false.
        return
    end if
    do j = 1, n
        do i = 1, m
            if (abs(x(i,j) - y(i,j)) > tol) then
                check = .false.
                return
            end if
        end do
    end do
end function

! ------------------------------------------------------------------------------
function is_vec_equal_double(x, y, tol) result(check)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, y
    real(real64), intent(in) :: tol
    logical :: check

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x)
    check = .true.

    ! Process
    if (size(y) /= n) then
        check = .false.
        return
    end if
    do i = 1, n
        if (abs(x(i) - y(i)) > tol) then
            check = .false.
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
function is_vec_equal_i32(x, y) result(check)
    ! Arguments
    integer(int32), intent(in), dimension(:) :: x, y
    logical :: check

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x)
    check = .true.

    ! Process
    if (size(y) /= n) then
        check = .false.
        return
    end if
    do i = 1, n
        if (x(i) /= y(i)) then
            check = .false.
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
function is_mtx_equal_complex(x, y, tol) result(check)
    ! Arguments
    complex(real64), intent(in), dimension(:,:) :: x, y
    real(real64), intent(in) :: tol
    logical :: check

    ! Local Variables
    integer(int32) :: i, j, m, n

    ! Initialization
    m = size(x, 1)
    n = size(x, 2)
    check = .true.

    ! Process
    if (size(y, 1) /= m .or. size(y, 2) /= n) then
        check = .false.
        return
    end if
    do j = 1, n
        do i = 1, m
            if (abs(real(x(i,j), real64) - real(y(i,j), real64)) > tol .or. &
                abs(aimag(x(i,j)) - aimag(y(i,j))) > tol) then
                check = .false.
                return
            end if
        end do
    end do
end function

! ------------------------------------------------------------------------------
function is_vec_equal_complex(x, y, tol) result(check)
    ! Arguments
    complex(real64), intent(in), dimension(:) :: x, y
    real(real64), intent(in) :: tol
    logical :: check

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x)
    check = .true.

    ! Process
    if (size(y) /= n) then
        check = .false.
        return
    end if
    do i = 1, n
        if (abs(real(x(i), real64) - real(y(i), real64)) > tol .or. &
            abs(aimag(x(i)) - aimag(y(i))) > tol) then
            check = .false.
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
end module
