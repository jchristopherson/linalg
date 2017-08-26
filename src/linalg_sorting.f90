! linalg_sorting.f90

!> @brief \b linalg_sorting
!!
!! @par Purpose
!! Provides sorting routines.
module linalg_sorting
    use ferror, only : errors
    use linalg_constants
    implicit none
    private
    public :: sort

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Sorts an array.
    interface sort
        module procedure :: sort_dbl_array
        module procedure :: sort_dbl_array_ind
    end interface


contains
! ******************************************************************************
! SORTING ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !!
    !! @par Remarks
    !! The routine utilizes a quick sort algorithm unless the size of the array
    !! is less than or equal to 20.  For such small arrays an insertion sort
    !! algorithm is utilized.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DLASRT.
    subroutine sort_dbl_array(x, ascend)
        ! Arguments
        real(dp), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: ascend

        ! Local Variables
        character :: id
        integer(i32) :: n, info

        ! Initialization
        if (present(ascend)) then
            if (ascend) then
                id = 'I'
            else
                id = 'D'
            end if
        else
            id = 'I'
        end if
        n = size(x)

        ! Process
        call DLASRT(id, n, x, info)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Sorts an array.
    !!
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in,out] ind On input, an integer array.  On output, the contents
    !!  of this array are shifted in the same order as that of @p x as a means
    !!  of tracking the sorting operation.  It is often useful to set this
    !!  array to an ascending group of values (1, 2, ... n) such that this
    !!  array tracks the original positions of the sorted array.  Such an array
    !!  can then be used to align other arrays.  This array must be the same
    !!  size as @p x.
    !! @param[in] ascend An optional input that, if specified, controls if the
    !!  the array is sorted in an ascending order (default), or a descending
    !!  order.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p ind is not sized to match @p x.
    !!
    !! @par Remarks
    !! This routine utilizes a quick sort algorithm explained at 
    !! http://www.fortran.com/qsort_c.f95.
    subroutine sort_dbl_array_ind(x, ind, ascend, err)
        ! Arguments
        real(dp), intent(inout), dimension(:) :: x
        integer(i32), intent(inout), dimension(:) :: ind
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(i32), parameter :: select = 20

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg
        integer(i32) :: n
        logical :: dir

        ! Initialization
        n = size(x)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (present(ascend)) then
            dir = ascend
        else
            dir = .true. ! Ascend == true
        end if

        ! Input Check
        if (size(ind) /= n) then
            write(errmsg, "(AI0AI0A)") &
                "Expected the tracking array to be of size ", n, &
                ", but found an array of size ", size(ind), "."
            call errmgr%report_error("sort_dbl_array_ind", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if
        if (n <= 1) return

        ! Process
        call qsort_dbl_ind(dir, x, ind)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! PRIVATE HELPER ROUTINES
! ------------------------------------------------------------------------------
    ! REF: http://www.fortran.com/qsort_c.f95
    recursive subroutine qsort_dbl_ind(ascend, x, ind)
        ! Arguments
        logical, intent(in) :: ascend
        real(dp), intent(inout), dimension(:) :: x
        integer(i32), intent(inout), dimension(:) :: ind

        ! Local Variables
        integer(i32) :: iq

        ! Process
        if (size(x) > 1) then
            call dbl_partition_ind(ascend, x, ind, iq)
            call qsort_dbl_ind(ascend, x(:iq-1), ind(:iq-1))
            call qsort_dbl_ind(ascend, x(iq:), ind(iq:))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !
    subroutine dbl_partition_ind(ascend, x, ind, marker)
        ! Arguments
        logical, intent(in) :: ascend
        real(dp), intent(inout), dimension(:) :: x
        integer(i32), intent(inout), dimension(:) :: ind
        integer(i32), intent(out) :: marker

        ! Local Variables
        integer(i32) :: i, j, itemp
        real(dp) :: temp, pivot

        ! Process
        pivot = x(1)
        i = 0
        j = size(x) + 1
        if (ascend) then
            ! Ascending Sort
            do
                j = j - 1
                do
                    if (x(j) <= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (x(i) >= pivot) exit
                    i = i + 1
                end do
                if (i < j) then
                    ! Exchage X(I) and X(J)
                    temp = x(i)
                    x(i) = x(j)
                    x(j) = temp

                    itemp = ind(i)
                    ind(i) = ind(j)
                    ind(j) = itemp
                else if (i == j) then
                    marker = i + 1
                    return
                else
                    marker = i
                    return
                end if
            end do
        else
            ! Descending Sort
            do
                j = j - 1
                do
                    if (x(j) >= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (x(i) <= pivot) exit
                    i = i + 1
                end do
                if (i < j) then
                    ! Exchage X(I) and X(J)
                    temp = x(i)
                    x(i) = x(j)
                    x(j) = temp

                    itemp = ind(i)
                    ind(i) = ind(j)
                    ind(j) = itemp
                else if (i == j) then
                    marker = i + 1
                    return
                else
                    marker = i
                    return
                end if
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
