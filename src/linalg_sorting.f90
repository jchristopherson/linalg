! linalg_sorting.f90

!> @brief \b linalg_sorting
!!
!! @par Purpose
!! Provides sorting routines.
submodule (linalg_core) linalg_sorting
contains
! ******************************************************************************
! SORTING ROUTINES
! ------------------------------------------------------------------------------
    module subroutine sort_dbl_array(x, ascend)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: ascend

        ! Local Variables
        character :: id
        integer(int32) :: n, info

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
    module subroutine sort_dbl_array_ind(x, ind, ascend, err)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg
        integer(int32) :: n
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
    module subroutine sort_cmplx_array(x, ascend)
        ! Arguments
        complex(real64), intent(inout), dimension(:) :: x
        logical, intent(in), optional :: ascend

        ! Local Variables
        logical :: dir

        ! Initialization
        if (present(ascend)) then
            dir = ascend
        else
            dir = .true.
        end if

        ! Process
        call qsort_cmplx(dir, x)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine sort_cmplx_array_ind(x, ind, ascend, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg
        integer(int32) :: n
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
            call errmgr%report_error("sort_cmplx_array_ind", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if
        if (n <= 1) return

        ! Process
        call qsort_cmplx_ind(dir, x, ind)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine sort_eigen_cmplx(vals, vecs, ascend, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:) :: vals
        complex(real64), intent(inout), dimension(:,:) :: vecs
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg
        integer(int32) :: i, n, flag
        logical :: dir
        integer(int32), allocatable, dimension(:) :: ind

        ! Initialization
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

        ! Ensure the eigenvector matrix is sized appropriately
        n = size(vals)
        if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
            ! ARRAY SIZE ERROR
            write(errmsg, '(AI0AI0AI0AI0A)') &
                "Expected the eigenvector matrix to be of size ", n, &
                "-by-", n, ", but found a matrix of size ", size(vecs, 1), &
                "-by-", size(vecs, 2), "."
            call errmgr%report_error("sort_eigen_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
        end if

        ! Allocate memory for the tracking array
        allocate(ind(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("sort_eigen_cmplx", &
                "Insufficient memory available.", LA_OUT_OF_MEMORY_ERROR)
            return
        end if
        do i = 1, n
            ind(i) = i
        end do

        ! Sort
        call qsort_cmplx_ind(dir, vals, ind)

        ! Shift the eigenvectors around to keep them associated with the
        ! appropriate eigenvalue
        vecs = vecs(:,ind)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine sort_eigen_dbl(vals, vecs, ascend, err)
        ! Arguments
        real(real64), intent(inout), dimension(:) :: vals
        real(real64), intent(inout), dimension(:,:) :: vecs
        logical, intent(in), optional :: ascend
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg
        integer(int32) :: i, n, flag
        logical :: dir
        integer(int32), allocatable, dimension(:) :: ind

        ! Initialization
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

        ! Ensure the eigenvector matrix is sized appropriately
        n = size(vals)
        if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
            ! ARRAY SIZE ERROR
            write(errmsg, '(AI0AI0AI0AI0A)') &
                "Expected the eigenvector matrix to be of size ", n, &
                "-by-", n, ", but found a matrix of size ", size(vecs, 1), &
                "-by-", size(vecs, 2), "."
            call errmgr%report_error("sort_eigen_dbl", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
        end if

        ! Allocate memory for the tracking array
        allocate(ind(n), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("sort_eigen_dbl", &
                "Insufficient memory available.", LA_OUT_OF_MEMORY_ERROR)
            return
        end if
        do i = 1, n
            ind(i) = i
        end do

        ! Sort
        call qsort_dbl_ind(dir, vals, ind)

        ! Shift the eigenvectors around to keep them associated with the
        ! appropriate eigenvalue
        vecs = vecs(:,ind)
    end subroutine

! ******************************************************************************
! PRIVATE HELPER ROUTINES
! ------------------------------------------------------------------------------
    !> @brief A recursive quick sort algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in,out] ind On input, a tracking array of the same length as @p x.
    !!  On output, the same array, but shuffled to match the sorting order of
    !!  @p x.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95.
    recursive subroutine qsort_dbl_ind(ascend, x, ind)
        ! Arguments
        logical, intent(in) :: ascend
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind

        ! Local Variables
        integer(int32) :: iq

        ! Process
        if (size(x) > 1) then
            call dbl_partition_ind(ascend, x, ind, iq)
            call qsort_dbl_ind(ascend, x(:iq-1), ind(:iq-1))
            call qsort_dbl_ind(ascend, x(iq:), ind(iq:))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief A routine to perform the partioning necessary for the quick sort
    !! algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in,out] ind On input, a tracking array of the same length as @p x.
    !!  On output, the same array, but shuffled to match the sorting order of
    !!  @p x.
    !! @param[out] marker The partioning marker.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95
    subroutine dbl_partition_ind(ascend, x, ind, marker)
        ! Arguments
        logical, intent(in) :: ascend
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        integer(int32), intent(out) :: marker

        ! Local Variables
        integer(int32) :: i, j, itemp
        real(real64) :: temp, pivot

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
    !> @brief A recursive quick sort algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !!
    !! @par Remarks
    !! As this routine operates on complex valued items, the complex values are
    !! sorted based upon the real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95
    recursive subroutine qsort_cmplx(ascend, x)
        ! Arguments
        logical, intent(in) :: ascend
        complex(real64), intent(inout), dimension(:) :: x

        ! Local Variables
        integer(int32) :: iq

        ! Process
        if (size(x) > 1) then
            call cmplx_partition(ascend, x, iq)
            call qsort_cmplx(ascend, x(:iq-1))
            call qsort_cmplx(ascend, x(iq:))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief A routine to perform the partioning necessary for the quick sort
    !! algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[out] marker The partioning marker.
    !!
    !! @par Remarks
    !! As this routine operates on complex valued items, the complex values are
    !! sorted based upon the real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95.
    subroutine cmplx_partition(ascend, x, marker)
        ! Arguments
        logical, intent(in) :: ascend
        complex(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(out) :: marker

        ! Local Variables
        integer(int32) :: i, j
        complex(real64) :: temp
        real(real64) :: pivot

        ! Process
        pivot = real(x(1), real64)
        i = 0
        j = size(x) + 1
        if (ascend) then
            ! Ascending Sort
            do
                j = j - 1
                do
                    if (real(x(j), real64) <= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (real(x(i), real64) >= pivot) exit
                    i = i + 1
                end do
                if (i < j) then
                    ! Exchage X(I) and X(J)
                    temp = x(i)
                    x(i) = x(j)
                    x(j) = temp
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
                    if (real(x(j), real64) >= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (real(x(i), real64) <= pivot) exit
                    i = i + 1
                end do
                if (i < j) then
                    ! Exchage X(I) and X(J)
                    temp = x(i)
                    x(i) = x(j)
                    x(j) = temp
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
    !> @brief A recursive quick sort algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in,out] ind On input, a tracking array of the same length as @p x.
    !!  On output, the same array, but shuffled to match the sorting order of
    !!  @p x.
    !!
    !! @par Remarks
    !! As this routine operates on complex valued items, the complex values are
    !! sorted based upon the real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95
    recursive subroutine qsort_cmplx_ind(ascend, x, ind)
        ! Arguments
        logical, intent(in) :: ascend
        complex(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind

        ! Local Variables
        integer(int32) :: iq

        ! Process
        if (size(x) > 1) then
            call cmplx_partition_ind(ascend, x, ind, iq)
            call qsort_cmplx_ind(ascend, x(:iq-1), ind(:iq-1))
            call qsort_cmplx_ind(ascend, x(iq:), ind(iq:))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief A routine to perform the partioning necessary for the quick sort
    !! algorithm.
    !!
    !! @param[in] ascend Set to true to sort in ascending order; else, false
    !!  to sort in descending order.
    !! @param[in,out] x On input, the array to sort.  On output, the sorted 
    !!  array.
    !! @param[in,out] ind On input, a tracking array of the same length as @p x.
    !!  On output, the same array, but shuffled to match the sorting order of
    !!  @p x.
    !! @param[out] marker The partioning marker.
    !!
    !! @par Remarks
    !! As this routine operates on complex valued items, the complex values are
    !! sorted based upon the real component of the number.
    !!
    !! @par Notes
    !! This implementation is a slight modification of the code presented at
    !! http://www.fortran.com/qsort_c.f95.
    subroutine cmplx_partition_ind(ascend, x, ind, marker)
        ! Arguments
        logical, intent(in) :: ascend
        complex(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(inout), dimension(:) :: ind
        integer(int32), intent(out) :: marker

        ! Local Variables
        integer(int32) :: i, j, itemp
        complex(real64) :: temp
        real(real64) :: pivot

        ! Process
        pivot = real(x(1), real64)
        i = 0
        j = size(x) + 1
        if (ascend) then
            ! Ascending Sort
            do
                j = j - 1
                do
                    if (real(x(j), real64) <= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (real(x(i), real64) >= pivot) exit
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
                    if (real(x(j), real64) >= pivot) exit
                    j = j - 1
                end do
                i = i + 1
                do
                    if (real(x(i), real64) <= pivot) exit
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
end submodule
