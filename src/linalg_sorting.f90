! linalg_sorting.f90

module linalg_sorting
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    use ferror
    implicit none
    private
    public :: sort
    
    interface sort
        !! An interface to the sorting routines.
        module procedure :: sort_dbl_array
        module procedure :: sort_dbl_array_ind
        module procedure :: sort_cmplx_array
        module procedure :: sort_cmplx_array_ind
        module procedure :: sort_eigen_cmplx
        module procedure :: sort_eigen_dbl
        module procedure :: sort_int32_array
        module procedure :: sort_int32_array_ind
    end interface

contains
! ******************************************************************************
! SORTING ROUTINES
! ------------------------------------------------------------------------------
subroutine sort_dbl_array(x, ascend)
    !! Sorts an array.
    real(real64), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.

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
subroutine sort_dbl_array_ind(x, ind, ascend, err)
    !! Sorts an array.
    real(real64), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    integer(int32), intent(inout), dimension(:) :: ind
        !! An array, the same size as x, that is sorted along with x.  This is
        !! often useful as a tracking array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
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
        call report_array_size_error("sort_dbl_array_ind", errmgr, "ind", &
            n, size(ind))
        return
    end if
    if (n <= 1) return

    ! Process
    call qsort_dbl_ind(dir, x, ind)
end subroutine

! ------------------------------------------------------------------------------
subroutine sort_cmplx_array(x, ascend)
    !! Sorts an array.
    complex(real64), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.

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
subroutine sort_cmplx_array_ind(x, ind, ascend, err)
    !! Sorts an array.
    complex(real64), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    integer(int32), intent(inout), dimension(:) :: ind
        !! An array, the same size as x, that is sorted along with x.  This is
        !! often useful as a tracking array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
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
        call report_array_size_error("sort_cmplx_array_ind", errmgr, "ind", &
            n, size(ind))
        return
    end if
    if (n <= 1) return

    ! Process
    call qsort_cmplx_ind(dir, x, ind)
end subroutine

! ------------------------------------------------------------------------------
subroutine sort_eigen_cmplx(vals, vecs, ascend, err)
    !! Sorts eigenvalues and their associated eigenvectors.
    complex(real64), intent(inout), dimension(:) :: vals
        !! On input, an N-element array containing the eigenvalues.  On output,
        !! the sored eigenvalues.
    complex(real64), intent(inout), dimension(:,:) :: vecs
        !! On input, the N-by-N matrix containing the eigenvectors (one vector
        !! per column) associated with vals.  On output, the sorted eigenvector
        !! matrix.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
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
        call report_matrix_size_error("sort_eigen_cmplx", errmgr, "vecs", &
            n, n, size(vecs, 1), size(vecs, 2))
    end if

    ! Allocate memory for the tracking array
    allocate(ind(n), stat = flag)
    if (flag /= 0) then
        call report_memory_error("sort_eigen_cmplx", errmgr, flag)
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
subroutine sort_eigen_dbl(vals, vecs, ascend, err)
    !! Sorts eigenvalues and their associated eigenvectors.
    real(real64), intent(inout), dimension(:) :: vals
        !! On input, an N-element array containing the eigenvalues.  On output,
        !! the sored eigenvalues.
    real(real64), intent(inout), dimension(:,:) :: vecs
        !! On input, the N-by-N matrix containing the eigenvectors (one vector
        !! per column) associated with vals.  On output, the sorted eigenvector
        !! matrix.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
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
        call report_matrix_size_error("sort_eigen_dbl", errmgr, "vecs", &
            n, n, size(vecs, 1), size(vecs, 2))
        return
    end if

    ! Allocate memory for the tracking array
    allocate(ind(n), stat = flag)
    if (flag /= 0) then
        call report_memory_error("sort_eigen_dbl", errmgr, flag)
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
    
! ------------------------------------------------------------------------------
subroutine sort_int32_array(x, ascend)
    !! Sorts an array.
    integer(int32), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.

    ! Local Variables
    logical :: dir

    ! Initialization
    if (present(ascend)) then
        dir = ascend
    else
        dir = .true.
    end if

    ! Process
    call qsort_int32(dir, x)
end subroutine

! ------------------------------------------------------------------------------
subroutine sort_int32_array_ind(x, ind, ascend, err)
    !! Sorts an array.
    integer(int32), intent(inout), dimension(:) :: x
        !! On input, the array to sort.  On output, the sorted array.
    integer(int32), intent(inout), dimension(:) :: ind
        !! An array, the same size as x, that is sorted along with x.  This is
        !! often useful as a tracking array.
    logical, intent(in), optional :: ascend
        !! An optional input that, if specified, controls if the array is 
        !! sorted in an ascending order (default), or a descending order.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
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
        call report_array_size_error("sort_int32_array_ind", errmgr, "ind", &
            n, size(ind))
        return
    end if
    if (n <= 1) return

    ! Process
    call qsort_int32_ind(dir, x, ind)
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
!> @brief A recursive quick sort algorithm.
!!
!! @param[in] ascend Set to true to sort in ascending order; else, false
!!  to sort in descending order.
!! @param[in,out] x On input, the array to sort.  On output, the sorted 
!!  array.
!!
!! @par Notes
!! This implementation is a slight modification of the code presented at
!! http://www.fortran.com/qsort_c.f95.
recursive subroutine qsort_int32(ascend, x)
    ! Arguments
    logical, intent(in) :: ascend
    integer(int32), intent(inout), dimension(:) :: x

    ! Local Variables
    integer(int32) :: iq

    ! Process
    if (size(x) > 1) then
        call int32_partition(ascend, x, iq)
        call qsort_int32(ascend, x(:iq-1))
        call qsort_int32(ascend, x(iq:))
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
!! @par Notes
!! This implementation is a slight modification of the code presented at
!! http://www.fortran.com/qsort_c.f95
subroutine int32_partition(ascend, x, marker)
    ! Arguments
    logical, intent(in) :: ascend
    integer(int32), intent(inout), dimension(:) :: x
    integer(int32), intent(out) :: marker

    ! Local Variables
    integer(int32) :: i, j, temp, pivot

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
!! @par Notes
!! This implementation is a slight modification of the code presented at
!! http://www.fortran.com/qsort_c.f95.
recursive subroutine qsort_int32_ind(ascend, x, ind)
    ! Arguments
    logical, intent(in) :: ascend
    integer(int32), intent(inout), dimension(:) :: x
    integer(int32), intent(inout), dimension(:) :: ind

    ! Local Variables
    integer(int32) :: iq

    ! Process
    if (size(x) > 1) then
        call int32_partition_ind(ascend, x, ind, iq)
        call qsort_int32_ind(ascend, x(:iq-1), ind(:iq-1))
        call qsort_int32_ind(ascend, x(iq:), ind(iq:))
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
subroutine int32_partition_ind(ascend, x, ind, marker)
    ! Arguments
    logical, intent(in) :: ascend
    integer(int32), intent(inout), dimension(:) :: x
    integer(int32), intent(inout), dimension(:) :: ind
    integer(int32), intent(out) :: marker

    ! Local Variables
    integer(int32) :: i, j, itemp, temp, pivot

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
end module
