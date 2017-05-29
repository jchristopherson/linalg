! rz.f90

!> @brief \b rz
!!
!! @par Purpose
!! Provides a set of routines supporting RZ factorization of trapezoidal 
!! systems.
module rz
    use ferror, only : errors
    use linalg_constants
    implicit none
    private
    public :: rz_factor
    public :: mult_rz

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Z from an
    !! RZ factorization.
    interface mult_rz
        module procedure :: mult_rz_mtx
        module procedure :: mult_rz_vec
    end interface


contains
! ******************************************************************************
! RZ FACTORIZATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Factors an upper trapezoidal matrix by means of orthogonal
    !! transformations such that A = R * Z = (R 0) * Z.  Z is an orthogonal
    !! matrix of dimension (M+L)-by-(M+L), and R is an M-by-M upper triangular
    !! matrix.
    !!
    !! @param[in,out] a On input, the M-by-N upper trapezoidal matrix to factor.
    !!  On output, the leading M-by-M upper triangular part of the matrix
    !!  contains the upper triangular matrix R, and elements N-L+1 to N of the
    !!  first M rows of A, with the array @p tau, represent the orthogonal
    !!  matrix Z as a product of M elementary reflectors.
    !! @param[out] tau
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Further Details
    !! @verbatim
    !!  The factorization is obtained by Householder's method.  The kth
    !!  transformation matrix, Z( k ), which is used to introduce zeros into
    !!  the ( m - k + 1 )th row of A, is given in the form
    !!
    !!     Z( k ) = ( I     0   ),
    !!              ( 0  T( k ) )
    !!
    !!  where
    !!
    !!     T( k ) = I - tau*u( k )*u( k )**T,   u( k ) = (   1    ),
    !!                                                   (   0    )
    !!                                                   ( z( k ) )
    !!
    !!  tau is a scalar and z( k ) is an l element vector. tau and z( k )
    !!  are chosen to annihilate the elements of the kth row of A2.
    !!
    !!  The scalar tau is returned in the kth element of TAU and the vector
    !!  u( k ) in the kth row of A2, such that the elements of z( k ) are
    !!  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in
    !!  the upper triangular part of A1.
    !!
    !!  Z is given by
    !!
    !!     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
    !! @endverbatim
    !!
    !! @par Notes
    !! This routine is based upon the LAPACK routine DTZRZF.
    !!
    !! @par See Also
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node44.html)
    subroutine rz_factor(a, tau, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, lwork, flag, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("rz_factor", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DTZRZF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(work(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("rz_factor", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("rz_factor", &
                        "Incorrectly sized input array WORK, argument 3.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
        else
            allocate(work(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("rz_factor", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DTZRZF
        call DTZRZF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Z from an
    !! RZ factorization such that: C = op(Z) * C, or C = C * op(Z).
    !!
    !! @param[in] lside Set to true to apply Z or Z**T from the left; else, set
    !!  to false to apply Z or Z**T from the right.
    !! @param[in] trans Set to true to apply Z**T; else, set to false.
    !! @param[in] l The number of columns in matrix @p a containing the
    !!  meaningful part of the Householder vectors.  If @p lside is true,
    !!  M >= L >= 0; else, if @p lside is false, N >= L >= 0.
    !! @param[in,out] a On input the K-by-LTA matrix C, where LTA = M if
    !!  @p lside is true; else, LTA = N if @p lside is false.  The I-th row must
    !!  contain the Householder vector in the last k rows. Notice, the contents
    !!  of this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of the
    !!  elementary reflectors, where M >= K >= 0 if @p lside is true; else,
    !!  N >= K >= 0 if @p lside is false.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Z and the original matrix C.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DORMRZ.
    subroutine mult_rz_mtx(lside, trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        integer(i32), intent(in) :: l
        real(dp), intent(inout), dimension(:,:) :: a, c
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(i32) :: m, n, k, lwork, flag, istat, lda
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        lda = size(a, 1)
        if (lside) then
            side = 'L'
        else
            side = 'R'
        end if
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (lside) then
            if (l > m .or. l < 0) then
               flag = 3
            else if (k > m) then
                flag = 5
            else if (size(a, 1) < k .or. size(a, 2) /= m) then
                flag = 4
            end if
        else
            if (l > n .or. l < 0) then
                flag = 3
            else if (k > n) then
                flag = 5
            else if (size(a, 1) < k .or. size(a, 2) /= n) then
                flag = 4
            end if
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(work(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("mult_rz_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("mult_rz_mtx", &
                        "Incorrectly sized input array WORK, argument 7.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
        else
            allocate(work(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMRZ
        call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a vector by the orthogonal matrix Z from an
    !! RZ factorization such that: C = op(Z) * C.
    !!
    !! @param[in] trans Set to true to apply Z**T; else, set to false.
    !! @param[in] l The number of columns in matrix @p a containing the
    !!  meaningful part of the Householder vectors.  If @p lside is true,
    !!  M >= L >= 0; else, if @p lside is false, N >= L >= 0.
    !! @param[in,out] a On input the K-by-LTA matrix C, where LTA = M if
    !!  @p lside is true; else, LTA = N if @p lside is false.  The I-th row must
    !!  contain the Householder vector in the last k rows. Notice, the contents
    !!  of this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of the
    !!  elementary reflectors, where M >= K >= 0 if @p lside is true; else,
    !!  N >= K >= 0 if @p lside is false.
    !! @param[in,out] c On input, the M-element array C.  On output, the product
    !!  of the orthogonal matrix Z and the original array C.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DORMRZ.
    subroutine mult_rz_vec(trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        integer(i32), intent(in) :: l
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: c
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(i32) :: m, k, lwork, flag, istat, lda
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        lda = size(a, 1)
        side = 'L'
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (l > m .or. l < 0) then
            flag = 2
        else if (k > m) then
            flag = 4
        else if (size(a, 1) < k .or. size(a, 2) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(work(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("mult_rz_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("mult_rz_vec", &
                        "Incorrectly sized input array WORK, argument 6.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
        else
            allocate(work(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMRZ
        call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
end module
