! linalg_factor.f90

!> @brief \b linalg_factor
!!
!! @par Purpose
!! Provides a set of matrix factorization routines.
submodule (linalg) linalg_factor
    use lapack
    use qrupdate
    implicit none
contains
! ******************************************************************************
! LU FACTORIZATION
! ------------------------------------------------------------------------------
    module subroutine lu_factor_dbl(a, ipvt, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(ipvt) /= mn) then
            ! ERROR: IPVT not sized correctly
            call errmgr%report_error("lu_factor_dbl", &
                "Incorrectly sized input array IPVT, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Compute the LU factorization by calling the LAPACK routine DGETRF
        call DGETRF(m, n, a, m, ipvt, flag)

        ! If flag > 0, the matrix is singular.  Notice, flag should not be
        ! able to be < 0 as we've already verrified inputs prior to making the
        ! call to LAPACK
        if (flag > 0) then
            ! WARNING: Singular matrix
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Singular matrix encountered (row ", flag, ")"
            call errmgr%report_warning("lu_factor_dbl", trim(errmsg), &
                LA_SINGULAR_MATRIX_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine lu_factor_cmplx(a, ipvt, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), dimension(:) :: ipvt
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(ipvt) /= mn) then
            ! ERROR: IPVT not sized correctly
            call errmgr%report_error("lu_factor_cmplx", &
                "Incorrectly sized input array IPVT, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Compute the LU factorization by calling the LAPACK routine ZGETRF
        call ZGETRF(m, n, a, m, ipvt, flag)

        ! If flag > 0, the matrix is singular.  Notice, flag should not be
        ! able to be < 0 as we've already verrified inputs prior to making the
        ! call to LAPACK
        if (flag > 0) then
            ! WARNING: Singular matrix
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) &
                "Singular matrix encountered (row ", flag, ")"
            call errmgr%report_warning("lu_factor_cmplx", trim(errmsg), &
                LA_SINGULAR_MATRIX_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lu_all(lu, ipvt, u, p, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: lu
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(out), dimension(:,:) :: u, p
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: j, jp, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_all", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Ensure P starts off as an identity matrix
        call DLASET('A', n, n, zero, one, p, n)

        ! Process
        do j = 1, n
            ! Define the pivot matrix
            jp = ipvt(j)
            if (j /= jp) call swap(p(j,1:n), p(jp,1:n))

            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = zero

            if (j > 1) lu(1:j-1,j) = zero
            lu(j,j) = one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lu_all_cmplx(lu, ipvt, u, p, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: lu
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(out), dimension(:,:) :: u
        real(real64), intent(out), dimension(:,:) :: p
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: j, jp, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        complex(real64), parameter :: c_zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: c_one = (1.0d0, 0.0d0)

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_all_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Ensure P starts off as an identity matrix
        call DLASET('A', n, n, zero, one, p, n)

        ! Process
        do j = 1, n
            ! Define the pivot matrix
            jp = ipvt(j)
            if (j /= jp) call swap(p(j,1:n), p(jp,1:n))

            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = c_zero

            if (j > 1) lu(1:j-1,j) = c_zero
            lu(j,j) = c_one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lu_only(lu, u, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: lu
        real(real64), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: j, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_only", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = zero

            if (j > 1) lu(1:j-1,j) = zero
            lu(j,j) = one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lu_only_cmplx(lu, u, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: lu
        complex(real64), intent(out), dimension(:,:) :: u
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: j, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Initialization
        n = size(lu, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(lu, 2) /= n) then
            flag = 2
        else if (size(u, 1) /= n .or. size(u, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lu_only_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        do j = 1, n
            ! Build L and U
            u(1:j,j) = lu(1:j,j)
            u(j+1:n,j) = zero

            if (j > 1) lu(1:j-1,j) = zero
            lu(j,j) = one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ******************************************************************************
! QR FACTORIZATION
! ------------------------------------------------------------------------------
    module subroutine qr_factor_no_pivot(a, tau, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        real(real64), dimension(1) :: temp
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(tau) /= mn) then
            ! ERROR: TAU not sized correctly
            call errmgr%report_error("qr_factor_no_pivot", &
                "Incorrectly sized input array TAU, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEQRF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_no_pivot", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGEQRF
        call DGEQRF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine qr_factor_no_pivot_cmplx(a, tau, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: tau
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        complex(real64), dimension(1) :: temp
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(tau) /= mn) then
            ! ERROR: TAU not sized correctly
            call errmgr%report_error("qr_factor_no_pivot_cmplx", &
                "Incorrectly sized input array TAU, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZGEQRF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_no_pivot_cmplx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_no_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZGEQRF
        call ZGEQRF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine qr_factor_pivot(a, tau, jpvt, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        integer(int32), intent(inout), dimension(:) :: jpvt
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        real(real64), dimension(1) :: temp
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_factor_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEQP3(m, n, a, m, jpvt, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_pivot", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGEQP3
        call DGEQP3(m, n, a, m, jpvt, tau, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine qr_factor_pivot_cmplx(a, tau, jpvt, work, olwork, rwork, &
            err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: tau
        integer(int32), intent(inout), dimension(:) :: jpvt
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        real(real64), intent(out), target, dimension(:), optional :: rwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        complex(real64), dimension(1) :: temp
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        real(real64), pointer, dimension(:) :: rptr
        real(real64), allocatable, target, dimension(:) :: rwrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(jpvt) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_factor_pivot_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if
        if (present(rwork)) then
            if (size(rwork) < 2 * n) then
                call errmgr%report_error("qr_factor_pivot_cmplx", &
                    "Incorrectly sized input array RWORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            rptr => rwork(1:2*n)
        else
            allocate(rwrk(2 * n), stat = flag)
            if (flag /= 0) then
                call errmgr%report_error("qr_factor_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            rptr => rwrk
        end if

        ! Workspace Query
        call ZGEQP3(m, n, a, m, jpvt, tau, temp, -1, rptr, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_factor_pivot_cmplx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_factor_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZGEQP3
        call ZGEQP3(m, n, a, m, jpvt, tau, wptr, lwork, rptr, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_qr_no_pivot(r, tau, q, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), dimension(:,:) :: q
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, m, n, mn, qcol, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        qcol = size(q, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
            flag = 3
        else if (qcol == n .and. m < n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_no_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("form_qr_no_pivot", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("form_qr_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Copy the sub-diagonal portion of R to Q, and then zero out the
        ! sub-diagonal portion of R
        do j = 1, mn
            q(j+1:m,j) = r(j+1:m,j)
            r(j+1:m,j) = zero
        end do

        ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
        call DORGQR(m, qcol, mn, q, m, tau, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_qr_no_pivot_cmplx(r, tau, q, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: r
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(out), dimension(:,:) :: q
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: j, m, n, mn, qcol, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        qcol = size(q, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
            flag = 3
        else if (qcol == n .and. m < n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_no_pivot_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("form_qr_no_pivot_cmplx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("form_qr_no_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Copy the sub-diagonal portion of R to Q, and then zero out the
        ! sub-diagonal portion of R
        do j = 1, mn
            q(j+1:m,j) = r(j+1:m,j)
            r(j+1:m,j) = zero
        end do

        ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
        call ZUNGQR(m, qcol, mn, q, m, tau, wptr, lwork, flag)

        ! End
        if (allocated(wrk)) deallocate(wrk)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_qr_pivot(r, tau, pvt, q, p, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: pvt
        real(real64), intent(out), dimension(:,:) :: q, p
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: j, jp, m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(pvt) /= n) then
            flag = 3
        else if (size(q, 1) /= m .or. &
            (size(q, 2) /= m .and. size(q, 2) /= n)) then
            flag = 4
        else if (size(q, 2) == n .and. m < n) then
            flag = 4
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 5
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Generate Q and R
        call form_qr_no_pivot(r, tau, q, work = work, olwork = olwork, &
            err = errmgr)
        if (present(olwork)) return ! Just a workspace query
        if (errmgr%has_error_occurred()) return

        ! Form P
        do j = 1, n
            jp = pvt(j)
            p(:,j) = zero
            p(jp,j) = one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_qr_pivot_cmplx(r, tau, pvt, q, p, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: r
        complex(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: pvt
        complex(real64), intent(out), dimension(:,:) :: q, p
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: j, jp, m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(r, 1)
        n = size(r, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(tau) /= mn) then
            flag = 2
        else if (size(pvt) /= n) then
            flag = 3
        else if (size(q, 1) /= m .or. &
            (size(q, 2) /= m .and. size(q, 2) /= n)) then
            flag = 4
        else if (size(q, 2) == n .and. m < n) then
            flag = 4
        else if (size(p, 1) /= n .or. size(p, 2) /= n) then
            flag = 5
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_pivot_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Generate Q and R
        call form_qr_no_pivot_cmplx(r, tau, q, work = work, olwork = olwork, &
            err = errmgr)
        if (present(olwork)) return ! Just a workspace query
        if (errmgr%has_error_occurred()) return

        ! Form P
        do j = 1, n
            jp = pvt(j)
            p(:,j) = zero
            p(jp,j) = one
        end do

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_qr_mtx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:,:) :: a, c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, nrowa, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        if (lside) then
            side = 'L'
            nrowa = m
        else
            side = 'R'
            nrowa = n
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
            ! A is M-by-K, M >= K >= 0
            if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        else
            ! A is N-by-K, N >= K >= 0
            if (size(a, 1) /= n .or. size(a, 2) < k) flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_mtx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMQR
        call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_qr_mtx_cmplx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(inout), dimension(:,:) :: a, c
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, nrowa, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        if (lside) then
            side = 'L'
            nrowa = m
        else
            side = 'R'
            nrowa = n
        end if
        if (trans) then
            t = 'C'
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
            ! A is M-by-K, M >= K >= 0
            if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        else
            ! A is N-by-K, N >= K >= 0
            if (size(a, 1) /= n .or. size(a, 2) < k) flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_mtx_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_mtx_cmplx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_mtx_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMQR
        call ZUNMQR(side, t, m, n, k, a, nrowa, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_qr_vec(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(int32) :: m, k, nrowa, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        side = 'L'
        nrowa = m
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
        if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMQR
        call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_qr_vec_cmplx(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(inout), dimension(:) :: c
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        character :: side, t
        integer(int32) :: m, k, nrowa, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        side = 'L'
        nrowa = m
        if (trans) then
            t = 'C'
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
        if (size(a, 1) /= m .or. size(a, 2) < k) flag = 3
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_qr_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_qr_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMQR
        call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine qr_rank1_update_dbl(q, r, u, v, work, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: q, r
        real(real64), intent(inout), dimension(:) :: u, v
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        logical :: full
        integer(int32) :: m, n, k, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(u, 1)
        n = size(r, 2)
        k = min(m, n)
        full = size(q, 2) == m
        lwork = 2 * k
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (.not.full .and. size(q, 2) /= k) then
            flag = 1
        else if (size(r, 1) /= m) then
            flag = 2
        else if (size(u) /= m) then
            flag = 3
        else if (size(v) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_rank1_update", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DQR1UP(m, n, k, q, m, r, m, u, v, wptr)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine qr_rank1_update_cmplx(q, r, u, v, work, rwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: q, r
        complex(real64), intent(inout), dimension(:) :: u, v
        complex(real64), intent(out), target, optional, dimension(:) :: work
        real(real64), intent(out), target, optional, dimension(:) :: rwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        logical :: full
        integer(int32) :: m, n, k, lwork, istat, flag, lrwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        real(real64), pointer, dimension(:) :: rwptr
        real(real64), allocatable, target, dimension(:) :: rwrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(u, 1)
        n = size(r, 2)
        k = min(m, n)
        full = size(q, 2) == m
        lwork = k
        lrwork = k
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (.not.full .and. size(q, 2) /= k) then
            flag = 1
        else if (size(r, 1) /= m) then
            flag = 2
        else if (size(u) /= m) then
            flag = 3
        else if (size(v) /= n) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_rank1_update_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_rank1_update_cmplx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_rank1_update_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        if (present(rwork)) then
            if (size(rwork) < lrwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("qr_rank1_update_cmplx", &
                    "Incorrectly sized input array RWORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lrwork)
        else
            allocate(rwrk(lrwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("qr_rank1_update_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            rwptr => rwrk
        end if

        ! Process
        call ZQR1UP(m, n, k, q, m, r, m, u, v, wptr, rwptr)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ******************************************************************************
! CHOLESKY FACTORIZATION
! ------------------------------------------------------------------------------
    module subroutine cholesky_factor_dbl(a, upper, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo
        integer(int32) :: i, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(upper)) then
            if (upper) then
                uplo = 'U'
            else
                uplo = 'L'
            end if
        else
            uplo = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("cholesky_factor", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRF(uplo, n, a, n, flag)
        if (flag > 0) then
            ! ERROR: Matrix is not positive definite
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "The leading minor of order ", flag, &
                " is not positive definite."
            call errmgr%report_error("cholesky_factor", trim(errmsg), &
                LA_MATRIX_FORMAT_ERROR)
        end if

        ! Zero out the non-used upper or lower diagonal
        if (uplo == 'U') then
            ! Zero out the lower
            do i = 1, n - 1
                a(i+1:n,i) = zero
            end do
        else
            ! Zero out the upper
            do i = 2, n
                a(1:i-1,i) = zero
            end do
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cholesky_factor_cmplx(a, upper, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        logical, intent(in), optional :: upper
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        character :: uplo
        integer(int32) :: i, n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(upper)) then
            if (upper) then
                uplo = 'U'
            else
                uplo = 'L'
            end if
        else
            uplo = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("cholesky_factor_cmplx", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call ZPOTRF(uplo, n, a, n, flag)
        if (flag > 0) then
            ! ERROR: Matrix is not positive definite
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "The leading minor of order ", flag, &
                " is not positive definite."
            call errmgr%report_error("cholesky_factor_cmplx", trim(errmsg), &
                LA_MATRIX_FORMAT_ERROR)
        end if

        ! Zero out the non-used upper or lower diagonal
        if (uplo == 'U') then
            ! Zero out the lower
            do i = 1, n - 1
                a(i+1:n,i) = zero
            end do
        else
            ! Zero out the upper
            do i = 2, n
                a(1:i-1,i) = zero
            end do
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cholesky_rank1_update_dbl(r, u, work, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_update", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DCH1UP(n, r, n, u, wptr)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cholesky_rank1_update_cmplx(r, u, work, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: r
        complex(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_update_cmplx", &
                trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_update_cmplx", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_update", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call ZCH1UP(n, r, n, u, wptr)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cholesky_rank1_downdate_dbl(r, u, work, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        real(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_downdate", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_downdate", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_downdate", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DCH1DN(n, r, n, u, wptr, flag)
        if (flag == 1) then
            ! ERROR: The matrix is not positive definite
            call errmgr%report_error("cholesky_rank1_downdate", &
                "The downdated matrix is not positive definite.", &
                LA_MATRIX_FORMAT_ERROR)
        else if (flag == 2) then
            ! ERROR: The matrix is singular
            call errmgr%report_error("cholesky_rank1_downdate", &
                "The input matrix is singular.", LA_SINGULAR_MATRIX_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine cholesky_rank1_downdate_cmplx(r, u, work, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: r
        complex(real64), intent(inout), dimension(:) :: u
        real(real64), intent(out), target, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        n = size(r, 1)
        lwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(r, 2) /= n) then
            flag = 1
        else if (size(u) /= n) then
            flag = 2
        end if
        if (flag /= 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
                trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: Workspace array is not sized correctly
                call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
                    "The workspace array is too short.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call ZCH1DN(n, r, n, u, wptr, flag)
        if (flag == 1) then
            ! ERROR: The matrix is not positive definite
            call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
                "The downdated matrix is not positive definite.", &
                LA_MATRIX_FORMAT_ERROR)
        else if (flag == 2) then
            ! ERROR: The matrix is singular
            call errmgr%report_error("cholesky_rank1_downdate_cmplx", &
                "The input matrix is singular.", LA_SINGULAR_MATRIX_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ******************************************************************************
! RZ FACTORIZATION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine rz_factor_dbl(a, tau, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, lwork, flag, istat
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("rz_factor", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DTZRZF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("rz_factor", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine rz_factor_cmplx(a, tau, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: tau
        complex(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, lwork, flag, istat
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("rz_factor_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZTZRZF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("rz_factor_cmplx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("rz_factor_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZTZRZF
        call ZTZRZF(m, n, a, m, tau, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_rz_mtx(lside, trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        integer(int32), intent(in) :: l
        real(real64), intent(inout), dimension(:,:) :: a, c
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, lwork, flag, istat, lda
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_mtx", &
                    "Incorrectly sized input array WORK, argument 7.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_rz_mtx_cmplx(lside, trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        integer(int32), intent(in) :: l
        complex(real64), intent(inout), dimension(:,:) :: a, c
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, lwork, flag, istat, lda
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            t = 'C'
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_mtx_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMRZ(side, t, m, n, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_mtx_cmplx", &
                    "Incorrectly sized input array WORK, argument 7.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_mtx_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMRZ
        call ZUNMRZ(side, t, m, n, k, l, a, lda, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_rz_vec(trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        integer(int32), intent(in) :: l
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, k, lwork, flag, istat, lda
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
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

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_rz_vec_cmplx(trans, l, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        integer(int32), intent(in) :: l
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(inout), dimension(:) :: c
        complex(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, k, lwork, flag, istat, lda
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c)
        k = size(tau)
        lda = size(a, 1)
        side = 'L'
        if (trans) then
            t = 'C'
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
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_rz_vec_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_rz_vec_cmplx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_rz_vec_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMRZ
        call ZUNMRZ(side, t, m, 1, k, l, a, lda, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ******************************************************************************
! SVD ROUTINES
! ------------------------------------------------------------------------------
    module subroutine svd_dbl(a, s, u, vt, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: s
        real(real64), intent(out), optional, dimension(:,:) :: u, vt
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: jobu, jobvt
        integer(int32) :: m, n, mn, istat, lwork, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(u)) then
            if (size(u, 2) == m) then
                jobu = 'A'
            else if (size(u, 2) == mn) then
                jobu = 'S'
            end if
        else
            jobu = 'N'
        end if
        if (present(vt)) then
            jobvt = 'A'
        else
            jobvt = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(s) /= mn) then
            flag = 2
        else if (present(u)) then
            if (size(u, 1) /= m) flag = 3
            if (size(u, 2) /= m .and. size(u, 2) /= mn) flag = 3
        else if (present(vt)) then
            if (size(vt, 1) /= n .or. size(vt, 2) /= n) flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("svd", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, temp, -1, &
            flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGESVD
        if (present(u) .and. present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, vt, n, wptr, lwork, &
                flag)
        else if (present(u) .and. .not.present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, u, m, temp, n, wptr, &
                lwork, flag)
        else if (.not.present(u) .and. present(vt)) then
            call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, vt, n, wptr, &
                lwork, flag)
        else
            call DGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, wptr, &
                lwork, flag)
        end if

        ! Check for convergence
        if (flag > 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 101) flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("svd", errmsg, LA_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
101     format(I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine svd_cmplx(a, s, u, vt, work, olwork, rwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: s
        complex(real64), intent(out), optional, dimension(:,:) :: u, vt
        complex(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        real(real64), intent(out), target, optional, dimension(:) :: rwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: jobu, jobvt
        integer(int32) :: m, n, mn, istat, lwork, flag, lrwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        real(real64), dimension(1) :: rtemp
        real(real64), pointer, dimension(:) :: rwptr
        real(real64), allocatable, target, dimension(:) :: rwrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        lrwork = 5 * mn
        if (present(u)) then
            if (size(u, 2) == m) then
                jobu = 'A'
            else if (size(u, 2) == mn) then
                jobu = 'S'
            end if
        else
            jobu = 'N'
        end if
        if (present(vt)) then
            jobvt = 'A'
        else
            jobvt = 'N'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(s) /= mn) then
            flag = 2
        else if (present(u)) then
            if (size(u, 1) /= m) flag = 3
            if (size(u, 2) /= m .and. size(u, 2) /= mn) flag = 3
        else if (present(vt)) then
            if (size(vt, 1) /= n .or. size(vt, 2) /= n) flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("svd_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, temp, -1, &
            rtemp, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd_cmplx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("svd_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        if (present(rwork)) then
            if (size(rwork) < lrwork) then
                ! ERROR: RWORK not sized correctly
                call errmgr%report_error("svd_cmplx", &
                    "Incorrectly sized input array RWORK, argument 7.", &
                    LA_ARRAY_SIZE_ERROR)
            end if
            rwptr => rwork(1:lrwork)
        else
            allocate(rwrk(lrwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("svd_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            rwptr => rwrk
        end if

        ! Call ZGESVD
        if (present(u) .and. present(vt)) then
            call ZGESVD(jobu, jobvt, m, n, a, m, s, u, m, vt, n, wptr, lwork, &
                rwptr, flag)
        else if (present(u) .and. .not.present(vt)) then
            call ZGESVD(jobu, jobvt, m, n, a, m, s, u, m, temp, n, wptr, &
                lwork, rwptr, flag)
        else if (.not.present(u) .and. present(vt)) then
            call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, vt, n, wptr, &
                lwork, rwptr, flag)
        else
            call ZGESVD(jobu, jobvt, m, n, a, m, s, temp, m, temp, n, wptr, &
                lwork, rwptr, flag)
        end if

        ! Check for convergence
        if (flag > 0) then
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 101) flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("svd_cmplx", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
101     format(I0, A)
    end subroutine

! ******************************************************************************
! LQ FACTORIZATION
! ------------------------------------------------------------------------------
    module subroutine lq_factor_no_pivot(a, tau, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:) :: tau
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        real(real64), dimension(1) :: temp
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(tau) /= mn) then
            ! ERROR: TAU not sized correctly
            call errmgr%report_error("lq_factor_no_pivot", &
                "Incorrectly sized input array TAU, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELQF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("lq_factor_no_pivot", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("lq_factor_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DGELQF
        call DGELQF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine lq_factor_no_pivot_cmplx(a, tau, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: a
        complex(real64), intent(out), dimension(:) :: tau
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, istat, lwork, flag
        complex(real64), dimension(1) :: temp
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(tau) /= mn) then
            ! ERROR: TAU not sized correctly
            call errmgr%report_error("lq_factor_no_pivot_cmplx", &
                "Incorrectly sized input array TAU, argument 2.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZGELQF(m, n, a, m, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("lq_factor_no_pivot_cmplx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("lq_factor_no_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZGELQF
        call ZGELQF(m, n, a, m, tau, wptr, lwork, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lq_no_pivot(l, tau, q, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: l
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), dimension(:,:) :: q
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: i, j, m, n, mn, k, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(l, 1)
        n = size(l, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m > n) then
            flag = 1
        else if (size(tau) /= mn) then
            flag = 2
        else if (size(q, 1) /= n .or. size(q, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lq_no_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORGLQ(n, n, mn, q, n, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("form_lq_no_pivot", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("form_lq_no_pivot", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Copy the upper triangular portion of L to Q, and then zero it out in L
        do j = 2, n
            k = min(j - 1, m)
            q(1:k,j) = l(1:k,j)
            l(1:k,j) = zero
        end do

        ! Build Q
        call DORGLQ(n, n, mn, q, n, tau, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine form_lq_no_pivot_cmplx(l, tau, q, work, olwork, err)
        ! Arguments
        complex(real64), intent(inout), dimension(:,:) :: l
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(out), dimension(:,:) :: q
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, j, m, n, mn, k, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(l, 1)
        n = size(l, 2)
        mn = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m > n) then
            flag = 1
        else if (size(tau) /= mn) then
            flag = 2
        else if (size(q, 1) /= n .or. size(q, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_lq_no_pivot_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNGLQ(n, n, mn, q, n, tau, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("form_lq_no_pivot_cmplx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("form_lq_no_pivot_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Copy the upper triangular portion of L to Q, and then zero it out in L
        do j = 2, n
            k = min(j - 1, m)
            q(1:k,j) = l(1:k,j)
            l(1:k,j) = zero
        end do

        ! Build Q
        call ZUNGLQ(n, n, mn, q, n, tau, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_lq_mtx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:,:) :: c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, ncola, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        if (lside) then
            side = 'L'
            ncola = m
        else
            side = 'R'
            ncola = n
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
        if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_lq_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_lq_mtx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_lq_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMLQ
        call DORMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_lq_mtx_cmplx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(inout), dimension(:,:) :: c
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, ncola, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c, 1)
        n = size(c, 2)
        k = size(tau)
        if (lside) then
            side = 'L'
            ncola = m
        else
            side = 'R'
            ncola = n
        end if
        if (trans) then
            t = 'C'
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
        if (size(a, 1) /= k .or. size(a, 2) /= ncola) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_lq_mtx_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_lq_mtx_cmplx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_lq_mtx_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMLQ
        call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_lq_vec(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: c
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, istat, flag, lwork
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c)
        n = 1
        k = size(tau)
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
        if (size(a, 1) /= k .or. size(a, 2) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_lq_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_lq_vec", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_lq_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call DORMLQ
        call DORMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mult_lq_vec_cmplx(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(in), dimension(:) :: tau
        complex(real64), intent(inout), dimension(:) :: c
        complex(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: side, t
        integer(int32) :: m, n, k, istat, flag, lwork
        complex(real64), pointer, dimension(:) :: wptr
        complex(real64), allocatable, target, dimension(:) :: wrk
        complex(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg

        ! Initialization
        m = size(c)
        n = 1
        k = size(tau)
        side = 'L'
        if (trans) then
            t = 'C'
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
        if (size(a, 1) /= k .or. size(a, 2) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            allocate(character(len = 256) :: errmsg)
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_lq_vec_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mult_lq_vec_cmplx", &
                    "Incorrectly sized input array WORK, argument 6.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mult_lq_vec_cmplx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Call ZUNMLQ
        call ZUNMLQ(side, t, m, n, k, a, k, tau, c, m, wptr, lwork, flag)

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
end submodule
