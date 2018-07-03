! linalg_solve.f90

!> @brief \b linalg_solve
!!
!! @par Purpose
!! Provides a set of routines for solving systems of linear equations.
submodule (linalg_core) linalg_solve
contains
! ******************************************************************************
! TRIANGULAR MATRIX SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine solve_tri_mtx(lside, upper, trans, nounit, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside, upper, trans, nounit
        real(real64), intent(in) :: alpha
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        character :: side, uplo, transa, diag

        ! Local Variables
        integer(int32) :: m, n, nrowa
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(b, 1)
        n = size(b, 2)
        if (lside) then
            nrowa = m
            side = 'L'
        else
            nrowa = n
            side = 'R'
        end if
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (trans) then
            transa = 'T'
        else
            transa = 'N'
        end if
        if (nounit) then
            diag = 'N'
        else
            diag = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check - matrix A must be square
        if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_mtx", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DTRSM
        call DTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, b, m)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_tri_mtx_cmplx(lside, upper, trans, nounit, alpha, a, b, err)
        ! Arguments
        logical, intent(in) :: lside, upper, trans, nounit
        complex(real64), intent(in) :: alpha
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        character :: side, uplo, transa, diag

        ! Local Variables
        integer(int32) :: m, n, nrowa
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(b, 1)
        n = size(b, 2)
        if (lside) then
            nrowa = m
            side = 'L'
        else
            nrowa = n
            side = 'R'
        end if
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (trans) then
            transa = 'C'
        else
            transa = 'N'
        end if
        if (nounit) then
            diag = 'N'
        else
            diag = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check - matrix A must be square
        if (size(a, 1) /= nrowa .or. size(a, 2) /= nrowa) then
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_mtx_cmplx", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZTRSM
        call ZTRSM(side, uplo, transa, diag, m, n, alpha, a, nrowa, b, m)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_tri_vec(upper, trans, nounit, a, x, err)
        ! Arguments
        logical, intent(in) :: upper, trans, nounit
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo, t, diag
        integer(int32) :: n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (trans) then
            t = 'T'
        else
            t = 'N'
        end if
        if (nounit) then
            diag = 'N'
        else
            diag = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_vec", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        else if (size(x) /= n) then
            ! ERROR: Inner matrix dimensions must agree
            call errmgr%report_error("solve_tri_vec", &
                "The inner matrix dimensions must be equal.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DTRSV
        call DTRSV(uplo, t, diag, n, a, n, x, 1)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_tri_vec_cmplx(upper, trans, nounit, a, x, err)
        ! Arguments
        logical, intent(in) :: upper, trans, nounit
        complex(real64), intent(in), dimension(:,:) :: a
        complex(real64), intent(inout), dimension(:) :: x
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        character :: uplo, t, diag
        integer(int32) :: n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (trans) then
            t = 'C'
        else
            t = 'N'
        end if
        if (nounit) then
            diag = 'N'
        else
            diag = 'U'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            ! ERROR: A must be square
            call errmgr%report_error("solve_tri_vec_cmplx", &
                "The input matrix must be square.", LA_ARRAY_SIZE_ERROR)
            return
        else if (size(x) /= n) then
            ! ERROR: Inner matrix dimensions must agree
            call errmgr%report_error("solve_tri_vec_cmplx", &
                "The inner matrix dimensions must be equal.", &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZTRSV
        call ZTRSV(uplo, t, diag, n, a, n, x, 1)
    end subroutine

! ******************************************************************************
! LU SOLUTION
! ------------------------------------------------------------------------------
    module subroutine solve_lu_mtx(a, ipvt, b, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_lu_mtx_cmplx(a, ipvt, b, err)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_mtx_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZGETRS
        call ZGETRS("N", n, nrhs, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_lu_vec(a, ipvt, b, err)
        ! Arguments
        real(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call DGETRS
        call DGETRS("N", n, 1, a, n, ipvt, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_lu_vec_cmplx(a, ipvt, b, err)
        ! Arguments
        complex(real64), intent(in), dimension(:,:) :: a
        integer(int32), intent(in), dimension(:) :: ipvt
        complex(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 1
        else if (size(ipvt) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_lu_vec_cmplx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Call ZGETRS
        call ZGETRS("N", n, 1, a, n, ipvt, b, n, flag)
    end subroutine

! ******************************************************************************
! QR SOLUTION
! ------------------------------------------------------------------------------
    module subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: m, n, nrhs, k, lwork, flag, istat
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b, 1) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., .true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., .true., a, tau, b, wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:n,1:n), b(1:n,:))
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, k, flag, lwork, istat
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        k = min(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (m < n) then
            flag = 1
        else if (size(tau) /= k) then
            flag = 2
        else if (size(b) /= m) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_no_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call mult_qr(.true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_no_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Compute Q**T * B, and store in B
        call mult_qr(.true., a, tau, b, work = wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
        call solve_triangular_system(.true., .false., .true., a(1:n,1:n), b)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:,:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(int32), parameter :: imin = 2
        integer(int32), parameter :: imax = 1
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
            rnk, maxmn, flag, istat, lwork1, lwork2, lwork3
        real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(real64), pointer, dimension(:) :: wptr, w, tau2
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
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
        else if (size(b, 1) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_qr(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
        call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b(1:n,:), &
            olwork = lwork3)
        lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(1:maxmn,:) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., .true., a, tau, b(1:m,:), w)

        ! Solve the triangular system T11 * B(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
        call solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), b(1:rnk,:))
        if (n > rnk) b(rnk+1:n,:) = zero

        ! Compute B(1:n,1:nrhs) = Y**T * B(1:n,1:nrhs)
        if (rnk < n) then
            call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b(1:n,:), w)
        end if

        ! Apply the pivoting: B(1:N,1:NRHS) = P * B(1:N,1:NRHS)
        do j = 1, nrhs
            do i = 1, n
                wptr(jpvt(i)) = b(i,j)
            end do
            b(:,j) = wptr(1:n)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(in), dimension(:) :: tau
        integer(int32), intent(in), dimension(:) :: jpvt
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(int32), parameter :: imin = 2
        integer(int32), parameter :: imax = 1
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, flag, &
            istat, lwork1, lwork2
        real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(real64), pointer, dimension(:) :: wptr, w, tau2
        real(real64), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        ismin = mn + 1
        ismax = 2 * mn + 1
        rcond = epsilon(rcond)
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
        else if (size(b) /= maxmn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_qr_pivot_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call rz_factor(a(1:mn,:), a(1:mn,1), olwork = lwork1)
        call mult_rz(.true., n, a(1:mn,:), a(1:mn,1), b(1:n), olwork = lwork2)
        lwork = max(lwork1, lwork2, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_qr_no_pivot_mtx", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_qr_pivot_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Determine the rank of R11 using an incremental condition estimation
        wptr(ismin) = one
        wptr(ismax) = one
        smax = abs(a(1,1))
        smin = smax
        if (abs(a(1,1)) == zero) then
            rnk = 0
            b(maxmn) = zero
            return
        else
            rnk = 1
        end if
        do
            if (rnk < mn) then
                i = rnk + 1
                call DLAIC1(imin, rnk, wptr(ismin:ismin+rnk-1), smin, &
                    a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
                call DLAIC1(imax, rnk, wptr(ismax:ismax+rnk-1), smax, &
                    a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
                if (smaxpr * rcond <= sminpr) then
                    do i = 1, rnk
                        wptr(ismin+i-1) = s1 * wptr(ismin+i-1)
                        wptr(ismax+i-1) = s2 * wptr(ismax+i-1)
                    end do
                    wptr(ismin+rnk) = c1
                    wptr(ismax+rnk) = c2
                    smin = sminpr
                    smax = smaxpr
                    rnk = rnk + 1
                    cycle
                end if
            end if
            exit
        end do

        ! Partition R = [R11 R12]
        !               [ 0  R22]
        tau2 => wptr(1:rnk)
        w => wptr(rnk+1:lwork)
        if (rnk < n) call rz_factor(a(1:rnk,:), tau2, w)

        ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
        call mult_qr(.true., a, tau, b(1:m))

        ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
        call solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            b(1:rnk))
        if (n > rnk) b(rnk+1:n) = zero

        ! Compute B(1:n) = Y**T * B(1:n)
        if (rnk < n) then
            call mult_rz(.true., n - rnk, a(1:rnk,:), tau2, b(1:n), w)
        end if

        ! Apply the pivoting: B(1:N) = P * B(1:N)
        do i = 1, n
            wptr(jpvt(i)) = b(i)
        end do
        b = wptr(1:n)
    end subroutine

! ******************************************************************************
! CHOLESKY SOLVE
! ------------------------------------------------------------------------------
    module subroutine solve_cholesky_mtx(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(int32) :: n, nrhs, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        nrhs = size(b, 2)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b, 1) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, nrhs, a, n, b, n, flag)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_cholesky_vec(upper, a, b, err)
        ! Arguments
        logical, intent(in) :: upper
        real(real64), intent(in), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        character :: uplo
        integer(int32) :: n, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        n = size(a, 1)
        if (upper) then
            uplo = 'U'
        else
            uplo = 'L'
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(a, 2) /= n) then
            flag = 2
        else if (size(b) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_cholesky_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Process
        call DPOTRS(uplo, n, 1, a, n, b, n, flag)
    end subroutine

! ******************************************************************************
! MATRIX INVERSION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine mtx_inverse_dbl(a, iwork, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        integer(int32), intent(out), target, optional, dimension(:) :: iwork
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: n, liwork, lwork, istat, flag
        integer(int32), pointer, dimension(:) :: iptr
        integer(int32), allocatable, target, dimension(:) :: iwrk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(a, 1)
        liwork = n
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(a, 2) /= n) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix must be squre to invert.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGETRI(n, a, n, istat, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Workspace Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Integer Workspace Allocation
        if (present(iwork)) then
            if (size(iwork) < liwork) then
                ! ERROR: IWORK not sized correctly
                call errmgr%report_error("svd", &
                    "Incorrectly sized input array IWORK, argument 2.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => iwork(1:liwork)
        else
            allocate(iwrk(liwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_inverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
        end if

        ! Compute the LU factorization of A
        call DGETRF(n, n, a, n, iptr, flag)

        ! Compute the inverse of the LU factored matrix
        call DGETRI(n, a, n, iptr, wptr, lwork, flag)

        ! Check for a singular matrix
        if (flag > 0) then
            call errmgr%report_error("mtx_inverse", &
                "The matrix is singular; therefore, the inverse could " // &
                "not be computed.", LA_SINGULAR_MATRIX_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine mtx_pinverse_dbl(a, ainv, tol, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(out), dimension(:,:) :: ainv
        real(real64), intent(in), optional :: tol
        real(real64), intent(out), target, dimension(:), optional :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! External Function Interfaces
        interface
            function DLAMCH(cmach) result(x)
                use, intrinsic :: iso_fortran_env, only : real64
                character, intent(in) :: cmach
                real(real64) :: x
            end function
        end interface

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, m, n, mn, lwork, istat, flag, i1, i2a, i2b, i3a, &
            i3b, i4
        real(real64), pointer, dimension(:) :: s, wptr, w
        real(real64), pointer, dimension(:,:) :: u, vt
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        real(real64) :: t, tref, tolcheck
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        i1 = m * mn
        i2a = i1 + 1
        i2b = i2a + n * n - 1
        i3a = i2b + 1
        i3b = i3a + mn - 1
        i4 = i3b + 1
        tolcheck = dlamch('s')
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(ainv, 1) /= n .or. size(ainv, 2) /= m) then
            write(errmsg, '(AI0AI0A)') &
                "The output matrix AINV is not sized appropriately.  " // &
                "It is expected to be ", n, "-by-", m, "."
            call errmgr%report_error("mtx_pinverse", errmsg, &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGESVD('S', 'A', m, n, a, m, a(1:mn,:), a, m, a, n, temp, -1, flag)
        lwork = int(temp(1), int32)
        lwork = lwork + m * mn + n * n + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("mtx_pinverse", &
                    "Incorrectly sized input array WORK, argument 4.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("mtx_pinverse", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if
        u(1:m,1:mn) => wptr(1:i1)
        vt(1:n,1:n) => wptr(i2a:i2b)
        s => wptr(i3a:i3b)
        w => wptr(i4:lwork)

        ! Compute the SVD of A
        call DGESVD('S', 'A', m, n, a, m, s, u, m, vt, n, w, size(w), flag)

        ! Check for convergence
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("mtx_pinverse", errmsg, &
                LA_CONVERGENCE_ERROR)
            return
        end if

        ! Determine the threshold tolerance for the singular values such that
        ! singular values less than the threshold result in zero when inverted.
        tref = max(m, n) * epsilon(t) * s(1)
        if (present(tol)) then
            t = tol
        else
            t = tref
        end if
        !if (t < safe_denom(t)) then
        if (t < tolcheck) then
            ! The supplied tolerance is too small, simply fall back to the
            ! default, but issue a warning to the user
            t = tref
            ! call errmgr%report_warning("pinverse_1", "The supplied tolerance was " // &
            !     "smaller than a value that would result in an overflow " // &
            !     "condition, or is negative; therefore, the tolerance has " // &
            !     "been reset to its default value.")
        end if

        ! Compute the pseudoinverse such that pinv(A) = V * inv(S) * U**T by
        ! first computing V * inv(S) (result is N-by-M), and store in the first
        ! MN rows of VT in a transposed manner.
        do i = 1, mn
            ! Apply 1 / S(I) to VT(I,:)
            if (s(i) < t) then
                vt(i,:) = zero
            else
                call recip_mult_array(s(i), vt(i,1:n))
            end if
        end do

        ! Compute (VT**T * inv(S)) * U**T
        call mtx_mult(.true., .true., one, vt(1:mn,:), u, zero, ainv)
    end subroutine

! ******************************************************************************
! LEAST SQUARES SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_mtx(a, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b, 1) /= maxmn) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, nrhs, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_vec(a, b, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, lwork, istat, flag
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        real(real64), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(b) /= maxmn) then
            call errmgr%report_error("solve_least_squares_vec", &
                "Input 2 is not sized correctly.", LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELS('N', m, n, 1, a, m, b, maxmn, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec", &
                    "Incorrectly sized input array WORK, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELS('N', m, n, 1, a, m, b, maxmn, wptr, lwork, flag)
        if (flag > 0) then
            call errmgr%report_error("solve_least_squares_mtx", &
                "The supplied matrix is not of full rank; therefore, " // &
                "the solution could not be computed via this routine.  " // &
                "Try a routine that utilizes column pivoting.", &
                LA_INVALID_OPERATION_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_mtx_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, nrhs, lwork, istat, flag, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        integer(int32), allocatable, target, dimension(:) :: iwrk
        integer(int32), pointer, dimension(:) :: iptr
        real(real64), dimension(1) :: temp
        integer(int32), dimension(1) :: itemp
        real(real64) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        nrhs = size(b, 2)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_mtx_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, nrhs, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, nrhs, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_vec_pvt(a, b, ipvt, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        integer(int32), intent(inout), target, optional, dimension(:) :: ipvt
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, maxmn, lwork, istat, flag, rnk
        real(real64), pointer, dimension(:) :: wptr
        real(real64), allocatable, target, dimension(:) :: wrk
        integer(int32), allocatable, target, dimension(:) :: iwrk
        integer(int32), pointer, dimension(:) :: iptr
        real(real64), dimension(1) :: temp
        integer(int32), dimension(1) :: itemp
        real(real64) :: rc
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        maxmn = max(m, n)
        rc = epsilon(rc)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        end if
        if (flag /= 0) then
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_vec_pvt", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSY(m, n, 1, a, m, b, maxmn, itemp, rc, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(ipvt)) then
            if (size(ipvt) < n) then
                ! ERROR: IPVT is not big enough
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Incorrectly sized pivot array, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            iptr => ipvt(1:n)
        else
            allocate(iwrk(n), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            iptr => iwrk
            iptr = 0
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec_pvt", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec_pvt", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSY(m, n, 1, a, m, b, maxmn, iptr, rc, rnk, wptr, lwork, flag)
        if (present(arnk)) arnk = rnk
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_mtx_svd(a, b, s, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a, b
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work, s
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, nrhs, mn, maxmn, istat, flag, lwork, rnk
        real(real64), pointer, dimension(:) :: wptr, sptr
        real(real64), allocatable, target, dimension(:) :: wrk, sing
        real(real64), dimension(1) :: temp
        real(real64) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        nrhs = size(b, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b, 1) /= maxmn) then
            flag = 2
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_mtx_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, nrhs, a, m, b, maxmn, temp, rcond, rnk, temp, -1, &
            flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(s)) then
            if (size(s) < mn) then
                ! ERROR: S not sized correctly
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Incorrectly sized input array S, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            sptr => s(1:mn)
        else
            allocate(sing(mn), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            sptr => sing
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_mtx_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, nrhs, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("solve_least_squares_mtx_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine solve_least_squares_vec_svd(a, b, s, arnk, work, olwork, err)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:) :: b
        integer(int32), intent(out), optional :: arnk
        real(real64), intent(out), target, optional, dimension(:) :: work, s
        integer(int32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: m, n, mn, maxmn, istat, flag, lwork, rnk
        real(real64), pointer, dimension(:) :: wptr, sptr
        real(real64), allocatable, target, dimension(:) :: wrk, sing
        real(real64), dimension(1) :: temp
        real(real64) :: rcond
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        maxmn = max(m, n)
        rcond = epsilon(rcond)
        if (present(arnk)) arnk = 0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        flag = 0
        if (size(b) /= maxmn) then
            flag = 2
        end if
        if (flag /= 0) then
            ! ERROR: One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("solve_least_squares_vec_svd", &
                trim(errmsg), LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGELSS(m, n, 1, a, m, b, maxmn, temp, rcond, rnk, temp, -1, flag)
        lwork = int(temp(1), int32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(s)) then
            if (size(s) < mn) then
                ! ERROR: S not sized correctly
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Incorrectly sized input array S, argument 3.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            sptr => s(1:mn)
        else
            allocate(sing(mn), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            sptr => sing
        end if

        if (present(work)) then
            if (size(work) < lwork) then
                ! ERROR: WORK not sized correctly
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Incorrectly sized input array WORK, argument 5.", &
                    LA_ARRAY_SIZE_ERROR)
                return
            end if
            wptr => work(1:lwork)
        else
            allocate(wrk(lwork), stat = istat)
            if (istat /= 0) then
                ! ERROR: Out of memory
                call errmgr%report_error("solve_least_squares_vec_svd", &
                    "Insufficient memory available.", &
                    LA_OUT_OF_MEMORY_ERROR)
                return
            end if
            wptr => wrk
        end if

        ! Process
        call DGELSS(m, n, 1, a, m, b, maxmn, sptr, rcond, rnk, wptr, lwork, &
            flag)
        if (present(arnk)) arnk = rnk
        if (flag > 0) then
            write(errmsg, '(I0A)') flag, " superdiagonals could not " // &
                "converge to zero as part of the QR iteration process."
            call errmgr%report_warning("solve_least_squares_vec_svd", errmsg, &
                LA_CONVERGENCE_ERROR)
        end if
    end subroutine


end submodule
