! linalg_eigen.f90

module linalg_eigen
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    use ferror
    implicit none
    private
    public :: eigen

    interface eigen
        !! An interface to the eigenvalue and eigenvector routines.
        module procedure :: eigen_symm
        module procedure :: eigen_asymm
        module procedure :: eigen_gen
        module procedure :: eigen_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
subroutine eigen_symm(vecs, a, vals, work, olwork, err)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is a symmetric matrix.
    logical, intent(in) :: vecs
        !! Set to true to compute the eigenvectors as well as the eigenvalues; 
        !! else, set to false to just compute the eigenvalues.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N symmetric matrix on which to operate.  On 
        !! output, and if vecs is set to true, the matrix will contain the 
        !! eigenvectors (one per column) corresponding to each eigenvalue in 
        !! vals.  If vecs is set to false, the lower triangular portion of the 
        !! matrix is overwritten.
    real(real64), intent(out), dimension(:) :: vals
        !! An N-element array that will contain the eigenvalues sorted into 
        !! ascending order.
    real(real64), intent(out), pointer, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: jobz
    integer(int32) :: n, istat, flag, lwork
    real(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64), dimension(1) :: temp
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    n = size(a, 1)
    if (vecs) then
        jobz = 'V'
    else
        jobz = 'N'
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("eigen_symm", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(vals) /= n) then
        call report_array_size_error("eigen_symm", errmgr, "vals", n, &
            size(vals))
        return
    end if

    ! Workspace Query
    call DSYEV(jobz, 'L', n, a, n, vals, temp, -1, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("eigen_symm", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("eigen_symm", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Process
    call DSYEV(jobz, 'L', n, a, n, vals, wptr, lwork, flag)
    if (flag > 0) then
        call errmgr%report_error("eigen_symm", &
            "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_asymm(a, vals, vecs, work, olwork, err)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is square, but not necessarily symmetric.
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix on which to operate.  On output, the 
        !! contents of this matrix are overwritten.
    complex(real64), intent(out), dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, dimension(:,:) :: vecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).  If not provided, only the 
        !! eigenvalues will be computed.
    real(real64), intent(out), pointer, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns without
        !! performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: i, j, jp1, n, n1, n2a, n2b, n3a, n3b, istat, flag, &
        lwork, lwork1
    real(real64) :: eps
    real(real64), dimension(1) :: dummy, temp
    real(real64), dimension(1,1) :: dummy_mtx
    real(real64), pointer, dimension(:) :: wr, wi, wptr, w
    real(real64), pointer, dimension(:,:) :: vr
    real(real64), allocatable, target, dimension(:) :: wrk
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    jobvl = 'N'
    if (present(vecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    n = size(a, 1)
    eps = two * epsilon(eps)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("eigen_asymm", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(vals) /= n) then
        call report_array_size_error("eigen_asymm", errmgr, "vals", n, &
            size(vals))
        return
    else if (present(vecs)) then
        if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
            call report_matrix_size_error("eigen_asymm", errmgr, "vecs", &
                n, n, size(vecs, 1), size(vecs, 2))
            return
        end if
    end if

    ! Workspace Query
    call DGEEV(jobvl, jobvr, n, a, n, dummy, dummy, dummy_mtx, n, &
        dummy_mtx, n, temp, -1, flag)
    lwork1 = int(temp(1), int32)
    if (present(vecs)) then
        lwork = lwork1 + 2 * n + n * n
    else
        lwork = lwork1 + 2 * n
    end if
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("eigen_asymm", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("eigen_asymm", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Locate each array within the workspace array
    n1 = n
    n2a = n1 + 1
    n2b = n2a + n - 1
    n3a = n2b + 1
    n3b = n3a + lwork1 - 1

    ! Assign pointers
    wr => wptr(1:n1)
    wi => wptr(n2a:n2b)
    w => wptr(n3a:n3b)

    ! Process
    if (present(vecs)) then
        ! Assign a pointer to the eigenvector matrix
        vr(1:n,1:n) => wptr(n3b+1:lwork)

        ! Compute the eigenvectors and eigenvalues
        call DGEEV(jobvl, jobvr, n, a, n, wr, wi, dummy_mtx, n, vr, n, &
            w, lwork1, flag)

        ! Check for convergence
        if (flag > 0) then
            call errmgr%report_error("eigen_asymm", &
                "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
            return
        end if

        ! Store the eigenvalues and eigenvectors
        j = 1
        do while (j <= n)
            if (abs(wi(j)) < eps) then
                ! We've got a real-valued eigenvalue
                vals(j) = cmplx(wr(j), zero, real64)
                do i = 1, n
                    vecs(i,j) = cmplx(vr(i,j), zero, real64)
                end do
            else
                ! We've got a complex cojugate pair of eigenvalues
                jp1 = j + 1
                vals(j) = cmplx(wr(j), wi(j), real64)
                vals(jp1) = conjg(vals(j))
                do i = 1, n
                    vecs(i,j) = cmplx(vr(i,j), vr(i,jp1), real64)
                    vecs(i,jp1) = conjg(vecs(i,j))
                end do

                ! Increment j and continue the loop
                j = j + 2
                cycle
            end if

            ! Increment j
            j = j + 1
        end do
    else
        ! Compute just the eigenvalues
        call DGEEV(jobvl, jobvr, n, a, n, wr, wi, dummy_mtx, n, &
            dummy_mtx, n, w, lwork1, flag)

        ! Check for convergence
        if (flag > 0) then
            call errmgr%report_error("eigen_asymm", &
                "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
            return
        end if

        ! Store the eigenvalues
        do i = 1, n
            vals(i) = cmplx(wr(i), wi(i), real64)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_gen(a, b, alpha, beta, vecs, work, olwork, err)
    !! Computes the eigenvalues, and optionally the eigenvectors, by solving
    !! the eigenvalue problem: \(A X = \lambda B X\).
    real(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix \(A\).  On output, the contents of this 
        !! matrix are overwritten.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the N-by-N matrix \(B\).  On output, the contents of this 
        !! matrix are overwritten.
    complex(real64), intent(out), dimension(:) :: alpha
        !! An N-element array that, if beta is not supplied, contains the 
        !! eigenvalues.  If beta is supplied however, the eigenvalues must be 
        !! computed as \(\lambda = \alpha / \beta\).  This however, is not as
        !! trivial as it seems as it is entirely possible, and likely, that
        !! \(\alpha / \beta\) can overflow or underflow.  With that said, the 
        !! values in \(\alpha\) will always be less than and usually comparable 
        !! with the NORM(\(A\)).
    real(real64), intent(out), optional, dimension(:) :: beta
        !! An optional N-element array that if provided forces alpha to return 
        !! the numerator, and this array contains the denominator used to 
        !! determine the eigenvalues as \(\lambda = \alpha / \beta\).  If used,
        !! the values in this array will always be less than and usually 
        !! comparable with the NORM(\(B\)).
    complex(real64), intent(out), optional, dimension(:,:) :: vecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).  If not provided, only the 
        !! eigenvalues will be computed.
    real(real64), intent(out), optional, pointer, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least 
        !! olwork.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns 
        !! without performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: two = 2.0d0

    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: i, j, jp1, n, n1, n2a, n2b, n3a, n3b, n4a, n4b, &
        istat, flag, lwork, lwork1
    real(real64), dimension(1) :: temp
    real(real64), dimension(1,1) :: dummy
    real(real64), pointer, dimension(:) :: wptr, w, alphar, alphai, bptr
    real(real64), pointer, dimension(:,:) :: vr
    real(real64), allocatable, target, dimension(:) :: wrk
    real(real64) :: eps
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    jobvl = 'N'
    jobvr = 'N'
    if (present(vecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    n = size(a, 1)
    eps = two * epsilon(eps)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("eigen_gen", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(b, 1) /= n .or. size(b, 2) /= n) then
        call report_matrix_size_error("eigen_gen", errmgr, "b", n, n, &
            size(b, 1), size(b, 2))
        return
    else if (size(alpha) /= n) then
        call report_array_size_error("eigen_gen", errmgr, "alpha", n, &
            size(alpha))
        return
    else if (present(beta)) then
        if (size(beta) /= n) then
            call report_array_size_error("eigen_gen", errmgr, "beta", n, &
                size(beta))
            return
        end if
    else if (present(vecs)) then
        if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
            call report_matrix_size_error("eigen_gen", errmgr, "vecs", n, n, &
                size(vecs, 1), size(vecs, 2))
            return
        end if
    end if

    ! Workspace Query
    call DGGEV(jobvl, jobvr, n, a, n, b, n, temp, temp, temp, dummy, n, &
        dummy, n, temp, -1, flag)
    lwork1 = int(temp(1), int32)
    lwork = lwork1 + 2 * n
    if (.not.present(beta)) then
        lwork = lwork + n
    end if
    if (present(vecs)) then
        lwork = lwork + n * n
    end if
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("eigen_gen", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work(1:lwork)
    else
        allocate(wrk(lwork), stat = istat)
        if (istat /= 0) then
            call report_memory_error("eigen_gen", errmgr, istat)
            return
        end if
        wptr => wrk
    end if

    ! Locate each array within the workspace array & assign pointers
    n1 = n
    n2a = n1 + 1
    n2b = n2a + n - 1
    n3a = n2b + 1
    n3b = n3a + lwork1 - 1
    n4b = n3b
    alphar => wptr(1:n1)
    alphai => wptr(n2a:n2b)
    w => wptr(n3a:n3b)
    if (.not.present(beta)) then
        n4a = n3b + 1
        n4b = n4a + n - 1
        bptr => wptr(n4a:n4b)
    end if

    ! Process
    if (present(vecs)) then
        ! Assign a pointer to the eigenvector matrix
        vr(1:n,1:n) => wptr(n4b+1:lwork)

        ! Compute the eigenvalues and eigenvectors
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                beta, dummy, n, vr, n, w, lwork1, flag)
        else
            call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                bptr, dummy, n, vr, n, w, lwork1, flag)
        end if

        ! Check for convergence
        if (flag > 0) then
            call errmgr%report_error("eigen_gen", &
                "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
            return
        end if

        ! Store the eigenvalues and eigenvectors
        j = 1
        do while (j <= n)
            if (abs(alphai(j)) < eps) then
                ! Real-Valued
                alpha(j) = cmplx(alphar(j), zero, real64)
                do i = 1, n
                    vecs(i,j) = cmplx(vr(i,j), zero, real64)
                end do
            else
                ! Complex-Valued
                jp1 = j + 1
                alpha(j) = cmplx(alphar(j), alphai(j), real64)
                alpha(jp1) = cmplx(alphar(jp1), alphai(jp1), real64)
                do i = 1, n
                    vecs(i,j) = cmplx(vr(i,j), vr(i,jp1), real64)
                    vecs(i,jp1) = conjg(vecs(i,j))
                end do

                ! Increment j and continue
                j = j + 2
                cycle
            end if

            ! Increment j
            j = j + 1
        end do
        if (.not.present(beta)) alpha = alpha / bptr
    else
        ! Compute just the eigenvalues
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                beta, dummy, n, dummy, n, w, lwork1, flag)
        else
            call DGGEV(jobvl, jobvr, n, a, n, b, n, alphar, alphai, &
                bptr, dummy, n, dummy, n, w, lwork1, flag)
        end if

        ! Check for convergence
        if (flag > 0) then
            call errmgr%report_error("eigen_gen", &
                "The algorithm failed to converge.", LA_CONVERGENCE_ERROR)
            return
        end if

        ! Store the eigenvalues
        do i = 1, n
            alpha(i) = cmplx(alphar(i), alphai(i), real64)
        end do
        if (.not.present(beta)) alpha = alpha / bptr
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_cmplx(a, vals, vecs, work, olwork, rwork, err)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is square, but not necessarily symmetric.
    complex(real64), intent(inout), dimension(:,:) :: a
        !! On input, the N-by-N matrix on which to operate.  On output, the 
        !! contents of this matrix are overwritten.
    complex(real64), intent(out), dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, dimension(:,:) :: vecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).  If not provided, only the 
        !! eigenvalues will be computed.
    complex(real64), intent(out), target, optional, dimension(:) :: work
        !! An optional input, that if provided, prevents any local memory 
        !! allocation.  If not provided, the memory required is allocated
        !! within.  If provided, the length of the array must be at least
        !! olwork.
    real(real64), intent(out), target, optional, dimension(:) :: rwork
        !! An optional input, that if provided, prevents any local memory 
        !! allocation for real-valued workspaces.  If not provided, the 
        !! memory required is allocated within.  If provided, the length of the 
        !! array must be at least 2 * N.
    integer(int32), intent(out), optional :: olwork
        !! An optional output used to determine workspace size.  If supplied, 
        !! the routine determines the optimal size for work, and returns without
        !! performing any actual calculations.
    class(errors), intent(inout), optional, target :: err
        !! An error object to report any errors that occur.

    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: n, flag, lwork, lrwork
    real(real64) :: rdummy(1)
    complex(real64) :: temp(1), dummy(1), dummy_mtx(1,1)
    complex(real64), allocatable, target, dimension(:) :: wrk
    complex(real64), pointer, dimension(:) :: wptr
    real(real64), allocatable, target, dimension(:) :: rwrk
    real(real64), pointer, dimension(:) :: rwptr
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    jobvl = 'N'
    if (present(vecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    n = size(a, 1)
    lrwork = 2 * n

    ! Input Check
    if (size(a, 2) /= n) then
        call report_square_matrix_error("eigen_cmplx", errmgr, "a", n, &
            size(a, 1), size(a, 2))
        return
    else if (size(vals) /= n) then
        call report_array_size_error("eigen_cmplx", errmgr, "vals", n, &
            size(vals))
        return
    else if (present(vecs)) then
        if (size(vecs, 1) /= n .or. size(vecs, 2) /= n) then
            call report_matrix_size_error("eigen_cmplx", errmgr, "vecs", n, n, &
                size(vecs, 1), size(vecs, 2))
            return
        end if
    end if

    ! Workspace Query
    call ZGEEV(jobvl, jobvr, n, a, n, dummy, dummy_mtx, n, dummy_mtx, n, temp, &
        -1, rdummy, flag)
    lwork = int(temp(1), int32)
    if (present(olwork)) then
        olwork = lwork
        return
    end if

    ! Local Memory Allocation
    if (present(work)) then
        if (size(work) < lwork) then
            call report_array_size_error("eigen_cmplx", errmgr, "work", lwork, &
                size(work))
            return
        end if
        wptr => work
    else
        allocate(wrk(lwork), stat = flag)
        if (flag /= 0) then
            call report_memory_error("eigen_cmplx", errmgr, flag)
            return
        end if
        wptr => wrk
    end if
    
    if (present(rwork)) then
        if (size(rwork) < lrwork) then
            call report_array_size_error("eigen_cmplx", errmgr, "rwork", &
                lrwork, size(rwork))
            return
        end if
        rwptr => rwork
    else
        allocate(rwrk(lrwork), stat = flag)
        if (flag /= 0) then
            call report_memory_error("eigen_cmplx", errmgr, flag)
            return
        end if
        rwptr => rwrk
    end if

    ! Process
    if (present(vecs)) then
        call ZGEEV(jobvl, jobvr, n, a, n, vals, dummy_mtx, n, vecs, n, &
            wptr, lwork, rwptr, flag)
    else
        call ZGEEV(jobvl, jobvr, n, a, n, vals, dummy_mtx, n, dummy_mtx, n, &
            wptr, lwork, rwptr, flag)
    end if

    if (flag > 0) then
        call errmgr%report_error("eigen_cmplx", &
            "The algorithm failed to converge.", &
            LA_CONVERGENCE_ERROR)
        return
    end if
end subroutine

! ------------------------------------------------------------------------------
end module
