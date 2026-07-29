! linalg_eigen.f90

module linalg_eigen
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
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
pure subroutine eigen_symm(a, vals, vecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is a symmetric matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N symmetric matrix on which to operate.
    real(real64), intent(out), dimension(:) :: vals
        !! An N-element array that will contain the eigenvalues sorted into 
        !! ascending order.
    real(real64), intent(out), optional, dimension(:,:) :: vecs
        !! If present, the eigenvectors will be computed and this matrix will 
        !! contain the eigenvectors (one per column) corresponding to each
        !! eigenvalue in vals.

    ! Local Variables
    character :: jobz
    integer(int32) :: n, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), allocatable, dimension(:,:) :: ac
    real(real64), dimension(1) :: temp

    ! Initialization
    n = size(a, 1)
    if (present(vecs)) then
        jobz = 'V'
    else
        jobz = 'N'
    end if

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(vals) /= n) then
        error stop 2
    end if
    if (present(vecs) .and. (size(vecs, 1) /= n .or. size(vecs, 2) /= n)) then
        error stop 3
    end if

    ! Workspace Query
    call DSYEV(jobz, 'L', n, temp, n, vals, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    if (present(vecs)) then
        vecs = a
        call DSYEV('V', 'L', n, vecs, n, vals, w, lwork, flag)
    else
        allocate(ac(n, n), source = a)
        call DSYEV('N', 'L', n, ac, n, vals, w, lwork, flag)
    end if
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine eigen_asymm(a, vals, rvecs, lvecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is square, but not necessarily symmetric.
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, the N-by-N matrix on which to operate.  On output, the 
        !! contents of this matrix are overwritten.
    complex(real64), intent(out), dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, dimension(:,:) :: lvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the
        !! left eigenvectors (one per column).
    
    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: n, flag, lwork
    real(real64), dimension(1) :: dummy, temp
    real(real64), allocatable, dimension(:) :: w, wr, wi
    real(real64), allocatable, dimension(:,:) :: ac, vr, vl

    ! Initialization
    if (present(rvecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    if (present(lvecs)) then
        jobvl = 'V'
    else
        jobvl = 'N'
    end if
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(vals) /= n) then
        error stop 2
    end if
    if (present(rvecs) .and. (size(rvecs, 1) /= n .or. size(rvecs, 2) /= n)) then
        error stop 3
    end if
    if (present(lvecs) .and. (size(lvecs, 1) /= n .or. size(lvecs, 2) /= n)) then
        error stop 4
    end if

    ! Workspace Query
    ! call DGEEV(jobvl, jobvr, n, dummy, n, dummy, dummy, dummy, n, &
    !     dummy, n, temp, -1, flag)
    ! lwork = int(temp(1), int32)
    ! allocate(w(lwork), wr(n), wi(n))

    ! Process
    allocate(ac(n, n), source = a)
    ! if (present(rvecs) .and. present(lvecs)) then
    !     ! Compute both the right and left eigenvectors
    !     allocate(vr(n, n), vl(n, n))
    !     call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, vl, n, vr, n, w, lwork, flag)
    !     if (flag > 0) error stop LA_CONVERGENCE_ERROR
    !     call extract_eigenvectors(wr, wi, vr, rvecs, vals)
    !     call extract_eigenvectors(wr, wi, vl, lvecs)
    ! else if (present(rvecs) .and. .not.present(lvecs)) then
    !     ! Compute the right eigenvectors
    !     allocate(vr(n, n))
    !     call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, dummy, n, vr, n, w, lwork, flag)
    !     if (flag > 0) error stop LA_CONVERGENCE_ERROR
    !     call extract_eigenvectors(wr, wi, vr, rvecs, vals)
    ! else if (.not.present(rvecs) .and. present(lvecs)) then
    !     ! Compute the left eigenvectors
    !     allocate(vl(n, n))
    !     call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, vl, n, dummy, n, w, lwork, flag)
    !     if (flag > 0) error stop LA_CONVERGENCE_ERROR
    !     call extract_eigenvectors(wr, wi, vl, lvecs, vals)
    ! else
    !     ! Only compute the eigenvalues
    !     call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, dummy, n, dummy, n, &
    !         w, lwork, flag)
    !     if (flag > 0) error stop LA_CONVERGENCE_ERROR
    !     vals = cmplx(wr, wi, real64)
    ! end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine eigen_gen(a, b, alpha, beta, rvecs, lvecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, by solving
    !! the eigenvalue problem: \(A X = \lambda B X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-N matrix \(B\).
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
    complex(real64), intent(out), optional, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, dimension(:,:) :: lvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the
        !! left eigenvectors (one per column).

    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: i, n, flag, lwork
    real(real64) :: eps, p, q
    real(real64), dimension(1) :: temp
    real(real64), dimension(1,1) :: dummy
    real(real64), allocatable, dimension(:) :: w, ar, ai, bt
    real(real64), allocatable, dimension(:,:) :: ac, bc, vr, vl

    ! Initialization
    if (present(lvecs)) then
        jobvl = 'V'
    else
        jobvl = 'N'
    end if
    if (present(rvecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    n = size(a, 1)
    eps = epsilon(eps)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop 2
    end if
    if (size(alpha) /= n) then
        error stop 3
    end if
    if (present(beta) .and. size(beta) /= n) then
        error stop 4
    end if
    if (present(rvecs) .and. (size(rvecs, 1) /= n .or. size(rvecs, 2) /= n)) then
        error stop 5
    end if
    if (present(lvecs) .and. (size(lvecs, 1) /= n .or. size(lvecs, 2) /= n)) then
        error stop 6
    end if

    ! Workspace Query
    call DGGEV(jobvl, jobvr, n, dummy, n, dummy, n, temp, temp, temp, dummy, n, &
        dummy, n, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), ar(n), ai(n))

    ! Process
    allocate(ac(n, n), source = a)
    allocate(bc(n, n), source = b)
    if (.not.present(beta)) allocate(bt(n))
    if (present(rvecs) .and. present(lvecs)) then
        ! Compute both the right and left eigenvectors
        allocate(vl(n, n), vr(n, n))
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, beta, vl, n, &
                vr, n, w, lwork, flag)
        else
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, bt, vl, n, &
                vr, n, w, lwork, flag)
        end if
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(ar, ai, vr, rvecs, alpha, .false.)
        call extract_eigenvectors(ar, ai, vl, lvecs)
    else if (present(rvecs) .and. .not.present(lvecs)) then
        ! Compute the right eigenvectors
        allocate(vr(n, n))
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, beta, dummy, n, &
                vr, n, w, lwork, flag)
        else
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, bt, dummy, n, &
                vr, n, w, lwork, flag)
        end if
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(ar, ai, vr, rvecs, alpha, .false.)
    else if (.not.present(rvecs) .and. present(lvecs)) then
        ! Compute the left eigenvectors
        allocate(vl(n, n))
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, beta, vl, n, &
                dummy, n, w, lwork, flag)
        else
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, bt, vl, n, &
                dummy, n, w, lwork, flag)
        end if
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(ar, ai, vl, lvecs, alpha, .false.)
    else
        ! Compute only the eigenvalues
        if (present(beta)) then
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, beta, dummy, n, &
                dummy, n, w, lwork, flag)
        else
            call DGGEV(jobvl, jobvr, n, ac, n, bc, n, ar, ai, bt, dummy, n, &
                dummy, n, w, lwork, flag)
        end if
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        alpha = cmplx(ar, ai, real64)
    end if

    ! Finish the eigenvalue calculation, if necessary
    if (.not.present(beta)) then
        ! Compute: alpha / bt
        do i = 1, n
            call DLADIV(real(alpha(i)), aimag(alpha(i)), bt(i), 0.0d0, p, q)
            alpha(i) = cmplx(p, q, real64)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine eigen_cmplx(a, vals, rvecs, lvecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is square, but not necessarily symmetric.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix on which to operate.
    complex(real64), intent(out), dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, dimension(:,:) :: lvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the
        !! left eigenvectors (one per column).

    ! Local Variables
    character :: jobvl, jobvr
    integer(int32) :: n, flag, lwork, lrwork
    real(real64) :: rdummy(1)
    complex(real64) :: temp(1), dummy(1)
    complex(real64), allocatable, dimension(:) :: w
    real(real64), allocatable, dimension(:) :: rw
    complex(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    if (present(lvecs)) then
        jobvl = 'V'
    else
        jobvl = 'N'
    end if
    if (present(rvecs)) then
        jobvr = 'V'
    else
        jobvr = 'N'
    end if
    n = size(a, 1)
    lrwork = 2 * n

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(vals) /= n) then
        error stop 2
    end if
    if (present(rvecs) .and. (size(rvecs, 1) /= n .or. size(rvecs, 2) /= n)) then
        error stop 3
    end if
    if (present(lvecs) .and. (size(lvecs, 1) /= n .or. size(lvecs, 2) /= n)) then
        error stop 4
    end if

    ! Workspace Query
    call ZGEEV(jobvl, jobvr, n, dummy, n, dummy, dummy, n, dummy, n, temp, &
        -1, rdummy, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(lrwork))

    ! Process
    allocate(ac(n, n), source = a)
    if (present(rvecs) .and. present(lvecs)) then
        ! Compute the right and left eigenvectors
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, lvecs, n, rvecs, n, w, lwork, &
            rw, flag)
    else if (present(rvecs) .and. .not.present(lvecs)) then
        ! Compute the right eigenvectors
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, dummy, n, rvecs, n, w, lwork, &
            rw, flag)
    else if (.not.present(rvecs) .and. present(lvecs)) then
        ! Compute the left eigenvectors
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, lvecs, n, dummy, n, w, lwork, &
            rw, flag)
    else
        ! Only compute the eigenvalues
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, dummy, n, dummy, n, w, lwork, &
            rw, flag)
    end if
    if (flag > 0) error stop LA_CONVERGENCE_ERROR
end subroutine

! ------------------------------------------------------------------------------
pure subroutine extract_eigenvectors(wr, wi, v, vecs, vals, conjgvals)
    !! Extracts the eigenvalues and eigenvectors from the compact form used by
    !! LAPACK into a complex-valued, full form.
    real(real64), intent(in), dimension(:) :: wr
        !! The real components of the eigenvalues.
    real(real64), intent(in), dimension(:) :: wi
        !! The imaginary components of the eigenvalues.
    real(real64), intent(in), dimension(:,:) :: v
        !! The eigenvectors, in compact form.
    complex(real64), intent(out), dimension(:,:) :: vecs
        !! The full form of the eigenvector matrix.
    complex(real64), intent(out), optional, dimension(:) :: vals
        !! The eigenvalues.
    logical, intent(in), optional :: conjgvals
        !! Default is true.  If true, use the conjugate to compute 
        !! conjugate-pair eigenvalues; else, false to use the direct inputs.

    ! Local Variables
    integer(int32) :: j, jp1, n
    real(real64) :: eps
    logical :: cv

    ! Initialization
    n = size(wr)
    eps = 2.0d0 * epsilon(eps)
    cv = .true.
    if (present(conjgvals)) cv = conjgvals

    ! Process
    j = 1
    do while (j <= n)
        if (abs(wi(j)) < eps) then
            ! We've got a real-valued eigenvalue
            if (present(vals)) then
                vals(j) = cmplx(wr(j), 0.0d0, real64)
            end if
            vecs(:,j) = cmplx(v(:,j), 0.0d0, real64)
        else
            ! We've got a complex conjugate pair of eigenvalues
            jp1 = j + 1
            if (present(vals)) then
                vals(j) = cmplx(wr(j), wi(j), real64)
                if (cv) then
                    vals(jp1) = conjg(vals(j))
                else
                    vals(jp1) = cmplx(wr(jp1), wi(jp1), real64)
                end if
            end if
            vecs(:,j) = cmplx(v(:,j), v(:,jp1), real64)
            vecs(:,jp1) = conjg(vecs(:,j))
            j = j + 2
            cycle
        end if
        j = j + 1
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
