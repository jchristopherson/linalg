! linalg_eigen.f90

module linalg_eigen
    !! Provides eigenvalue and eigenvector computations for dense and sparse matrices.
    use iso_fortran_env, only : int32, real64
    use lapack
    use arpack
    use linalg_errors
    use linalg_sparse
    use linalg_lu, only : lu_factor
    use ieee_arithmetic, only : ieee_is_nan
    implicit none
    private
    public :: eigen

    interface eigen
        !! An interface to the eigenvalue and eigenvector routines.
        module procedure :: eigen_symm
        module procedure :: eigen_asymm
        module procedure :: eigen_gen
        module procedure :: eigen_cmplx
        module procedure :: eigen_sparse_symm
        module procedure :: eigen_sparse_asymm
        module procedure :: eigen_sparse_gen_symm
        module procedure :: eigen_sparse_gen_asymm
    end interface
contains
! ------------------------------------------------------------------------------
pure subroutine eigen_symm(a, vals, vecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, of a matrix
    !! by solving the eigenvalue problem \(A \vec{v} = \lambda \vec{v}\) when
    !! \(A\) is a symmetric matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N symmetric matrix on which to operate.
    real(real64), intent(out), allocatable, dimension(:) :: vals
        !! An N-element array that will contain the eigenvalues sorted into 
        !! ascending order.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: vecs
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
    allocate(vals(n))

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Workspace Query
    call DSYEV(jobz, 'L', n, temp, n, vals, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    if (present(vecs)) then
        allocate(vecs(n, n), source = a)
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
    complex(real64), intent(out), allocatable, dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: lvecs
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
    allocate(vals(n))

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Workspace Query
    call DGEEV(jobvl, jobvr, n, dummy, n, dummy, dummy, dummy, n, &
        dummy, n, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), wr(n), wi(n))

    ! Process
    allocate(ac(n, n), source = a)
    if (present(rvecs) .and. present(lvecs)) then
        ! Compute both the right and left eigenvectors
        allocate(vr(n, n), vl(n, n), rvecs(n, n), lvecs(n, n))
        call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, vl, n, vr, n, w, lwork, flag)
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(wr, wi, vr, rvecs, vals)
        call extract_eigenvectors(wr, wi, vl, lvecs)
    else if (present(rvecs) .and. .not.present(lvecs)) then
        ! Compute the right eigenvectors
        allocate(vr(n, n), rvecs(n, n))
        call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, dummy, n, vr, n, w, lwork, flag)
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(wr, wi, vr, rvecs, vals)
    else if (.not.present(rvecs) .and. present(lvecs)) then
        ! Compute the left eigenvectors
        allocate(vl(n, n), lvecs(n, n))
        call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, vl, n, dummy, n, w, lwork, flag)
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        call extract_eigenvectors(wr, wi, vl, lvecs, vals)
    else
        ! Only compute the eigenvalues
        call DGEEV(jobvl, jobvr, n, ac, n, wr, wi, dummy, n, dummy, n, &
            w, lwork, flag)
        if (flag > 0) error stop LA_CONVERGENCE_ERROR
        vals = cmplx(wr, wi, real64)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine eigen_gen(a, b, alpha, beta, rvecs, lvecs)
    !! Computes the eigenvalues, and optionally the eigenvectors, by solving
    !! the eigenvalue problem: \(A X = \lambda B X\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-N matrix \(B\).
    complex(real64), intent(out), allocatable, dimension(:) :: alpha
        !! An N-element array that, if beta is not supplied, contains the 
        !! eigenvalues.  If beta is supplied however, the eigenvalues must be 
        !! computed as \(\lambda = \alpha / \beta\).  This however, is not as
        !! trivial as it seems as it is entirely possible, and likely, that
        !! \(\alpha / \beta\) can overflow or underflow.  With that said, the 
        !! values in \(\alpha\) will always be less than and usually comparable 
        !! with the NORM(\(A\)).
    real(real64), intent(out), optional, allocatable, dimension(:) :: beta
        !! An optional N-element array that if provided forces alpha to return 
        !! the numerator, and this array contains the denominator used to 
        !! determine the eigenvalues as \(\lambda = \alpha / \beta\).  If used,
        !! the values in this array will always be less than and usually 
        !! comparable with the NORM(\(B\)).
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: lvecs
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
    allocate(alpha(n))

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop 2
    end if

    ! Workspace Query
    call DGGEV(jobvl, jobvr, n, dummy, n, dummy, n, temp, temp, temp, dummy, n, &
        dummy, n, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), ar(n), ai(n))

    ! Process
    allocate(ac(n, n), source = a)
    allocate(bc(n, n), source = b)
    if (present(beta)) then
        allocate(beta(n))
    else
        allocate(bt(n))
    end if
    if (present(rvecs) .and. present(lvecs)) then
        ! Compute both the right and left eigenvectors
        allocate(vl(n, n), vr(n, n), rvecs(n, n), lvecs(n, n))
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
        allocate(vr(n, n), rvecs(n, n))
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
        allocate(vl(n, n), lvecs(n, n))
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
    complex(real64), intent(out), allocatable, dimension(:) :: vals
        !! An N-element array containing the eigenvalues of the matrix.  The 
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: rvecs
        !! An optional N-by-N matrix, that if supplied, signals to compute the 
        !! right eigenvectors (one per column).
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: lvecs
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
    allocate(vals(n))

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
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
        allocate(rvecs(n, n), lvecs(n, n))
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, lvecs, n, rvecs, n, w, lwork, &
            rw, flag)
    else if (present(rvecs) .and. .not.present(lvecs)) then
        ! Compute the right eigenvectors
        allocate(rvecs(n, n))
        call ZGEEV(jobvl, jobvr, n, ac, n, vals, dummy, n, rvecs, n, w, lwork, &
            rw, flag)
    else if (.not.present(rvecs) .and. present(lvecs)) then
        ! Compute the left eigenvectors
        allocate(lvecs(n, n))
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
subroutine eigen_sparse_symm(a, k, vals, vecs, which, ncv, maxiter, tol)
    !! Computes the k largest, or smallest, eigenvalues, and optionally the
    !! corresponding eigenvectors, of a sparse matrix stored in compressed
    !! sparse row (CSR) format by solving the eigenvalue problem 
    !! \(A \vec{v} = \lambda \vec{v}\) when \(A\) is a symmetric matrix.  The
    !! calculations are performed by way of the implicitly restarted Lanczos
    !! iteration provided by ARPACK.
    class(csr_matrix), intent(in) :: a
        !! The N-by-N symmetric CSR matrix on which to operate.
    integer(int32), intent(in) :: k
        !! The number of eigenvalues to compute.  This value must be greater
        !! than zero, but less than N.
    real(real64), intent(out), allocatable, dimension(:) :: vals
        !! A K-element array containing the requested eigenvalues.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: vecs
        !! If present, the eigenvectors will be computed and this N-by-K matrix
        !! will contain the eigenvectors (one per column) corresponding to each
        !! eigenvalue in vals.
    character(len = 2), intent(in), optional :: which
        !! An optional parameter identifying which eigenvalues to compute.  The
        !! default is 'LM'.  The available options are as follows.
        !!
        !!  - 'LA': The K eigenvalues of largest algebraic value.
        !!
        !!  - 'SA': The K eigenvalues of smallest algebraic value.
        !!
        !!  - 'LM': The K eigenvalues of largest magnitude.
        !!
        !!  - 'SM': The K eigenvalues of smallest magnitude.
        !!
        !!  - 'BE': K eigenvalues taken from both ends of the spectrum.
    integer(int32), intent(in), optional :: ncv
        !! An optional parameter controlling the number of Lanczos basis
        !! vectors to utilize.  This value must be greater than k, but no
        !! greater than N.  The default is MIN(N, MAX(2 * k + 1, 20)).
    integer(int32), intent(in), optional :: maxiter
        !! An optional parameter controlling the maximum number of iterations
        !! allowed.  The default is 300.
    real(real64), intent(in), optional :: tol
        !! An optional parameter controlling the convergence tolerance.  If
        !! zero, or less than zero, machine precision is utilized.  The default
        !! is zero.

    ! Local Variables
    character(len = 1), parameter :: bmat = 'I'
    character(len = 2) :: w
    integer(int32) :: n, nc, ni, nev, lworkl, ido, info, nconv
    integer(int32) :: iparam(11), ipntr(11)
    real(real64) :: t
    real(real64), allocatable, dimension(:) :: resid, workd, workl, d
    real(real64), allocatable, dimension(:,:) :: v, z
    logical, allocatable, dimension(:) :: sel

    ! Initialization
    n = size(a, 1)
    w = 'LM'
    if (present(which)) w = which
    nc = min(n, max(2 * k + 1, 20))
    if (present(ncv)) nc = ncv
    ni = 300
    if (present(maxiter)) ni = maxiter
    t = 0.0d0
    if (present(tol)) t = tol

    ! Input Check
    if (size(a, 2) /= n) then
        error stop LA_ARRAY_SIZE_ERROR
    end if
    if (k < 1 .or. k >= n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (nc <= k .or. nc > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (ni < 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if

    ! Memory Allocation
    lworkl = nc * (nc + 8)
    allocate(resid(n), workd(3 * n), workl(lworkl), v(n, nc), d(k), &
        z(n, k), sel(nc))

    ! Set up the iteration parameters - regular mode
    iparam = 0
    iparam(1) = 1   ! use the exact shift strategy
    iparam(3) = ni  ! the iteration limit
    iparam(4) = 1   ! the block size - must be 1
    iparam(7) = 1   ! mode 1: A * x = lambda * x

    ! The reverse communication loop
    ido = 0
    info = 0
    nev = k
    do
        call dsaupd(ido, bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, &
            workd, workl, lworkl, info)
        if (ido /= -1 .and. ido /= 1) exit
        workd(ipntr(2):ipntr(2)+n-1) = &
            matmul(a, workd(ipntr(1):ipntr(1)+n-1))
    end do
    if (info == 1) error stop LA_CONVERGENCE_ERROR
    if (info < 0) error stop LA_INVALID_OPERATION_ERROR

    ! Extract the results
    call dseupd(present(vecs), 'A', sel, d, z, n, 0.0d0, bmat, n, w, nev, t, &
        resid, nc, v, n, iparam, ipntr, workd, workl, lworkl, info)
    if (info /= 0) error stop LA_CONVERGENCE_ERROR

    nconv = iparam(5)
    allocate(vals(nconv), source = d(1:nconv))
    if (present(vecs)) then
        allocate(vecs(n, nconv), source = z(:,1:nconv))
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_sparse_asymm(a, k, vals, vecs, which, ncv, maxiter, tol)
    !! Computes the k largest, or smallest, eigenvalues, and optionally the
    !! corresponding eigenvectors, of a sparse matrix stored in compressed
    !! sparse row (CSR) format by solving the eigenvalue problem 
    !! \(A \vec{v} = \lambda \vec{v}\) when \(A\) is square, but not 
    !! necessarily symmetric.  The calculations are performed by way of the
    !! implicitly restarted Arnoldi iteration provided by ARPACK.
    class(csr_matrix), intent(in) :: a
        !! The N-by-N CSR matrix on which to operate.
    integer(int32), intent(in) :: k
        !! The number of eigenvalues to compute.  This value must be greater
        !! than zero, but less than N - 1.
    complex(real64), intent(out), allocatable, dimension(:) :: vals
        !! A K-element array containing the requested eigenvalues.  The
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: vecs
        !! If present, the right eigenvectors will be computed and this N-by-K
        !! matrix will contain the eigenvectors (one per column) corresponding
        !! to each eigenvalue in vals.
    character(len = 2), intent(in), optional :: which
        !! An optional parameter identifying which eigenvalues to compute.  The
        !! default is 'LM'.  The available options are as follows.
        !!
        !!  - 'LM': The K eigenvalues of largest magnitude.
        !!
        !!  - 'SM': The K eigenvalues of smallest magnitude.
        !!
        !!  - 'LR': The K eigenvalues of largest real part.
        !!
        !!  - 'SR': The K eigenvalues of smallest real part.
        !!
        !!  - 'LI': The K eigenvalues of largest imaginary part.
        !!
        !!  - 'SI': The K eigenvalues of smallest imaginary part.
    integer(int32), intent(in), optional :: ncv
        !! An optional parameter controlling the number of Arnoldi basis
        !! vectors to utilize.  This value must be at least k + 2, but no
        !! greater than N.  The default is MIN(N, MAX(2 * k + 1, 20)).
    integer(int32), intent(in), optional :: maxiter
        !! An optional parameter controlling the maximum number of iterations
        !! allowed.  The default is 300.
    real(real64), intent(in), optional :: tol
        !! An optional parameter controlling the convergence tolerance.  If
        !! zero, or less than zero, machine precision is utilized.  The default
        !! is zero.

    ! Local Variables
    character(len = 1), parameter :: bmat = 'I'
    character(len = 2) :: w
    integer(int32) :: n, nc, ni, nev, lworkl, ido, info, nconv
    integer(int32) :: iparam(11), ipntr(14)
    real(real64) :: t
    real(real64), allocatable, dimension(:) :: resid, workd, workl, workev, &
        dr, di
    real(real64), allocatable, dimension(:,:) :: v, z
    logical, allocatable, dimension(:) :: sel

    ! Initialization
    n = size(a, 1)
    w = 'LM'
    if (present(which)) w = which
    nc = min(n, max(2 * k + 1, 20))
    if (present(ncv)) nc = ncv
    ni = 300
    if (present(maxiter)) ni = maxiter
    t = 0.0d0
    if (present(tol)) t = tol

    ! Input Check
    if (size(a, 2) /= n) then
        error stop LA_ARRAY_SIZE_ERROR
    end if
    if (k < 1 .or. k >= n - 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (nc < k + 2 .or. nc > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (ni < 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if

    ! Memory Allocation
    lworkl = 3 * nc**2 + 6 * nc
    allocate(resid(n), workd(3 * n), workl(lworkl), workev(3 * nc), &
        v(n, nc), dr(k + 1), di(k + 1), z(n, k + 1), sel(nc))

    ! Set up the iteration parameters - regular mode
    iparam = 0
    iparam(1) = 1   ! use the exact shift strategy
    iparam(3) = ni  ! the iteration limit
    iparam(4) = 1   ! the block size - must be 1
    iparam(7) = 1   ! mode 1: A * x = lambda * x

    ! The reverse communication loop
    ido = 0
    info = 0
    nev = k
    do
        call dnaupd(ido, bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, &
            workd, workl, lworkl, info)
        if (ido /= -1 .and. ido /= 1) exit
        workd(ipntr(2):ipntr(2)+n-1) = &
            matmul(a, workd(ipntr(1):ipntr(1)+n-1))
    end do
    if (info == 1) error stop LA_CONVERGENCE_ERROR
    if (info < 0) error stop LA_INVALID_OPERATION_ERROR

    ! Extract the results
    call dneupd(present(vecs), 'A', sel, dr, di, z, n, 0.0d0, 0.0d0, workev, &
        bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, workd, workl, &
        lworkl, info)
    if (info /= 0) error stop LA_CONVERGENCE_ERROR

    nconv = iparam(5)
    allocate(vals(nconv))
    if (present(vecs)) then
        allocate(vecs(n, nconv))
        call extract_eigenvectors(dr(1:nconv), di(1:nconv), z(:,1:nconv), &
            vecs, vals, .false.)
    else
        vals = cmplx(dr(1:nconv), di(1:nconv), real64)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_sparse_gen_symm(a, b, k, vals, vecs, sigma, which, ncv, &
    maxiter, tol)
    !! Computes k eigenvalues, and optionally the corresponding eigenvectors,
    !! of a pair of sparse matrices stored in compressed sparse row (CSR)
    !! format by solving the generalized eigenvalue problem
    !! \(A \vec{v} = \lambda B \vec{v}\) when \(A\) is symmetric and \(B\) is
    !! symmetric positive definite.  The calculations are performed by way of
    !! the implicitly restarted Lanczos iteration provided by ARPACK.
    class(csr_matrix), intent(in) :: a
        !! The N-by-N symmetric CSR matrix \(A\).
    class(csr_matrix), intent(in) :: b
        !! The N-by-N symmetric positive definite CSR matrix \(B\).
    integer(int32), intent(in) :: k
        !! The number of eigenvalues to compute.  This value must be greater
        !! than zero, but less than N.
    real(real64), intent(out), allocatable, dimension(:) :: vals
        !! A K-element array containing the requested eigenvalues.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: vecs
        !! If present, the eigenvectors will be computed and this N-by-K matrix
        !! will contain the eigenvectors (one per column) corresponding to each
        !! eigenvalue in vals.
    real(real64), intent(in), optional :: sigma
        !! An optional shift.  If supplied, the shift-invert transformation
        !! \((A - \sigma B)^{-1} B\) is utilized, which is the preferred means
        !! of finding the eigenvalues nearest to \(\sigma\).  If not supplied,
        !! the regular inverse transformation \(B^{-1} A\) is utilized.  Notice,
        !! \(\sigma\) must not be an eigenvalue of the pencil as the shifted
        !! matrix must be factorable.
    character(len = 2), intent(in), optional :: which
        !! An optional parameter identifying which eigenvalues to compute.  The
        !! default is 'LM'.  The available options are 'LA', 'SA', 'LM', 'SM',
        !! and 'BE'.  Note that when a shift is supplied, the selection applies
        !! to the transformed spectrum such that 'LM' returns the eigenvalues
        !! nearest to \(\sigma\).
    integer(int32), intent(in), optional :: ncv
        !! An optional parameter controlling the number of Lanczos basis
        !! vectors to utilize.  This value must be greater than k, but no
        !! greater than N.  The default is MIN(N, MAX(2 * k + 1, 20)).
    integer(int32), intent(in), optional :: maxiter
        !! An optional parameter controlling the maximum number of iterations
        !! allowed.  The default is 300.
    real(real64), intent(in), optional :: tol
        !! An optional parameter controlling the convergence tolerance.  If
        !! zero, or less than zero, machine precision is utilized.  The default
        !! is zero.

    ! Local Variables
    character(len = 1), parameter :: bmat = 'G'
    character(len = 2) :: w
    integer(int32) :: n, nc, ni, nev, lworkl, ido, info, nconv, mode, i1, i2, i3
    integer(int32) :: iparam(11), ipntr(11)
    integer(int32), allocatable, dimension(:) :: ju
    real(real64) :: t, s
    real(real64), allocatable, dimension(:) :: resid, workd, workl, d, y
    real(real64), allocatable, dimension(:,:) :: v, z
    logical, allocatable, dimension(:) :: sel
    type(csr_matrix) :: c
    type(msr_matrix) :: lu

    ! Initialization
    n = size(a, 1)
    w = 'LM'
    if (present(which)) w = which
    nc = min(n, max(2 * k + 1, 20))
    if (present(ncv)) nc = ncv
    ni = 300
    if (present(maxiter)) ni = maxiter
    t = 0.0d0
    if (present(tol)) t = tol
    s = 0.0d0
    mode = 2
    if (present(sigma)) then
        s = sigma
        mode = 3
    end if

    ! Input Check
    if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop LA_ARRAY_SIZE_ERROR
    end if
    if (k < 1 .or. k >= n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (nc <= k .or. nc > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (ni < 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if

    ! Factor the matrix against which the inner solutions are performed
    if (mode == 3) then
        c = a - s * b
    else
        c = b
    end if
    allocate(ju(n))
    call lu_factor(c, lu, ju)
    if (any(ieee_is_nan(lu%values))) error stop LA_SINGULAR_MATRIX_ERROR

    ! Memory Allocation
    lworkl = nc * (nc + 8)
    allocate(resid(n), workd(3 * n), workl(lworkl), v(n, nc), d(k), &
        z(n, k), sel(nc))

    ! Set up the iteration parameters
    iparam = 0
    iparam(1) = 1       ! use the exact shift strategy
    iparam(3) = ni      ! the iteration limit
    iparam(4) = 1       ! the block size - must be 1
    iparam(7) = mode    ! 2: inv(B) * A, or 3: inv(A - sigma * B) * B

    ! The reverse communication loop
    ido = 0
    info = 0
    nev = k
    do
        call dsaupd(ido, bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, &
            workd, workl, lworkl, info)
        i1 = ipntr(1)
        i2 = ipntr(2)
        i3 = ipntr(3)
        select case (ido)
        case (-1, 1)
            if (mode == 2) then
                ! ARPACK requires A * x be left in the input slot as well
                y = matmul(a, workd(i1:i1+n-1))
                workd(i1:i1+n-1) = y
            else if (ido == -1) then
                y = matmul(b, workd(i1:i1+n-1))
            else
                ! B * x has already been computed
                y = workd(i3:i3+n-1)
            end if
            workd(i2:i2+n-1) = pgmres_solver(c, lu, ju, y)
        case (2)
            workd(i2:i2+n-1) = matmul(b, workd(i1:i1+n-1))
        case default
            exit
        end select
    end do
    if (info == 1) error stop LA_CONVERGENCE_ERROR
    if (info < 0) error stop LA_INVALID_OPERATION_ERROR

    ! Extract the results
    call dseupd(present(vecs), 'A', sel, d, z, n, s, bmat, n, w, nev, t, &
        resid, nc, v, n, iparam, ipntr, workd, workl, lworkl, info)
    if (info /= 0) error stop LA_CONVERGENCE_ERROR

    nconv = iparam(5)
    allocate(vals(nconv), source = d(1:nconv))
    if (present(vecs)) then
        allocate(vecs(n, nconv), source = z(:,1:nconv))
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine eigen_sparse_gen_asymm(a, b, k, vals, vecs, sigma, which, ncv, &
    maxiter, tol)
    !! Computes k eigenvalues, and optionally the corresponding eigenvectors,
    !! of a pair of sparse matrices stored in compressed sparse row (CSR)
    !! format by solving the generalized eigenvalue problem
    !! \(A \vec{v} = \lambda B \vec{v}\) when \(A\) is square, but not
    !! necessarily symmetric, and \(B\) is symmetric positive semi-definite.
    !! The calculations are performed by way of the implicitly restarted
    !! Arnoldi iteration provided by ARPACK.
    class(csr_matrix), intent(in) :: a
        !! The N-by-N CSR matrix \(A\).
    class(csr_matrix), intent(in) :: b
        !! The N-by-N symmetric positive semi-definite CSR matrix \(B\).
    integer(int32), intent(in) :: k
        !! The number of eigenvalues to compute.  This value must be greater
        !! than zero, but less than N - 1.
    complex(real64), intent(out), allocatable, dimension(:) :: vals
        !! A K-element array containing the requested eigenvalues.  The
        !! eigenvalues are not sorted.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: vecs
        !! If present, the right eigenvectors will be computed and this N-by-K
        !! matrix will contain the eigenvectors (one per column) corresponding
        !! to each eigenvalue in vals.
    real(real64), intent(in), optional :: sigma
        !! An optional real-valued shift.  If supplied, the shift-invert
        !! transformation \((A - \sigma B)^{-1} B\) is utilized, which is the
        !! preferred means of finding the eigenvalues nearest to \(\sigma\).
        !! If not supplied, the regular inverse transformation \(B^{-1} A\) is
        !! utilized.  Notice, \(\sigma\) must not be an eigenvalue of the pencil
        !! as the shifted matrix must be factorable.
    character(len = 2), intent(in), optional :: which
        !! An optional parameter identifying which eigenvalues to compute.  The
        !! default is 'LM'.  The available options are 'LM', 'SM', 'LR', 'SR',
        !! 'LI', and 'SI'.  Note that when a shift is supplied, the selection
        !! applies to the transformed spectrum such that 'LM' returns the
        !! eigenvalues nearest to \(\sigma\).
    integer(int32), intent(in), optional :: ncv
        !! An optional parameter controlling the number of Arnoldi basis
        !! vectors to utilize.  This value must be at least k + 2, but no
        !! greater than N.  The default is MIN(N, MAX(2 * k + 1, 20)).
    integer(int32), intent(in), optional :: maxiter
        !! An optional parameter controlling the maximum number of iterations
        !! allowed.  The default is 300.
    real(real64), intent(in), optional :: tol
        !! An optional parameter controlling the convergence tolerance.  If
        !! zero, or less than zero, machine precision is utilized.  The default
        !! is zero.

    ! Local Variables
    character(len = 1), parameter :: bmat = 'G'
    character(len = 2) :: w
    integer(int32) :: n, nc, ni, nev, lworkl, ido, info, nconv, mode, i1, i2, i3
    integer(int32) :: iparam(11), ipntr(14)
    integer(int32), allocatable, dimension(:) :: ju
    real(real64) :: t, s
    real(real64), allocatable, dimension(:) :: resid, workd, workl, workev, &
        dr, di, y
    real(real64), allocatable, dimension(:,:) :: v, z
    logical, allocatable, dimension(:) :: sel
    type(csr_matrix) :: c
    type(msr_matrix) :: lu

    ! Initialization
    n = size(a, 1)
    w = 'LM'
    if (present(which)) w = which
    nc = min(n, max(2 * k + 1, 20))
    if (present(ncv)) nc = ncv
    ni = 300
    if (present(maxiter)) ni = maxiter
    t = 0.0d0
    if (present(tol)) t = tol
    s = 0.0d0
    mode = 2
    if (present(sigma)) then
        s = sigma
        mode = 3
    end if

    ! Input Check
    if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
        error stop LA_ARRAY_SIZE_ERROR
    end if
    if (k < 1 .or. k >= n - 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (nc < k + 2 .or. nc > n) then
        error stop LA_INVALID_INPUT_ERROR
    end if
    if (ni < 1) then
        error stop LA_INVALID_INPUT_ERROR
    end if

    ! Factor the matrix against which the inner solutions are performed
    if (mode == 3) then
        c = a - s * b
    else
        c = b
    end if
    allocate(ju(n))
    call lu_factor(c, lu, ju)
    if (any(ieee_is_nan(lu%values))) error stop LA_SINGULAR_MATRIX_ERROR

    ! Memory Allocation
    lworkl = 3 * nc**2 + 6 * nc
    allocate(resid(n), workd(3 * n), workl(lworkl), workev(3 * nc), &
        v(n, nc), dr(k + 1), di(k + 1), z(n, k + 1), sel(nc))

    ! Set up the iteration parameters
    iparam = 0
    iparam(1) = 1       ! use the exact shift strategy
    iparam(3) = ni      ! the iteration limit
    iparam(4) = 1       ! the block size - must be 1
    iparam(7) = mode    ! 2: inv(B) * A, or 3: inv(A - sigma * B) * B

    ! The reverse communication loop
    ido = 0
    info = 0
    nev = k
    do
        call dnaupd(ido, bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, &
            workd, workl, lworkl, info)
        i1 = ipntr(1)
        i2 = ipntr(2)
        i3 = ipntr(3)
        select case (ido)
        case (-1, 1)
            if (mode == 2) then
                ! ARPACK requires A * x be left in the input slot as well
                y = matmul(a, workd(i1:i1+n-1))
                workd(i1:i1+n-1) = y
            else if (ido == -1) then
                y = matmul(b, workd(i1:i1+n-1))
            else
                ! B * x has already been computed
                y = workd(i3:i3+n-1)
            end if
            workd(i2:i2+n-1) = pgmres_solver(c, lu, ju, y)
        case (2)
            workd(i2:i2+n-1) = matmul(b, workd(i1:i1+n-1))
        case default
            exit
        end select
    end do
    if (info == 1) error stop LA_CONVERGENCE_ERROR
    if (info < 0) error stop LA_INVALID_OPERATION_ERROR

    ! Extract the results
    call dneupd(present(vecs), 'A', sel, dr, di, z, n, s, 0.0d0, workev, &
        bmat, n, w, nev, t, resid, nc, v, n, iparam, ipntr, workd, workl, &
        lworkl, info)
    if (info /= 0) error stop LA_CONVERGENCE_ERROR

    nconv = iparam(5)
    allocate(vals(nconv))
    if (present(vecs)) then
        allocate(vecs(n, nconv))
        call extract_eigenvectors(dr(1:nconv), di(1:nconv), z(:,1:nconv), &
            vecs, vals, .false.)
    else
        vals = cmplx(dr(1:nconv), di(1:nconv), real64)
    end if
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
        ! The j == n test traps a conjugate pair whose second half lies
        ! outside of the requested subset of eigenvalues.
        if (abs(wi(j)) < eps .or. j == n) then
            ! We've got a real-valued eigenvalue
            if (present(vals)) then
                vals(j) = cmplx(wr(j), wi(j), real64)
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
