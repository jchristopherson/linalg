module linalg_svd
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    implicit none
    private
    public :: svd

    interface svd
        module procedure :: svd_dbl
        module procedure :: svd_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
pure subroutine svd_dbl(a, s, u, vt)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^T\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), optional, allocatable, target, dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors.
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! transpose of the N-by-N orthogonal matrix \(V\).  The matrix \(V\) 
        !! contains the right singular vectors.

    ! Local Variables
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), allocatable, target, dimension(:) :: w, sc
    real(real64), pointer, dimension(:) :: sptr
    real(real64), allocatable, dimension(:,:) :: ac
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(ac(m, n), source = a)
    if (present(u)) then
        jobu = 'A'
        allocate(u(m,m))
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
        allocate(vt(n,n))
    else
        jobvt = 'N'
    end if
    if (present(s)) then
        allocate(s(mn))
        sptr => s
    else
        allocate(sc(mn))
        sptr => sc
    end if

    ! Workspace Query
    call DGESVD(jobu, jobvt, m, n, ac, m, temp, temp, m, temp, n, temp, -1, &
        flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGESVD
    if (present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, vt, n, w, lwork, &
            flag)
    else if (present(u) .and. .not.present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, temp, n, w, &
            lwork, flag)
    else if (.not.present(u) .and. present(vt)) then
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, vt, n, w, &
            lwork, flag)
    else
        call DGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, temp, n, w, &
            lwork, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine svd_cmplx(a, s, u, vt)
    !! Computes the singular value decomposition of an M-by-N matrix \(A\) such 
    !! that \(A = U S V^H\) where \(U\) is an M-by-M orthogonal matrix, \(S\)
    !! is an M-by-N diagonal matrix containing the singular values, and \(V\)
    !! is an N-by-N orthogonal matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), optional, allocatable, target, dimension(:) :: s
        !! A MIN(M, N)-element array containing the singular values of a sorted 
        !! in descending order.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: u
        !! An optional argument, that if supplied, is used to contain the 
        !! orthogonal matrix \(U\) from the decomposition.  The matrix \(U\) 
        !! contains the left singular vectors.
    complex(real64), intent(out), optional, allocatable, dimension(:,:) :: vt
        !! An optional argument, that if supplied, is used to contain the 
        !! conjugate transpose of the N-by-N orthogonal matrix \(V\).  The 
        !! matrix \(V\) contains the right singular vectors.

    ! Local Variables
    character :: jobu, jobvt
    integer(int32) :: m, n, mn, lwork, flag, lrwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), allocatable, dimension(:,:) :: ac
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64), allocatable, target, dimension(:) :: rw, sc
    real(real64), pointer, dimension(:) :: sptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    allocate(ac(m, n), source = a)
    lrwork = 5 * mn
    if (present(u)) then
        jobu = 'A'
        allocate(u(m,m))
    else
        jobu = 'N'
    end if
    if (present(vt)) then
        jobvt = 'A'
        allocate(vt(n,n))
    else
        jobvt = 'N'
    end if
    if (present(s)) then
        allocate(s(mn))
        sptr => s
    else
        allocate(sc(mn))
        sptr => sc
    end if

    ! Workspace Query
    call ZGESVD(jobu, jobvt, m, n, ac, m, rtemp, temp, m, temp, n, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(lrwork))

    ! Call ZGESVD
    if (present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, vt, n, w, lwork, &
            rw, flag)
    else if (present(u) .and. .not.present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, u, m, temp, n, w, &
            lwork, rw, flag)
    else if (.not.present(u) .and. present(vt)) then
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, vt, n, w, &
            lwork, rw, flag)
    else
        call ZGESVD(jobu, jobvt, m, n, ac, m, sptr, temp, m, temp, n, w, &
            lwork, rw, flag)
    end if

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
end module