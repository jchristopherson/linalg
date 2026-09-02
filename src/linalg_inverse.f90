module linalg_inverse
    use iso_fortran_env, only : int32, real64
    use lapack
    use blas
    use linalg_errors
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: mtx_inverse
    public :: mtx_pinverse

    interface mtx_inverse
        module procedure :: mtx_inverse_dbl
        module procedure :: mtx_inverse_cmplx
    end interface

    interface mtx_pinverse
        module procedure :: mtx_pinverse_dbl
        module procedure :: mtx_pinverse_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
pure function mtx_inverse_dbl(a) result(rst)
    !! Computes the inverse of a square matrix.  In the event of a singular
    !! matrix, the resulting matrix is filled with NaN's.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to invert.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-N inverted matrix.

    ! Local Variables
    integer(int32) :: n, lwork, flag
    integer(int32), allocatable, dimension(:) :: ipvt
    real(real64) :: nan
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    n = size(a, 1)
    allocate(ipvt(n))
    allocate(rst(n, n), source = a)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Workspace Query
    call DGETRI(n, rst, n, ipvt, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Compute the LU factorization of A
    call DGETRF(n, n, rst, n, ipvt, flag)

    ! Compute the inverse of the LU factored matrix
    call DGETRI(n, rst, n, ipvt, w, lwork, flag)

    ! Check for a singular matrix
    if (flag > 0) then
        ! Singular Matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        rst = nan
    end if
end function

! ------------------------------------------------------------------------------
pure function mtx_inverse_cmplx(a) result(rst)
    !! Computes the inverse of a square matrix.  In the event of a singular
    !! matrix, the resulting matrix is filled with NaN's.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to invert.
    complex(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: n, lwork, flag
    integer(int32), allocatable, dimension(:) :: ipvt
    real(real64) :: nan
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    n = size(a, 1)
    allocate(ipvt(n))
    allocate(rst(n, n), source = a)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if

    ! Workspace Query
    call ZGETRI(n, rst, n, ipvt, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Compute the LU factorization of A
    call ZGETRF(n, n, rst, n, ipvt, flag)

    ! Compute the inverse of the LU factored matrix
    call ZGETRI(n, rst, n, ipvt, w, lwork, flag)

    ! Check for a singular matrix
    if (flag > 0) then
        ! Singular Matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        rst = cmplx(nan, nan)
    end if
end function

! ------------------------------------------------------------------------------
pure function mtx_pinverse_dbl(a, tol) result(ainv)
    !! Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix using the
    !! singular value decomposition of the matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to invert.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default tolerance 
        !! on singular values such that singular values less than this
        !! tolerance are forced to have a reciprocal of zero, as opposed to 
        !! 1/S(I).  The default tolerance is: MAX(M, N) * EPS * MAX(S).
    real(real64), allocatable, dimension(:,:) :: ainv
        !! The N-by-M inverted matrix.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, flag
    real(real64), allocatable, dimension(:,:) :: u, vt, ac
    real(real64), allocatable, dimension(:) :: w, s
    real(real64), dimension(1) :: temp
    real(real64) :: t, tref, tolcheck, ss

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    tolcheck = dlamch('s')
    allocate(ac(m, n), source = a)

    ! Workspace Query
    call DGESVD('S', 'S', m, n, ac, m, temp, ac, m, ac, mn, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), u(m, mn), vt(mn, n), s(mn), ainv(n, m))

    ! Compute the SVD of A
    call DGESVD('S', 'S', m, n, ac, m, s, u, m, vt, mn, w, lwork, flag)

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if

    ! Determine the threshold tolerance for the singular values such that
    ! singular values less than the threshold result in zero when inverted.
    tref = max(m, n) * epsilon(t) * s(1)
    if (present(tol)) then
        t = tol
    else
        t = tref
    end if
    if (t < tolcheck) then
        ! The supplied tolerance is too small, simply fall back to the
        ! default
        t = tref
    end if

    ! Compute the pseudoinverse such that pinv(A) = V * inv(S) * U**T by
    ! first computing inv(S) * U**T
    do i = 1, mn
        if (s(i) < t) then
            ss = s(i)
        else
            ss = one / s(i)
        end if
        call DSCAL(m, ss, u(:,i), 1)
    end do
    call DGEMM("T", "T", n, m, mn, one, vt, mn, u, m, zero, ainv, n)
end function

! ------------------------------------------------------------------------------
pure function mtx_pinverse_cmplx(a, tol) result(ainv)
    !! Computes the Moore-Penrose pseudo-inverse of a M-by-N matrix using the
    !! singular value decomposition of the matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to invert.
    real(real64), intent(in), optional :: tol
        !! An optional input, that if supplied, overrides the default tolerance 
        !! on singular values such that singular values less than this
        !! tolerance are forced to have a reciprocal of zero, as opposed to 
        !! 1/S(I).  The default tolerance is: MAX(M, N) * EPS * MAX(S).
    complex(real64), allocatable, dimension(:,:) :: ainv
        !! The N-by-M inverted matrix.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, flag, lrwork
    real(real64), allocatable, dimension(:) :: rw, s
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), allocatable, dimension(:,:) :: u, vt, ac
    complex(real64) :: temp(1), val
    real(real64) :: ss, t, tref, tolcheck, rtemp(1)

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    tolcheck = dlamch('s')
    lrwork = 5 * mn
    allocate(ac(m, n), source = a)

    ! Workspace Query
    call ZGESVD('S', 'S', m, n, ac, m, rtemp, ac, m, ac, mn, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(lrwork), s(mn), u(m, mn), vt(mn, n), ainv(n, m))

    ! Compute the SVD of A
    call ZGESVD('S', 'S', m, n, ac, m, s, u, m, vt, mn, w, lwork, rw, flag)

    ! Check for convergence
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if

    ! Determine the threshold tolerance for the singular values such that
    ! singular values less than the threshold result in zero when inverted.
    tref = max(m, n) * epsilon(t) * s(1)
    if (present(tol)) then
        t = tol
    else
        t = tref
    end if
    if (t < tolcheck) then
        ! The supplied tolerance is too small, simply fall back to the
        ! default
        t = tref
    end if

    ! Compute the pseudoinverse such that pinv(A) = V * inv(S) * U**H by
    ! first computing inv(S) * U**H
    do i = 1, mn
        if (s(i) < t) then
            ss = s(i)
        else
            ss = 1.0d0 / s(i)
        end if
        call ZDSCAL(m, ss, u(:,i), 1)
    end do
    call ZGEMM("C", "C", n, m, mn, one, vt, mn, u, m, zero, ainv, n)
end function

! ------------------------------------------------------------------------------
end module