module linalg_least_squares
    use iso_fortran_env, only : int32, real64
    use lapack
    use linalg_errors
    implicit none
    private
    public :: solve_least_squares
    public :: solve_least_squares_full
    public :: solve_least_squares_svd

    interface solve_least_squares
        module procedure :: solve_least_squares_mtx
        module procedure :: solve_least_squares_mtx_cmplx
        module procedure :: solve_least_squares_vec
        module procedure :: solve_least_squares_vec_cmplx
    end interface

    interface solve_least_squares_full
        module procedure :: solve_least_squares_mtx_pvt
        module procedure :: solve_least_squares_mtx_pvt_cmplx
        module procedure :: solve_least_squares_vec_pvt
        module procedure :: solve_least_squares_vec_pvt_cmplx
    end interface

    interface solve_least_squares_svd
        module procedure :: solve_least_squares_mtx_svd
        module procedure :: solve_least_squares_mtx_svd_cmplx
        module procedure :: solve_least_squares_vec_svd
        module procedure :: solve_least_squares_vec_svd_cmplx
    end interface
contains
! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx(a, b) result(x)
    !! Solves the system of equations \(A X = B\) assuming matrix \(A\) is of 
    !! full rank.
    real(real64), intent(in), dimension(:,:) :: a 
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, flag
    real(real64), allocatable, dimension(:) :: w
    real(real64), allocatable, dimension(:,:) :: ac, bc
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELS('N', m, n, nrhs, temp, m, temp, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    allocate(x(n, nrhs))
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call DGELS('N', m, n, nrhs, ac, m, x, n, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call DGELS('N', m, n, nrhs, ac, m, bc, m, w, lwork, flag)
        x = bc(1:n,:)
    end if
    if (flag > 0) then
        error stop LA_INVALID_OPERATION_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx_cmplx(a, b) result(x)
    !! Solves the system of equations \(A X = B\) assuming matrix \(A\) is of 
    !! full rank.
    complex(real64), intent(in), dimension(:,:) :: a 
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, flag
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), allocatable, dimension(:,:) :: ac, bc
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELS('N', m, n, nrhs, temp, m, temp, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    allocate(x(n, nrhs))
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call ZGELS('N', m, n, nrhs, ac, m, x, n, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call ZGELS('N', m, n, nrhs, ac, m, bc, m, w, lwork, flag)
        x = bc(1:n,:)
    end if
    if (flag > 0) then
        error stop LA_INVALID_OPERATION_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) assuming matrix 
    !! \(A\) is of full rank.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, lwork, flag
    real(real64), allocatable, dimension(:) :: w, bc
    real(real64), allocatable, dimension(:,:) :: ac
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELS('N', m, n, 1, temp, m, temp, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    allocate(x(n))
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call DGELS('N', m, n, 1, ac, m, x, n, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call DGELS('N', m, n, 1, ac, m, bc, m, w, lwork, flag)
        x = bc(1:n)
    end if
    if (flag > 0) then
        error stop LA_INVALID_OPERATION_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec_cmplx(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) assuming matrix 
    !! \(A\) is of full rank.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, lwork, flag
    complex(real64), allocatable, dimension(:) :: w, bc
    complex(real64), allocatable, dimension(:,:) :: ac
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELS('N', m, n, 1, temp, m, temp, maxmn, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Process
    allocate(x(n))
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call ZGELS('N', m, n, 1, ac, m, x, n, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call ZGELS('N', m, n, 1, ac, m, bc, m, w, lwork, flag)
        x = bc(1:n)
    end if
    if (flag > 0) then
        error stop LA_INVALID_OPERATION_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx_pvt(a, b) result(x)
    !! Solves the system of equations \(A X = B\) using a full orthogonal
    !! factorization of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, flag, rnk
    real(real64), allocatable, dimension(:) :: w
    integer(int32), allocatable, dimension(:) :: iw
    real(real64), allocatable, dimension(:,:) :: ac, bc
    real(real64), dimension(1) :: temp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELSY(m, n, nrhs, temp, m, temp, maxmn, itemp, rc, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), x(n, nrhs))
    allocate(iw(n), source = 0)

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call DGELSY(m, n, nrhs, ac, m, x, n, iw, rc, rnk, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call DGELSY(m, n, nrhs, ac, m, bc, m, iw, rc, rnk, w, lwork, flag)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx_pvt_cmplx(a, b) result(x)
    !! Solves the system of equations \(A X = B\) using a full orthogonal
    !! factorization of \(A\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, nrhs, lwork, flag, rnk, lrwork
    complex(real64), allocatable, dimension(:) :: w
    real(real64), allocatable, dimension(:) :: rw
    integer(int32), allocatable, dimension(:) :: iw
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    complex(real64), allocatable, dimension(:,:) :: ac, bc

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    lrwork = 2 * n
    rc = epsilon(rc)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELSY(m, n, nrhs, temp, m, temp, maxmn, itemp, rc, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(lwork), x(n, nrhs))
    allocate(iw(n), source = 0)

    ! Process
    allocate(ac(m, n), source = a)
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call ZGELSY(m, n, nrhs, ac, m, x, n, iw, rc, rnk, w, lwork, rw, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call ZGELSY(m, n, nrhs, ac, m, bc,m, iw, rc, rnk, w, lwork, rw, flag)
        x = bc(1:n,:)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec_pvt(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a full 
    !! orthogonal factorization of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, lwork, flag, rnk
    real(real64), allocatable, dimension(:) :: w, bc
    integer(int32), allocatable, dimension(:) :: iw
    real(real64), dimension(1) :: temp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    real(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    flag = 0
    if (size(b) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELSY(m, n, 1, temp, m, temp, maxmn, itemp, rc, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), x(n))
    allocate(iw(n), source = 0)

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call DGELSY(m, n, 1, ac, m, x, n, iw, rc, rnk, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call DGELSY(m, n, 1, ac, m, bc, m, iw, rc, rnk, w, lwork, flag)
        x = b(1:n)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec_pvt_cmplx(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a full 
    !! orthogonal factorization of \(A\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, maxmn, lwork, lrwork, flag, rnk
    complex(real64), allocatable, dimension(:) :: w, bc
    real(real64), allocatable, dimension(:) :: rw
    integer(int32), allocatable, dimension(:) :: iw
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    integer(int32), dimension(1) :: itemp
    real(real64) :: rc
    complex(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    maxmn = max(m, n)
    lrwork = 2 * n
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELSY(m, n, 1, temp, m, temp, maxmn, itemp, rc, rnk, temp, -1, rtemp, &
        flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), x(n))
    allocate(iw(n), source = 0)

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call ZGELSY(m, n, 1, ac, m, x, n, iw, rc, rnk, w, lwork, rw, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call ZGELSY(m, n, 1, ac, m, bc, m, iw, rc, rnk, w, lwork, rw, flag)
        x = b(1:n)
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx_svd(a, b) result(x)
    !! Solves the system of equations \(A X = B\) using a singular value
    !! decomposition of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, nrhs, mn, maxmn, flag, lwork, rnk
    real(real64), allocatable, dimension(:) :: w, s
    real(real64), dimension(1) :: temp
    real(real64) :: rc
    real(real64), allocatable, dimension(:,:) :: ac, bc

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    flag = 0
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELSS(m, n, nrhs, temp, m, temp, maxmn, temp, rc, rnk, temp, -1, &
        flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), s(mn), x(n, nrhs))

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call DGELSS(m, n, nrhs, ac, m, x, n, s, rc, rnk, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call DGELSS(m, n, nrhs, ac, m, bc, m, s, rc, rnk, w, lwork, flag)
    end if
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_mtx_svd_cmplx(a, b) result(x)
    !! Solves the system of equations \(A X = B\) using a singular value
    !! decomposition of \(A\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS matrix \(B\).
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS matrix \(X\).

    ! Local Variables
    integer(int32) :: m, n, nrhs, mn, maxmn, flag, lwork, rnk, lrwork
    complex(real64), allocatable, target, dimension(:) :: w
    real(real64), allocatable, target, dimension(:) :: rw, s
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64) :: rc
    complex(real64), allocatable, dimension(:,:) :: ac, bc

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    mn = min(m, n)
    lrwork = 5 * mn
    maxmn = max(m, n)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b, 1) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELSS(m, n, nrhs, temp, m, temp, maxmn, rtemp, rc, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), s(mn), x(n, nrhs))

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m,:) = b
        call ZGELSS(m, n, nrhs, ac, m, x, n, s, rc, rnk, w, lwork, rw, flag)
    else
        ! Create a copy of B
        allocate(bc(m, nrhs), source = b)
        call ZGELSS(m, n, nrhs, ac, m, bc, m, s, rc, rnk, w, lwork, rw, flag)
        x = bc(1:n,:)
    end if
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec_svd(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a singular 
    !! value decomposition of \(A\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    real(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    real(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, mn, maxmn, flag, lwork, rnk
    real(real64), allocatable, target, dimension(:) :: w, s, bc
    real(real64), dimension(1) :: temp
    real(real64) :: rc
    real(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call DGELSS(m, n, 1, temp, m, temp, maxmn, temp, rc, rnk, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), s(mn), x(n))

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call DGELSS(m, n, 1, ac, m, x, n, s, rc, rnk, w, lwork, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call DGELSS(m, n, 1, ac, m, bc, m, s, rc, rnk, w, lwork, flag)
        x = bc(1:m)
    end if
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end function

! ------------------------------------------------------------------------------
pure function solve_least_squares_vec_svd_cmplx(a, b) result(x)
    !! Solves the system of equations \(A \vec{x} = \vec{b}\) using a singular 
    !! value decomposition of \(A\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix \(A\).
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element vector \(\vec{b}\).
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element vector \(\vec{x}\).

    ! Local Variables
    integer(int32) :: m, n, mn, maxmn, flag, lwork, rnk, lrwork
    real(real64), allocatable, dimension(:) :: rw, s
    complex(real64), allocatable, dimension(:) :: w, bc
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    real(real64) :: rc
    complex(real64), allocatable, dimension(:,:) :: ac

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    lrwork = 5 * mn
    maxmn = max(m, n)
    rc = epsilon(rc)
    allocate(ac(m, n), source = a)

    ! Input Check
    if (size(b) /= m) then
        error stop 2
    end if

    ! Workspace Query
    call ZGELSS(m, n, 1, temp, m, temp, maxmn, rtemp, rc, rnk, temp, -1, &
        rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), s(mn), x(n))

    ! Process
    if (m <= n) then
        ! Use X to store B
        x(1:m) = b
        call ZGELSS(m, n, 1, ac, m, x, n, s, rc, rnk, w, lwork, rw, flag)
    else
        ! Create a copy of B
        allocate(bc(m), source = b)
        call ZGELSS(m, n, 1, ac, m, bc, m, s, rc, rnk, w, lwork, rw, flag)
        x = bc(1:n)
    end if
    if (flag > 0) then
        error stop LA_CONVERGENCE_ERROR
    end if
end function

! ------------------------------------------------------------------------------
end module