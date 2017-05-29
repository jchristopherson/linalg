! qr.f90

!> @brief \b qr
!!
!! @par Purpose
!! Provides a set of routines for solving systems of equations using QR
!! factorization.
module qr
    use ferror, only : errors
    use linalg_constants
    use linalg_core, only : solve_triangular_system
    use rz, only : rz_factor, mult_rz
    implicit none
    private
    public :: qr_factor
    public :: form_qr
    public :: mult_qr
    public :: solve_qr
    public :: qr_rank1_update

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/QR_decomposition)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/QRDecomposition.html)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface qr_factor
        module procedure :: qr_factor_no_pivot
        module procedure :: qr_factor_pivot
    end interface

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/QR_decomposition)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface form_qr
        module procedure :: form_qr_no_pivot
        module procedure :: form_qr_pivot
    end interface

! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization.
    interface mult_qr
        module procedure :: mult_qr_mtx
        module procedure :: mult_qr_vec
    end interface

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/QR_decomposition)
    !! - [LAPACK Users Manual](http://netlib.org/lapack/lug/node39.html)
    interface solve_qr
        module procedure :: solve_qr_no_pivot_mtx
        module procedure :: solve_qr_no_pivot_vec
        module procedure :: solve_qr_pivot_mtx
        module procedure :: solve_qr_pivot_vec
    end interface


contains
! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix without
    !! pivoting.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
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
    !!  - LA_ARRAY_SIZE_ERROR: Occurs if @p tau or @p work are not sized 
    !!      appropriately.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!
    !! @remarks
    !! QR factorization without pivoting is best suited to solving an
    !! overdetermined system in least-squares terms, or to solve a normally
    !! defined system.  To solve an underdetermined system, it is recommended to
    !! use either LQ factorization, or a column-pivoting based QR factorization.
    !!
    !! @par Usage
    !! To solve a system of M equations of N unknowns using QR factorization,
    !! the following code will suffice assuming M >= N.
    !! @code {.f90}
    !! ! Solve the system: A*X = B in a least-squares sense, where A is an
    !! ! M-by-N matrix, B is an M-by-NRHS matrix, and X is an N-by-NRHS matrix.
    !!
    !! ! Variables
    !! real(dp), dimension(m, n) :: a
    !! real(dp), dimension(m, nrhs) :: b, qtb
    !! real(dp), dimension(n, nrhs) :: x
    !! real(dp), dimension(n) :: tau
    !! real(dp), dimension(m, m) :: q
    !!
    !! ! Initialize A and B...
    !!
    !! ! Compute the QR factorization. We're intentionally not forming the full
    !! ! Q matrix, but instead storing it in terms of its elementary reflector
    !! ! components in the sub-diagonal portions of A, and the corresponding
    !! ! scalar factors in TAU.  Additionally, we'll let the algorithm allocate
    !! ! it's own workspace array; therefore, the call to factor A is:
    !! call qr_factor(a, tau)
    !!
    !! ! Solve A*X = B for X.  The first N rows of B are used to store X.
    !! call solve_qr(a, tau, b)
    !!
    !! ! Also note, we could form Q and R explicitly.  Then solution of the
    !! ! system of equations can be found.  First we form Q and R.
    !! call form_qr(a, tau, q) ! Forms Q, and R is stored in A
    !!
    !! ! Since we now have Q and R, we seek a solution to the equation:
    !! ! Q*R*X = B, but Q is an orthogonal matrix (i.e. Q**T = inv(Q)).
    !! ! Then: R*X = Q**T * B, and R is upper triangular; therefore, back
    !! ! substitution will suffice for a solution procedure.
    !! !
    !! ! Next, compute Q**T * B, and store in QTB.
    !! call mtx_mult(.true., .false., 1.0d0, q, b, 0.0d0, qtb)
    !!
    !! ! Copy the first N rows of Q**T * B into X for the solution process.
    !! ! Notice, only the first N rows are needed as rows N+1:M are all zero in
    !! ! matrix R.
    !! x = qtb(1:n,nrhs)
    !!
    !! ! Compute the solution and store in X
    !! call solve_triangular_system(.true., .true., .false., .true., 1.0d0, &
    !!  a(1:n,1:n), x)
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGEQRF.
    subroutine qr_factor_no_pivot(a, tau, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, mn, istat, lwork, flag
        real(dp), dimension(1) :: temp
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
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
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
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
    !> @brief Computes the QR factorization of an M-by-N matrix with column
    !! pivoting such that A * P = Q * R.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !!  elements on and above the diagonal contain the MIN(M, N)-by-N upper
    !!  trapezoidal matrix R (R is upper triangular if M >= N).  The elements
    !!  below the diagonal, along with the array @p tau, represent the
    !!  orthogonal matrix Q as a product of elementary reflectors.
    !! @param[out] tau A MIN(M, N)-element array used to store the scalar
    !!  factors of the elementary reflectors.
    !! @param[in,out] jpvt On input, an N-element array that if JPVT(I) .ne. 0,
    !!  the I-th column of A is permuted to the front of A * P; if JPVT(I) = 0,
    !!  the I-th column of A is a free column.  On output, if JPVT(I) = K, then
    !!  the I-th column of A * P was the K-th column of A.
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
    !! @par Usage
    !! To solve a system of M equations of N unknowns using QR factorization,
    !! the following code will suffice for any M and N.
    !! @code {.f90}
    !! ! Solve the least-squares (M >= N), or the underdetermined (M < N)
    !! ! problem A*X = B, where A is an M-by-N matrix, B is an M-by-NRHS matrix,
    !! ! and X is an N-by-NRHS matrix.  In the underdetermined case, or the
    !! ! case where the rank of matrix A is less than N, the solution obtained
    !! ! contains the fewest possible non-zero entries.
    !!
    !! ! Variables
    !! real(dp), dimension(m, n) :: a
    !! real(dp), dimension(n, nrhs) :: b
    !! real(dp), dimension(k) :: tau ! k = min(m, n)
    !! real(dp), dimension(m, m) :: q
    !! real(dp), dimension(n, n) :: p
    !! integer(i32), dimension(n) :: pvt
    !!
    !! ! Initialize A and B...
    !!
    !! ! Allow all columns to be free.
    !! pvt = 0
    !!
    !! ! Compute the QR factorization. We're intentionally not forming the full
    !! ! Q matrix, but instead storing it in terms of its elementary reflector
    !! ! components in the sub-diagonal portions of A, and the corresponding
    !! ! scalar factors in TAU.  Additionally, we'll let the algorithm allocate
    !! ! it's own workspace array; therefore, the call to factor A is:
    !! call qr_factor(a, tau, pvt)
    !!
    !! ! Solve A*X = B for X.  If M > N, the first N rows of B are used to store
    !! ! X.  If M < N, the input matrix B must be N-by-NRHS, and only the first
    !! ! M rows are used for the actual matrix B.  The remaining N-M rows
    !! ! can contain whatever as they are not referenced until they are
    !! ! overwritten by the N-by-NRHS solution matrix X.
    !! call solve_qr(a, tau, pvt, b)
    !!
    !! ! Notice, if the explicit Q matrix from the factorization is desired,
    !! ! the form_qr routine works similarly as in the no-pivot case;
    !! ! however, the permutation matrix P is also constructed.  The call would
    !! ! be as follows.  Also, as with the no-pivot algorithm, the matrix R is
    !! ! stored in matrix A.
    !! call form_qr(a, tau, pvt, q, p)
    !!
    !! ! Solution can proceed as per typical, but with a full Q matrix.  Also
    !! ! note, the problem is of the form: A*P = Q*R.  Solution is straight
    !! ! forward, as with the no-pivot case; however, if M < N, then R is upper
    !! ! trapezoidal, and must be appropriately partitioned to solve.  The rank
    !! ! of matrix R should be considered when applying the partition.
    !! @endcode
    !!
    !! @par Notes
    !! This routine utilizes the LAPACK routine DGEQP3.
    subroutine qr_factor_pivot(a, tau, jpvt, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(out), dimension(:) :: tau
        integer(i32), intent(inout), dimension(:) :: jpvt
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, mn, istat, lwork, flag
        real(dp), dimension(1) :: temp
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_factor_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DGEQP3(m, n, a, m, jpvt, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
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
    end subroutine

! ******************************************************************************
! FORM Q AND R MATRICES
! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p h.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
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
    !! This routine utilizes the LAPACK routine DORGQR.
    subroutine form_qr_no_pivot(r, tau, q, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), dimension(:,:) :: q
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n, mn, qcol, istat, flag, lwork
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("form_qr_no_pivot", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("form_qr_no_pivot", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("form_qr_no_pivot", &
                        "Incorrectly sized input array WORK, argument 4.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Forms the full M-by-M orthogonal matrix Q from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    !!
    !! @param[in,out] r On input, an M-by-N matrix where the elements below the
    !!  diagonal contain the elementary reflectors generated from the QR
    !!  factorization.  On and above the diagonal, the matrix contains the
    !!  matrix R.  On output, the elements below the diagonal are zeroed such
    !!  that the remaining matrix is simply the M-by-N matrix R.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  each elementary reflector defined in @p h.
    !! @param[in] pvt An N-element column pivot array as returned by the QR
    !!  factorization.
    !! @param[out] q An M-by-M matrix where the full orthogonal matrix Q will be
    !!  written.  In the event that M > N, Q may be supplied as M-by-N, and
    !!  therefore only return the useful submatrix Q1 (Q = [Q1, Q2]) as the
    !!  factorization can be written as Q * R = [Q1, Q2] * [R1; 0].
    !! @param[out] p An N-by-N matrix where the pivot matrix will be written.
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
    !! This routine utilizes the LAPACK routine DORGQR.
    subroutine form_qr_pivot(r, tau, pvt, q, p, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: r
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: pvt
        real(dp), intent(out), dimension(:,:) :: q, p
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: j, jp, m, n, mn, flag
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
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
    end subroutine

! ******************************************************************************
! Q MULTIPLICATION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Multiplies a general matrix by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C, or C = C * op(Q).
    !!
    !! @param[in] lside Set to true to apply Q or Q**T from the left; else, set
    !!  to false to apply Q or Q**T from the right.
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] a On input, an LDA-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  If @p lside is set to
    !!  true, LDA = M, and M >= K >= 0; else, if @p lside is set to false,
    !!  LDA = N, and N >= K >= 0.  Notice, the contents of this matrix are
    !!  restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-by-N matrix C.  On output, the product
    !!  of the orthogonal matrix Q and the original matrix C.
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
    !! This routine utilizes the LAPACK routine DORMQR.
    subroutine mult_qr_mtx(lside, trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: lside, trans
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:,:) :: a, c
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(i32) :: m, n, k, nrowa, istat, flag, lwork
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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_mtx", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, n, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("mult_qr_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("mult_qr_mtx", &
                        "Incorrectly sized input array WORK, argument 6.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Multiplies a vector by the orthogonal matrix Q from a QR
    !! factorization such that: C = op(Q) * C.
    !!
    !! @param[in] trans Set to true to apply Q**T; else, set to false.
    !! @param[in] a On input, an M-by-K matrix containing the elementary
    !!  reflectors output from the QR factorization.  Notice, the contents of
    !!  this matrix are restored on exit.
    !! @param[in] tau A K-element array containing the scalar factors of each
    !!  elementary reflector defined in @p a.
    !! @param[in,out] c On input, the M-element vector C.  On output, the
    !!  product of the orthogonal matrix Q and the original vector C.
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
    !! This routine is based upon the LAPACK routine DORM2R.
    subroutine mult_qr_vec(trans, a, tau, c, work, olwork, err)
        ! Arguments
        logical, intent(in) :: trans
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: c
        real(dp), intent(out), pointer, dimension(:), optional :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        character :: side, t
        integer(i32) :: m, k, nrowa, istat, flag, lwork
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        real(dp), dimension(1) :: temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("mult_qr_vec", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Workspace Query
        call DORMQR(side, t, m, 1, k, a, nrowa, tau, c, m, temp, -1, flag)
        lwork = int(temp(1), i32)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("mult_qr_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("mult_qr_vec", &
                        "Incorrectly sized input array WORK, argument 6.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
    end subroutine

! ******************************************************************************
! QR SOLUTION ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] b On input, the M-by-NRHS right-hand-side matrix.  On output,
    !!  the first N columns are overwritten by the solution matrix X.
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
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_mtx(a, tau, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a, b
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, nrhs, k, lwork, flag, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
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
        call mult_qr_mtx(.true., .true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 4.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
    !> @brief Solves a system of M QR-factored equations of N unknowns where
    !! M >= N.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are restored.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] b On input, the M-element right-hand-side vector.  On output,
    !!  the first N elements are overwritten by the solution vector X.
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
    !! This routine is based upon a subset of the LAPACK routine DGELS.
    subroutine solve_qr_no_pivot_vec(a, tau, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(i32) :: m, n, k, flag, lwork, istat
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
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
        call mult_qr_vec(.true., a, tau, b, olwork = lwork)
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_no_pivot_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_vec", &
                        "Incorrectly sized input array WORK, argument 4.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
        call mult_qr_vec(.true., a, tau, b, work = wptr)

        ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
        call solve_triangular_system(.true., .false., .true., a(1:n,1:n), b)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are altered.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] jpvt An N-element array, as output by @ref qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-by-NRHS matrix where the first M
    !!  rows contain the right-hand-side matrix B.  On output, the first N rows
    !!  are overwritten by the solution matrix X.
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
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_mtx(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: jpvt
        real(dp), intent(inout), dimension(:,:) :: b
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(i32), parameter :: imin = 2
        integer(i32), parameter :: imax = 1
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
            rnk, maxmn, flag, istat, lwork1, lwork2, lwork3
        real(dp) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(dp), pointer, dimension(:) :: wptr, w, tau2
        real(dp), allocatable, target, dimension(:) :: wrk
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
        call mult_qr_mtx(.true., .true., a, tau, b(1:m,:), olwork = lwork2)
        call mult_rz(.true., .true., n, a(1:mn,:), a(1:mn,1), b, &
            olwork = lwork3)
        lwork = max(lwork1, lwork2, lwork3, 2 * mn + 1) + mn
        if (present(olwork)) then
            olwork = lwork
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_pivot_mtx", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 5.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
            b(maxmn,nrhs) = zero
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
            call mult_rz(.true., .true., n - rnk, a(1:rnk,:), tau2, b, w)
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
    !> @brief Solves a system of M QR-factored equations of N unknowns where the
    !! QR factorization made use of column pivoting.
    !!
    !! @param[in] a On input, the M-by-N QR factored matrix as returned by
    !!  @ref qr_factor.  On output, the contents of this matrix are altered.
    !!  Notice, M must be greater than or equal to N.
    !! @param[in] tau A MIN(M, N)-element array containing the scalar factors of
    !!  the elementary reflectors as returned by @ref qr_factor.
    !! @param[in] jpvt An N-element array, as output by @ref qr_factor, used to
    !!  track the column pivots.
    !! @param[in] b On input, the MAX(M, N)-element array where the first M
    !!  elements contain the right-hand-side vector B.  On output, the first N
    !!  elements are overwritten by the solution vector X.
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
    !! This routine is based upon a subset of the LAPACK routine DGELSY.
    subroutine solve_qr_pivot_vec(a, tau, jpvt, b, work, olwork, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: a
        real(dp), intent(in), dimension(:) :: tau
        integer(i32), intent(in), dimension(:) :: jpvt
        real(dp), intent(inout), dimension(:) :: b
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        integer(i32), intent(out), optional :: olwork
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        integer(i32), parameter :: imin = 2
        integer(i32), parameter :: imax = 1
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn, flag, &
            istat, lwork1, lwork2
        real(dp) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
        real(dp), pointer, dimension(:) :: wptr, w, tau2
        real(dp), allocatable, target, dimension(:) :: wrk
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
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("solve_qr_pivot_vec", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("solve_qr_no_pivot_mtx", &
                        "Incorrectly sized input array WORK, argument 5.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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
! RANK 1 UPDATE
! ------------------------------------------------------------------------------
    !> @brief Computes the rank 1 update to an M-by-N QR factored matrix A
    !! (M >= N) where A = Q * R, and A1 = A + U * V**T such that A1 = Q1 * R1.
    !!
    !! @param[in,out] q On input, the original M-by-K orthogonal matrix Q.  On
    !!  output, the updated matrix Q1.
    !! @param[in,out] r On input, the M-by-N matrix R.  On output, the updated
    !!  matrix R1.
    !! @param[in,out] u On input, the M-element U update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[in,out] v On input, the N-element V update vector.  On output,
    !!  the original content of the array is overwritten.
    !! @param[out] work An optional argument that if supplied prevents local
    !!  memory allocation.  If provided, the array must have at least 2*K
    !!  elements.
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
    !! @par Remarks
    !! @verbatim
    !! Notice, K must either be equal to M, or to N.  In the event that K = N,
    !! only the submatrix Qa is updated.  This is appropriate as the QR
    !! factorization for an overdetermined system can be written as follows:
    !!  A = Q * R = [Qa, Qb] * [Ra]
    !!                         [0 ]
    !!
    !! Note: Ra is upper triangular of dimension N-by-N.
    !! @endverbatim
    !!
    !! @par Notes
    !! This routine utilizes the QRUPDATE routine DQR1UP.
    !!
    !! @par See Also
    !! [Source](https://sourceforge.net/projects/qrupdate/)
    subroutine qr_rank1_update(q, r, u, v, work, err)
        ! Arguments
        real(dp), intent(inout), dimension(:,:) :: q, r
        real(dp), intent(inout), dimension(:) :: u, v
        real(dp), intent(out), pointer, optional, dimension(:) :: work
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        logical :: full
        integer(i32) :: m, n, k, lwork, istat, flag
        real(dp), pointer, dimension(:) :: wptr
        real(dp), allocatable, target, dimension(:) :: wrk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

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
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qr_rank1_update", trim(errmsg), &
                LA_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        if (present(work)) then
            if (.not.associated(work)) then
                allocate(wrk(lwork), stat = istat)
                if (istat /= 0) then
                    ! ERROR: Out of memory
                    call errmgr%report_error("qr_rank1_update", &
                        "Insufficient memory available.", &
                        LA_OUT_OF_MEMORY_ERROR)
                    return
                end if
                wptr => wrk
            else
                if (size(work) < lwork) then
                    ! ERROR: WORK not sized correctly
                    call errmgr%report_error("qr_rank1_update", &
                        "Incorrectly sized input array WORK, argument 5.", &
                        LA_ARRAY_SIZE_ERROR)
                    return
                end if
                wptr => work(1:lwork)
            end if
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

        ! End
        if (allocated(wrk)) deallocate(wrk)
    end subroutine

! ------------------------------------------------------------------------------
end module
