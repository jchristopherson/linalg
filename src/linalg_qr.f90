! linalg_qr.f90
! SPDX-License-Identifier: MIT
module linalg_qr
    !! Provides QR factorization, updates, factor formation, and least-squares 
    !! solvers.
    use iso_fortran_env
    use linalg_errors
    use linalg_rz
    use linalg_tri
    use lapack
    use blas
    use qrupdate
    use linalg_types
    use linalg_basic, only : identity
    implicit none
    private
    public :: qr_factorization
    public :: qr_factor
    public :: form_qr
    public :: mult_qr
    public :: qr_rank1_update
    public :: solve_qr

    type, extends(matrix_factorization) :: qr_factorization
        !! A container for QR factorization of a matrix and associated
        !! operations.
        real(real64), private, allocatable, dimension(:,:) :: a
            !! The matrix, and after a call to factor, the elements on and above
            !! the diagonal contain the upper trapezoidal matrix R.  The 
            !! remainder of the matrix contains the elementray reflectors used
            !! to define Q along with the array tau.
        real(real64), private, allocatable, dimension(:) :: tau
            !! The scalar factors of the elementary reflectors.
        logical, private :: factored = .false.
            !! True if factored; else, false.
        logical, private :: use_pivoting = .false.
            !! True if column pivoting should be used; else, false.
        integer, private, allocatable, dimension(:) :: pivot
            !! The column pivot tracking array.
        real(real64), private, allocatable, dimension(:,:) :: solve_matrix
            !! Scratch matrix for the RZ factorization used by pivoted solves.
        real(real64), private, allocatable, dimension(:) :: solve_tau
            !! Reflector coefficients for the scratch RZ factorization.
        real(real64), private, allocatable, dimension(:) :: workspace
            !! Workspace array that is large enough for the factorization,
            !! multiplication, or solve calls.  This one workspace can be used
            !! for all operations
    contains
        procedure, public :: is_factored => qrf_is_factored
            !! Checks to see if the matrix has been factored.
        procedure, public :: size => qrf_size
            !! Gets the size of the factored matrix.
        procedure, public :: factor => qrf_factor
            !! Factors the matrix.
        procedure, public :: solve_multiple_rhs => qrf_solve_multiple_rhs
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure, public :: solve_single_rhs => qrf_solve_single_rhs
            !! Solves the factored system with a single right-hand-side vector.
        procedure, public :: solve_multiple_rhs_in_place => qrf_solve_multiple_rhs_in_place
            !! Solves the factored system with multiple right-hand-side vectors in place.
        procedure, public :: solve_single_rhs_in_place => qrf_solve_single_rhs_in_place
            !! Solves the factored system with a single right-hand-side vector.
        procedure, private :: qrf_multiply_mtx
            !! Multiplies the Q matrix with another matrix without explicitly
            !! forming the full Q matrix.
        procedure, private :: qrf_multiply_vec
            !! Multiplies the Q matrix with a vector without explicitly
            !! forming the full Q matrix.
        generic, public :: q_multiply => qrf_multiply_mtx, qrf_multiply_vec
            !! Multiplies a matrix or a vector with the orthogonal matrix 
            !! \(Q\) or its transpose.
        procedure, private :: estimate_rank => qrf_estimate_rank
            !! Estimates the rank of the pivoted triangular factor.
        procedure, public :: get_factored_form => qrf_get_factored_matrices
            !! Gets the factored matrices.
        procedure, public :: update => qrf_update
            !! Updates the factorization.
    end type

    interface qr_factorization
        module procedure :: qrf_init
    end interface

    interface qr_factor
        module procedure :: qr_factor_no_pivot
        module procedure :: qr_factor_no_pivot_cmplx
        module procedure :: qr_factor_pivot
        module procedure :: qr_factor_pivot_cmplx
    end interface

    interface form_qr
        module procedure :: form_qr_no_pivot
        module procedure :: form_qr_no_pivot_cmplx
        module procedure :: form_qr_pivot
        module procedure :: form_qr_pivot_cmplx
    end interface

    interface mult_qr
        module procedure :: mult_qr_mtx
        module procedure :: mult_qr_mtx_cmplx
        module procedure :: mult_qr_vec
        module procedure :: mult_qr_vec_cmplx
    end interface

    interface qr_rank1_update
        module procedure :: qr_rank1_update_dbl
        module procedure :: qr_rank1_update_cmplx
    end interface

    interface solve_qr
        module procedure :: solve_qr_no_pivot_mtx
        module procedure :: solve_qr_no_pivot_mtx_cmplx
        module procedure :: solve_qr_no_pivot_vec
        module procedure :: solve_qr_no_pivot_vec_cmplx
        module procedure :: solve_qr_pivot_mtx
        module procedure :: solve_qr_pivot_mtx_cmplx
        module procedure :: solve_qr_pivot_vec
        module procedure :: solve_qr_pivot_vec_cmplx
    end interface
contains
! ******************************************************************************
! QR_FACTORIZATION MEMBERS
! ------------------------------------------------------------------------------
function qrf_init(a, pvt) result(rst)
    !! Initializes a new [[qr_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to factor.  A copy of this matrix is made and stored by 
        !! this object; therefore, this matrix is not altered.
    logical, intent(in), optional :: pvt
        !! An optional parameter that, if set to true, enables column pivoting.
        !! The default is false such that no column pivoting will be utilized.
    type(qr_factorization) :: rst
        !! The new [[qr_factorization]] object.

    call rst%update(a, pvt = pvt)
end function

! ------------------------------------------------------------------------------
subroutine qrf_update(this, a, pvt)
    !! Updates the [[qr_factorization]] object with a new matrix.  The matrix
    !! is only stored, not factored; therefore a call to factor is still 
    !! required.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to factor.  A copy of this matrix is made and stored by
        !! the object; therefore, this matrix is not altered.  Also, any
        !! existing factorization is destroyed by this routine.
    logical, intent(in), optional :: pvt
        !! An optional parameter that, if set to true, enables column pivoting.
        !! The default is false such that no column pivoting will be utilized.

    integer :: m, n, mn, lwork, info, idummy(1)
    real(real64) :: dummy(1), w(1)

    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)

    allocate(this%a(m, n), source = a)
    allocate(this%tau(mn))

    this%factored = .false.
    this%use_pivoting = .false.
    if (present(pvt)) this%use_pivoting = pvt

    if (this%use_pivoting) then
        call DGEQP3(m, n, dummy, m, idummy, dummy, w, -1, info)
        allocate(this%pivot(n), source = 0)
        allocate(this%solve_matrix(max(1, m), n), source = 0.0d0)
        allocate(this%solve_tau(mn))
    else
        call DGEQRF(m, n, dummy, m, dummy, w, -1, info)
    end if
    lwork = max(int(w(1)), 2 * max(m, n) + 1)
    allocate(this%workspace(lwork))
end subroutine

! ------------------------------------------------------------------------------
pure function qrf_is_factored(this) result(rst)
    !! Gets a value querying whether the matrix has been factored.
    class(qr_factorization), intent(in) :: this
        !! The [[qr_factorization]] object.
    logical :: rst
        !! True if the matrix has been factored; else, false.
    rst = this%factored
end function

! ------------------------------------------------------------------------------
pure function qrf_size(this, dim) result(rst)
    !! Retrieves the size of teh requested matrix dimension.
    class(qr_factorization), intent(in) :: this
        !! The [[qr_factorization]] object.
    integer, intent(in) :: dim
        !! The dimension to query.
    integer :: rst
        !! The requested dimension size.
    if (allocated(this%a)) then
        rst = size(this%a, dim)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine qrf_factor(this, flag)
    !! Factors the matrix.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    integer, intent(out), optional :: flag
        !! An error flag

    ! Local Variables
    integer :: m, n, info, lwork

    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    m = this%size(1)
    n = this%size(2)
    lwork = size(this%workspace)

    ! If not allocated, return with an error
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if

    ! If factored, do nothing
    if (this%is_factored()) return

    ! Factor the matrix
    if (this%use_pivoting) then
        call DGEQP3(m, n, this%a, m, this%pivot, this%tau, this%workspace, &
            lwork, info)
    else
        call DGEQRF(m, n, this%a, m, this%tau, this%workspace, lwork, info)
    end if
    this%factored = .true.
end subroutine

! ------------------------------------------------------------------------------
function qrf_solve_multiple_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining multiple right-hand-side vectors.
    !! The solution has N rows.  Underdetermined systems require column
    !! pivoting; otherwise, an LA_INVALID_OPERATION_ERROR is returned.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: b
        !! The right-hand-side matrix.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The solution matrix.
    real(real64), allocatable, dimension(:,:) :: rhs_work

    ! Local Variables
    integer :: m, n, nrhs, info

    ! Initialization
    m = this%size(1)
    n = this%size(2)
    nrhs = size(b, 2)
    if (present(flag)) flag = LA_NO_ERROR

    ! Validation
    if (size(b, 1) /= m) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (.not.this%use_pivoting) then
        if (m < n) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end if

    ! Process
    if (m >= n) then
        allocate(rhs_work(m, nrhs), source = b)
    else
        allocate(rhs_work(n, nrhs), source = 0.0d0)
        rhs_work(1:m,:) = b
    end if
    call this%solve_multiple_rhs_in_place(rhs_work, info)
    if (info /= LA_NO_ERROR) then
        if (present(flag)) flag = info
        return
    end if
    allocate(rst(n, nrhs), source = rhs_work(1:n,:))
end function

! ------------------------------------------------------------------------------
function qrf_solve_single_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a single right-hand-side vector.
    !! The solution has N elements.  Underdetermined systems require column
    !! pivoting; otherwise, an LA_INVALID_OPERATION_ERROR is returned.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand-side vector.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:) :: rst
        !! The solution vector.
    real(real64), allocatable, dimension(:) :: rhs_work

    ! Local Variables
    integer :: m, n, info

    ! Initialization
    m = this%size(1)
    n = this%size(2)
    if (present(flag)) flag = LA_NO_ERROR

    ! Validation
    if (size(b) /= m) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (.not.this%use_pivoting) then
        if (m < n) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end if

    ! Process
    if (m >= n) then
        allocate(rhs_work(m), source = b)
    else
        allocate(rhs_work(n), source = 0.0d0)
        rhs_work(1:m) = b
    end if
    call this%solve_single_rhs_in_place(rhs_work, info)
    if (info /= LA_NO_ERROR) then
        if (present(flag)) flag = info
        return
    end if
    allocate(rst(n), source = rhs_work(1:n))
end function

! ------------------------------------------------------------------------------
subroutine qrf_solve_multiple_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining multiple right-hand-side vectors.
    !! The solution is returned in the first N rows of B, which must have
    !! MAX(M,N) rows.  Underdetermined systems require column pivoting.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the M-by-NRHS right-hand side occupies the first M rows.
        !! On output, the N-by-NRHS solution occupies the first N rows.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: m, n, nrhs, mn, info, lwork, block_size, first, last, rank
    integer :: i, j, ldrz
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    m = this%size(1)
    n = this%size(2)
    nrhs = size(b, 2)
    if (present(flag)) flag = LA_NO_ERROR

    ! Validation
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (.not.this%use_pivoting) then
        if (m < n) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end if
    if (size(b, 1) /= max(m, n)) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if

    ! Ensure the matrix is factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Process
    if (n == 0 .or. nrhs == 0) return
    mn = min(m, n)
    lwork = size(this%workspace)
    block_size = min(nrhs, lwork)
    if (this%use_pivoting) then
        call this%estimate_rank(rank)
    else
        rank = n
    end if
    if (rank == 0) then
        b = 0.0d0
        return
    end if

    ! Apply Q^T to each RHS block.  The block width is bounded by the stored
    ! workspace, so the in-place solve never allocates a temporary work array.
    do first = 1, nrhs, block_size
        last = min(nrhs, first + block_size - 1)
        call DORMQR('L', 'T', m, last - first + 1, mn, this%a, m, this%tau, &
            b(:,first:last), size(b,1), this%workspace, lwork, info)
        if (info /= 0) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end do

    if (this%use_pivoting .and. rank < n) then
        ! RZ-factor the numerical-rank portion of R.  Clear the lower triangle
        ! because it contains the QR reflectors, not entries of R.
        ldrz = size(this%solve_matrix, 1)
        this%solve_matrix(1:rank,1:n) = this%a(1:rank,1:n)
        do j = 1, min(rank, n)
            if (j < rank) this%solve_matrix(j+1:rank,j) = 0.0d0
        end do
        call DTZRZF(rank, n, this%solve_matrix, ldrz, this%solve_tau, &
            this%workspace, lwork, info)
        if (info /= 0) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
        call DTRSM('L', 'U', 'N', 'N', rank, nrhs, one, this%solve_matrix, &
            ldrz, b, size(b,1))
        b(rank+1:n,:) = 0.0d0
        do first = 1, nrhs, block_size
            last = min(nrhs, first + block_size - 1)
            call DORMRZ('L', 'T', n, last - first + 1, rank, n - rank, &
                this%solve_matrix, ldrz, this%solve_tau, b(:,first:last), &
                size(b,1), this%workspace, lwork, info)
            if (info /= 0) then
                if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
                return
            end if
        end do
    else
        call DTRSM('L', 'U', 'N', 'N', n, nrhs, one, this%a, m, b, size(b,1))
    end if

    if (this%use_pivoting) then
        do j = 1, nrhs
            this%workspace(1:n) = b(1:n,j)
            do i = 1, n
                b(this%pivot(i),j) = this%workspace(i)
            end do
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine qrf_solve_single_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining a single right-hand-side vector.
    !! The solution is returned in the first N elements of B, which must have
    !! MAX(M,N) elements.  Underdetermined systems require column pivoting.
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(inout), dimension(:) :: b
        !! On input, the M-element right-hand side.  On output, the N-element
        !! solution occupies the first N elements.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: m, n, mn, info, lwork, rank, i, ldrz
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    m = this%size(1)
    n = this%size(2)
    if (present(flag)) flag = LA_NO_ERROR

    ! Validation
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (.not.this%use_pivoting) then
        if (m < n) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end if
    if (size(b) /= max(m, n)) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if

    ! Ensure the matrix is factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Process
    if (n == 0) return
    mn = min(m, n)
    lwork = size(this%workspace)
    if (this%use_pivoting) then
        call this%estimate_rank(rank)
    else
        rank = n
    end if
    if (rank == 0) then
        b = 0.0d0
        return
    end if

    call DORMQR('L', 'T', m, 1, mn, this%a, m, this%tau, b, size(b), &
        this%workspace, lwork, info)
    if (info /= 0) then
        if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
        return
    end if

    if (this%use_pivoting .and. rank < n) then
        ldrz = size(this%solve_matrix, 1)
        this%solve_matrix(1:rank,1:n) = this%a(1:rank,1:n)
        do i = 1, min(rank, n)
            if (i < rank) this%solve_matrix(i+1:rank,i) = 0.0d0
        end do
        call DTZRZF(rank, n, this%solve_matrix, ldrz, this%solve_tau, &
            this%workspace, lwork, info)
        if (info /= 0) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
        call DTRSV('U', 'N', 'N', rank, this%solve_matrix, ldrz, b, 1)
        b(rank+1:n) = 0.0d0
        call DORMRZ('L', 'T', n, 1, rank, n - rank, this%solve_matrix, &
            ldrz, this%solve_tau, b, size(b), this%workspace, lwork, info)
        if (info /= 0) then
            if (present(flag)) flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    else
        call DTRSV('U', 'N', 'N', n, this%a, m, b, 1)
    end if

    if (this%use_pivoting) then
        this%workspace(1:n) = b(1:n)
        do i = 1, n
            b(this%pivot(i)) = this%workspace(i)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine qrf_estimate_rank(this, rank)
    !! Estimates the numerical rank of the pivoted triangular factor using the
    !! incremental condition estimator used by the direct pivoted QR solvers.
    class(qr_factorization), intent(inout) :: this
        !! The factored QR object whose workspace is used for temporary vectors.
    integer, intent(out) :: rank
        !! Estimated numerical rank.

    ! Parameters
    integer, parameter :: imin = 2
    integer, parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer :: m, n, mn, ismin, ismax, i
    real(real64) :: rcond, smax, smin, smax_next, smin_next
    real(real64) :: s_min, c_min, s_max, c_max

    m = size(this%a, 1)
    n = size(this%a, 2)
    mn = min(m, n)
    rank = 0
    if (mn == 0) return

    ismin = 1
    ismax = mn + 1
    rcond = epsilon(one)
    smax = abs(this%a(1,1))
    smin = smax
    if (smax == zero) return

    rank = 1
    this%workspace(ismin) = one
    this%workspace(ismax) = one
    do while (rank < mn)
        i = rank + 1
        call DLAIC1(imin, rank, this%workspace(ismin:ismin+rank-1), smin, &
            this%a(1:rank,i), this%a(i,i), smin_next, s_min, c_min)
        call DLAIC1(imax, rank, this%workspace(ismax:ismax+rank-1), smax, &
            this%a(1:rank,i), this%a(i,i), smax_next, s_max, c_max)
        if (smax_next * rcond > smin_next) exit

        this%workspace(ismin:ismin+rank-1) = &
            s_min * this%workspace(ismin:ismin+rank-1)
        this%workspace(ismax:ismax+rank-1) = &
            s_max * this%workspace(ismax:ismax+rank-1)
        this%workspace(ismin+rank) = c_min
        this%workspace(ismax+rank) = c_max
        smin = smin_next
        smax = smax_next
        rank = rank + 1
    end do
end subroutine

! ------------------------------------------------------------------------------
function qrf_multiply_mtx(this, lside, trans, c, flag) result(rst)
    !! Multiplies the orthogonal matrix \(Q\) with another matrix
    !! \(C\) such that the operation is either \(op(Q) C \) or \(C op(Q)\)
    !! where \(op(Q)\) is either \(Q**T\) or \(Q\).
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^T\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^T\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:,:) :: c
        !! The matrix
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    character :: side, trns
    integer :: m, n, k, lda, lwork, info

    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    m = size(c, 1)
    n = size(c, 2)
    k = size(this%tau)
    if (lside) then
        lda = m
    else
        lda = n
    end if
    if (lside) then
        side = "L"
    else
        side = "R"
    end if
    if (trans) then
        trns = "T"
    else
        trns = "N"
    end if
    lwork = size(this%workspace)

    ! Validation
    if (this%size(1) /= lda) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if

    ! Ensure the matrix has been factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Process
    allocate(rst(m, n), source = c)
    call DORMQR(side, trns, m, n, k, this%a, lda, this%tau, rst, m, &
        this%workspace, lwork, info)
end function

! ------------------------------------------------------------------------------
function qrf_multiply_vec(this, trans, c, flag) result(rst)
    !! Multiplies tthe orthogonal matrix \(Q\) with a vector \(c\) such that 
    !! the operation is either \(op(Q) c \) where \(op(Q)\) is either \(Q**T\)
    !! or \(Q\).
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:) :: c
        !! The vector
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    character :: side, trns
    integer :: m, k, lwork, info

    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    m = size(c)
    k = size(this%tau)
    side = "L"
    if (trans) then
        trns = "T"
    else
        trns = "N"
    end if
    lwork = size(this%workspace)

    ! Validation
    if (this%size(1) /= m) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if

    ! Ensure the matrix has been factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Process
    allocate(rst(m), source = c)
    call DORMQR(side, trns, m, 1, k, this%a, m, this%tau, rst, m, &
        this%workspace, lwork, info)
end function

! ------------------------------------------------------------------------------
subroutine qrf_get_factored_matrices(this, q, r, p, flag)
    !! Forms the full M-by-M orthogonal matrix \(Q\), the upper triangular
    !! matrix \(R\), and if pivoting was used, the pivot matrix \(P\), where
    !! either \(A = Q R\) or \(A P = Q R\).
    class(qr_factorization), intent(inout) :: this
        !! The [[qr_factorization]] object.
    real(real64), intent(out), allocatable, dimension(:,:) :: q
        !! The M-by-M matrix \(Q\).
    real(real64), intent(out), allocatable, dimension(:,:) :: r
        !! The M-by-N matrix \(R\).
    real(real64), intent(out), optional, allocatable, dimension(:,:) :: p
        !! The N-by-N column pivot tracking matrix \(P\)
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: m, n, info

    ! Initialization
    m = this%size(1)
    n = this%size(2)

    ! Ensure the matrix has been factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Form the matrices
    allocate(q(m, m))
    allocate(r(m, n), source = this%a)
    if (present(p)) then
        if (this%use_pivoting) then
            allocate(p(n, n))
            call form_qr(r, this%tau, this%pivot, q, p)
        else
            ! P is requested, but no pivoting was utilized
            p = identity(n)
            call form_qr(r, this%tau, q)
        end if
    else
        call form_qr(r, this%tau, q)
    end if
end subroutine

! ******************************************************************************
! PROCEDURAL INTERFACE
! ------------------------------------------------------------------------------
pure subroutine qr_factor_no_pivot(a, tau, qr, q, r)
    !! Computes the QR factorization of an M-by-N matrix.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    real(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, target, dimension(:,:) :: ac, qc, rc
    real(real64), pointer, dimension(:) :: tptr
    real(real64), pointer, dimension(:,:) :: aptr, qptr, rptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r)

    ! Workspace Query
    call DGEQRF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGEQRF
    call DGEQRF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if
        call form_qr(rptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_no_pivot_cmplx(a, tau, qr, q, r)
    !! Computes the QR factorization of an M-by-N matrix.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    complex(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    complex(real64), dimension(1) :: temp
    complex(real64), allocatable, target, dimension(:) :: w, tc
    complex(real64), allocatable, target, dimension(:,:) :: ac, qc, rc
    complex(real64), pointer, dimension(:) :: tptr
    complex(real64), pointer, dimension(:,:) :: aptr, qptr, rptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r)
    
    ! Workspace Query
    call ZGEQRF(m, n, temp, m, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZGEQRF
    call ZGEQRF(m, n, aptr, m, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if
        call form_qr(rptr, tptr, qptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_pivot(a, jpvt, tau, qr, q, r, p)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
    real(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N column-pivot tracking matrix \(P\) where \(A P = Q R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    real(real64), dimension(1) :: temp
    real(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, target, dimension(:,:) :: ac, qc, rc, pc
    real(real64), pointer, dimension(:) :: tptr
    real(real64), pointer, dimension(:,:) :: aptr, qptr, rptr, pptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r) .or. present(p)

    ! Input Check
    if (size(jpvt) /= n) then
        error stop 2
    end if

    ! Workspace Query
    call DGEQP3(m, n, temp, m, jpvt, temp, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DGEQP3
    call DGEQP3(m, n, aptr, m, jpvt, tptr, w, lwork, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Form the matrices
        call form_qr(rptr, tptr, jpvt, qptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_factor_pivot_cmplx(a, jpvt, tau, qr, q, r, p)
    !! Computes the QR factorization of an M-by-N matrix using column pivoting
    !! such that \(A P = Q R\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N matrix to factor.
    integer(int32), intent(inout), dimension(:) :: jpvt
        !! On input, an N-element array that if JPVT(I) .ne. 0, the I-th column 
        !! of A is permuted to the front of A * P; if JPVT(I) = 0, the I-th 
        !! column of A is a free column.  On output, if JPVT(I) = K, then the 
        !! I-th column of A * P was the K-th column of A.
    complex(real64), intent(out), allocatable, optional, target, dimension(:) :: tau
        !! A MIN(M, N)-element array used to store the scalar factors of the 
        !! elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: qr
        !! The M-by-N factored matrix stored such that the elements on and above
        !! the diagonal contain the MIN(M, N)-by-N upper trapezoidal matrix 
        !! \(R\) (\(R\) is upper triangular if M >= N).  The elements below the
        !! diagonal, along with the array tau, represent the orthogonal matrix
        !! \(Q\) as a product of elementary reflectors.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: q
        !! The M-by-M orthogonal matrix \(Q\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: r
        !! The M-by-N upper trapezoidal matrix \(R\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N column-pivot tracking matrix \(P\) where \(A P = Q R\).

    ! Local Variables
    logical :: buildqr
    integer(int32) :: m, n, mn, lwork, flag
    complex(real64), dimension(1) :: temp
    real(real64), dimension(1) :: rtemp
    complex(real64), allocatable, target, dimension(:) :: w, tc
    real(real64), allocatable, dimension(:) :: rw
    complex(real64), allocatable, target, dimension(:,:) :: ac, qc, rc, pc
    complex(real64), pointer, dimension(:) :: tptr
    complex(real64), pointer, dimension(:,:) :: aptr, qptr, rptr, pptr

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    if (present(tau)) then
        allocate(tau(mn))
        tptr => tau
    else
        allocate(tc(mn))
        tptr => tc
    end if
    if (present(qr)) then
        allocate(qr(m, n), source = a)
        aptr => qr
    else
        allocate(ac(m, n), source = a)
        aptr => ac
    end if
    buildqr = present(q) .or. present(r) .or. present(p)

    ! Input Check
    if (size(jpvt) /= n) then
        error stop 2
    end if

    ! Workspace Query
    call ZGEQP3(m, n, temp, m, jpvt, temp, temp, -1, rtemp, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork), rw(2 * n))

    ! Call ZGEQP3
    call ZGEQP3(m, n, aptr, m, jpvt, tptr, w, lwork, rw, flag)

    ! Build Q & R
    if (buildqr) then
        ! Q
        if (present(q)) then
            allocate(q(m, m))
            qptr => q
        else
            allocate(qc(m, m))
            qptr => qc
        end if

        ! R
        if (present(r)) then
            allocate(r(m, n), source = aptr)
            rptr => r
        else
            if (allocated(ac)) then
                rptr => ac
            else
                allocate(rc(m, n), source = aptr)
                rptr => rc
            end if
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Form the matrices
        call form_qr(rptr, tptr, jpvt, qptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_no_pivot(r, tau, q)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    real(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).

    ! Parameters
    real(real64), parameter :: zero = 0.0d0

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        error stop 3
    end if
    if (qcol == n .and. m < n) then
        error stop 3
    end if

    ! Workspace Query
    call DORGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call DORGQR(m, qcol, mn, q, m, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_no_pivot_cmplx(r, tau, q)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    complex(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, m, n, mn, qcol, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)
    qcol = size(q, 2)

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(q, 1) /= m .or. (qcol /= m .and. qcol /= n)) then
        error stop 3
    end if
    if (qcol == n .and. m < n) then
        error stop 3
    end if

    ! Workspace Query
    call ZUNGQR(m, qcol, mn, q, m, tau, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Copy the sub-diagonal portion of R to Q, and then zero out the
    ! sub-diagonal portion of R
    do j = 1, mn
        q(j+1:m,j) = r(j+1:m,j)
        r(j+1:m,j) = zero
    end do

    ! Build Q - Build M-by-M or M-by-N, but M-by-N only for M >= N
    call ZUNGQR(m, qcol, mn, q, m, tau, w, lwork, flag)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_pivot(r, tau, pvt, q, p)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    integer(int32), intent(in), dimension(:) :: pvt
        !! An N-element column pivot array as returned by the QR factorization.
    real(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the pivot matrix will be written.

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: j, jp, m, n, mn

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 5
    end if

    ! Generate Q and R
    call form_qr_no_pivot(r, tau, q)

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_qr_pivot_cmplx(r, tau, pvt, q, p)
    !! Forms the full M-by-M orthogonal matrix \(Q\) from the elementary
    !! reflectors returned by the base QR factorization algorithm.
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, an M-by-N matrix where the elements below the diagonal 
        !! contain the elementary reflectors generated from the QR 
        !! factorization.  On and above the diagonal, the matrix contains the
        !! matrix \(R\).  On output, the elements below the diagonal are zeroed 
        !! such that the remaining matrix is simply the M-by-N matrix \(R\).
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of each 
        !! elementary reflector defined in \(R\).
    integer(int32), intent(in), dimension(:) :: pvt
        !! An N-element column pivot array as returned by the QR factorization.
    complex(real64), intent(out), dimension(:,:) :: q
        !! An M-by-M matrix where the full orthogonal matrix \(Q\) will be
        !! written.  In the event that M > N, \(Q\) may be supplied as M-by-N, 
        !! and therefore only return the useful submatrix \(Q_1\)
        !! \(Q = [Q_1 Q_2]\) as the factorization can be written as 
        !! \(Q R = [Q_1, Q_2] [R1 0]^T\).
    complex(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the pivot matrix will be written.

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: j, jp, m, n, mn, flag

    ! Initialization
    m = size(r, 1)
    n = size(r, 2)
    mn = min(m, n)

    ! Input Check
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 5
    end if

    ! Generate Q and R
    call form_qr_no_pivot_cmplx(r, tau, q)

    ! Form P
    do j = 1, n
        jp = pvt(j)
        p(:,j) = zero
        p(jp,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function mult_qr_mtx(lside, trans, a, tau, c) result(qc)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^T\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^T\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).
    real(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of \(C\) and \(Q\).

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m, n), source = c)
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

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            error stop 3
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            error stop 3
        end if
    end if

    ! Workspace Query
    call DORMQR(side, t, m, n, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMQR
    call DORMQR(side, t, m, n, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_mtx_cmplx(lside, trans, a, tau, c) result(qc)
    !! Multiplies a general matrix by the orthogonal matrix \(Q\) from a QR
    !! factorization such that \(C = op(Q) C\) or \(C = C op(Q)\).
    logical, intent(in) :: lside
        !! Set to true to apply \(Q\) or \(Q^H\) from the left; else, set to 
        !! false to apply \(Q\) or \(Q^H\) from the right.
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an LDA-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.  If lside is set to true, LDA = M, 
        !! and M >= K >= 0; else, if lside is set to false, LDA = N, and 
        !! N >= K >= 0.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(in), dimension(:,:) :: c
        !! The M-by-N matrix \(C\).
    complex(real64), allocatable, dimension(:,:) :: qc
        !! The M-by-N product of \(C\) and \(Q\).

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, n, k, nrowa, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c, 1)
    n = size(c, 2)
    k = size(tau)
    allocate(qc(m, n), source = c)
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

    ! Input Check
    if (lside) then
        ! A is M-by-K, M >= K >= 0
        if (size(a, 1) /= m .or. size(a, 2) < k) then
            error stop 3
        end if
    else
        ! A is N-by-K, N >= K >= 0
        if (size(a, 1) /= n .or. size(a, 2) < k) then
            error stop 3
        end if
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, n, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_vec(trans, a, tau, c) result(qc)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^T\); else, set to false to apply \(Q\).
    real(real64), intent(in), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization.
    real(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    real(real64), intent(in), dimension(:) :: c
        !! The M-element vector \(\vec{c}\).
    real(real64), allocatable, dimension(:) :: qc
        !! The product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, flag, lwork
    real(real64), allocatable, dimension(:) :: w
    real(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    allocate(qc(m), source = c)
    if (trans) then
        t = 'T'
    else
        t = 'N'
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        error stop 2
    end if

    ! Workspace Query
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call DORMQR
    call DORMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure function mult_qr_vec_cmplx(trans, a, tau, c) result(qc)
    !! Multiplies a vector by the orthogonal matrix \(Q\) from a QR 
    !! factorization such that \(\vec{c} = op(Q) \vec{c}\).
    logical, intent(in) :: trans
        !! Set to true to apply \(Q^H\); else, set to false to apply \(Q\).
    complex(real64), intent(in), dimension(:,:) :: a
        !! On input, an M-by-K matrix containing the elementary reflectors 
        !! output from the QR factorization. Notice, the contents of this matrix
        !! are restored on exit.
    complex(real64), intent(in), dimension(:) :: tau
        !! A K-element array containing the scalar factors of each elementary 
        !! reflector defined in\(A\).
    complex(real64), intent(in), dimension(:) :: c
        !! The M-element vector \(\vec{c}\).
    complex(real64), allocatable, dimension(:) :: qc
        !! The product of the orthogonal matrix \(Q\) and the original vector 
        !! \(\vec{c}\).

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    character :: side, t
    integer(int32) :: m, k, nrowa, flag, lwork
    complex(real64), allocatable, dimension(:) :: w
    complex(real64), dimension(1) :: temp

    ! Initialization
    m = size(c)
    k = size(tau)
    side = 'L'
    nrowa = m
    allocate(qc(m), source = c)
    if (trans) then
        t = 'C'
    else
        t = 'N'
    end if

    ! Input Check
    flag = 0
    if (size(a, 1) /= m .or. size(a, 2) < k) then
        error stop 2
    end if

    ! Workspace Query
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, temp, -1, flag)
    lwork = int(temp(1), int32)
    allocate(w(lwork))

    ! Call ZUNMQR
    call ZUNMQR(side, t, m, 1, k, a, nrowa, tau, qc, m, w, lwork, flag)
end function

! ------------------------------------------------------------------------------
pure subroutine qr_rank1_update_dbl(q, r, u, v)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^T\) such that
    !! \(A_1 = Q_1 R_1\).
    real(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    real(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    real(real64), intent(in), dimension(:) :: u
        !! The M-element \(\vec{u}\) update vector.
    real(real64), intent(in), dimension(:) :: v
        !! The N-element \(\vec{v}\) update vector.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork
    real(real64), allocatable, dimension(:) :: w, uc, vc

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = 2 * k
    allocate(w(lwork))

    ! Input Check
    if (m < n) then
        error stop 2
    end if
    if (.not.full .and. size(q, 2) /= k) then
        error stop 1
    end if
    if (size(r, 1) /= m) then
        error stop 2
    end if
    if (size(u) /= m) then
        error stop 3
    end if
    if (size(v) /= n) then
        error stop 4
    end if

    ! Local Memory Allocation
    allocate(uc(m), source = u)
    allocate(vc(n), source = v)

    ! Process
    call DQR1UP(m, n, k, q, m, r, m, uc, vc, w)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine qr_rank1_update_cmplx(q, r, u, v)
    !! Computes the rank-1 update to an M-by-N QR factored matrix \(A\) where
    !! \(M \ge N\), \(A = Q R\), and \(A_1 = A + \vec{u} \vec{v}^H\) such that
    !! \(A_1 = Q_1 R_1\).
    complex(real64), intent(inout), dimension(:,:) :: q
        !! On input, the original M-by-K orthogonal matrix \(Q\).  On output, 
        !! the updated matrix \(Q_1\).
    complex(real64), intent(inout), dimension(:,:) :: r
        !! On input, the M-by-N matrix \(R\).  On output, the updated matrix 
        !! \(R_1\).
    complex(real64), intent(in), dimension(:) :: u
        !! On input, the M-element \(\vec{u}\) update vector.  On output, the 
        !! original content of the array is overwritten.
    complex(real64), intent(in), dimension(:) :: v
        !! On input, the N-element \(\vec{v}\) update vector.  On output, the 
        !! original content of the array is overwritten.

    ! Local Variables
    logical :: full
    integer(int32) :: m, n, k, lwork, lrwork
    complex(real64), allocatable, dimension(:) :: w, uc, vc
    real(real64), allocatable, dimension(:) :: rw

    ! Initialization
    m = size(u, 1)
    n = size(r, 2)
    k = min(m, n)
    full = size(q, 2) == m
    lwork = k
    lrwork = k
    allocate(w(lwork), rw(lrwork))

    ! Input Check
    if (m < n) then
        error stop 2
    end if
    if (.not.full .and. size(q, 2) /= k) then
        error stop 1
    end if
    if (size(r, 1) /= m) then
        error stop 2
    end if
    if (size(u) /= m) then
        error stop 3
    end if
    if (size(v) /= n) then
        error stop 4
    end if

    ! Local Memory Allocation
    allocate(uc(m), source = u)
    allocate(vc(n), source = v)

    ! Process
    call ZQR1UP(m, n, k, q, m, r, m, uc, vc, w, rw)
end subroutine

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_mtx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: m, n, nrhs, k
    real(real64), allocatable, dimension(:,:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
    x = solve_triangular_system(.true., .true., .false., .true., one, &
        a(1:n,1:n), qtb(1:n,:))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_mtx_cmplx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS right-hand-side matrix.

    ! Parameters
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: m, n, nrhs, k
    complex(real64), allocatable, dimension(:,:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b, 1) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N,:)
    x = solve_triangular_system(.true., .true., .false., .true., one, &
        a(1:n,1:n), qtb(1:n,:))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_vec(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: m, n, k
    real(real64), allocatable, dimension(:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
    x = solve_triangular_system(.true., .false., .true., a(1:n,1:n), qtb(1:n))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_no_pivot_vec_cmplx(a, tau, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.  M must be
    !! greater than or equal to N.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: m, n, k
    complex(real64), allocatable, dimension(:) :: qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    ! Input Check
    if (m < n) then
        error stop 1
    end if
    if (size(tau) /= k) then
        error stop 2
    end if
    if (size(b) /= m) then
        error stop 3
    end if

    ! Compute Q**T * B
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system: A(1:N,1:N)*X = B(1:N)
    x = solve_triangular_system(.true., .false., .true., a(1:n,1:n), qtb(1:n))
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_mtx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(in), dimension(:,:) :: b
        !! The MAX(M, N)-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), allocatable, dimension(:) :: w, t2
    real(real64), allocatable, dimension(:,:) :: rz, qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n, nrhs))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call DLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call DLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system T11 * X(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
    if (rnk < n) then
        x(1:rnk,:) = solve_triangular_system(.true., .true., .false., .true., &
            one, rz(1:rnk,1:rnk), qtb(1:rnk,:))
        x(rnk+1:n,:) = zero

        ! Compute X(1:n,1:nrhs) = Y**T * X(1:n,1:nrhs)
        x = mult_rz(.true., .true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), qtb(1:rnk,:))
    end if

    ! Apply the pivoting: X(1:N,1:NRHS) = P * X(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            w(jpvt(i)) = x(i,j)
        end do
        x(1:n,j) = w(1:n)
    end do
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_mtx_cmplx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The M-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, j, m, n, mn, nrhs, lwork, ismin, ismax, &
        rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), allocatable, dimension(:) :: w, t2
    complex(real64), allocatable, dimension(:,:) :: rz, qtb

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    nrhs = size(b, 2)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n, nrhs))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call ZLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., .true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk,1:nrhs) = B(1:rnk,1:nrhs)
    if (rnk < n) then
        x(1:rnk,:) = solve_triangular_system(.true., .true., .false., .true., &
            one, rz(1:rnk,1:rnk), qtb(1:rnk,:))
        x(rnk+1:n,:) = zero

        ! Compute X(1:n,1:nrhs) = Y**T * X(1:n,1:nrhs)
        x = mult_rz(.true., .true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .true., .false., .true., one, &
            a(1:rnk,1:rnk), qtb(1:rnk,:))
    end if

    ! Apply the pivoting: X(1:N,1:NRHS) = P * X(1:N,1:NRHS)
    do j = 1, nrhs
        do i = 1, n
            w(jpvt(i)) = x(i,j)
        end do
        x(1:n,j) = w(1:n)
    end do
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_vec(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    real(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    real(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    real(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr, s1, c1, s2, c2
    real(real64), allocatable, dimension(:) :: w, t2, qtb
    real(real64), allocatable, dimension(:,:) :: rz

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call DLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call DLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
    if (rnk < n) then
        x(1:rnk) = solve_triangular_system(.true., .false., .true., &
            rz(1:rnk,1:rnk), qtb(1:rnk))
        x(rnk+1:n) = zero

        ! Compute X(1:n) = Y**T * X(1:n)
        x = mult_rz(.true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            qtb(1:rnk))
    end if

    ! Apply the pivoting: X(1:N) = P * X(1:N)
    do i = 1, n
        w(jpvt(i)) = x(i)
    end do
    x(1:n) = w(1:n)
end function

! ------------------------------------------------------------------------------
pure function solve_qr_pivot_vec_cmplx(a, tau, jpvt, b) result(x)
    !! Solves a system of M QR-factored equations of N unknowns.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The M-by-N QR factored matrix as returned by qr_factor.
    complex(real64), intent(in), dimension(:) :: tau
        !! A MIN(M, N)-element array containing the scalar factors of the 
        !! elementary reflectors as returned by qr_factor.
    integer(int32), intent(in), dimension(:) :: jpvt
        !! An N-element array, as output by qr_factor, used to track the 
        !! column pivots.
    complex(real64), intent(in), dimension(:) :: b
        !! The M-element right-hand-side vector.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution vector.

    ! Parameters
    integer(int32), parameter :: imin = 2
    integer(int32), parameter :: imax = 1
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: i, m, n, mn, lwork, ismin, ismax, rnk, maxmn
    real(real64) :: rcond, smax, smin, smaxpr, sminpr
    complex(real64) :: s1, c1, s2, c2
    complex(real64), allocatable, dimension(:) :: w, t2, qtb
    complex(real64), allocatable, dimension(:,:) :: rz

    ! Initialization
    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)
    maxmn = max(m, n)
    ismin = 1
    ismax = mn + 1
    rcond = epsilon(rcond)
    lwork = 2 * maxmn + 1
    allocate(w(lwork), source = one)
    allocate(x(n))

    ! Input Check
    if (size(tau) /= mn) then
        error stop 2
    end if
    if (size(jpvt) /= n) then
        error stop 3
    end if
    if (size(b, 1) /= m) then
        error stop 4
    end if

    ! Determine the rank of R11 using an incremental condition estimation
    smax = abs(a(1,1))
    smin = smax
    if (abs(a(1,1)) == zero) then
        rnk = 0
        x = zero
        return
    else
        rnk = 1
    end if
    do
        if (rnk < mn) then
            i = rnk + 1
            call ZLAIC1(imin, rnk, w(ismin:ismin+rnk-1), smin, &
                a(1:rnk-1,i), a(i,i), sminpr, s1, c1)
            call ZLAIC1(imax, rnk, w(ismax:ismax+rnk-1), smax, &
                a(1:rnk-1,i), a(i,i), smaxpr, s2, c2)
            if (smaxpr * rcond <= sminpr) then
                do i = 1, rnk
                    w(ismin+i-1) = s1 * w(ismin+i-1)
                    w(ismax+i-1) = s2 * w(ismax+i-1)
                end do
                w(ismin+rnk) = c1
                w(ismax+rnk) = c2
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
    ! [R11 R12] = [T11 0] * Y
    if (rnk < n) then
        call rz_factor(a(1:rnk,:), t2, rz)
    end if

    ! Compute B(1:m,1:NRHS) = Q**T * B(1:M,1:NRHS)
    qtb = mult_qr(.true., a, tau, b)

    ! Solve the triangular system T11 * B(1:rnk) = B(1:rnk)
    if (rnk < n) then
        x(1:rnk) = solve_triangular_system(.true., .false., .true., &
            rz(1:rnk,1:rnk), qtb(1:rnk))
        x(rnk+1:n) = zero

        ! Compute X(1:n) = Y**T * X(1:n)
        x = mult_rz(.true., n - rnk, rz, t2, x)
    else
        x = solve_triangular_system(.true., .false., .true., a(1:rnk,1:rnk), &
            qtb(1:rnk))
    end if

    ! Apply the pivoting: X(1:N) = P * X(1:N)
    do i = 1, n
        w(jpvt(i)) = x(i)
    end do
    x(1:n) = w(1:n)
end function

! ------------------------------------------------------------------------------
end module
