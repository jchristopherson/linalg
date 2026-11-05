! linalg_lu.f90
! SPDX-License-Identifier: MIT
module linalg_lu
    !! Provides LU factorization, factor extraction, and linear-system solvers.
    use iso_fortran_env
    use linalg_errors
    use linalg_sparse
    use linalg_basic
    use lapack
    use sparskit
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use linalg_types
    implicit none
    private
    public :: lu_factorization
    public :: sparse_lu_factorization
    public :: lu_factor
    public :: form_lu
    public :: solve_lu

    type, extends(matrix_factorization) :: lu_factorization
        !! A container for LU factorization of a matrix and associated 
        !! operations.
        real(real64), private, allocatable, dimension(:,:) :: a
            !! The matrix, and after a call to factor, the factored matrix.
        logical, private :: factored = .false.
            !! True if factored; else, false.
        integer, private, allocatable, dimension(:) :: pivot
            !! The pivot array.
        real(real64), private, allocatable, dimension(:) :: workspace
            !! A workspace array for matrix inversion.
    contains
        procedure, public :: is_factored => luf_is_factored
            !! Checks to see if the matrix has been factored.
        procedure, public :: size => luf_size
            !! Gets the size of the factored matrix.
        procedure, public :: factor => luf_factor
            !! Factors the matrix.
        procedure, public :: solve_multiple_rhs => luf_solve_multiple_rhs
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure, public :: solve_single_rhs => luf_solve_single_rhs
            !! Solves the factored system with a single right-hand-side vector.
        procedure, public :: solve_multiple_rhs_in_place => luf_solve_multiple_rhs_in_place
            !! Solves the factored system with multiple right-hand-side vectors.
        procedure, public :: solve_single_rhs_in_place => luf_solve_single_rhs_in_place
            !! Solves the factored system with a single right-hand-side vector.
        procedure, public :: get_factored_form => luf_get_factored_matrices
            !! Gets the factored matrices.
        procedure, public :: determinant => luf_det
            !! Gets the determinant of the original matrix.
        procedure, public :: inverse => luf_inverse
            !! Computes the inverse of the matrix.
        procedure, public :: update => luf_update
            !! Updates the factorization.
    end type

    type, extends(matrix_factorization) :: sparse_lu_factorization
        !! Owns a real CSR system and an ILUT preconditioner. SOLVE and
        !! SOLVE_IN_PLACE solve the original A * X = B using restarted GMRES;
        !! they do not merely apply the incomplete factors. This is an
        !! iterative solver, not a pivoted sparse direct LU factorization.
        !! FACTOR prepares the preconditioner explicitly, or solves prepare
        !! it lazily. UPDATE invalidates it. No dense matrix copy is formed.
        !! Errors return through optional FLAG arguments; if omitted, errors
        !! terminate with ERROR STOP. No inverse or determinant is provided.
        type(csr_matrix), private :: a
            !! Owned copy of the original coefficient matrix.
        type(msr_matrix), private :: lu
            !! Cached incomplete factors, with inverse upper diagonal.
        integer, private, allocatable :: ju(:)
            !! Pointers to the starts of the upper-factor rows in LU.
        logical, private :: initialized = .false.
            !! Whether UPDATE or the constructor has stored a matrix.
        logical, private :: factored = .false.
            !! Whether the cached preconditioner is ready for use.
        real(real64), private :: drop_tolerance = sqrt(epsilon(0.0d0))
            !! ILUT threshold, independent of the GMRES convergence tolerance.
        real(real64), private :: tolerance = sqrt(epsilon(0.0d0))
            !! Required relative true-residual tolerance for every RHS.
        integer, private :: max_iterations = 100
            !! GMRES iteration budget, checked at restart boundaries.
        integer, private :: restart = 50
            !! GMRES restart length, bounded by the matrix order and 50.
    contains
        procedure, public :: is_factored => sluf_is_factored
            !! Queries whether the cached ILUT preconditioner is ready.
        procedure, public :: size => sluf_size
            !! Queries dimensions of the owned CSR matrix.
        procedure, public :: factor => sluf_factor
            !! Validates settings and prepares the ILUT preconditioner once.
        procedure, public :: update => sluf_update
            !! Copies a new CSR matrix, resets controls, and invalidates factors.
        procedure, public :: solve_multiple_rhs => sluf_solve_multiple_rhs
            !! Allocates complete solutions for columns of a RHS matrix.
        procedure, public :: solve_single_rhs => sluf_solve_single_rhs
            !! Allocates the complete solution for one RHS vector.
        procedure, public :: solve_multiple_rhs_in_place => sluf_solve_multiple_rhs_in_place
            !! Replaces a RHS matrix only if all column solves converge.
        procedure, public :: solve_single_rhs_in_place => sluf_solve_single_rhs_in_place
            !! Replaces a RHS vector only if its solve converges.
        procedure, public :: solve_prepared => sluf_solve_prepared
            !! Read-only solve into caller storage; requires explicit FACTOR.
        procedure, public :: get_factored_form => sluf_get_factored_form
            !! Copies the incomplete preconditioner, not exact factors of A.
    end type

    interface sparse_lu_factorization
        !! Constructs an owned CSR system without performing factorization.
        module procedure :: sluf_init
    end interface

    interface lu_factorization
        module procedure :: luf_init
    end interface

    interface lu_factor
        module procedure :: lu_factor_dbl
        module procedure :: lu_factor_cmplx
        module procedure :: csr_lu_factor
    end interface

    interface form_lu
        module procedure :: form_lu_all
        module procedure :: form_lu_all_cmplx
        module procedure :: form_lu_only
        module procedure :: form_lu_only_cmplx
    end interface

    interface solve_lu
        module procedure :: solve_lu_mtx
        module procedure :: solve_lu_mtx_cmplx
        module procedure :: solve_lu_vec
        module procedure :: solve_lu_vec_cmplx
        module procedure :: csr_lu_solve
    end interface
contains
! ******************************************************************************
! SPARSE_LU_FACTORIZATION MEMBERS
! ------------------------------------------------------------------------------
function sluf_init(a, droptol, tol, maxits, im) result(rst)
    !! Constructs a sparse solver. The matrix and controls are stored without
    !! factorization; FACTOR or the first solve validates them.
    class(csr_matrix), intent(in) :: a
        !! Coefficient matrix to copy. Solving requires valid nonempty square
        !! CSR storage with unique column indices within each row.
    real(real64), intent(in), optional :: droptol
        !! Finite nonnegative ILUT drop threshold; default SQRT(EPSILON).
    real(real64), intent(in), optional :: tol
        !! Finite relative true-residual tolerance, at least EPSILON;
        !! default SQRT(EPSILON). Independent of DROPTOL.
    integer, intent(in), optional :: maxits
        !! Positive GMRES iteration budget, default 100. Checked at restart
        !! boundaries, so up to IM-1 additional iterations may be performed.
    integer, intent(in), optional :: im
        !! Restart length in [1, MIN(N, 50)]; default MIN(N, 50).
    type(sparse_lu_factorization) :: rst
        !! New unfactored solver owning its own copy of A and the controls.

    call rst%update(a, droptol, tol, maxits, im)
end function

subroutine sluf_update(this, a, droptol, tol, maxits, im)
    !! Replaces the owned system and discards all cached factors. Omitted
    !! controls reset to defaults, rather than retaining previous settings.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose matrix, controls, and factorization state are replaced.
    class(csr_matrix), intent(in) :: a
        !! New coefficient matrix, copied without changing the caller's data.
    real(real64), intent(in), optional :: droptol
        !! Finite nonnegative ILUT threshold; default SQRT(EPSILON).
    real(real64), intent(in), optional :: tol
        !! Finite relative residual tolerance at least EPSILON;
        !! default SQRT(EPSILON). Applies to the original system, not its factors.
    integer, intent(in), optional :: maxits
        !! Positive GMRES budget checked at restart boundaries; default 100.
    integer, intent(in), optional :: im
        !! Restart length in [1, MIN(N, 50)]; default MIN(N, 50).

    this%a = a
    this%lu = msr_matrix()
    if (allocated(this%ju)) deallocate(this%ju)
    this%initialized = .true.
    this%factored = .false.
    this%drop_tolerance = sqrt(epsilon(0.0d0))
    this%tolerance = sqrt(epsilon(0.0d0))
    this%max_iterations = 100
    this%restart = min(size(a, 1), 50)
    if (present(droptol)) this%drop_tolerance = droptol
    if (present(tol)) this%tolerance = tol
    if (present(maxits)) this%max_iterations = maxits
    if (present(im)) this%restart = im
end subroutine

pure function sluf_is_factored(this) result(rst)
    !! Queries preconditioner readiness without modifying the solver.
    class(sparse_lu_factorization), intent(in) :: this
        !! Solver to query.
    logical :: rst
        !! True after successful FACTOR, until UPDATE invalidates the factors.
    rst = this%factored
end function

pure function sluf_size(this, dim) result(rst)
    !! Queries the dimensions of the original CSR system.
    class(sparse_lu_factorization), intent(in) :: this
        !! Solver to query; need not be factored or initialized.
    integer, intent(in) :: dim
        !! Dimension: 1 for rows, 2 for columns; other values return zero.
    integer :: rst
        !! Requested extent, or zero when no matrix has been stored.
    rst = 0
    if (this%initialized) rst = size(this%a, dim)
end function

pure subroutine report_sparse_lu_status(info, flag)
    !! Applies the sparse solver's optional-status/error-stop policy.
    integer, intent(in) :: info
        !! Library status: LA_NO_ERROR for success, an LA_* error otherwise.
    integer, intent(out), optional :: flag
        !! Receives INFO when present. A nonzero INFO terminates execution
        !! when FLAG is absent; successful calls never terminate.
    if (present(flag)) then
        flag = info
    else if (info /= LA_NO_ERROR) then
        error stop info
    end if
end subroutine

subroutine sluf_factor(this, flag)
    !! Validates the matrix and controls and prepares its ILUT preconditioner.
    !! Repeated successful calls reuse the existing factors. ILUT does not
    !! pivot and may regularize zero pivots; readiness alone does not guarantee
    !! nonsingularity of A or convergence of a subsequent GMRES solve.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose incomplete factors are created and cached.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR without a stored
        !! matrix; LA_INVALID_INPUT_ERROR for invalid controls; or CSR/ILUT
        !! validation errors from [[lu_factor]]. Omit to use ERROR STOP.
    integer :: info

    info = LA_NO_ERROR
    if (.not.this%initialized) then
        info = LA_NOT_INITIALIZED_ERROR
    else if (.not.this%factored) then
        if (.not.ieee_is_finite(this%tolerance)) then
            info = LA_INVALID_INPUT_ERROR
        else if (this%tolerance < epsilon(0.0d0) .or. &
            this%max_iterations < 1 .or. this%restart < 1 .or. &
            this%restart > min(this%size(1), 50)) then
            info = LA_INVALID_INPUT_ERROR
        else
            if (allocated(this%ju)) deallocate(this%ju)
            allocate(this%ju(this%size(1)))
            call csr_lu_factor(this%a, this%lu, this%ju, &
                this%drop_tolerance, info)
            if (info == LA_NO_ERROR) then
                if (.not.all(ieee_is_finite(this%lu%values))) &
                    info = LA_SINGULAR_MATRIX_ERROR
            end if
        end if
        this%factored = info == LA_NO_ERROR
    end if
    call report_sparse_lu_status(info, flag)
end subroutine

subroutine sluf_prepare_solve(this, rows, flag)
    !! Validates RHS row count and lazily prepares the preconditioner.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose factors may be prepared if not yet available.
    integer, intent(in) :: rows
        !! Number of elements in a RHS vector or rows of a RHS matrix.
    integer, intent(out) :: flag
        !! LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR,
        !! LA_MATRIX_NOT_SQUARE_ERROR, LA_MATRIX_DIMENSION_ERROR, or FACTOR's
        !! status otherwise. This internal routine never intentionally stops.
    flag = LA_NO_ERROR
    if (.not.this%initialized) then
        flag = LA_NOT_INITIALIZED_ERROR
    else if (this%size(1) /= this%size(2)) then
        flag = LA_MATRIX_NOT_SQUARE_ERROR
    else if (rows /= this%size(1)) then
        flag = LA_MATRIX_DIMENSION_ERROR
    else
        call this%factor(flag)
    end if
end subroutine

subroutine sluf_solve_prepared(this, b, x, flag)
    !! Solves the original CSR system using already prepared factors. No
    !! factorization, caching, or mutation of THIS occurs. Workspace is local
    !! to this call; external-library reentrancy remains a separate requirement.
    class(sparse_lu_factorization), intent(in) :: this
        !! Read-only solver after successful FACTOR and before any UPDATE.
    real(real64), intent(in) :: b(:)
        !! N-element RHS vector, not modified; must not overlap X.
    real(real64), intent(out) :: x(:)
        !! Caller-owned N-element solution. Set to zero on any reported error.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR if NORM2(B-A*X) <= TOL*NORM2(B);
        !! LA_INVALID_OPERATION_ERROR if not factored;
        !! LA_MATRIX_DIMENSION_ERROR for incorrect vector lengths;
        !! LA_INVALID_INPUT_ERROR for nonfinite B; or LA_CONVERGENCE_ERROR
        !! on iteration/residual failure. Omit to use ERROR STOP.
    real(real64), allocatable :: solution(:)
    integer :: info

    x = 0.0d0
    info = LA_NO_ERROR
    if (.not.this%factored) then
        info = LA_INVALID_OPERATION_ERROR
    else if (size(b) /= this%size(1) .or. size(x) /= this%size(2)) then
        info = LA_MATRIX_DIMENSION_ERROR
    else
        call pgmres_solve(this%a, this%lu, this%ju, b, solution, &
            im = this%restart, tol = this%tolerance, &
            maxits = this%max_iterations, flag = info)
        if (info == LA_NO_ERROR) x = solution
    end if
    call report_sparse_lu_status(info, flag)
end subroutine

function sluf_solve_single_rhs(this, b, flag) result(rst)
    !! Allocates a complete GMRES solution for one RHS, preparing factors lazily.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose preconditioner is prepared if necessary.
    real(real64), intent(in) :: b(:)
        !! RHS with N elements, not modified.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on convergence; otherwise the preparation or prepared
        !! solve error described by SOLVE_SINGLE_RHS_IN_PLACE. Omit to stop.
    real(real64), allocatable :: rst(:)
        !! N-element solution on success; an unchanged copy of B on error.
        !! Always has the input shape, even if an incorrect RHS size is given.
    rst = b
    call this%solve_single_rhs_in_place(rst, flag)
end function

function sluf_solve_multiple_rhs(this, b, flag) result(rst)
    !! Solves each RHS column with GMRES, reusing one cached preconditioner.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose preconditioner is prepared if necessary.
    real(real64), intent(in) :: b(:,:)
        !! N-by-NRHS matrix with a RHS in each column, not modified.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR if every column converges; otherwise the first
        !! preparation or solve error. Omit to use ERROR STOP.
    real(real64), allocatable :: rst(:,:)
        !! N-by-NRHS solutions on success; an unchanged copy of B on error.
        !! Always has the input shape; no partial column results are returned.
    rst = b
    call this%solve_multiple_rhs_in_place(rst, flag)
end function

subroutine sluf_solve_single_rhs_in_place(this, b, flag)
    !! Replaces a RHS with the solution only after successful preparation and
    !! convergence against the true residual of the original system.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose incomplete factors are prepared lazily and then reused.
    real(real64), intent(inout) :: b(:)
        !! N-element RHS on entry; complete solution on success. Unchanged
        !! on any reported error, including nonconvergence.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on success; LA_NOT_INITIALIZED_ERROR,
        !! LA_MATRIX_NOT_SQUARE_ERROR, LA_MATRIX_DIMENSION_ERROR,
        !! LA_MATRIX_FORMAT_ERROR, LA_INVALID_INPUT_ERROR,
        !! LA_SINGULAR_MATRIX_ERROR, LA_INVALID_OPERATION_ERROR, or
        !! LA_CONVERGENCE_ERROR otherwise. Omit to use ERROR STOP.
    real(real64), allocatable :: solution(:)
    integer :: info

    call sluf_prepare_solve(this, size(b), info)
    if (info == LA_NO_ERROR) then
        allocate(solution(size(b)))
        call this%solve_prepared(b, solution, info)
        if (info == LA_NO_ERROR) b = solution
    end if
    call report_sparse_lu_status(info, flag)
end subroutine

subroutine sluf_solve_multiple_rhs_in_place(this, b, flag)
    !! Replaces all RHS columns only if every independent GMRES solve succeeds.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose incomplete factors are prepared lazily and reused.
    real(real64), intent(inout) :: b(:,:)
        !! N-by-NRHS right-hand sides on entry; complete solutions on success.
        !! The entire matrix is unchanged if any column fails.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on success; otherwise the first preparation or solve
        !! status listed for SOLVE_SINGLE_RHS_IN_PLACE. Omit to use ERROR STOP.
    real(real64), allocatable :: solution(:,:)
    integer :: info, column

    call sluf_prepare_solve(this, size(b, 1), info)
    if (info == LA_NO_ERROR) then
        allocate(solution(size(b, 1), size(b, 2)))
        do column = 1, size(b, 2)
            call this%solve_prepared(b(:,column), solution(:,column), info)
            if (info /= LA_NO_ERROR) exit
        end do
        if (info == LA_NO_ERROR) b = solution
    end if
    call report_sparse_lu_status(info, flag)
end subroutine

subroutine sluf_get_factored_form(this, lu, ju, flag)
    !! Copies the ILUT preconditioner, preparing it if needed. These are not
    !! exact factors of A; using LUSOL alone is not the complete system solve.
    class(sparse_lu_factorization), intent(inout) :: this
        !! Solver whose preconditioner may be prepared and then copied.
    type(msr_matrix), intent(out) :: lu
        !! Owned copy of the incomplete MSR factors, including inverse U
        !! diagonal. Reset to an empty matrix on error.
    integer, intent(out), allocatable :: ju(:)
        !! Owned copy of the N upper-row pointers; unallocated on error.
    integer, intent(out), optional :: flag
        !! FACTOR's status; LA_NO_ERROR on success. Omit to use ERROR STOP.
    integer :: info

    call this%factor(info)
    if (info == LA_NO_ERROR) then
        lu = this%lu
        ju = this%ju
    end if
    call report_sparse_lu_status(info, flag)
end subroutine

! ******************************************************************************
! LU_FACTORIZATION MEMBERS
! ------------------------------------------------------------------------------
function luf_init(a) result(rst)
    !! Initializes a new [[lu_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to factor.  A copy of this matrix is
        !! made and stored by this object; therefore, this
        !! matrix is not altered.
    type(lu_factorization) :: rst
        !! The new [[lu_factorization]] object.

    call rst%update(a)
end function

! ------------------------------------------------------------------------------
subroutine luf_update(this, a)
    !! Updates the [[lu_factorization]] object with a new matrix.  The matrix
    !! is only stored, not factored; therefore a call to factor is still 
    !! required.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix to factor.  A copy of this matrix is made and stored by
        !! the object; therefore, this matrix is not altered.  Also, any
        !! existing factorization is destroyed by this routine.

    integer :: m, n, mn

    m = size(a, 1)
    n = size(a, 2)
    mn = min(m, n)

    if (allocated(this%a)) deallocate(this%a)
    if (allocated(this%pivot)) deallocate(this%pivot)
    allocate(this%a(m, n), source = a)
    allocate(this%pivot(mn))
    this%factored = .false.
end subroutine

! ------------------------------------------------------------------------------
pure function luf_is_factored(this) result(rst)
    !! Gets a value querying whether the matrix has been factored.
    class(lu_factorization), intent(in) :: this
        !! The [[lu_factorization]] object.
    logical :: rst
        !! True if the matrix has been factored; else, false.
    rst = this%factored
end function

! ------------------------------------------------------------------------------
pure function luf_size(this, dim) result(rst)
    !! Retrieves the size of the requested matrix dimension.
    class(lu_factorization), intent(in) :: this
        !! The [[lu_factorization]] object.
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
subroutine luf_factor(this, flag)
    !! Factors the matrix.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: m, n, info

    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR

    ! If not allocated, return with an error
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if

    ! If factored, do nothing
    if (this%is_factored()) return

    ! Factor the matrix
    m = this%size(1)
    n = this%size(2)
    call DGETRF(m, n, this%a, m, this%pivot, info)
    if (info > 0) then
        if (present(flag)) flag = LA_SINGULAR_MATRIX_ERROR
    else
        this%factored = .true.
    end if
end subroutine

! ------------------------------------------------------------------------------
function luf_solve_multiple_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining multiple right-hand-side vectors.
    !! The encapsulated matrix must be square; else, an 
    !! LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(in), dimension(:,:) :: b
        !! The right-hand-side matrix.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The solution matrix.

    ! Local Variables
    integer :: n, nrhs, info

    ! Initialization
    n = this%size(1)
    nrhs = size(b, 2)
    if (present(flag)) flag = LA_NO_ERROR
    
    ! Solve
    allocate(rst(n, nrhs), source = b)
    call this%solve_in_place(rst, flag = flag)
end function

! ------------------------------------------------------------------------------
function luf_solve_single_rhs(this, b, flag) result(rst)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining a single right-hand-side vector.
    !! The encapsulated matrix must be square; else, an 
    !! LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(in), dimension(:) :: b
        !! The right-hand-side vector.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64), allocatable, dimension(:) :: rst
        !! The solution vector.

    ! Local Variables
    integer :: n, info

    ! Initialization
    n = this%size(1)
    if (present(flag)) flag = LA_NO_ERROR
    
    ! Solve
    allocate(rst(n), source = b)
    call this%solve_in_place(rst, flag = flag)
end function

! ------------------------------------------------------------------------------
subroutine luf_solve_multiple_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining multiple right-hand-side vectors.
    !! The encapsulated matrix must be square; else, an 
    !! LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(inout), dimension(:,:) :: b
        !! On input, the right-hand-side matrix \(B\).  On output, the solution
        !! matrix \(X\).
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: n, nrhs, info

    ! Validation
    n = this%size(1)
    nrhs = size(b, 2)
    if (present(flag)) flag = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (this%size(2) /= n) then
        if (present(flag)) flag = LA_MATRIX_NOT_SQUARE_ERROR
        return
    end if
    if (size(b, 1) /= n) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if

    ! If the matrix is not factored, be sure to factor
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Solve
    call DGETRS("N", n, nrhs, this%a, n, this%pivot, b, n, info)
end subroutine

! ------------------------------------------------------------------------------
subroutine luf_solve_single_rhs_in_place(this, b, flag)
    !! Solves the linear system \(A X = B \) where \(B\)
    !! is a matrix defining a single right-hand-side vector.
    !! The encapsulated matrix must be square; else, an 
    !! LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(inout), dimension(:) :: b
        !! The right-hand-side vector.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: n, info

    ! Validation
    n = this%size(1)
    if (present(flag)) flag = LA_NO_ERROR
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (this%size(2) /= n) then
        if (present(flag)) flag = LA_MATRIX_NOT_SQUARE_ERROR
        return
    end if
    if (size(b) /= n) then
        if (present(flag)) flag = LA_MATRIX_DIMENSION_ERROR
        return
    end if

    ! If the matrix is not factored, be sure to factor
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Solve
    call DGETRS("N", n, 1, this%a, n, this%pivot, b, n, info)
end subroutine

! ------------------------------------------------------------------------------
subroutine luf_get_factored_matrices(this, l, u, p, flag)
    !! Extracts the factored matrices \(L\), \(U\), and the row pivot 
    !! matrix \(P\), where \(P A = L U\), and \(A\) is M-by-N.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(out), allocatable, dimension(:,:) :: l
        !! The M-by-M \(L\) matrix.
    real(real64), intent(out), allocatable, dimension(:,:) :: u
        !! The M-by-N \(U\) matrix.
    real(real64), intent(out), allocatable, dimension(:,:) :: p
        !! The M-by-M \(P\) matrix.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: j, jp, m, n, mn, info

    ! Initialization
    m = this%size(1)
    n = this%size(2)
    mn = min(m, n)
    if (present(flag)) flag = LA_NO_ERROR

    ! If not allocated, return with an error
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if

    ! Is the matrix factored?
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= LA_NO_ERROR) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Form L
    allocate(l(m, m), u(m, n), p(m, m), source = 0.0d0)
    do j = 1, m
        l(j,j) = 1.0d0
        l(j+1:,j) = this%a(j+1:,j)
    end do

    ! Form U
    do j = 1, n
        u(1:j,j) = this%a(1:j,j)
    end do

    ! Form P
    call DLASET('A', mn, mn, 0.0d0, 1.0d0, p, m)
    do j = 1, mn
        jp = this%pivot(j)
        if (j /= jp) call swap(p(j,1:mn), p(jp,1:mn))
    end do
end subroutine

! ------------------------------------------------------------------------------
function luf_det(this, flag) result(rst)
    !! Returns the determinant of the matrix.  The matrix must be square; else,
    !! an LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    integer, intent(out), optional :: flag
        !! An error flag.
    real(real64) :: rst
        !! The determinant of the matrix.

    ! Local Variables
    integer :: i, ep, n, info
    real(real64) :: temp

    ! Validation
    n = this%size(1)
    if (present(flag)) flag = LA_NO_ERROR
    rst = ieee_value(0.0d0, ieee_quiet_nan)
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
        return
    end if
    if (this%size(2) /= n) then
        if (present(flag)) flag = LA_MATRIX_NOT_SQUARE_ERROR
        return
    end if

    ! Ensure the matrix is factored
    if (.not.this%is_factored()) then
        call this%factor(flag = info)
        if (info /= 0) then
            if (present(flag)) flag = info
            return
        end if
    end if

    ! Compute the determinant of the factored form
    temp = 1.0d0
    ep = 0
    do i = 1, n
        if (this%pivot(i) /= i) temp = -temp

        temp = this%a(i, i) * temp
        if (temp == 0.0d0) then
            rst = 0.0d0
            exit
        end if

        do while (abs(temp) < 1.0d0)
            temp = 1.0d1 * temp
            ep = ep - 1
        end do

        do while (abs(temp) > 1.0d1)
            temp = 1.0d-1 * temp
            ep = ep + 1
        end do
    end do
    rst = temp * (1.0d1)**ep
end function

! ------------------------------------------------------------------------------
subroutine luf_inverse(this, ainv, flag)
    !! Computes the inverse of the matrix.  The matrix must be square; else,
    !! an LA_MATRIX_NOT_SQUARE_ERROR is returned.
    class(lu_factorization), intent(inout) :: this
        !! The [[lu_factorization]] object.
    real(real64), intent(out), allocatable, dimension(:,:) :: ainv
        !! The inverted matrix.
    integer, intent(out), optional :: flag
        !! An error flag.

    ! Local Variables
    integer :: n, lwork, info, idummy(1)
    real(real64) :: dummy(1), w(1)

    ! Validation
    n = this%size(1)
    if (present(flag)) flag = LA_NO_ERROR
    if (this%size(2) /= n) then
        if (present(flag)) flag = LA_MATRIX_NOT_SQUARE_ERROR
        return
    end if
    if (.not.allocated(this%a)) then
        if (present(flag)) flag = LA_NOT_INITIALIZED_ERROR
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

    ! Determine workspace requirements
    call DGETRI(n, dummy, n, idummy, w, -1, info)
    lwork = int(w(1))
    if (.not.allocated(this%workspace)) then
        allocate(this%workspace(lwork))
    else
        if (size(this%workspace) < lwork) then
            deallocate(this%workspace)
            allocate(this%workspace(lwork))
        end if
    end if

    ! Compute the inverse
    allocate(ainv(n, n), source = this%a)
    call DGETRI(n, ainv, n, this%pivot, this%workspace, lwork, info)
    if (info > 0) then
        if (present(flag)) flag = LA_SINGULAR_MATRIX_ERROR
    end if
end subroutine

! ******************************************************************************
! PROCEDURAL INTERFACE
! ------------------------------------------------------------------------------
pure subroutine lu_factor_dbl(a, ipvt, lu, l, u, p)
    !! Computes the LU factorization of an M-by-N matrix.  In the event of a 
    !! singular matrix, the output matrices are populated with NaN's.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    integer, intent(out), allocatable, optional, target, dimension(:) :: ipvt
        !! An N-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lu
        !! The N-by-N factored matrix in the form [\(L\)\\ \(U\)] where the unit
        !! diagonal elements of \(L\) are not stored.
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The N-by-N lower triangular matrix \(L\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: u
        !! The N-by-N upper triangular matrix \(U\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N row permutation matrix.

    ! Local Variables
    logical :: buildlu
    integer :: n, flag
    integer, allocatable, target, dimension(:) :: ic
    integer, pointer, dimension(:) :: iptr
    real(real64) :: nan
    real(real64), allocatable, target, dimension(:,:) :: ac, lc, uc, pc
    real(real64), pointer, dimension(:,:) :: aptr, lptr, uptr, pptr

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) error stop 1
    if (present(ipvt)) then
        allocate(ipvt(n))
        iptr => ipvt
    else
        allocate(ic(n))
        iptr => ic
    end if
    if (present(lu)) then
        allocate(lu(n, n), source = a)
        aptr => lu
    else
        allocate(ac(n, n), source = a)
        aptr => ac
    end if
    buildlu = present(l) .or. present(u) .or. present(p)

    ! Compute the LU factorization by calling the LAPACK routine DGETRF
    call DGETRF(n, n, aptr, n, iptr, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        if (present(lu)) lu = nan
        if (present(l)) allocate(l(n,n), source = nan)
        if (present(u)) allocate(u(n,n), source = nan)
        if (present(p)) allocate(p(n,n), source = nan)
        if (present(ipvt)) ipvt = 0
        return
    end if

    ! Build L & U?
    if (buildlu) then
        ! L
        if (present(l)) then
            allocate(l(n, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(n, n), source = aptr)
                lptr => lc
            end if
        end if

        ! U
        if (present(u)) then
            allocate(u(n, n))
            uptr => u
        else
            allocate(uc(n, n))
            uptr => uc
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Build the matrices
        call form_lu(lptr, iptr, uptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine lu_factor_cmplx(a, ipvt, lu, l, u, p)
    !! Computes the LU factorization of an M-by-N matrix.  In the event of a 
    !! singular matrix, the output matrices are populated with NaN's.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N matrix to factor.
    integer, intent(out), allocatable, optional, target, dimension(:) :: ipvt
        !! An N-element array used to track row-pivot operations.  The 
        !! array stored pivot information such that row I is interchanged with 
        !! row IPVT(I).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: lu
        !! The N-by-N factored matrix in the form [\(L\)\\ \(U\)] where the unit
        !! diagonal elements of \(L\) are not stored.
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: l
        !! The N-by-N lower triangular matrix \(L\).
    complex(real64), intent(out), allocatable, optional, target, dimension(:,:) :: u
        !! The N-by-N upper triangular matrix \(U\).
    real(real64), intent(out), allocatable, optional, target, dimension(:,:) :: p
        !! The N-by-N row permutation matrix.

    ! Local Variables
    logical :: buildlu
    integer :: n, flag
    integer, allocatable, target, dimension(:) :: ic
    integer, pointer, dimension(:) :: iptr
    real(real64), allocatable, target, dimension(:,:) :: pc
    real(real64), pointer, dimension(:,:) :: pptr
    real(real64) :: nan
    complex(real64) :: cnan
    complex(real64), allocatable, target, dimension(:,:) :: ac, lc, uc
    complex(real64), pointer, dimension(:,:) :: aptr, lptr, uptr

    ! Initialization
    n = size(a, 1)
    if (size(a, 2) /= n) error stop 1
    if (present(ipvt)) then
        allocate(ipvt(n))
        iptr => ipvt
    else
        allocate(ic(n))
        iptr => ic
    end if
    if (present(lu)) then
        allocate(lu(n, n), source = a)
        aptr => lu
    else
        allocate(ac(n, n), source = a)
        aptr => ac
    end if
    buildlu = present(l) .or. present(u) .or. present(p)

    ! Compute the LU factorization by calling the LAPACK routine ZGETRF
    call ZGETRF(n, n, aptr, n, iptr, flag)

    ! If flag > 0, the matrix is singular.  Notice, flag should not be
    ! able to be < 0 as we've already verrified inputs prior to making the
    ! call to LAPACK
    if (flag > 0) then
        ! WARNING: Singular matrix
        nan = ieee_value(nan, ieee_quiet_nan)
        cnan = cmplx(nan, nan)
        if (present(lu)) lu = cnan
        if (present(l)) allocate(l(n,n), source = cnan)
        if (present(u)) allocate(u(n,n), source = cnan)
        if (present(p)) allocate(p(n,n), source = nan)
        if (present(ipvt)) ipvt = 0
        return
    end if

    ! Build L & U?
    if (buildlu) then
        ! L
        if (present(l)) then
            allocate(l(n, n), source = aptr)
            lptr => l
        else
            if (allocated(ac)) then
                lptr => ac
            else
                allocate(lc(n, n), source = aptr)
                lptr => lc
            end if
        end if

        ! U
        if (present(u)) then
            allocate(u(n, n))
            uptr => u
        else
            allocate(uc(n, n))
            uptr => uc
        end if

        ! P
        if (present(p)) then
            allocate(p(n, n))
            pptr => p
        else
            allocate(pc(n, n))
            pptr => pc
        end if

        ! Build the matrices
        call form_lu(lptr, iptr, uptr, pptr)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function validate_csr_factor_input(a) result(flag)
    !! Validates CSR storage before passing raw arrays to SPARSKIT ILUT.
    class(csr_matrix), intent(in) :: a
        !! Matrix to inspect without modification. Rows may be unsorted, but
        !! column indices within each row must be unique and in range.
    integer :: flag
        !! LA_NO_ERROR for valid nonempty square storage; otherwise
        !! LA_NOT_INITIALIZED_ERROR, LA_MATRIX_NOT_SQUARE_ERROR,
        !! LA_INVALID_INPUT_ERROR, LA_MATRIX_FORMAT_ERROR, or
        !! LA_SINGULAR_MATRIX_ERROR for a numerically zero row.
    integer :: order, row, first, last

    flag = validate_csr_storage(a)
    if (flag /= LA_NO_ERROR) return
    order = size(a, 1)
    flag = LA_MATRIX_NOT_SQUARE_ERROR
    if (size(a, 2) /= order) return
    do row = 1, order
        first = a%row_indices(row)
        last = a%row_indices(row+1) - 1
        if (all(a%values(first:last) == 0.0d0)) then
            flag = LA_SINGULAR_MATRIX_ERROR
            return
        end if
    end do
    flag = LA_NO_ERROR
end function

pure subroutine csr_lu_factor(a, lu, ju, droptol, flag)
    !! Builds incomplete LU factors using SPARSKIT ILUT. These approximate
    !! factors are intended for preconditioning, not an exact direct solve.
    !! A detected zero pivot retains the legacy NaN-valued factor output.
    class(csr_matrix), intent(in) :: a
        !! Nonempty square CSR matrix, copied into neither dense storage nor A.
        !! Storage must be valid with unique column indices within each row.
    type(msr_matrix), intent(out) :: lu
        !! Combined incomplete factors in MSR format, with inverse U diagonal.
        !! Not usable when FLAG reports an error.
    integer, intent(out), dimension(:) :: ju
        !! Caller-owned N-element array receiving upper-factor row pointers.
    real(real64), intent(in), optional :: droptol
        !! Finite nonnegative ILUT drop tolerance; default SQRT(EPSILON).
        !! Zero disables threshold dropping, but the per-row fill cap remains.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on success. Validation returns the statuses documented
        !! by VALIDATE_CSR_FACTOR_INPUT; incorrect JU size returns
        !! LA_MATRIX_DIMENSION_ERROR. ILUT breakdown returns
        !! LA_SINGULAR_MATRIX_ERROR, LA_MATRIX_FORMAT_ERROR, or
        !! LA_INVALID_OPERATION_ERROR. When omitted, validation and format
        !! errors use ERROR STOP; a zero pivot retains the legacy NaN output.

    ! Local Variables
    integer :: i, m, n, nn, nnz, lfil, iwk, ierr
    integer, allocatable, dimension(:) :: jlu, jw
    real(real64), allocatable, dimension(:) :: alu, w
    real(real64) :: dt, nan
    
    ! Initialization
    if (present(flag)) flag = LA_NO_ERROR
    if (present(droptol)) then
        dt = droptol
    else
        dt = sqrt(epsilon(dt))
    end if
    m = size(a, 1)
    n = size(a, 2)
    nnz = nonzero_count(a)

    ! Input Check
    ierr = validate_csr_factor_input(a)
    if (ierr == LA_NO_ERROR .and. size(ju) /= m) &
        ierr = LA_MATRIX_DIMENSION_ERROR
    if (.not.ieee_is_finite(dt)) then
        ierr = LA_INVALID_INPUT_ERROR
    else if (dt < 0.0d0) then
        ierr = LA_INVALID_INPUT_ERROR
    end if
    if (ierr /= LA_NO_ERROR) then
        if (present(flag)) then
            flag = ierr
            return
        end if
        error stop ierr
    end if

    ! Parameter Determination
    lfil = 1
    do i = 1, m
        lfil = max(lfil, a%row_indices(i+1) - a%row_indices(i))
    end do
    iwk = max(n + 1, lfil * m + 1, nnz)

    ! Local Memory Allocation
    allocate(alu(iwk), w(n+1), jlu(iwk), jw(2 * n))

    ! Factorization
    do
        ! Factor the matrix
        call ilut(n, a%values, a%column_indices, a%row_indices, lfil, dt, &
            alu, jlu, ju, iwk, w, jw, ierr)

        ! Check the error flag
        if (ierr == 0) then
            ! Success
            exit
        else if (ierr > 0) then
            ! WARNING: Singular matrix - a zero pivot was found at step IERR.
            ! Return a structurally valid, but NaN-valued, diagonal factor so
            ! that the condition propagates into any subsequent solution.
            nan = ieee_value(nan, ieee_quiet_nan)
            lu%m = m
            lu%n = n
            lu%nnz = 0
            allocate(lu%values(m + 1), source = nan)
            allocate(lu%indices(m + 1), source = m + 2)
            ju = m + 2
            if (present(flag)) flag = LA_SINGULAR_MATRIX_ERROR
            return
        else if (ierr == -1) then
            ! The input matrix is not formatted correctly
            if (present(flag)) then
                flag = LA_MATRIX_FORMAT_ERROR
                return
            end if
            error stop LA_MATRIX_FORMAT_ERROR
        else if (ierr == -2 .or. ierr == -3) then
            ! ALU and JLU are too small - try something larger
            ! This is the main reason for the loop - to offload worrying about
            ! workspace size from the user
            iwk = iwk + max(m + n, iwk / 2)
            deallocate(alu)
            deallocate(jlu)
            allocate(alu(iwk), jlu(iwk))
        else if (ierr == -4) then
            ! Illegal value for LFIL - reset and try again
            lfil = n
        else if (ierr == -5) then
            ! Zero row encountered
            if (present(flag)) then
                flag = LA_SINGULAR_MATRIX_ERROR
                return
            end if
            error stop LA_MATRIX_FORMAT_ERROR
        else
            ! We should never get here, but just in case
            if (present(flag)) then
                flag = LA_INVALID_OPERATION_ERROR
                return
            end if
            error stop LA_INVALID_OPERATION_ERROR
        end if
    end do

    ! Determine the actual number of non-zero elements
    nnz = jlu(m+1) - 1

    ! Copy the contents to the output arrays
    lu%m = m
    lu%n = n
    lu%nnz = nnz
    nn = nnz
    allocate(lu%values(nn), source = alu(:nn))
    allocate(lu%indices(nn), source = jlu(:nn))
    ! The unused MSR separator must be defined for whole-array finite checks.
    lu%values(n + 1) = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lu_all(lu, ipvt, u, p)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    real(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the row permutation matrix will be written.

    ! Local Variables
    integer :: j, jp, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 3
    end if
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 4
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
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lu_all_cmplx(lu, ipvt, u, p)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.
    real(real64), intent(out), dimension(:,:) :: p
        !! An N-by-N matrix where the row permutation matrix will be written.

    ! Local Variables
    integer :: j, jp, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0
    complex(real64), parameter :: c_zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: c_one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 3
    end if
    if (size(p, 1) /= n .or. size(p, 2) /= n) then
        error stop 4
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
end subroutine
! ------------------------------------------------------------------------------
pure subroutine form_lu_only(lu, u)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    real(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    real(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.

    ! Local Variables
    integer :: j, n

    ! Parameters
    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: one = 1.0d0

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 2
    end if

    ! Process
    do j = 1, n
        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = zero

        if (j > 1) lu(1:j-1,j) = zero
        lu(j,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine form_lu_only_cmplx(lu, u)
    !! Extracts the L and U matrices from the condensed [L\\U] storage format 
    !! used by the lu_factor.
    complex(real64), intent(inout), dimension(:,:) :: lu
        !! On input, the N-by-N matrix as output by lu_factor.  On output, the 
        !! N-by-N lower triangular matrix L.
    complex(real64), intent(out), dimension(:,:) :: u
        !! An N-by-N matrix where the U matrix will be written.

    ! Local Variables
    integer :: j, n

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    ! Initialization
    n = size(lu, 1)

    ! Input Check
    if (size(lu, 2) /= n) then
        error stop 1
    end if
    if (size(u, 1) /= n .or. size(u, 2) /= n) then
        error stop 2
    end if

    ! Process
    do j = 1, n
        ! Build L and U
        u(1:j,j) = lu(1:j,j)
        u(j+1:n,j) = zero

        if (j > 1) lu(1:j-1,j) = zero
        lu(j,j) = one
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function solve_lu_mtx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS right-hand-side matrix.
    real(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Local Variables
    integer :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b, 1) /= n) then
        error stop 3
    end if

    ! Call DGETRS
    allocate(x(n, nrhs), source = b)
    call DGETRS("N", n, nrhs, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_mtx_cmplx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(in), dimension(:,:) :: b
        !! The N-by-NRHS right-hand-side matrix.
    complex(real64), allocatable, dimension(:,:) :: x
        !! The N-by-NRHS solution matrix.

    ! Local Variables
    integer :: n, nrhs, flag

    ! Initialization
    n = size(a, 1)
    nrhs = size(b, 2)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b, 1) /= n) then
        error stop 3
    end if

    ! Call ZGETRS
    allocate(x(n, nrhs), source = b)
    call ZGETRS("N", n, nrhs, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_vec(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    real(real64), intent(in), dimension(:) :: b
        !! The N-element right-hand-side array.
    real(real64), allocatable, dimension(:) :: x
        !! The N-element solution array.

    ! Local Variables
    integer :: n, flag

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b) /= n) then
        error stop 3
    end if

    ! Call DGETRS
    allocate(x(n), source = b)
    call DGETRS("N", n, 1, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function solve_lu_vec_cmplx(a, ipvt, b) result(x)
    !! Solves a system of LU-factored equations.
    complex(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU factored matrix as output by lu_factor.
    integer, intent(in), dimension(:) :: ipvt
        !! The N-element pivot array as output by lu_factor.
    complex(real64), intent(in), dimension(:) :: b
        !! The N-element right-hand-side array.
    complex(real64), allocatable, dimension(:) :: x
        !! The N-element solution array.

    ! Local Variables
    integer :: n, flag

    ! Initialization
    n = size(a, 1)

    ! Input Check
    if (size(a, 2) /= n) then
        error stop 1
    end if
    if (size(ipvt) /= n) then
        error stop 2
    end if
    if (size(b) /= n) then
        error stop 3
    end if

    ! Call ZGETRS
    allocate(x(n), source = b)
    call ZGETRS("N", n, 1, a, n, ipvt, x, n, flag)
end function

! ------------------------------------------------------------------------------
pure function csr_lu_solve(lu, ju, b) result(x)
    !! Applies the combined triangular factors with SPARSKIT LUSOL. If the
    !! factors came from ILUT, this solves the incomplete-factor system only;
    !! use PGMRES_SOLVE or SPARSE_LU_FACTORIZATION for the original CSR system.
    !! Malformed factors use ERROR STOP; singular/nonfinite factors yield NaNs.
    class(msr_matrix), intent(in) :: lu
        !! Nonempty square MSR factors, with inverse U pivots on the diagonal.
    integer, intent(in), dimension(:) :: ju
        !! N pointers splitting each row into strict lower and upper entries.
    real(real64), intent(in), dimension(:) :: b
        !! Finite N-element right-hand side, not modified.
    real(real64), allocatable, dimension(:) :: x
        !! N-element solution of the factor system. Singular/nonfinite factors
        !! yield NaNs, preserving legacy propagation. No GMRES refinement occurs.

    ! Local Variables
    integer :: m, n, info
    
    ! Initialization
    m = size(lu, 1)
    n = size(lu, 2)

    ! Input Check
    if (m /= n) then
        error stop LA_MATRIX_NOT_SQUARE_ERROR
    end if
    if (size(ju) /= m) then
        error stop LA_MATRIX_DIMENSION_ERROR
    end if
    if (size(b) /= m) then
        error stop LA_MATRIX_DIMENSION_ERROR
    end if
    if (.not.all(ieee_is_finite(b))) error stop LA_INVALID_INPUT_ERROR
    info = validate_msr_lu(lu, ju)
    if (info == LA_SINGULAR_MATRIX_ERROR) then
        allocate(x(m), source = ieee_value(0.0d0, ieee_quiet_nan))
        return
    else if (info /= LA_NO_ERROR) then
        error stop info
    end if

    ! Process
    allocate(x(m), source = b)
    call lusol(m, b, x, lu%values, lu%indices, ju)
end function

! ------------------------------------------------------------------------------
end module
