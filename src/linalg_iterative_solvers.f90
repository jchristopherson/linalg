module linalg_iterative_solvers
    !! Storage-independent iterative solves using SPARSKIT reverse communication.
    !! The bundled GMRES backend uses SAVE state. Entire solves must be
    !! serialized; nested solves using this backend from callbacks are forbidden.
    use iso_fortran_env, only : int32, int64, real64
    use ieee_arithmetic, only : ieee_is_finite
    use linalg_linear_operators, only : linear_operator
    use linalg_sparse, only : csr_matrix, msr_matrix, size, matmul, &
        validate_csr_storage, validate_msr_lu
    use linalg_errors
    use sparskit, only : gmres, lusol
    implicit none
    private
    public :: gmres_solve, iterative_solver_info
    public :: csr_operator, ilu_preconditioner

    type, extends(linear_operator) :: csr_operator
        !! Read-only operator owning a validated CSR copy. Construction copies
        !! storage once; subsequent applications neither copy nor densify A.
        type(csr_matrix), private :: matrix
            !! Owned coefficients; independent of the source matrix's lifetime.
    contains
        procedure :: size => csr_op_size
            !! Queries dimensions of the owned CSR matrix.
        procedure :: apply => csr_op_apply
            !! Computes a CSR matrix-vector product into caller storage.
    end type

    type, extends(linear_operator) :: ilu_preconditioner
        !! Fixed approximate inverse backed by prepared ILUT triangular factors.
        !! This applies M**(-1), not the original matrix or a full GMRES solve.
        !! Constructor copies factors and pointers; applications use local output.
        type(msr_matrix), private :: factors
            !! Owned combined triangular factors with inverse U pivots.
        integer(int32), allocatable, private :: upper_rows(:)
            !! Owned pointers separating strict lower and upper entries.
    contains
        procedure :: size => ilu_op_size
            !! Queries the square dimensions of the prepared preconditioner.
        procedure :: apply => ilu_op_apply
            !! Applies the fixed approximate inverse with SPARSKIT LUSOL.
    end type

    interface csr_operator
        !! Constructs an owned operator from validated finite CSR storage.
        module procedure :: csr_op_init
    end interface

    interface ilu_preconditioner
        !! Constructs an owned approximate-inverse operator from ILUT factors.
        module procedure :: ilu_op_init
    end interface

    type :: iterative_solver_info
        !! Diagnostics for one solve, reset at entry and returned on all paths.
        integer :: status = LA_NO_ERROR
            !! Library status returned through FLAG; zero denotes success.
        integer :: backend_status = 0
            !! Last SPARSKIT termination code, or -1 if the wrapper enforces
            !! its product budget. Zero also means the backend was not invoked;
            !! inspect STATUS and CONVERGED to distinguish those cases.
        integer :: matvecs = 0
            !! Actual A applications, including the final true-residual check.
        integer :: preconditioner_applications = 0
            !! Actual right-preconditioner applications, including corrections.
        real(real64) :: residual_norm = huge(0.0d0)
            !! Final NORM2(B-A*X), or HUGE if it could not be evaluated safely.
        real(real64) :: target_residual = 0.0d0
            !! Required bound TOL*NORM2(B)+ATOL after successful validation.
        logical :: converged = .false.
            !! True only when the original-system residual meets its bound.
    end type

contains

pure function csr_op_init(a) result(rst)
    !! Copies a CSR matrix after validating its storage. Does not factor it.
    class(csr_matrix), intent(in) :: a
        !! Nonempty finite CSR matrix with unique columns within each row.
        !! Rectangular matrices are allowed here; GMRES requires square ones.
    type(csr_operator) :: rst
        !! New operator owning a deep copy of A. Invalid storage uses ERROR STOP
        !! with the status returned by VALIDATE_CSR_STORAGE.
    integer :: status
    status = validate_csr_storage(a)
    if (status /= LA_NO_ERROR) error stop status
    rst%matrix = a
end function

pure function csr_op_size(this, dim) result(rst)
    !! Queries one dimension of the owned coefficient matrix.
    class(csr_operator), intent(in) :: this
        !! Operator to inspect without modification.
    integer, intent(in) :: dim
        !! Dimension: 1 for rows, 2 for columns; otherwise returns zero.
    integer :: rst
        !! Requested extent, or zero for an uninitialized operator.
    rst = size(this%matrix, dim)
end function

subroutine csr_op_apply(this, x, y)
    !! Computes Y=A*X using validated CSR storage, without allocating Y.
    class(csr_operator), intent(in) :: this
        !! Initialized read-only coefficient matrix.
    real(real64), intent(in) :: x(:)
        !! Finite input with THIS%SIZE(2) elements; must not overlap Y.
    real(real64), intent(out) :: y(:)
        !! Caller-owned product with THIS%SIZE(1) elements. Incorrect dimensions,
        !! uninitialized storage, or nonfinite X use ERROR STOP.
    if (this%size(1) < 1 .or. this%size(2) < 1) error stop LA_NOT_INITIALIZED_ERROR
    if (size(x) /= this%size(2) .or. size(y) /= this%size(1)) &
        error stop LA_MATRIX_DIMENSION_ERROR
    if (.not.all(ieee_is_finite(x))) error stop LA_INVALID_INPUT_ERROR
    y = matmul(this%matrix, x)
end subroutine

pure function ilu_op_init(lu, ju) result(rst)
    !! Copies precomputed ILUT factors and upper-row pointers; no setup occurs
    !! during APPLY. Factors can come from LU_FACTOR or GET_FACTORED_FORM.
    class(msr_matrix), intent(in) :: lu
        !! Nonempty finite square MSR factors, with inverse upper pivots.
    integer(int32), intent(in) :: ju(:)
        !! N pointers to the upper parts of the factor rows.
    type(ilu_preconditioner) :: rst
        !! New fixed approximate-inverse operator owning independent copies.
        !! Invalid factors use ERROR STOP with VALIDATE_MSR_LU's status.
    integer :: status
    status = validate_msr_lu(lu, ju)
    if (status /= LA_NO_ERROR) error stop status
    rst%factors = lu
    rst%upper_rows = ju
end function

pure function ilu_op_size(this, dim) result(rst)
    !! Queries the prepared factor dimensions.
    class(ilu_preconditioner), intent(in) :: this
        !! Preconditioner to inspect without modification.
    integer, intent(in) :: dim
        !! Dimension: 1 for rows, 2 for columns; otherwise returns zero.
    integer :: rst
        !! Requested extent, or zero for an uninitialized preconditioner.
    rst = size(this%factors, dim)
end function

subroutine ilu_op_apply(this, x, y)
    !! Computes Y=M**(-1)*X by triangular solves, not by an inner Krylov solve.
    class(ilu_preconditioner), intent(in) :: this
        !! Prepared fixed ILUT factors and row pointers, never modified.
    real(real64), intent(in) :: x(:)
        !! Finite N-element input; must not overlap Y.
    real(real64), intent(out) :: y(:)
        !! Caller-owned N-element approximate-inverse result. Uninitialized
        !! factors, invalid vector sizes, or nonfinite X use ERROR STOP.
    if (.not.allocated(this%upper_rows)) error stop LA_NOT_INITIALIZED_ERROR
    if (size(x) /= this%size(2) .or. size(y) /= this%size(1)) &
        error stop LA_MATRIX_DIMENSION_ERROR
    if (.not.all(ieee_is_finite(x))) error stop LA_INVALID_INPUT_ERROR
    call lusol(this%size(1), x, y, this%factors%values, &
        this%factors%indices, this%upper_rows)
end subroutine

subroutine gmres_solve(a, b, x, preconditioner, tol, atol, im, maxits, x0, flag, info)
    !! Solves A*X=B with restarted GMRES and optional fixed right preconditioning.
    !! Only operator applications are requested; no entries or transpose are
    !! needed. All callbacks must represent fixed linear maps during a solve.
    !! Nonfinite callback outputs are rejected, but callbacks must manage their
    !! own failures because LINEAR_OPERATOR%APPLY has no status argument.
    !! The original residual is checked independently of SPARSKIT's estimates,
    !! including after a reported Arnoldi breakdown (which may be convergence).
    !! This routine is NOT reentrant because the bundled backend uses SAVE.
    class(linear_operator), intent(in) :: a
        !! Fixed real N-by-N operator, with N positive. APPLY and SIZE must not
        !! modify shared data or invoke this GMRES backend recursively.
    real(real64), intent(in) :: b(:)
        !! Finite N-element RHS, preserved. Must not overlap X.
    real(real64), intent(out), allocatable :: x(:)
        !! Allocated N-element solution on success, latest backend iterate on
        !! iteration/callback failure. Unallocated on input/allocation errors.
        !! An iterate returned with a nonzero FLAG must not be assumed accurate.
    class(linear_operator), intent(in), optional :: preconditioner
        !! Fixed N-by-N action Z=M**(-1)*R for RIGHT preconditioning. It applies
        !! the approximate inverse, not M itself. Omit for no preconditioning.
        !! Variable/adaptive preconditioning requires a flexible method instead.
    real(real64), intent(in), optional :: tol
        !! Finite RHS-relative tolerance in [0,1); default SQRT(EPSILON).
        !! Success requires NORM2(B-A*X) <= TOL*NORM2(B)+ATOL.
    real(real64), intent(in), optional :: atol
        !! Finite nonnegative absolute residual tolerance, default zero.
        !! TOL and ATOL must not both be zero.
    integer, intent(in), optional :: im
        !! Restart length in [2,MAX(2,N)], default MAX(2,MIN(N,30)). There is
        !! no legacy PGMRES limit of 50. One is rejected because this backend
        !! silently replaces lengths <=1 with 15. Scalar systems use length 2.
    integer, intent(in), optional :: maxits
        !! Positive budget for backend-requested A products, default 300.
        !! Includes initial/restart residual products, not preconditioner calls
        !! or the wrapper's final residual check. Not a Krylov-step count.
    real(real64), intent(in), optional :: x0(:)
        !! Finite N-element initial guess, preserved. Must not overlap X.
        !! Omit for zero.
    integer, intent(out), optional :: flag
        !! LA_NO_ERROR on verified convergence; LA_INVALID_INPUT_ERROR for
        !! invalid controls, empty dimensions, or nonfinite B/X0;
        !! LA_MATRIX_NOT_SQUARE_ERROR or LA_MATRIX_DIMENSION_ERROR for shapes;
        !! LA_OUT_OF_MEMORY_ERROR for allocation or 32-bit workspace overflow;
        !! LA_INVALID_OPERATION_ERROR for unexpected requests/nonfinite outputs;
        !! LA_CONVERGENCE_ERROR for budget, breakdown, or residual failure.
        !! Omit to use ERROR STOP on error, even if INFO is present.
    type(iterative_solver_info), intent(out), optional :: info
        !! Diagnostics, reset and populated even when a supplied FLAG reports
        !! failure. Counts include only callbacks that were actually performed.
    type(iterative_solver_info) :: details
    integer(int32) :: ipar(16)
    integer :: order, restart, budget, status, allocation_status
    integer :: input_start, output_start, backend_products
    integer(int64) :: basis, workspace_length
    real(real64) :: fpar(16), relative_tol, absolute_tol, rhs_norm
    real(real64), allocatable :: workspace(:), residual(:)

    order = a%size(1)
    restart = max(2, min(order, 30))
    budget = 300
    relative_tol = sqrt(epsilon(0.0d0))
    absolute_tol = 0.0d0
    if (present(im)) restart = im
    if (present(maxits)) budget = maxits
    if (present(tol)) relative_tol = tol
    if (present(atol)) absolute_tol = atol
    status = LA_NO_ERROR
    if (order < 1) then
        status = LA_INVALID_INPUT_ERROR
    else if (a%size(2) /= order) then
        status = LA_MATRIX_NOT_SQUARE_ERROR
    else if (size(b) /= order) then
        status = LA_MATRIX_DIMENSION_ERROR
    else if (.not.all(ieee_is_finite(b)) .or. &
        .not.ieee_is_finite(relative_tol) .or. .not.ieee_is_finite(absolute_tol)) then
        status = LA_INVALID_INPUT_ERROR
    else if (relative_tol < 0.0d0 .or. relative_tol >= 1.0d0 .or. &
        absolute_tol < 0.0d0 .or. &
        (relative_tol == 0.0d0 .and. absolute_tol == 0.0d0) .or. &
        restart < 2 .or. restart > max(2, order) .or. budget < 1) then
        status = LA_INVALID_INPUT_ERROR
    end if
    if (present(x0)) then
        if (size(x0) /= order) then
            status = LA_MATRIX_DIMENSION_ERROR
        else if (.not.all(ieee_is_finite(x0))) then
            status = LA_INVALID_INPUT_ERROR
        end if
    end if
    if (present(preconditioner)) then
        if (preconditioner%size(1) /= order .or. &
            preconditioner%size(2) /= order) status = LA_MATRIX_DIMENSION_ERROR
    end if
    if (status /= LA_NO_ERROR) then
        details%status = status
        call report_gmres_result(details, flag, info)
        return
    end if

    rhs_norm = norm2(b)
    details%target_residual = relative_tol * rhs_norm + absolute_tol
    basis = int(restart, int64)
    workspace_length = (int(order, int64)+3) * (basis+2) + basis*(basis+1)/2
    if (.not.ieee_is_finite(rhs_norm) .or. &
        .not.ieee_is_finite(details%target_residual)) then
        status = LA_INVALID_INPUT_ERROR
    else if (workspace_length > int(huge(0_int32), int64) .or. &
        basis*(basis+1) > int(huge(0_int32), int64)) then
        status = LA_OUT_OF_MEMORY_ERROR
    end if
    if (status /= LA_NO_ERROR) then
        details%status = status
        call report_gmres_result(details, flag, info)
        return
    end if
    allocate(x(order), workspace(int(workspace_length)), residual(order), &
        stat = allocation_status)
    if (allocation_status /= 0) then
        if (allocated(x)) deallocate(x)
        details%status = LA_OUT_OF_MEMORY_ERROR
        call report_gmres_result(details, flag, info)
        return
    end if
    x = 0.0d0
    if (present(x0)) x = x0
    workspace = 0.0d0
    ipar = 0
    fpar = 0.0d0
    ipar(2) = 0
    if (present(preconditioner)) ipar(2) = 2
    ipar(3) = 2
    ipar(4) = int(workspace_length, int32)
    ipar(5) = restart
    ipar(6) = budget
    fpar(1) = relative_tol
    fpar(2) = absolute_tol
    backend_products = 0

    do
        call gmres(order, b, x, ipar, fpar, workspace)
        if (ipar(1) <= 0) exit
        input_start = ipar(8)
        output_start = ipar(9)
        if (input_start < 1 .or. output_start < 1 .or. &
            input_start > size(workspace)-order+1 .or. &
            output_start > size(workspace)-order+1) then
            status = LA_INVALID_OPERATION_ERROR
            exit
        end if
        if (abs(input_start-output_start) < order) then
            status = LA_INVALID_OPERATION_ERROR
            exit
        end if
        if (.not.all(ieee_is_finite(workspace(input_start:input_start+order-1)))) then
            status = LA_INVALID_OPERATION_ERROR
            exit
        end if
        select case (ipar(1))
        case (1)
            if (backend_products >= budget) then
                ipar(1) = -1
                status = LA_CONVERGENCE_ERROR
                exit
            end if
            call a%apply(workspace(input_start:input_start+order-1), &
                workspace(output_start:output_start+order-1))
            backend_products = backend_products + 1
            details%matvecs = details%matvecs + 1
        case (5)
            if (.not.present(preconditioner)) then
                status = LA_INVALID_OPERATION_ERROR
                exit
            end if
            call preconditioner%apply(workspace(input_start:input_start+order-1), &
                workspace(output_start:output_start+order-1))
            details%preconditioner_applications = details%preconditioner_applications + 1
        case default
            status = LA_INVALID_OPERATION_ERROR
            exit
        end select
        if (.not.all(ieee_is_finite(workspace(output_start:output_start+order-1)))) then
            status = LA_INVALID_OPERATION_ERROR
            exit
        end if
    end do
    details%backend_status = ipar(1)

    if (status /= LA_INVALID_OPERATION_ERROR .and. all(ieee_is_finite(x))) then
        call a%apply(x, residual)
        details%matvecs = details%matvecs + 1
        residual = b - residual
        if (all(ieee_is_finite(residual))) then
            details%residual_norm = norm2(residual)
            details%converged = details%residual_norm <= details%target_residual
            if (details%converged) then
                status = LA_NO_ERROR
            else
                status = LA_CONVERGENCE_ERROR
            end if
        else
            status = LA_INVALID_OPERATION_ERROR
        end if
    else
        status = LA_INVALID_OPERATION_ERROR
    end if
    details%status = status
    call report_gmres_result(details, flag, info)
end subroutine

subroutine report_gmres_result(details, flag, info)
    !! Copies diagnostics and applies the optional-status/error-stop policy.
    type(iterative_solver_info), intent(in) :: details
        !! Completed diagnostics, including the final library status.
    integer, intent(out), optional :: flag
        !! Receives DETAILS%STATUS; omission makes nonzero statuses terminate.
    type(iterative_solver_info), intent(out), optional :: info
        !! Receives all diagnostics regardless of success or failure.
    if (present(info)) info = details
    if (present(flag)) then
        flag = details%status
    else if (details%status /= LA_NO_ERROR) then
        error stop details%status
    end if
end subroutine

end module