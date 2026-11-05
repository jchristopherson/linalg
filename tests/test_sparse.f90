! test_sparse.f90
! SPDX-License-Identifier: MIT
module test_sparse
    use iso_fortran_env
    use linalg
    use fortran_test_helper
    implicit none
    private :: stencil_test_operator

    type, extends(linear_operator) :: stencil_test_operator
        !! Matrix-free nonsymmetric three-diagonal stencil for solver tests.
        real(real64), allocatable :: diagonal(:)
            !! Diagonal coefficients; allocation determines the row count.
        real(real64) :: lower = -1.0d0
            !! Constant coefficient immediately below the diagonal.
        real(real64) :: upper = -0.5d0
            !! Constant coefficient immediately above the diagonal.
        integer :: extra_columns = 0
            !! Added column extent, used only to test shape rejection.
        logical :: invalid_output = .false.
            !! If true, APPLY deliberately emits NaNs to test callback checks.
    contains
        procedure :: size => stencil_test_size
            !! Queries row/column counts without materializing the matrix.
        procedure :: apply => stencil_test_apply
            !! Applies the stencil to caller-owned vector storage.
    end type

contains
pure function stencil_test_size(this, dim) result(rst)
    !! Queries dimensions of the test stencil.
    class(stencil_test_operator), intent(in) :: this
        !! Operator whose diagonal determines the row count.
    integer, intent(in) :: dim
        !! Dimension: 1 for rows, 2 for columns; otherwise returns zero.
    integer :: rst
        !! Requested extent, or zero when no diagonal has been allocated.
    rst = 0
    if (.not.allocated(this%diagonal)) return
    if (dim == 1) rst = size(this%diagonal)
    if (dim == 2) rst = size(this%diagonal) + this%extra_columns
end function

subroutine stencil_test_apply(this, x, y)
    !! Applies a three-diagonal stencil without assembled matrix storage.
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    class(stencil_test_operator), intent(in) :: this
        !! Fixed stencil coefficients and optional deliberate-failure switch.
    real(real64), intent(in) :: x(:)
        !! Input vector with THIS%SIZE(2) elements; must not overlap Y.
    real(real64), intent(out) :: y(:)
        !! Caller-owned product with THIS%SIZE(1) elements.
    integer :: order
    if (this%invalid_output) then
        y = ieee_value(0.0d0, ieee_quiet_nan)
        return
    end if
    order = size(this%diagonal)
    y = this%diagonal * x
    y(:order-1) = y(:order-1) + this%upper * x(2:)
    y(2:) = y(2:) + this%lower * x(:order-1)
end subroutine

function test_operator_gmres() result(rst)
    !! Tests matrix-free SPARSKIT communication, fixed right preconditioning,
    !! restarts, guesses, true residuals, input guards, breakdown, and budgets.
    use linalg_errors
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    logical :: rst
        !! True when all solutions, diagnostics, and expected failures agree.
    integer, parameter :: order = 12
    type(stencil_test_operator) :: op, preconditioner, invalid, scalar, large
    type(iterative_solver_info) :: details
    real(real64) :: matrix(order,order), expected(order), rhs(order), nan
    real(real64) :: product(order), invalid_rhs(order)
    real(real64), allocatable :: solution(:), large_rhs(:)
    integer :: index, flag

    rst = .true.
    nan = ieee_value(0.0d0, ieee_quiet_nan)
    op%diagonal = [(4.0d0 + real(index, real64)/10.0d0, index = 1, order)]
    matrix = 0.0d0
    do index = 1, order
        matrix(index,index) = op%diagonal(index)
        if (index < order) then
            matrix(index,index+1) = op%upper
            matrix(index+1,index) = op%lower
        end if
    end do
    expected = [(real(index, real64)/3.0d0, index = 1, order)]
    rhs = matmul(matrix, expected)
    call gmres_solve(op, rhs, solution, tol = 1.0d-11, im = 4, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = assert(solution, expected, tol = 1.0d-8) .and. rst
    rst = details%residual_norm <= details%target_residual .and. rst
    rst = details%matvecs > 2 .and. details%preconditioner_applications == 0 .and. rst
    rst = assert(details%residual_norm, norm2(rhs-matmul(matrix, solution)), &
        tol = 1.0d-12) .and. rst

    preconditioner%diagonal = 1.0d0 / op%diagonal
    preconditioner%lower = 0.0d0
    preconditioner%upper = 0.0d0
    call gmres_solve(op, rhs, solution, preconditioner = preconditioner, &
        tol = 1.0d-11, im = 4, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = details%preconditioner_applications > 0 .and. rst
    rst = assert(solution, expected, tol = 1.0d-8) .and. rst
    call gmres_solve(op, rhs, solution, x0 = expected, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = details%matvecs == 2 .and. rst
    rst = assert(solution, expected, tol = 1.0d-12) .and. rst
    call gmres_solve(op, rhs, solution, tol = 0.0d0, atol = 1.0d-9, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%residual_norm <= 1.0d-9 .and. rst
    call gmres_solve(op, 0.0d0 * rhs, solution, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. all(solution == 0.0d0) .and. rst
    rst = details%converged .and. details%residual_norm == 0.0d0 .and. rst

    scalar%diagonal = [2.0d0]
    call gmres_solve(scalar, [6.0d0], solution, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = assert(solution, [3.0d0], tol = 1.0d-12) .and. rst
    call gmres_solve(op, rhs, solution, im = 2, maxits = 2, &
        tol = 1.0d-12, flag = flag, info = details)
    rst = flag == LA_CONVERGENCE_ERROR .and. .not.details%converged .and. rst
    rst = allocated(solution) .and. details%status == flag .and. rst
    invalid = op
    invalid%diagonal = 0.0d0
    invalid%lower = 0.0d0
    invalid%upper = 0.0d0
    call gmres_solve(invalid, rhs, solution, flag = flag, info = details)
    rst = flag == LA_CONVERGENCE_ERROR .and. .not.details%converged .and. rst

    call gmres_solve(op, rhs(:order-1), solution, flag = flag)
    rst = flag == LA_MATRIX_DIMENSION_ERROR .and. .not.allocated(solution) .and. rst
    invalid = op
    invalid%extra_columns = 1
    call gmres_solve(invalid, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_NOT_SQUARE_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, preconditioner = invalid, flag = flag)
    rst = flag == LA_MATRIX_DIMENSION_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, x0 = expected(:order-1), flag = flag)
    rst = flag == LA_MATRIX_DIMENSION_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, im = 1, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, maxits = 0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, tol = -1.0d0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, tol = 1.0d0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, tol = 0.0d0, atol = 0.0d0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, atol = -1.0d0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, tol = nan, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_rhs = rhs
    invalid_rhs(1) = nan
    call gmres_solve(op, invalid_rhs, solution, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call gmres_solve(op, rhs, solution, x0 = invalid_rhs, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    invalid = op
    invalid%invalid_output = .true.
    call gmres_solve(invalid, rhs, solution, flag = flag, info = details)
    rst = flag == LA_INVALID_OPERATION_ERROR .and. .not.details%converged .and. rst
    call gmres_solve(op, rhs, solution, preconditioner = invalid, flag = flag, info = details)
    rst = flag == LA_INVALID_OPERATION_ERROR .and. .not.details%converged .and. rst
    call gmres_solve(op, rhs, solution, tol = 1.0d-11)
    rst = assert(solution, expected, tol = 1.0d-8) .and. rst
    call op%apply(expected, product)
    rst = assert(product, rhs) .and. rst

    allocate(large%diagonal(50000), large_rhs(50000), source = 1.0d0)
    call gmres_solve(large, large_rhs, solution, im = 50000, flag = flag, info = details)
    rst = flag == LA_OUT_OF_MEMORY_ERROR .and. .not.allocated(solution) .and. rst
    rst = details%matvecs == 0 .and. .not.details%converged .and. rst
    large%lower = 0.0d0
    large%upper = 0.0d0
    large%diagonal = large%diagonal(:64)
    call gmres_solve(large, large_rhs(:64), solution, im = 64, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = all(abs(solution - 1.0d0) < 1.0d-10) .and. rst
    if (.not.rst) print '(A)', 'Test Failed: Matrix-Free Reverse-Communication GMRES'
end function

! ------------------------------------------------------------------------------
function test_gmres_sparse_adapters() result(rst)
    !! Verifies owned adapter storage, fixed ILUT application, hybrid matrix-free
    !! solves, and agreement with the unchanged CSR GMRES entry point.
    use linalg_errors, only : LA_NO_ERROR
    logical :: rst
        !! True if adapter products, solutions, ownership, and diagnostics agree.
    integer, parameter :: order = 8
    type(csr_matrix) :: csr
    type(msr_matrix) :: factors
    type(csr_operator) :: op
    type(ilu_preconditioner) :: preconditioner
    type(stencil_test_operator) :: stencil
    type(iterative_solver_info) :: details
    real(real64) :: matrix(order,order), expected(order), rhs(order), product(order)
    real(real64), allocatable :: solution(:), legacy_solution(:)
    integer :: upper_rows(order)
    integer :: index, flag

    rst = .true.
    matrix = 0.0d0
    stencil%diagonal = [(4.0d0 + real(index, real64)/10.0d0, index = 1, order)]
    do index = 1, order
        matrix(index,index) = stencil%diagonal(index)
        if (index < order) then
            matrix(index,index+1) = stencil%upper
            matrix(index+1,index) = stencil%lower
        end if
    end do
    expected = [(real(index, real64)/3.0d0, index = 1, order)]
    rhs = matmul(matrix, expected)
    csr = dense_to_csr(matrix)
    op = csr_operator(csr)
    call lu_factor(csr, factors, upper_rows, droptol = 1.0d0)
    preconditioner = ilu_preconditioner(factors, upper_rows)
    call preconditioner%apply(rhs, product)
    rst = assert(product, solve_lu(factors, upper_rows, rhs), tol = 1.0d-12) .and. rst
    legacy_solution = pgmres_solver(csr, factors, upper_rows, rhs, tol = 1.0d-11)

    csr%values = 0.0d0
    factors%values = 0.0d0
    upper_rows = 1
    call op%apply(expected, product)
    rst = assert(product, rhs) .and. rst
    rst = op%size(1) == order .and. preconditioner%size(2) == order .and. rst
    call gmres_solve(op, rhs, solution, preconditioner = preconditioner, &
        tol = 1.0d-11, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = details%preconditioner_applications > 0 .and. rst
    rst = assert(solution, legacy_solution, tol = 1.0d-8) .and. rst
    rst = assert(solution, expected, tol = 1.0d-8) .and. rst
    call gmres_solve(stencil, rhs, solution, preconditioner = preconditioner, &
        tol = 1.0d-11, flag = flag, info = details)
    rst = flag == LA_NO_ERROR .and. details%converged .and. rst
    rst = assert(solution, expected, tol = 1.0d-8) .and. rst
    if (.not.rst) print '(A)', 'Test Failed: GMRES Sparse Operator Adapters'
end function

! ------------------------------------------------------------------------------
function test_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer, parameter :: m = 4
    integer, parameter :: n = 6
    integer, parameter :: nnz = 8
    real(real64) :: dense(m, n), v(nnz), check(m, n)
    integer :: ja(nnz), ia(m + 1), rows(nnz), cols(nnz)
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.

    ! Construct a small sparse matrix, but in dense form
    !
    ! | 1   2   0   0   0   0|
    ! | 0   3   0   4   0   0|
    ! | 0   0   5   6   7   0|
    ! | 0   0   0   0   0   8|
    dense = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [m, n])

    ! Define the sparse format - CSR
    v = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]
    ja = [1, 2, 2, 4, 3, 4, 5, 6]
    ia = [1, 3, 5, 8, 9]

    ! Form the sparse matrix
    sparse = dense_to_csr(dense)

    ! Test
    if (.not.assert(v, sparse%values)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -1"
    end if
    if (.not.assert(ja, sparse%column_indices)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -2"
    end if
    if (.not.assert(ia, sparse%row_indices)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -3"
    end if
    if (nonzero_count(sparse) /= nnz) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -4"
    end if

    ! Convert back to dense
    check = sparse

    ! Test
    if (.not.assert(dense, check)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -6"
    end if

    ! ----------
    ! Construct from coordinate format
    rows = [1, 1, 2, 2, 3, 3, 3, 4]
    cols = [1, 2, 2, 4, 3, 4, 5, 6]
    sparse = create_csr_matrix(m, n, rows, cols, v)

    ! Convert back to dense
    check = sparse

    ! Test
    if (.not.assert(dense, check)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_1 -7"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = matmul(sa, sa)
    ans = matmul(a, a)

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_mult_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_mult_2() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), b(4, 6), ans(4, 6), dense(4, 6)
    type(csr_matrix) :: sa, sb, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])
    b = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)
    sb = dense_to_csr(b)

    ! Multiply
    sc = matmul(sa, sb)
    ans = matmul(a, b)

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_mult_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_add_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Add
    sc = sa + sa
    ans = a + a

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_add_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_subtract_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), dense(4, 4)
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Add
    sc = sa + sa - sa

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_subtract_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_scalar_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4), s
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.
    call random_number(s)

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = s * sa
    ans = s * a

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_mult_1 -1"
    end if

    ! Test 2
    sc = sa * s
    dense = sc
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_mult_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_scalar_divide_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 4), ans(4, 4), dense(4, 4), s
    type(csr_matrix) :: sa, sc

    ! Initialization
    rst = .true.
    call random_number(s)

    ! Initialize the dense matrices
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0], [4, 4])

    ! Construct a sparse matrix from the dense matrix
    sa = dense_to_csr(a)

    ! Multiply
    sc = sa / s
    ans = a / s

    ! Convert the sparse matrix to a dense matrix for testing purposes
    dense = sc

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_scalar_divide_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_transpose_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4, 6), ans(6, 4), dense(6, 4)
    type(csr_matrix) :: sa, sat

    ! Initialization
    rst = .true.
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    ans = transpose(a)
    sa = dense_to_csr(a)

    ! Operation
    sat = transpose(sa)
    dense = sat

    ! Test 1
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_transpose_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_diag_mult_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: alpha, a(4, 6), diag(4, 4), d(4), ans(4, 6), dense(4, 6)
    type(csr_matrix) :: sa

    ! Initialization
    rst = .true.
    call random_number(alpha)
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    diag = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 6.0d0], [4, 4])
    d = [diag(1,1), diag(2,2), diag(3,3), diag(4,4)]
    sa = dense_to_csr(a)

    ! Test 1
    ans = alpha * matmul(diag, a)
    call diag_mtx_mult(.true., alpha, d, sa)
    dense = sa
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_1 -1"
    end if

    ! Test 2
    sa = dense_to_csr(a)
    ans = matmul(diag, a)
    call diag_mtx_mult(.true., 1.0d0, d, sa)
    dense = sa
    if (.not.assert(dense, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_diag_mult_2() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: alpha, a(4, 6), diag(6), dense(4, 6)
    type(csr_matrix) :: sa

    ! Initialization
    rst = .true.
    call random_number(alpha)
    a = reshape([ &
        1.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 5.0d0, 0.0d0, &
        0.0d0, 4.0d0, 6.0d0, 0.0d0, &
        0.0d0, 0.0d0, 7.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 8.0d0], [4, 6])
    call random_number(diag)
    sa = dense_to_csr(a)

    ! Test 1
    call diag_mtx_mult(.false., alpha, diag, a) ! solution
    call diag_mtx_mult(.false., alpha, diag, sa)
    dense = sa
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_2 -1"
    end if

    ! Test 2
    sa = dense_to_csr(a)
    call diag_mtx_mult(.false., 1.0d0, diag, a) ! solution
    call diag_mtx_mult(.false., 1.0d0, diag, sa)
    dense = sa
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_diag_mult_2 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_sparse_direct_solve_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer, allocatable :: ipiv(:)
    real(real64) :: a(4, 4), b(4), x(4), ans(4)
    real(real64), allocatable :: dlu(:,:)
    type(csr_matrix) :: sa

    ! Initialization
    rst = .true.
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])
    call random_number(b)
    sa = dense_to_csr(a)

    ! Compute the solution
    call lu_factor(a, ipvt = ipiv, lu = dlu)
    ans = solve_lu(dlu, ipiv, b)

    ! Test
    x = sparse_direct_solve(sa, b)
    if (.not.assert(x, ans, tol = 1.0d-6)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_sparse_direct_solve_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_sparse_direct_complete() result(rst)
    !! Verifies that the legacy convenience solver solves the original CSR
    !! system even when ILUT drops substantial off-diagonal information.
    logical :: rst
        !! True if the solution and the original-system residual are accurate.
    integer, parameter :: order = 8
    real(real64) :: matrix(order,order), expected(order), rhs(order)
    real(real64), allocatable :: solution(:)
    type(csr_matrix) :: csr
    integer :: index

    matrix = 0.0d0
    do index = 1, order
        matrix(index,index) = 4.0d0
        if (index < order) then
            matrix(index,index+1) = -0.5d0
            matrix(index+1,index) = -1.0d0
        end if
    end do
    expected = [(real(index, real64), index = 1, order)]
    rhs = matmul(matrix, expected)
    csr = dense_to_csr(matrix)
    solution = sparse_direct_solve(csr, rhs, droptol = 1.0d0)
    rst = assert(solution, expected, tol = 1.0d-6)
    rst = assert(matmul(matrix, solution), rhs, tol = 1.0d-6) .and. rst
    if (.not.rst) print '(A)', 'Test Failed: Complete CSR Convenience Solve'
end function

! ------------------------------------------------------------------------------
subroutine csr_solver_fixture(csr, factors, upper_rows, rhs)
    !! Builds an independent nonsymmetric system and its exact triangular
    !! factors without calling any of the factorization routines under test.
    type(csr_matrix), intent(out) :: csr
        !! Matrix A = L*U, with rows [2,1,-2], [4,5,0], and [-2,8,19].
    type(msr_matrix), intent(out) :: factors
        !! Combined L/U factors with inverse upper pivots in MSR storage.
    integer, intent(out) :: upper_rows(3)
        !! Pointers separating each row's lower and upper entries.
    real(real64), intent(out) :: rhs(3)
        !! A*[1,2,3], yielding [-2,14,71].

    csr = dense_to_csr(reshape([2.0d0, 4.0d0, -2.0d0, &
        1.0d0, 5.0d0, 8.0d0, -2.0d0, 0.0d0, 19.0d0], [3,3]))
    factors%m = 3
    factors%n = 3
    factors%nnz = 9
    factors%values = [0.5d0, 1.0d0/3.0d0, 0.2d0, 0.0d0, &
        1.0d0, -2.0d0, 2.0d0, 4.0d0, -1.0d0, 3.0d0]
    factors%indices = [5, 7, 9, 11, 2, 3, 1, 3, 1, 2]
    upper_rows = [5, 8, 11]
    rhs = [-2.0d0, 14.0d0, 71.0d0]
end subroutine

function test_csr_solver_validation() result(rst)
    !! Checks independent triangular solves, both GMRES APIs, scalar and zero
    !! systems, true residuals, malformed storage, invalid controls, and failures.
    use linalg_errors
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_is_nan
    logical :: rst
        !! True if every solution, residual, storage check, and status is correct.
    type(csr_matrix) :: csr, invalid_csr, scalar, uninitialized_csr
    type(msr_matrix) :: factors, invalid_factors, uninitialized_factors
    integer :: upper_rows(3), invalid_rows(3), scalar_rows(1), flag
    real(real64) :: rhs(3), saved_rhs(3), nonfinite_rhs(3), nan, infinity
    real(real64), allocatable :: solution(:), function_solution(:)

    rst = .true.
    call csr_solver_fixture(csr, factors, upper_rows, rhs)
    saved_rhs = rhs
    nan = ieee_value(0.0d0, ieee_quiet_nan)
    infinity = ieee_value(0.0d0, ieee_positive_inf)
    rst = validate_csr_storage(csr) == LA_NO_ERROR .and. rst
    rst = validate_msr_lu(factors, upper_rows) == LA_NO_ERROR .and. rst
    solution = solve_lu(factors, upper_rows, rhs)
    rst = assert(solution, [1.0d0, 2.0d0, 3.0d0], tol = 1.0d-12) .and. rst
    rst = assert(matmul(csr, solution), rhs, tol = 1.0d-12) .and. rst
    function_solution = pgmres_solver(csr, factors, upper_rows, rhs, tol = 1.0d-12)
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, &
        tol = 1.0d-12, flag = flag)
    rst = flag == LA_NO_ERROR .and. rst
    rst = assert(solution, function_solution, tol = 1.0d-12) .and. rst
    rst = norm2(rhs - matmul(csr, solution)) <= 1.0d-12 * norm2(rhs) .and. rst
    rst = assert(rhs, saved_rhs) .and. rst

    call pgmres_solve(csr, factors, upper_rows, [0.0d0, 0.0d0, 0.0d0], &
        solution, flag = flag)
    rst = flag == LA_NO_ERROR .and. rst
    rst = all(solution == 0.0d0) .and. rst
    function_solution = pgmres_solver(csr, factors, upper_rows, [0.0d0, 0.0d0, 0.0d0])
    rst = all(function_solution == 0.0d0) .and. rst
    solution = solve_lu(factors, upper_rows, [0.0d0, 0.0d0, 0.0d0])
    rst = all(solution == 0.0d0) .and. rst
    solution = sparse_direct_solve(csr, [0.0d0, 0.0d0, 0.0d0])
    rst = all(solution == 0.0d0) .and. rst

    scalar = dense_to_csr(reshape([2.0d0], [1,1]))
    solution = sparse_direct_solve(scalar, [6.0d0])
    rst = assert(solution, [3.0d0], tol = 1.0d-12) .and. rst
    call lu_factor(scalar, invalid_factors, scalar_rows)
    solution = solve_lu(invalid_factors, scalar_rows, [6.0d0])
    rst = assert(solution, [3.0d0], tol = 1.0d-12) .and. rst
    call pgmres_solve(scalar, invalid_factors, scalar_rows, [6.0d0], solution, flag = flag)
    rst = flag == LA_NO_ERROR .and. rst
    rst = assert(solution, [3.0d0], tol = 1.0d-12) .and. rst

    invalid_csr = csr
    invalid_csr%values = 0.0d0
    solution = sparse_direct_solve(invalid_csr, rhs)
    rst = all(ieee_is_nan(solution)) .and. rst
    invalid_factors = factors
    invalid_factors%values = nan
    solution = solve_lu(invalid_factors, upper_rows, rhs)
    rst = all(ieee_is_nan(solution)) .and. rst
    call pgmres_solve(csr, invalid_factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_SINGULAR_MATRIX_ERROR .and. .not.allocated(solution) .and. rst
    invalid_factors = factors
    invalid_factors%values(4) = nan
    rst = validate_msr_lu(invalid_factors, upper_rows) == LA_NO_ERROR .and. rst
    solution = solve_lu(invalid_factors, upper_rows, rhs)
    rst = assert(solution, [1.0d0, 2.0d0, 3.0d0], tol = 1.0d-12) .and. rst

    call pgmres_solve(csr, factors, upper_rows, rhs(:2), solution, flag = flag)
    rst = flag == LA_MATRIX_DIMENSION_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows(:2), rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_DIMENSION_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, im = 0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, im = 51, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, maxits = 0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, tol = 0.0d0, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, tol = nan, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    call pgmres_solve(csr, factors, upper_rows, rhs, solution, tol = infinity, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    nonfinite_rhs = rhs
    nonfinite_rhs(1) = nan
    call pgmres_solve(csr, factors, upper_rows, nonfinite_rhs, solution, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst

    rst = validate_csr_storage(uninitialized_csr) == LA_NOT_INITIALIZED_ERROR .and. rst
    rst = validate_msr_lu(uninitialized_factors, upper_rows(:0)) == &
        LA_NOT_INITIALIZED_ERROR .and. rst
    invalid_csr = csr
    invalid_csr%row_indices(1) = 0
    call pgmres_solve(invalid_csr, factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_FORMAT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_csr = csr
    invalid_csr%column_indices(1) = 4
    call pgmres_solve(invalid_csr, factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_FORMAT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_csr = csr
    invalid_csr%values(1) = infinity
    call pgmres_solve(invalid_csr, factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_INVALID_INPUT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_csr = dense_to_csr(reshape([1.0d0, 2.0d0, 3.0d0], [3,1]))
    call pgmres_solve(invalid_csr, factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_NOT_SQUARE_ERROR .and. .not.allocated(solution) .and. rst
    invalid_factors = factors
    invalid_factors%indices(4) = 100
    call pgmres_solve(csr, invalid_factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_FORMAT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_rows = upper_rows
    invalid_rows(2) = 6
    call pgmres_solve(csr, factors, invalid_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_FORMAT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_factors = factors
    invalid_factors%indices(5) = 1
    call pgmres_solve(csr, invalid_factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_MATRIX_FORMAT_ERROR .and. .not.allocated(solution) .and. rst
    invalid_factors = factors
    deallocate(invalid_factors%values)
    call pgmres_solve(csr, invalid_factors, upper_rows, rhs, solution, flag = flag)
    rst = flag == LA_NOT_INITIALIZED_ERROR .and. .not.allocated(solution) .and. rst

    invalid_factors%m = 3
    invalid_factors%n = 3
    invalid_factors%values = [1.0d0, 1.0d0, 1.0d0, 0.0d0]
    invalid_factors%indices = [5, 5, 5, 5]
    invalid_rows = [5, 5, 5]
    call pgmres_solve(csr, invalid_factors, invalid_rows, rhs, solution, &
        im = 1, maxits = 1, tol = 1.0d-12, flag = flag)
    rst = flag == LA_CONVERGENCE_ERROR .and. allocated(solution) .and. rst
    rst = norm2(rhs - matmul(csr, solution)) > 1.0d-12 * norm2(rhs) .and. rst
    rst = assert(rhs, saved_rhs) .and. rst
    if (.not.rst) print '(A)', 'Test Failed: CSR Solver Validation'
end function

subroutine run_csr_solver_error_case(case_name)
    !! Exercises ERROR STOP paths in separate expected-failure CTest processes.
    !! If the selected invalid call returns, the driver exits successfully so
    !! CTest's WILL_FAIL property marks that missing rejection as a failure.
    character(len = *), intent(in) :: case_name
        !! Registered invalid-input case name; an unknown name simply returns.
    type(csr_matrix) :: csr
    type(msr_matrix) :: factors
    type(csr_operator) :: op
    type(ilu_preconditioner) :: preconditioner
    integer :: upper_rows(3)
    real(real64) :: rhs(3)
    real(real64), allocatable :: solution(:)

    call csr_solver_fixture(csr, factors, upper_rows, rhs)
    select case (case_name)
    case ('csr_direct_bad_rhs')
        solution = sparse_direct_solve(csr, rhs(:2))
    case ('csr_direct_bad_storage')
        csr%row_indices(1) = 0
        solution = sparse_direct_solve(csr, rhs)
    case ('csr_lu_bad_rhs')
        solution = solve_lu(factors, upper_rows, rhs(:2))
    case ('csr_lu_bad_pointer')
        upper_rows(2) = 100
        solution = solve_lu(factors, upper_rows, rhs)
    case ('csr_gmres_bad_restart')
        solution = pgmres_solver(csr, factors, upper_rows, rhs, im = 51)
    case ('csr_gmres_bad_storage')
        factors%indices(5) = 1
        solution = pgmres_solver(csr, factors, upper_rows, rhs)
    case ('operator_gmres_bad_restart')
        op = csr_operator(csr)
        call gmres_solve(op, rhs, solution, im = 1)
    case ('csr_operator_bad_storage')
        csr%row_indices(1) = 0
        op = csr_operator(csr)
    case ('ilu_operator_bad_pointer')
        upper_rows(2) = 100
        preconditioner = ilu_preconditioner(factors, upper_rows)
    end select
end subroutine

! ------------------------------------------------------------------------------
function test_diag_to_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer, parameter :: n = 20
    integer :: i
    real(real64) :: d(n), diag(n,n), dense(n,n)
    type(csr_matrix) :: sd

    ! Initialization
    rst = .true.
    call random_number(d)
    diag = 0.0d0
    do i = 1, n
        diag(i,i) = d(i)
    end do

    ! Test 1
    sd = diag_to_csr(d)
    dense = sd
    if (.not.assert(dense, diag)) then
        rst = .false.
        print "(A)", "Test Failed: test_diag_to_csr_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_banded_to_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer, parameter :: kl = 4
    integer, parameter :: ku = 5
    integer, parameter :: m = 15
    integer, parameter :: n = 15
    real(real64) :: banded(kl+ku+1,n), dense(m,n), check(m,n)
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.
    call random_number(banded)

    ! Construct the dense matrix directly from the banded matrix
    dense = banded_to_dense(m, kl, ku, banded)

    ! Construct the sparse matrix directly from the banded matrix
    sparse = banded_to_csr(m, kl, ku, banded)

    ! Test
    check = sparse
    if (.not.assert(dense, check)) then
        rst = .false.
        print "(A)", "Test Failed: test_banded_to_csr_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_extract_diagonal_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer, parameter :: n = 100
    integer :: i
    real(real64) :: diag(n), ans(n), dense(n,n), densediag(n)
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.
    call random_number(ans)
    dense = 0.0d0
    do i = 1, n
        dense(i,i) = ans(i)
    end do
    sparse = dense_to_csr(dense)

    ! Extract the diagonal from the sparse matrix
    diag = extract_diagonal(sparse)

    ! Test 1
    if (.not.assert(diag, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_extract_diagonal_csr_1 -1"
    end if

    ! Ensure the dense extract routine works as well
    densediag = extract_diagonal(dense)

    ! Test 2
    if (.not.assert(densediag, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_extract_diagonal_csr_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_msr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: a(4,4), dense(4,4)
    type(msr_matrix) :: msr

    ! Initialization
    rst = .true.
    a = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])

    ! Convert from dense to MSR
    msr = dense_to_msr(a)

    ! Convert from MSR to dense
    dense = msr

    ! Check
    if (.not.assert(dense, a)) then
        rst = .false.
        print "(A)", "Test Failed: test_msr_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_lu_factor_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer :: i, ju(4)
    integer, allocatable :: ipiv(:)
    real(real64) :: dense(4, 4), check(4, 4), x(4), b(4), xd(4)
    real(real64), allocatable :: dlu(:,:)
    type(csr_matrix) :: sparse
    type(msr_matrix) :: slu

    ! Initialization
    rst = .true.
    dense = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])
    sparse = dense
    call random_number(b)

    ! Compute the factorization of the sparse matrix
    call lu_factor(sparse, slu, ju)

    ! Compute the factorization of the dense matrix
    call lu_factor(dense, ipvt = ipiv, lu = dlu)

    ! Test - the diagonal must be inverted
    check = slu
    do i = 1, size(check, 1)
        check(i,i) = 1.0d0 / check(i,i)
    end do
    if (.not.assert(check, dlu)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_lu_factor_1 -1"
    end if

    ! Solve the sparse system
    x = solve_lu(slu, ju, b)

    ! Now solve the dense system for comparison
    xd = solve_lu(dlu, ipiv, b)

    ! Test
    if (.not.assert(x, xd)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_lu_factor_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_pgmres_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer :: ju(4)
    integer, allocatable :: ipiv(:)
    real(real64) :: dense(4, 4), m(4, 4), b(4), bc(4), x(4), xd(4)
    real(real64), allocatable :: dlu(:,:)
    type(csr_matrix) :: a, am
    type(msr_matrix) :: lu

    ! Initialization
    rst = .true.
    dense = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])
    m = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])
    a = dense
    am = m
    b = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
    bc = b

    ! Compute the preconditioner
    call lu_factor(am, lu, ju)

    ! Solve the sparse system
    x = pgmres_solver(a, lu, ju, b, tol = 1.0d-6, maxits = 500)

    ! Solve the dense system directly
    call lu_factor(dense, ipvt = ipiv, lu = dlu)
    xd = solve_lu(dlu, ipiv, bc)

    ! Test
    if (.not.assert(x, xd)) then
        rst = .false.
        print "(A)", "Test Failed: test_pgmres_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_csr_eigen_symm_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: n = 50
    integer, parameter :: k = 4
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    integer :: i
    real(real64) :: dense(n, n)
    real(real64), allocatable :: vals(:), vecs(:,:), dvals(:), r(:)
    type(csr_matrix) :: a

    ! Initialization
    rst = .true.
    dense = 0.0d0
    do i = 1, n
        dense(i,i) = 2.0d0 + real(i, real64)
        if (i > 1) then
            dense(i,i-1) = -1.0d0
            dense(i-1,i) = -1.0d0
        end if
    end do
    a = dense

    ! Compute the k largest magnitude eigenvalues of the sparse matrix
    call eigen(a, k, vals, vecs)

    ! Compute all eigenvalues of the dense matrix for comparison
    call eigen(dense, dvals)

    ! The ARPACK results are returned in ascending order of algebraic value
    if (.not.assert(vals, dvals(n-k+1:n), tol)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_eigen_symm_1 -1"
    end if

    ! Verify the eigenvectors satisfy A * v = lambda * v
    do i = 1, k
        r = matmul(a, vecs(:,i)) - vals(i) * vecs(:,i)
        if (norm2(r) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_symm_1 -2"
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_csr_eigen_asymm_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: n = 50
    integer, parameter :: k = 4
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    integer :: i
    real(real64) :: dense(n, n)
    complex(real64), allocatable :: vals(:), vecs(:,:), dvals(:), r(:)
    complex(real64) :: cvec(n)
    type(csr_matrix) :: a

    ! Initialization
    rst = .true.
    dense = 0.0d0
    do i = 1, n
        dense(i,i) = 2.0d0 + real(i, real64)
        if (i > 1) dense(i,i-1) = -2.0d0
        if (i < n) dense(i,i+1) = 1.0d0
    end do
    a = dense

    ! Compute the k largest magnitude eigenvalues of the sparse matrix
    call eigen(a, k, vals, vecs)

    ! Compute all eigenvalues of the dense matrix for comparison
    call eigen(dense, dvals)

    ! Ensure each computed eigenvalue appears in the dense solution
    do i = 1, size(vals)
        if (minval(abs(dvals - vals(i))) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_asymm_1 -1"
        end if
    end do

    ! Verify the eigenvectors satisfy A * v = lambda * v
    do i = 1, size(vals)
        cvec = vecs(:,i)
        r = matmul(cmplx(dense, 0.0d0, real64), cvec) - vals(i) * cvec
        if (norm2(abs(r)) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_asymm_1 -2"
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_csr_eigen_gen_symm_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: n = 30
    integer, parameter :: k = 3
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    integer :: i
    real(real64) :: da(n, n), db(n, n)
    real(real64), allocatable :: vals(:), vecs(:,:), beta(:), dvals(:), r(:)
    complex(real64), allocatable :: alpha(:)
    type(csr_matrix) :: a, b

    ! Initialization
    rst = .true.
    da = 0.0d0
    db = 0.0d0
    do i = 1, n
        da(i,i) = 2.0d0 + real(i, real64)
        db(i,i) = 2.0d0
        if (i > 1) then
            da(i,i-1) = -1.0d0
            da(i-1,i) = -1.0d0
            db(i,i-1) = 0.5d0
            db(i-1,i) = 0.5d0
        end if
    end do
    a = da
    b = db

    ! A zero shift exercises the shift-invert path and returns the eigenvalues
    ! nearest zero
    call eigen(a, b, k, vals, vecs, sigma = 0.0d0)

    ! Compute the dense solution for comparison
    call eigen(da, db, alpha, beta)
    dvals = real(alpha, real64) / beta

    ! Test
    do i = 1, size(vals)
        if (minval(abs(dvals - vals(i))) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_gen_symm_1 -1"
        end if

        r = matmul(a, vecs(:,i)) - vals(i) * matmul(b, vecs(:,i))
        if (norm2(r) / norm2(matmul(a, vecs(:,i))) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_gen_symm_1 -2"
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_csr_eigen_gen_asymm_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: n = 30
    integer, parameter :: k = 3
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    integer :: i
    real(real64) :: da(n, n), db(n, n)
    real(real64), allocatable :: beta(:)
    complex(real64), allocatable :: vals(:), vecs(:,:), alpha(:), dvals(:), r(:)
    complex(real64) :: cvec(n)
    type(csr_matrix) :: a, b

    ! Initialization
    rst = .true.
    da = 0.0d0
    db = 0.0d0
    do i = 1, n
        da(i,i) = 2.0d0 + real(i, real64)
        db(i,i) = 2.0d0
        if (i > 1) then
            da(i,i-1) = -2.0d0
            db(i,i-1) = 0.5d0
            db(i-1,i) = 0.5d0
        end if
        if (i < n) da(i,i+1) = 1.0d0
    end do
    a = da
    b = db

    ! No shift exercises the regular inverse path
    call eigen(a, b, k, vals, vecs)

    ! Compute the dense solution for comparison
    call eigen(da, db, alpha, beta)
    dvals = alpha / beta

    ! Test
    do i = 1, size(vals)
        if (minval(abs(dvals - vals(i))) > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_gen_asymm_1 -1"
        end if

        cvec = vecs(:,i)
        r = matmul(cmplx(da, 0.0d0, real64), cvec) - &
            vals(i) * matmul(cmplx(db, 0.0d0, real64), cvec)
        if (norm2(abs(r)) / norm2(abs(matmul(cmplx(da, 0.0d0, real64), cvec))) &
            > tol) then
            rst = .false.
            print "(A)", "Test Failed: test_csr_eigen_gen_asymm_1 -2"
        end if
    end do
end function

! ------------------------------------------------------------------------------
end module
