module test_sparse
    use iso_fortran_env
    use linalg
    use fortran_test_helper
    implicit none

contains
! ------------------------------------------------------------------------------
function test_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer(int32), parameter :: m = 4
    integer(int32), parameter :: n = 6
    integer(int32), parameter :: nnz = 8
    real(real64) :: dense(m, n), v(nnz), check(m, n)
    integer(int32) :: ja(nnz), ia(m + 1)
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.

    ! Construct a small sparse matrix, but in dense form
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
    integer(int32) :: ipiv(4)
    real(real64) :: a(4, 4), b(4), x(4), ans(4)
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
    ans = b
    call lu_factor(a, ipiv)
    call solve_lu(a, ipiv, ans)

    ! Test
    call sparse_direct_solve(sa, b, x)
    if (.not.assert(x, ans, tol = 1.0d-6)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_sparse_direct_solve_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_diag_to_csr_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer(int32), parameter :: n = 20
    integer(int32) :: i
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
    integer(int32), parameter :: kl = 4
    integer(int32), parameter :: ku = 5
    integer(int32), parameter :: m = 15
    integer(int32), parameter :: n = 15
    real(real64) :: banded(kl+ku+1,n), dense(m,n), check(m,n)
    type(csr_matrix) :: sparse

    ! Initialization
    rst = .true.
    call random_number(banded)

    ! Construct the dense matrix directly from the banded matrix
    call banded_to_dense(m, kl, ku, banded, dense)

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
    integer(int32), parameter :: n = 100
    integer(int32) :: i
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
    call extract_diagonal(sparse, diag)

    ! Test 1
    if (.not.assert(diag, ans)) then
        rst = .false.
        print "(A)", "Test Failed: test_extract_diagonal_csr_1 -1"
    end if

    ! Ensure the dense extract routine works as well
    call extract_diagonal(dense, densediag)

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
    integer(int32) :: i, ipiv(4), ju(4)
    real(real64) :: dense(4, 4), check(4, 4), x(4), b(4)
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
    call lu_factor(dense, ipiv)

    ! Test - the diagonal must be inverted
    check = slu
    do i = 1, size(check, 1)
        check(i,i) = 1.0d0 / check(i,i)
    end do
    if (.not.assert(check, dense)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_lu_factor_1 -1"
    end if

    ! Solve the sparse system
    call solve_lu(slu, ju, b, x)

    ! Now solve the dense system for comparison
    call solve_lu(dense, ipiv, b)

    ! Test
    if (.not.assert(x, b)) then
        rst = .false.
        print "(A)", "Test Failed: test_csr_lu_factor_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
end module