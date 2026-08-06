program example
    use iso_fortran_env
    use linalg
    implicit none

    ! Local Variables
    integer(int32), allocatable :: ipiv(:)
    real(real64) :: dense(4, 4), b(4), x(4), bc(4)
    real(real64), allocatable :: lu(:,:)
    type(csr_matrix) :: sparse

    ! Build the matrices as dense matrices
    dense = reshape([ &
        5.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 8.0d0, 0.0d0, 6.0d0, &
        0.0d0, 0.0d0, 3.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 5.0d0], [4, 4])
    b = [2.0d0, -1.5d0, 8.0d0, 1.0d0]

    ! Convert to sparse (CSR format)
    ! Note, the assignment operator is overloaded to allow conversion.
    sparse = dense

    ! Compute the solution to the sparse equations
    x = sparse_direct_solve(sparse, b)

    ! Print the solution
    print "(A)", "Sparse Solution:"
    print *, x

    ! Perform a sanity check on the solution
    ! Note, matmul is overloaded to allow multiplication with sparse matrices
    bc = matmul(sparse, x)
    print "(A)", "Computed RHS:"
    print *, bc
    print "(A)", "Original RHS:"
    print *, b

    ! For comparison, solve the dense system via LU decomposition
    call lu_factor(dense, ipvt = ipiv, lu = lu)
    x = solve_lu(lu, ipiv, b)
    print "(A)", "Dense Solution:"
    print *, x
end program