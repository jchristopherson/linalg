! linalg_sparse_direct_example.f90
! SPDX-License-Identifier: MIT
program example
    use iso_fortran_env
    use linalg
    implicit none

    ! Local Variables
    real(real64) :: dense(4, 4), b(4), x(4), bc(4)
    type(csr_matrix) :: sparse
    type(sparse_lu_factorization) :: sparse_solver
    type(lu_factorization) :: dense_solver

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
    sparse_solver = sparse_lu_factorization(sparse)
    call sparse_solver%factor()
    x = sparse_solver%solve(b)

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
    dense_solver = lu_factorization(dense)
    call dense_solver%factor()
    x = dense_solver%solve(b)
    print "(A)", "Dense Solution:"
    print *, x
end program
