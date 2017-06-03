# linalg
A linear algebra library that provides a user-friendly interface to several BLAS and LAPACK routines.


## Example 1
This example solves a normally defined system of N equations of N unknowns.

```fortran
    use linalg_constants, only : dp, i32
    use linalg_factor, only : lu_factor
    use linalg_solve, only : solve_lu

    real(dp), dimension(20, 20) :: a
    real(dp), dimension(20) :: b
    integer(i32), dimension(20) :: ipvt

    ! Initialize A and B ...

    ! Solve the problem: A * X = B by first computing the LU factorization of matrix A
    call lu_factor(a, ipvt)

    ! Notice, the factored form of A overwrites the original matrix A

    ! Now, solve for X as follows (X overwrites B)
    call solve_lu(a, ipvt, b)

```

## Example 2
This example solves an overdefined system of M equations of N uknowns.

```fortran
    use linalg_constants, only : dp, i32
    use linalg_solve, only : least_squares_solve

    integer(i32), parameter :: m = 40
    integer(i32), parameter :: n = 30
    real(dp), dimension(m, n) :: a
    real(dp), dimension(m) :: b
    real(dp), dimension(n) :: x

    ! Initialize A and B ...
    
    ! Compute the solution using a least-squares approach
    call least_squares_solve(a, b)

    ! Notice, X overwrites the first N elements in B, so if you'd prefer:
    x = b(1:n)

```

## Example 3
This example computes the eigenvalues and eigenvectors of a mechanical system consisting of several masses connected by springs.

```fortran
    use linalg_constants, only : dp
    use linalg_eigen, only : eigen

    real(dp), parameter :: pi = 3.14159265359d0

    real(dp), dimension(20, 20) :: m, k
    real(dp), dimension(20) :: natFreq
    complex(dp), dimension(20, 20) :: modeShapes
    complex(dp), dimension(n) :: eigVals

    ! Initialize the mass matrix (M) and stiffness matrix (K) ...
    
    ! Solve the eigenvalue problem.
    call eigen(k, m, vals, vecs = modeShapes)

    ! The eigenvalues represent the square of the system natural frequencies (Hz).
    natFreq = sqrt(real(eigVals)) / (2.0d0 * pi)

    ! Note: It was assumed that the system was a normally constrained system such that all eigenvalues are positive and real-valued.

```