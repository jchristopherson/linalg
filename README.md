# linalg
A linear algebra library that provides a user-friendly interface to several BLAS and LAPACK routines.  The examples below provide an illustration of just how simple it is to perform a few common linear algebra operations.

## Status
![Build Status](https://github.com/jchristopherson/linalg/actions/workflows/cmake.yml/badge.svg)
[![Actions Status](https://github.com/jchristopherson/linalg/workflows/fpm/badge.svg)](https://github.com/jchristopherson/linalg/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/linalg/).

## Building LINALG
### CMake
This library can be built using [CMake](https://cmake.org/).  For general instructions see [Running CMake](https://cmake.org/runningcmake/).
```txt
cmake -B build -S .
cmake --build build
```
CMake will search for existing installations of BLAS, LAPACK, and ARPACK (via `find_package`).  If any of these are not found on the system, CMake will automatically download and build reference versions of the missing libraries as part of the build, so no additional setup is required to get a working build.

To install the library, run the following from the build directory.
```txt
cmake --install build
```

### FPM
[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```txt
fpm build
```
**Note:** Unlike the CMake build, FPM does not download or build missing dependencies.  BLAS, LAPACK, and [ARPACK](https://github.com/opencollab/arpack-ng) must already be installed on your system, with their libraries discoverable by the linker (e.g. `libopenblas-dev`, `liblapack-dev`, and `libarpack2-dev` on Debian/Ubuntu, or `mingw-w64-x86_64-openblas` and `mingw-w64-x86_64-arpack` via MSYS2 on Windows), before running `fpm build` or `fpm test`.

The LINALG library can be used within your FPM project by adding the following to your fpm.toml file.
```toml
[dependencies]
linalg = { git = "https://github.com/jchristopherson/linalg" }
```

## Standard Solution Example
This example solves a normally defined system of 3 equations of 3 unknowns.

```fortran
program example
    use iso_fortran_env, only : real64, int32
    use linalg
    implicit none

    ! Local Variables
    real(real64) :: a(3,3), b(3), x(3)
    integer(int32) :: i
    type(lu_factorization) :: lu

    ! Build the 3-by-3 matrix A.
    !     | 1   2   3 |
    ! A = | 4   5   6 |
    !     | 7   8   0 |
    a = reshape( &
        [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
        [3, 3])
    
    ! Build the right-hand-side vector B.
    !     | -1 |
    ! b = | -2 |
    !     | -3 |
    b = [-1.0d0, -2.0d0, -3.0d0]

    ! The solution is:
    !     |  1/3 |
    ! x = | -2/3 |
    !     |   0  |

    ! Compute the LU factorization and solve
    lu = lu_factorization(a)
    call lu%factor()
    x = lu%solve(b)

    ! Display the results.
    print '(A)', "LU Solution: X = "
    print '(F8.4)', (x(i), i = 1, size(x))
end program
```
The above program produces the following output.
```text
LU Solution: X =
  0.3333
 -0.6667
  0.0000
```

## Sparse System Solver
`sparse_lu_factorization` mirrors the dense type's construction, `update`,
`factor`, `solve`, and `solve_in_place` workflow. It owns a copy of the original
CSR matrix and prepares incomplete LU factors as a preconditioner. Its solve
methods then use restarted GMRES to solve the **original system**, not just the
incomplete triangular factors. This is an iterative solve to a requested
tolerance, not a pivoted sparse direct LU factorization.

```fortran
program sparse_solution_example
    use iso_fortran_env, only : real64
    use linalg
    use linalg_errors, only : LA_NO_ERROR
    implicit none

    real(real64) :: coefficients(3,3), rhs(3)
    real(real64), allocatable :: solution(:)
    type(csr_matrix) :: sparse_matrix
    type(sparse_lu_factorization) :: factorization
    integer :: status

    coefficients = reshape([4.0d0, -1.0d0, 0.0d0, &
        -1.0d0, 4.0d0, -1.0d0, 0.0d0, -1.0d0, 4.0d0], [3,3])
    rhs = [2.0d0, 4.0d0, 10.0d0]
    sparse_matrix = dense_to_csr(coefficients)
    factorization = sparse_lu_factorization(sparse_matrix, tol = 1.0d-10)

    call factorization%factor(status)
    if (status /= LA_NO_ERROR) error stop status
    solution = factorization%solve(rhs, status)
    if (status /= LA_NO_ERROR) error stop status
    print '(3F12.6)', solution
end program
```

The solution is `[1, 2, 3]`. Dense-to-CSR conversion is used only to keep this
small example readable; the constructor accepts an existing CSR matrix without
forming a dense copy.

The constructor and `update` accept the same optional controls:

| Argument | Meaning |
| --- | --- |
| `a` | Required real CSR coefficient matrix. Factorization requires nonempty square storage with unique column indices within each row. The caller's matrix is copied, never modified. |
| `droptol` | Finite nonnegative ILUT drop threshold, default `sqrt(epsilon(0.0d0))`. Controls preconditioner quality, not solution accuracy. ILUT also limits per-row fill; zero threshold does not guarantee exact factors. |
| `tol` | Finite relative residual tolerance at least machine epsilon, default `sqrt(epsilon(0.0d0))`. Success requires `norm2(b - matmul(a, x)) <= tol * norm2(b)` for each RHS. A zero RHS returns zero. |
| `maxits` | Positive GMRES iteration budget, default 100. The bundled solver checks the budget at restart boundaries and may perform up to `im-1` additional iterations. |
| `im` | Restart length from 1 through `min(N, 50)`, default `min(N, 50)`. The bundled GMRES implementation has a hard maximum of 50. |

Construction and `update` only store the matrix and controls; validation occurs
in `factor` or the first solve. `update` discards cached factors and resets
omitted controls to their defaults. `is_factored()` queries preconditioner
readiness, and `size(dim)` returns matrix dimensions, or zero before initialization.

`solve(b, flag)` allocates a vector or matrix with the same shape as `b`;
`solve_in_place(b, flag)` replaces caller-owned RHS storage. Both prepare factors
lazily and solve each RHS column independently. If any column fails, the whole
result remains an unchanged copy of the RHS. The optional `flag` receives
`LA_NO_ERROR` on success or an initialization, shape, CSR-format, control,
factorization, or convergence error. Omitting it makes errors use `error stop`.

`get_factored_form(lu, ju, flag)` returns owned copies of the incomplete MSR
factors and upper-row pointers. These are preconditioner data, not exact factors
of the original matrix. Determinant and inverse methods are intentionally absent.

The legacy `sparse_direct_solve(a, b, droptol)` name is retained for compatibility,
but it now refines the ILUT-preconditioned solution using GMRES against the
original matrix. It is not a true sparse direct solver. It uses GMRES's default
controls; use `pgmres_solve` or `sparse_lu_factorization` to configure convergence.
By contrast, `solve_lu(lu, ju, b)` only applies the supplied triangular factors.

`validate_csr_storage(a)` and `validate_msr_lu(lu, ju)` return library status
codes for malformed storage before raw SPARSKIT calls. GMRES validates both
formats, dimensions, finite inputs, and the 50-vector restart limit. Factor
application validates MSR row boundaries, upper-row splits, and triangular
column references; nonfinite or zero inverse pivots propagate as NaN results.
The unused MSR separator value does not affect validation or solution.

### Matrix-Free Integration
SPARSKIT's reverse-communication Krylov routines can request `A*x` from a
`linear_operator%apply` callback and preconditioner actions from another
operator representing `M**(-1)`. This supports a matrix-free original operator
with an assembled CSR approximation used only to build a preconditioner.
The current `pgmres_solver` function and status-returning `pgmres_solve`
subroutine are CSR-specific; accepting arbitrary operators requires a separate
reverse-communication driver, not a change to the CSR storage type.

For prepared full-system solves, `solve_prepared(b, x, flag)` takes the
factorization object as `intent(in)` and writes into caller-owned `x`. Call
`factor` first: this binding never initializes or changes the object. Its input
and output vectors must have N elements and must not overlap. On failure it
zeros `x` and reports a status, or stops if `flag` is omitted. It is suitable
for a `linear_solver` adapter used by ARPACK, with an inner tolerance chosen
for the eigenproblem. An incomplete-factor application alone is not a suitable
substitute for an accurate shift-invert solve.

Mutable lazy-setup calls must not run concurrently on the same object. Prepared
calls use local workspace and read-only coefficients and factors; concurrent
use still depends on the reentrancy of linked SPARSKIT/BLAS routines. The
`matrix_factorization` and `linear_solver` abstractions have different solve
contracts, so an adapter is preferable to making one inherit from the other.

## Overdetermined System Example
This example solves an overdetermined system of 3 equations of 2 uknowns.

```fortran
program example
    use iso_fortran_env
    use linalg
    implicit none

    ! Local Variables
    real(real64) :: a(3,2), b(3), x(2)
    integer(int32) :: i

    ! Build the 3-by-2 matrix A
    !     | 2   1 |
    ! A = |-3   1 |
    !     |-1   1 |
    a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])

    ! Build the right-hand-side vector B.
    !     |-1 |
    ! b = |-2 |
    !     | 1 |
    b = [-1.0d0, -2.0d0, 1.0d0]

    ! The solution is:
    ! x = [0.13158, -0.57895]**T

    ! Compute the solution via a least-squares approach.
    x = solve_least_squares(a, b)

    ! Display the results
    print '(A)', "Least Squares Solution: X = "
    print '(F9.5)', (x(i), i = 1, size(x))
end program
```
The above program produces the following output.
```text
Least Squares Solution: X =
  0.13158
 -0.57895
```

## Eigen Analysis Example
This example computes the eigenvalues and eigenvectors of a mechanical system consisting of several masses connected by springs.

```fortran
! This is an example illustrating the use of the eigenvalue and eigenvector
! routines to solve a free vibration problem of 3 masses connected by springs.
!
!     k1           k2           k3           k4
! |-\/\/\-| m1 |-\/\/\-| m2 |-\/\/\-| m3 |-\/\/\-|
!
! As illustrated above, the system consists of 3 masses connected by springs.
! Spring k1 and spring k4 connect the end masses to ground.  The equations of
! motion for this system are as follows.
!
! | m1  0   0 | |x1"|   | k1+k2  -k2      0  | |x1|   |0|
! | 0   m2  0 | |x2"| + |  -k2  k2+k3    -k3 | |x2| = |0|
! | 0   0   m3| |x3"|   |   0    -k3    k3+k4| |x3|   |0|
!
! Notice: x1" = the second time derivative of x1.
program example
    use iso_fortran_env
    use linalg
    implicit none

    ! Define the model parameters
    real(dp), parameter :: pi = 3.14159265359d0
    real(dp), parameter :: m1 = 0.5d0
    real(dp), parameter :: m2 = 2.5d0
    real(dp), parameter :: m3 = 0.75d0
    real(dp), parameter :: k1 = 5.0d6
    real(dp), parameter :: k2 = 10.0d6
    real(dp), parameter :: k3 = 10.0d6
    real(dp), parameter :: k4 = 5.0d6

    ! Local Variables
    integer(int32) :: i, j
    real(real64) :: m(3,3), k(3,3), natFreq(3)
    complex(real64), allocatable :: vals(:), modeShapes(:,:)

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the eigenvalues and eigenvectors.
    call eigen(k, m, vals, rvecs = modeShapes)

    ! Sort the eigenvalues and eigenvectors
    call sort(vals, modeShapes)

    ! Compute the natural frequency values, and return them with units of Hz.  
    ! Notice, all eigenvalues and eigenvectors are real for this example.
    natFreq = sqrt(real(vals)) / (2.0d0 * pi)

    ! Display the natural frequency and mode shape values.
    print '(A)', "Modal Information:"
    do i = 1, size(natFreq)
        print '(AI0AF8.4A)', "Mode ", i, ": (", natFreq(i), " Hz)"
        print '(F10.3)', (real(modeShapes(j,i)), j = 1, size(natFreq))
    end do
end program
```
The above program produces the following output.
```text
Modal Information:
Mode 1: (232.9225 Hz)
    -0.718
    -1.000
    -0.747
Mode 2: (749.6189 Hz)
    -0.419
    -0.164
     1.000
Mode 3: (923.5669 Hz)
     1.000
    -0.184
     0.179
```

## Sparse Matrix Example
The following example solves a sparse system of equations using a direct solver.  The solution is compared to the solution of the same system of equations but in dense format for comparison.
```fortran
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
    x = solve_lu(dense, ipiv, b)
    print "(A)", "Dense Solution:"
    print *, x
end program
```
The above program produces the following output.
```text
Sparse Solution:
  0.40000000000000002      -0.18750000000000000        2.6666666666666665       0.42500000000000004     
Computed RHS:
   2.0000000000000000       -1.5000000000000000        8.0000000000000000        1.0000000000000000
Original RHS:
   2.0000000000000000       -1.5000000000000000        8.0000000000000000        1.0000000000000000
Dense Solution:
  0.40000000000000002      -0.18750000000000000        2.6666666666666665       0.42499999999999999
```

## External Libraries
Here is a list of external code libraries utilized by this library.
- [BLAS](http://www.netlib.org/blas/)
- [LAPACK](http://www.netlib.org/lapack/)
- [ARPACK](https://github.com/opencollab/arpack-ng)
- [QRUpdate](https://sourceforge.net/projects/qrupdate/)
- [SPARSKIT](https://www-users.cse.umn.edu/~saad/software/SPARSKIT/)

When building with CMake, these dependencies do not necessarily have to be pre-installed.  The build will initially look for installed items, but if not found, will then download and build the latest version as part of the build process.  When building with FPM, however, BLAS, LAPACK, and ARPACK must already be installed on the system as FPM does not build missing dependencies (see [Building LINALG](#building-linalg)).

## License

This project is licensed under the **MIT License**.

See the [LICENSE](LICENSE) file for details.
