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
## Matrix-Free Iterative Solves
`gmres_solve` wraps SPARSKIT's reverse-communication restarted GMRES. It accepts
any real `linear_operator` implementing read-only `size(this, dim)` and
`apply(this, x, y)` bindings; it does not require assembled matrix entries or
a transpose.

```fortran
program operator_solution_example
    use iso_fortran_env, only : real64, int32
    use linalg
    use linalg_errors, only : LA_NO_ERROR
    implicit none

    type(csr_matrix) :: coefficients
    type(msr_matrix) :: factors
    type(csr_operator) :: op
    type(ilu_preconditioner) :: preconditioner
    type(iterative_solver_info) :: details
    integer(int32) :: upper_rows(3)
    integer :: status
    real(real64) :: rhs(3)
    real(real64), allocatable :: solution(:)

    coefficients = dense_to_csr(reshape([4.0d0, -1.0d0, 0.0d0, &
        -1.0d0, 4.0d0, -1.0d0, 0.0d0, -1.0d0, 4.0d0], [3,3]))
    rhs = [2.0d0, 4.0d0, 10.0d0]
    call lu_factor(coefficients, factors, upper_rows)
    op = csr_operator(coefficients)
    preconditioner = ilu_preconditioner(factors, upper_rows)
    call gmres_solve(op, rhs, solution, preconditioner = preconditioner, &
        tol = 1.0d-10, flag = status, info = details)
    if (status /= LA_NO_ERROR) error stop status
    print '(3F12.6)', solution
end program
```

The solution is `[1, 2, 3]`. The CSR conversion is only for this small example.
A custom matrix-free operator can replace `op`, including when the ILUT factors
come from a separately assembled approximation to the original operator.


| Argument | Contract |
| --- | --- |
| `a` | Required fixed real N-by-N `linear_operator`, N positive. `apply` assigns every output element and does not change the operator. |
| `b` | Required finite N-element `real(real64)` RHS, preserved. |
| `x` | Required allocatable output vector. N-element solution on success, latest backend iterate on iteration/callback failure, unallocated on input/allocation errors. Check `flag` before treating an iterate as a solution. |
| `preconditioner` | Optional fixed N-by-N `linear_operator` applying the right approximate inverse: `z = M**(-1)*r`, not `M*r`. Omit for no preconditioning. Adaptive/variable preconditioning requires a flexible method, not this driver. |
| `tol` | Optional finite relative residual tolerance in `[0,1)`, default `sqrt(epsilon(0.0d0))`. |
| `atol` | Optional finite nonnegative absolute residual tolerance, default zero. `tol` and `atol` cannot both be zero. Success requires `norm2(b-A*x) <= tol*norm2(b)+atol`. |
| `im` | Optional restart length from 2 through `max(2,N)`, default `max(2,min(N,30))`. This backend has no restart limit of 50. Length 1 is rejected because SPARSKIT silently replaces it with 15; scalar systems use length 2. |
| `maxits` | Optional positive budget for backend-requested matrix products, default 300. Includes initial/restart residual products, but not the final wrapper residual check or preconditioner applications. It is not a count of Krylov steps. |
| `x0` | Optional finite N-element initial guess, preserved; omitted means zero. |
| `flag` | Optional library status. `LA_NO_ERROR` requires a verified original-system residual. Invalid inputs, shape errors, workspace overflow/allocation failure, invalid callbacks, or nonconvergence return the documented `LA_*` errors. Omitting `flag` makes errors terminate with `error stop`, even when `info` is supplied. |
| `info` | Optional `iterative_solver_info` output, reset and populated on every return. |

Diagnostics contain `status` (the same library code as `flag`), `backend_status`
(SPARSKIT's raw termination code, or -1 for a wrapper-enforced budget),
`matvecs` (all actual A applications, including final verification),
`preconditioner_applications`, `residual_norm`, `target_residual`, and `converged`.
The residual is independently checked even after an Arnoldi breakdown, which
can indicate an exact solution. A raw nonzero backend code therefore does not
alone imply failure; use `status` and `converged`. When no residual can be
evaluated safely, `residual_norm` is initialized to `huge(0.0d0)`.

`csr_operator(a)` validates and owns a deep copy of a nonempty finite CSR matrix;
its `size` and `apply` bindings provide matrix products. Rectangular operators
can be constructed, but GMRES rejects them. `ilu_preconditioner(lu, ju)` validates
and owns copies of prepared MSR factors and upper-row pointers. Its `apply`
binding performs LUSOL only, not an inner iterative solve. Both adapters require
correctly sized, nonoverlapping input/output vectors; invalid construction or
direct binding calls use `error stop`. Default-constructed adapters have size
zero until initialized. Changing or releasing source arrays after construction
does not change the adapters.

Callbacks must remain fixed linear maps and must handle their own failures:
`linear_operator%apply` has no status argument. The wrapper rejects nonfinite
callback inputs and outputs, checks workspace slices for bounds and overlap,
and checks 32-bit backend workspace arithmetic before allocating. Krylov memory
scales as `O(N*im + im**2)`.

**Concurrency:** the bundled reverse-communication GMRES contains saved backend
state. Complete solves must be serialized, and callbacks must not recursively
invoke this backend. Read-only operators and separate caller workspaces do not
make this particular SPARSKIT implementation reentrant.

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
