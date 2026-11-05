# linalg
A linear algebra library that provides a user-friendly interface to several BLAS and LAPACK routines.  The examples below provide an illustration of just how simple it is to perform a few common linear algebra operations.

## Status
![Build Status](https://github.com/jchristopherson/linalg/actions/workflows/cmake.yml/badge.svg)
[![Actions Status](https://github.com/jchristopherson/linalg/workflows/fpm/badge.svg)](https://github.com/jchristopherson/linalg/actions)

## Capabilities

| Area | Capabilities |
| --- | --- |
| Dense factorizations | LU, Cholesky, QR, LQ, and RZ factorization routines |
| Linear systems | Direct solves, least-squares solutions, matrix inversion, and Moore-Penrose pseudoinverses |
| Spectral analysis | Singular value decomposition, and standard and generalized eigenvalue problems |
| Sparse linear algebra | CSR and MSR matrices, sparse direct solves, and sparse matrix operations |
| Iterative solves | Restarted GMRES for user-defined linear operators, with optional ILU preconditioning |
| Factorization updates | QR and Cholesky updates and downdates |
| Build options | CMake or FPM; CMake can fetch missing external dependencies |

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/linalg/).

## Getting Started

Build and install the library with [CMake](#cmake), or build it with
[FPM](#fpm). For FPM projects, add LINALG as a dependency in `fpm.toml` as
shown below. CMake can fetch missing external dependencies; with FPM, install
BLAS, LAPACK, and ARPACK separately.

Once LINALG is available to your project, import its public API with
`use linalg`. For a first solve, see the [standard solution example](#standard-solution-example)
or its [complete source](examples/linalg_lu_example.f90).

## API Overview

The `linalg` module re-exports the library's dense and sparse APIs. Many
operations are provided as procedures (for example, `solve_least_squares_svd`),
while factorizations can also be represented by types that expose factor and
solve methods (for example, `lu_factorization`). This lets you factor a matrix
once and reuse the factorization for multiple right-hand sides. The examples
below show both styles; see the [API documentation](https://jchristopherson.github.io/linalg/)
for argument details and additional routines.

### Procedural and type-bound factorization APIs

For dense LU, both styles can factor a matrix once and solve one or more
right-hand sides. The procedural API returns the factors and pivot array for
you to pass to `solve_lu`:

```fortran
call lu_factor(a, ipvt = ipvt, lu = factors)
x = solve_lu(factors, ipvt, b)
```

The type-bound API keeps the factorization state in an object and exposes
operations as methods:

```fortran
solver = lu_factorization(a)
call solver%factor()
x = solver%solve(b)
```

Choose the procedural form when you want explicit control over the factor
arrays or need to pass them between procedures. Choose the type-bound form
when you prefer the matrix and its factorization state to live together; it
also offers methods such as `update`, `determinant`, and `inverse`. Both forms
support solving multiple right-hand sides.

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
This example solves a normally defined system of N equations of N unknowns.

Full example: [examples/linalg_lu_example.f90](examples/linalg_lu_example.f90).

```fortran
type(lu_factorization) :: lu
! a is the coefficient matrix and b is the right-hand side.
lu = lu_factorization(a)
call lu%factor()
x = lu%solve(b)
```

## Overdetermined System Example
This example solves an overdetermined system of M equations of N unknowns.

Full example: [examples/linalg_od_example.f90](examples/linalg_od_example.f90).

```fortran
! a is an overdetermined matrix and b is the right-hand side.
x = solve_least_squares_svd(a, b)
```

## Eigen Analysis Example
This example computes the eigenvalues and eigenvectors of a mechanical system consisting of several masses connected by springs.

Full example: [examples/linalg_eigen_example.f90](examples/linalg_eigen_example.f90).

```fortran
! k and m are the stiffness and mass matrices.
call eigen(k, m, vals, rvecs = modeShapes)
call sort(vals, modeShapes)
natFreq = sqrt(real(vals)) / (2.0d0 * pi)
```

## Sparse Matrix Example
The following example solves a sparse system of equations using a direct solver.  The solution is compared to the solution of the same system of equations but in dense format for comparison.

Full example: [examples/linalg_sparse_example.f90](examples/linalg_sparse_example.f90).

```fortran
! sparse is a CSR matrix and b is the right-hand side.
sparse_solver = sparse_lu_factorization(sparse)
call sparse_solver%factor()
x = sparse_solver%solve(b)
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
    integer :: upper_rows(3)
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
