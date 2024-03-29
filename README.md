# linalg
A linear algebra library that provides a user-friendly interface to several BLAS and LAPACK routines.  The examples below provide an illustration of just how simple it is to perform a few common linear algebra operations.  There is also an optional C API that is available as part of this library.

## Status
![Build Status](https://github.com/jchristopherson/linalg/actions/workflows/cmake.yml/badge.svg)
[![Actions Status](https://github.com/jchristopherson/linalg/workflows/fpm/badge.svg)](https://github.com/jchristopherson/linalg/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/linalg/).

## Building LINALG
[CMake](https://cmake.org/)This library can be built using CMake.  For instructions see [Running CMake](https://cmake.org/runningcmake/).

[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```txt
fpm build
```
The LINALG library can be used within your FPM project by adding the following to your fpm.toml file.
```toml
[dependencies]
linalg = { git = "https://github.com/jchristopherson/linalg" }
```

## Standard Solution Example
This example solves a normally defined system of 3 equations of 3 unknowns.

```fortran
program example
    use iso_fortran_env
    use linalg
    implicit none

    ! Local Variables
    real(dp) :: a(3,3), b(3)
    integer(i32) :: i, pvt(3)

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

    ! Compute the LU factorization
    call lu_factor(a, pvt)

    ! Compute the solution.  The results overwrite b.
    call solve_lu(a, pvt, b)

    ! Display the results.
    print '(A)', "LU Solution: X = "
    print '(F8.4)', (b(i), i = 1, size(b))
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
    real(dp) :: a(3,2), b(3)
    integer(i32) :: i

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

    ! Compute the solution via a least-squares approach.  The results overwrite
    ! the first 2 elements in b.
    call solve_least_squares(a, b)

    ! Display the results
    print '(A)', "Least Squares Solution: X = "
    print '(F9.5)', (b(i), i = 1, size(a, 2))
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
    integer(i32) :: i, j
    real(dp) :: m(3,3), k(3,3), natFreq(3)
    complex(dp) :: vals(3), modeShapes(3,3)

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the eigenvalues and eigenvectors.
    call eigen(k, m, vals, vecs = modeShapes)

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
    integer(int32) :: ipiv(4)
    real(real64) :: dense(4, 4), b(4), x(4), bc(4)
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
    call sparse_direct_solve(sparse, b, x)  ! Results stored in x

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
    call lu_factor(dense, ipiv)
    call solve_lu(dense, ipiv, b)   ! Results stored in b
    print "(A)", "Dense Solution:"
    print *, b
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
- [QRUpdate](https://sourceforge.net/projects/qrupdate/)
- [FERROR](https://github.com/jchristopherson/ferror)
- [SPARSKIT](https://www-users.cse.umn.edu/~saad/software/SPARSKIT/)

The dependencies do not necessarily have to be installed to be used.  The build will initially look for installed items, but if not found, will then download and build the latest version as part of the build process.
