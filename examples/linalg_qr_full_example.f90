! linalg_qr_full_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg_factor, only : qr_factor, form_qr
    use linalg_solve, only : solve_triangular_system
    implicit none

    ! Variables
    real(real64) :: a(3,3), b(3), q(3,3), tau(3)
    integer(int32) :: i

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

    ! Compute the QR factorization without column pivoting
    call qr_factor(a, tau)

    ! Build Q and R.  A is overwritten with R
    call form_qr(a, tau, q)

    ! As this system is square, matrix R is upper triangular.  Also, Q is
    ! always orthogonal such that it's inverse and transpose are equal.  As the
    ! system is now factored, its form is: Q * R * X = B.  Solving this system
    ! is then as simple as solving the upper triangular system: 
    ! R * X = Q**T * B.

    ! Compute Q**T * B, and store the results in B
    b = matmul(transpose(q), b)

    ! Solve the upper triangular system R * X = Q**T * B for X
    call solve_triangular_system(.true., .false., .true., a, b)

    ! Display the results
    print '(A)', "QR Solution: X = "
    print '(F8.4)', (b(i), i = 1, size(b))

    ! Notice, QR factorization with column pivoting could be accomplished via
    ! a similar approach, but the column pivoting would need to be accounted
    ! for by noting that Q * R = A * P, where P is an N-by-N matrix describing
    ! the column pivoting operations.
end program
