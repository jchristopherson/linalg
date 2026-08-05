! linalg_qr_mult_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg
    implicit none

    ! Variables
    real(real64) :: a(3,3), b(3), qtb(3), x(3)
    real(real64), allocatable :: qr(:,:), tau(:)
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
    call qr_factor(a, tau = tau, qr = qr)

    ! As this system is square, matrix R is upper triangular.  Also, Q is
    ! always orthogonal such that it's inverse and transpose are equal.  As the
    ! system is now factored, its form is: Q * R * X = B.  Solving this system
    ! is then as simple as solving the upper triangular system: 
    ! R * X = Q**T * B.

    ! Compute Q**T * B, and store the results in B.  Notice, using mult_qr
    ! avoids direct construction of the full Q and R matrices.
    qtb = mult_qr(.true., qr, tau, b)

    ! Solve the upper triangular system R * X = Q**T * B for X
    x = solve_triangular_system(.true., .false., .true., qr, qtb)

    ! Display the results
    print '(A)', "QR Solution: X = "
    print '(F8.4)', (x(i), i = 1, size(x))

    ! Notice, QR factorization with column pivoting could be accomplished via
    ! a similar approach, but the column pivoting would need to be accounted
    ! for by noting that Q * R = A * P, where P is an N-by-N matrix describing
    ! the column pivoting operations.
end program
