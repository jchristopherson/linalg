! linalg_lu_full_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg
    implicit none

    ! Variables
    real(real64) :: a(3,3), b(3), y(3), x(3), pb(3)
    real(real64), allocatable :: l(:,:), u(:,:), p(:,:)
    integer(int32) :: i, pvt(3)

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
    call lu_factor(a, l = l, u = u, p = p)

    ! Solve the lower triangular system L * Y = P * B for Y
    pb = matmul(p, b)

    ! Now, compute the solution to the lower triangular system.  Store the
    ! result in B.  Remember, L is unit diagonal (ones on its diagonal)
    y = solve_triangular_system(.false., .false., .false., l, pb)
    
    ! Solve the upper triangular system U * X = Y for X.
    x = solve_triangular_system(.true., .false., .true., u, y)

    ! Display the results.
    print '(A)', "LU Solution: X = "
    print '(F8.4)', (x(i), i = 1, size(x))
end program
