! linalg_od_example.f90

! Example Source: https://en.wikipedia.org/wiki/Overdetermined_system
program example
    use iso_fortran_env, only : real64, int32
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
    x = solve_least_squares_svd(a, b)

    ! Display the results
    print '(A)', "Least Squares Solution: X = "
    print '(F9.5)', (x(i), i = 1, size(x))
end program