! linalg_od_example.f90

! Example Source: https://en.wikipedia.org/wiki/Overdetermined_system
program example
    use iso_fortran_env, only : real64, int32
    use linalg_solve, only : solve_least_squares_full
    implicit none

    ! Local Variables
    real(real64) :: a(3,2), b(3)
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

    ! Compute the solution via a least-squares approach.  The results overwrite
    ! the first 2 elements in b.
    call solve_least_squares_full(a, b)

    ! Display the results
    print '(A)', "Least Squares Solution: X = "
    print '(F9.5)', (b(i), i = 1, size(a, 2))
end program