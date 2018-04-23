! linalg_lu_example_immutable.f90

program example
    use iso_fortran_env
    use linalg_immutable

    ! Local Variables
    type(lu_results) :: lu
    real(real64) :: a(3,3), b(3), pb(3), y(3), x(3)

    ! Build the 3-by-3 matrix A
    !     | 1   2   3 |
    ! A = | 4   5   6 |
    !     | 7   8   0 |
    a = reshape( &
        [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
        [3, 3])

    ! Compute the LU factorization
    lu = mat_lu(a)

    ! Build the right-hand-side vector B.
    !     | -1 |
    ! b = | -2 |
    !     | -3 |
    b = [-1.0d0, -2.0d0, -3.0d0]

    ! Apply the row pivots (P * B)
    pb = matmul(lu%p, b)

    ! Compute the solution to the lower triangular problem L Y = P B, for Y
    y = mat_solve_lower_tri(lu%l, pb)

    ! Compute the solution to the upper triangular problem U X = Y
    x = mat_solve_upper_tri(lu%u, y)

    ! The solution is:
    !     |  1/3 |
    ! x = | -2/3 |
    !     |   0  |
    print '(A)', "LU Solution: X = "
    print '(F8.4)', (x(i), i = 1, size(x))
end program
