! linalg_cholesky_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg
    implicit none

    ! Variables
    real(real64) :: a(3, 3), r(3, 3), b(3), x(3), y(3)
    integer(int32) :: i

    ! Build the 3-by-3 positive-definite matrix A.
    !     | 4   12   -16 |
    ! A = | 12  37   -43 |
    !     |-16 -43    98 |
    a = reshape([4.0d0, 12.0d0, -16.0d0, 12.0d0, 37.0d0, -43.0d0, -16.0d0, &
        -43.0d0, 98.0d0], [3, 3])
    
    ! Build the 3-element array B
    !     | 5 |
    ! b = | 1 |
    !     | 3 |
    b = [5.0d0, 1.0d0, 3.0d0]

    ! Compute the Cholesky factorization of A considering only the upper 
    ! triangular portion of A (the default configuration).
    r = cholesky_factor(a)

    ! Compute the solution
    x = solve_cholesky(.true., r, b)

    ! Display the results
    print '(A)', "Cholesky Solution: X = "
    print '(F8.4)', (x(i), i = 1, size(x))

    ! The solution could also be computed manually noting the Cholesky 
    ! factorization causes A = R**T * R.  Then R**T * R * X = B.  
    
    ! Step 1 would then be to solve the problem R**T * Y = B, for Y.
    y = solve_triangular_system(.true., .true., .true., r, b)

    ! Now, solve the problem R * X = Y, for X
    x = solve_triangular_system(.true., .false., .true., r, y)

    ! Display the results
    print '(A)', "Cholesky Solution (Manual Approach): X = "
    print '(F8.4)', (x(i), i = 1, size(x))
end program