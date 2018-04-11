! linalg_cholesky_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg_factor, only : cholesky_factor
    use linalg_solve, only : solve_cholesky, solve_triangular_system
    implicit none

    ! Variables
    real(real64) :: a(3, 3), b(3), bu(3)
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

    ! Make a copy of B for later use - not necessary, but just for example to
    ! illustrate the long or manual method of solving a Cholesky factored system
    bu = b

    ! Compute the Cholesky factorization of A considering only the upper 
    ! triangular portion of A (the default configuration).
    call cholesky_factor(a)

    ! Compute the solution
    call solve_cholesky(.true., a, b)

    ! Display the results
    print '(A)', "Cholesky Solution: X = "
    print '(F8.4)', (b(i), i = 1, size(b))

    ! The solution could also be computed manually noting the Cholesky 
    ! factorization causes A = R**T * R.  Then R**T * R * X = B.  
    
    ! Step 1 would then be to solve the problem R**T * Y = B, for Y.
    call solve_triangular_system(.true., .true., .true., a, bu)

    ! Now, solve the problem R * X = Y, for X
    call solve_triangular_system(.true., .false., .true., a, bu)

    ! Display the results
    print '(A)', "Cholesky Solution (Manual Approach): X = "
    print '(F8.4)', (bu(i), i = 1, size(bu))
end program