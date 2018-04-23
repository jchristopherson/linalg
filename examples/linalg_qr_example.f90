! linalg_qr_example.f90

! Example Source:
! https://www.mathworks.com/help/matlab/ref/lu.html?s_tid=srchtitle
! 
! This is a repeat of the LU factorization example, but using QR factorization
! instead.
program example
    use iso_fortran_env, only : real64, int32
    use linalg_core
    implicit none

    ! Local Variables
    real(real64) :: a(3,3), tau(3), b(3)
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

    ! Compute the QR factorization, using pivoting
    pvt = 0     ! Zero every entry in order not to lock any column in place
    call qr_factor(a, tau, pvt)

    ! Compute the solution.  The results overwrite b.
    call solve_qr(a, tau, pvt, b)

    ! Display the results.
    print '(A)', "QR Solution: X = "
    print '(F8.4)', (b(i), i = 1, size(b))

    ! Notice, QR factorization without pivoting could be accomplished in the
    ! same manner.  The only difference is to omit the PVT array (column pivot
    ! tracking array).
end program