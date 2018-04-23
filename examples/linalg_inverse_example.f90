! linalg_inverse_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg_core
    implicit none

    ! Variables
    real(real64) :: a(3,3), ai(3,3), c(3,3)
    integer(int32) :: i

    ! Construct the 3-by-3 matrix A to invert
    !     | 1   2   3 |
    ! A = | 4   5   6 |
    !     | 7   8   0 |
    a = reshape([1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, &
        0.0d0], [3, 3])
    
    ! Compute the inverse of A.  Notice, the original matrix is overwritten
    ! with it's inverse.
    ai = a
    call mtx_inverse(ai)

    ! Show that A * inv(A) = I
    c = matmul(a, ai)

    ! Display the inverse
    print '(A)', "Inverse:"
    do i = 1, size(ai, 1)
        print *, ai(i,:)
    end do

    ! Display the result of A * inv(A)
    print '(A)', "A * A**-1:"
    do i = 1, size(c, 1)
        print *, c(i,:)
    end do
end program