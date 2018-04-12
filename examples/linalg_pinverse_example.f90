! linalg_pinverse_example.f90

program example
    use iso_fortran_env, only : int32, real64
    use linalg_solve, only : mtx_pinverse
    implicit none

    ! Variables
    real(real64) :: a(3,2), ai(2,3), ao(3,2), c(2,2)
    integer(int32) :: i

    ! Create the 3-by-2 matrix A
    !     | 1   0 |
    ! A = | 0   1 |
    !     | 0   1 |
    a = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0], [3, 2])
    ao = a  ! Just making a copy for later as mtx_pinverse will destroy the
            ! contents of the original matrix

    ! The Moore-Penrose pseudo-inverse of this matrix is:
    !         | 1   0    0  |
    ! A**-1 = |             |
    !         | 0  1/2  1/2 |
    call mtx_pinverse(a, ai)

    ! Notice, A**-1 * A is an identity matrix.
    c = matmul(ai, ao)

    ! Display the inverse
    print '(A)', "Inverse:"
    do i = 1, size(ai, 1)
        print *, ai(i,:)
    end do

    ! Display the result of inv(A) * A
    print '(A)', "A**-1 * A:"
    do i = 1, size(c, 1)
        print *, c(i,:)
    end do
end program