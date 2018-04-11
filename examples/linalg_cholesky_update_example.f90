! linalg_cholesky_update_example

program example
    use iso_fortran_env, only : real64, int32
    use linalg_factor, only : cholesky_factor, cholesky_rank1_update
    use linalg_core, only : rank1_update
    implicit none

    ! Variables
    real(real64) :: a(3,3), u(3), au(3,3)
    integer(int32) :: i

    ! Build the 3-by-3 positive-definite matrix A.
    !     | 4   12   -16 |
    ! A = | 12  37   -43 |
    !     |-16 -43    98 |
    a = reshape([4.0d0, 12.0d0, -16.0d0, 12.0d0, 37.0d0, -43.0d0, -16.0d0, &
        -43.0d0, 98.0d0], [3, 3])
    
    ! Build the update vector U
    u = [0.5d0, -1.5d0, 2.0d0]

    ! Compute the rank 1 update of A
    au = a
    call rank1_update(1.0d0, u, u, au)

    ! Compute the Cholesky factorization of the original matrix
    call cholesky_factor(a)

    ! Apply the rank 1 update to the factored matrix
    call cholesky_rank1_update(a, u)

    ! Compute the Cholesky factorization of the update of the original matrix
    call cholesky_factor(au)

    ! Display the matrices
    print '(A)', "Updating the Factored Form:"
    do i = 1, size(a, 1)
        print *, a(i,:)
    end do
    
    print '(A)', "Updating A Directly:"
    do i = 1, size(au, 1)
        print *, au(i,:)
    end do
end program