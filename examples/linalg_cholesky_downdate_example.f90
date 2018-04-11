! linalg_cholesky_downdate_example.f90

program example
    use iso_fortran_env, only : real64, int32
    use linalg_factor, only : cholesky_factor, cholesky_rank1_downdate
    use linalg_core, only : rank1_update
    implicit none

    ! Variables
    real(real64) :: a(3,3), u(3), ad(3,3)
    integer(int32) :: i

    ! Build the 3-by-3 matrix A.
    !     |  4.25   11.25   -15 |
    ! A = | 11.25   39.25   -46 |
    !     |  -15     -46    102 |
    a = reshape([4.25d0, 11.25d0, -15.0d0, 11.25d0, 39.25d0, -46.0d0, &
        -15.0d0, -46.0d0, 102.0d0], [3, 3])
    
    ! The downdate vector
    !     |  0.5 |
    ! u = | -1.5 |
    !     |   2  |
    u = [0.5d0, -1.5d0, 2.0d0]
    
    ! Compute the rank 1 downdate of A
    ad = a
    call rank1_update(-1.0d0, u, u, ad)

    ! Compute the Cholesky factorization of the original matrix
    call cholesky_factor(a)

    ! Apply the rank 1 downdate to the factored matrix
    call cholesky_rank1_downdate(a, u)

    ! Compute the Cholesky factorization of the downdate to the original matrix
    call cholesky_factor(ad)

    ! Display the matrices
    print '(A)', "Downdating the Factored Form:"
    do i = 1, size(a, 1)
        print *, a(i,:)
    end do
    
    print '(A)', "Downdating A Directly:"
    do i = 1, size(ad, 1)
        print *, ad(i,:)
    end do
end program