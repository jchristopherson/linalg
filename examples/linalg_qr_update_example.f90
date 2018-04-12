! linalg_qr_update_example.f90

program example
    use iso_fortran_env
    use linalg_core
    implicit none

    ! Variables
    real(real64) :: a(3,3), u(3), v(3), r(3,3), tau(3), q(3,3), qu(3,3)
    integer(int32) :: i

    ! Build the 3-by-3 matrix A.
    !     | 1   2   3 |
    ! A = | 4   5   6 |
    !     | 7   8   0 |
    a = reshape( &
        [1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 0.0d0], &
        [3, 3])

    ! Build the update vectors
    !     | 1/2 |      | 1 |
    ! u = | 3/2 |, v = | 5 |
    !     |  3  |      | 2 |
    u = [0.5d0, 1.5d0, 3.0d0]
    v = [1.0d0, 5.0d0, 2.0d0]

    ! Compute the QR factorization of the original matrix
    r = a   ! Making a copy as the matrix will be overwritten by qr_factor
    call qr_factor(r, tau)

    ! Form Q & R
    call form_qr(r, tau, q)

    ! Compute the rank 1 update to the original matrix such that: 
    ! A = A + u * v**T
    call rank1_update(1.0d0, u, v, a)

    ! Compute the rank 1 update to the factorization.  Notice, the contents 
    ! of U & V are destroyed as part of this process.
    call qr_rank1_update(q, r, u, v)

    ! As comparison, compute the QR factorization of the rank 1 updated matrix
    call qr_factor(a, tau)
    call form_qr(a, tau, qu)

    ! Display the matrices
    print '(A)', "Updating the Factored Form:"
    print '(A)', "Q = "
    do i = 1, size(q, 1)
        print *, q(i,:)
    end do
    print '(A)', "R = "
    do i = 1, size(r, 1)
        print *, r(i,:)
    end do
    
    print '(A)', "Updating A Directly:"
    print '(A)', "Q = "
    do i = 1, size(qu, 1)
        print *, qu(i,:)
    end do
    print '(A)', "R = "
    do i = 1, size(a, 1)
        print *, a(i,:)
    end do
end program 