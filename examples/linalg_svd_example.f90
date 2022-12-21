! linalg_svd_example.f90

program example
    use iso_fortran_env, only : int32, real64
    use linalg
    implicit none

    ! Variables
    real(real64) :: a(3,2), s(2), u(3,3), vt(2,2), ac(3,2)
    integer(int32) :: i

    ! Initialize the 3-by-2 matrix A
    !     | 2   1 |
    ! A = |-3   1 |
    !     |-1   1 |
    a = reshape([2.0d0, -3.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0], [3, 2])

    ! Compute the singular value decomposition of A.  Notice, V**T is returned
    ! instead of V.  Also note, A is overwritten.
    call svd(a, s, u, vt)

    ! Display the results
    print '(A)', "U ="
    do i = 1, size(u, 1)
        print *, u(i,:)
    end do

    print '(A)', "S ="
    print '(F9.5)', (s(i), i = 1, size(a, 2))

    print '(A)', "V**T ="
    do i = 1, size(vt, 1)
        print *, vt(i,:)
    end do

    ! Compute U * S * V**T, but first establish S in full form
    call diag_mtx_mult(.true., 1.0d0, s, vt) ! Compute: VT = S * V**T
    ac = matmul(u(:,1:2), vt)
    print '(A)', "U * S * V**T ="
    do i = 1, size(ac, 1)
        print *, ac(i,:)
    end do
end program
