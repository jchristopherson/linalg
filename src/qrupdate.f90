!> @brief A module providing explicit interfaces for the QRUPDATE library.
module qrupdate
    implicit none

    interface
        subroutine DQR1UP(m, n, k, q, ldq, r, ldr, u, v, w)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, ldq, ldr
            real(real64), intent(inout) :: q(ldq,*), r(ldr,*), u(*), v(*)
            real(real64), intent(out) :: w(*)
        end subroutine

        subroutine ZQR1UP(m, n, k, q, ldq, r, ldr, u, v, w, rw)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: m, n, k, ldq, ldr
            complex(real64), intent(inout) :: q(ldq,*), r(ldr,*), u(*), v(*)
            complex(real64), intent(out) :: w(*)
            real(real64), intent(out) :: rw(*)
        end subroutine

        subroutine DCH1UP(n, r, ldr, u, w)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ldr
            real(real64), intent(inout) :: r(ldr,*), u(*)
            real(real64), intent(out) :: w(*)
        end subroutine

        subroutine ZCH1UP(n, r, ldr, u, w)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ldr
            complex(real64), intent(inout) :: r(ldr,*), u(*)
            real(real64), intent(out) :: w(*)
        end subroutine

        subroutine DCH1DN(n, r, ldr, u, w, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ldr
            real(real64), intent(inout) :: r(ldr,*), u(*)
            real(real64), intent(out) :: w(*)
            integer(int32), intent(out) :: info
        end subroutine

        subroutine ZCH1DN(n, r, ldr, u, rw, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, ldr
            complex(real64), intent(inout) :: r(ldr,*), u(*)
            real(real64), intent(out) :: rw(*)
            integer(int32), intent(out) :: info
        end subroutine
    end interface
end module