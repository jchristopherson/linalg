C DISDOT - Re-implementation of DDOT required for ITERS.F solvers.
      double precision function distdot(n, x, ix, y, iy)
      implicit none
      integer n, ix, iy
      real*8 x(*), y(*)
      real*8 ddot
      external ddot
      distdot = ddot(n, x, ix, y, iy)
      end function