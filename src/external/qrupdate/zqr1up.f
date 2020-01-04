c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine zqr1up(m,n,k,Q,ldq,R,ldr,u,v,w,rw)
c purpose:      updates a QR factorization after rank-1 modification
c               i.e., given a m-by-k unitary Q and m-by-n upper
c               trapezoidal R, an m-vector u and n-vector v,
c               this subroutine updates Q -> Q1 and R -> R1 so that
c               Q1*R1 = Q*R + u*v', and Q1 is again unitary
c               and R1 upper trapezoidal.
c               (complex version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c k (in)        number of columns of Q, and rows of R. Must be
c               either k = m (full Q) or k = n < m (economical form).
c Q (io)        on entry, the unitary m-by-k matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      the leading dimension of Q. ldq >= m.
c R (io)        on entry, the upper trapezoidal m-by-n matrix R..
c               on exit, the updated matrix R1.
c ldr (in)      the leading dimension of R. ldr >= k.
c u (io)        the left m-vector. On exit, if k < m, u is destroyed.
c v (io)        the right n-vector. On exit, v is destroyed.
c w (out)       a workspace vector of size k.
c rw (out)      a real workspace vector of size k.
c
      integer m,n,k,ldq,ldr
      double complex Q(ldq,*),R(ldr,*),u(*),v(*),w(*)
      double precision rw(*)
      external zqrqh,zqhqr,zqrot,zqrtv1,zaxcpy
      external zdotc,dznrm2,dlamch,zdscal,zrot
      double complex zdotc
      double precision dznrm2,dlamch,ru,ruu
      integer info,i
      logical full
c quick return if possible.
      if (k == 0 .or. n == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (k /= m .and. (k /= n .or. n > m)) then
        info = 3
      else if (ldq < m) then
        info = 5
      else if (ldr < k) then
        info = 7
      endif
      if (info /= 0) then
        call xerbla('ZQR1UP',info)
        return
      end if

      full = k == m
c in the non-full case, we shall need the norm of u.
      if (.not.full) ru = dznrm2(m,u,1)
c form Q'*u. In the non-full case, form also u - Q*Q'u.
      do i = 1,k
        w(i) = zdotc(m,Q(1,i),1,u,1)
        if (.not.full) call zaxpy(m,-w(i),Q(1,i),1,u,1)
      end do
c generate rotations to eliminate Q'*u.
      call zqrtv1(k,w,rw)
c apply rotations to R.
      call zqrqh(k,n,R,ldr,rw,w(2))
c apply rotations to Q.
      call zqrot('B',m,k,Q,ldq,rw,w(2))
c update the first row of R.
      call zaxcpy(n,w(1),v,1,R(1,1),ldr)
c retriangularize R.
      call zqhqr(k,n,R,ldr,rw,w)
c apply rotations to Q.
      call zqrot('F',m,min(k,n+1),Q,ldq,rw,w)
c in the full case, we're finished
      if (full) return
c compute relative residual norm
      ruu = dznrm2(m,u,1)
      ru = ru * dlamch('e')
      if (ruu <= ru) return
c update the orthogonal basis.
      call zdscal(n,ruu,v,1)
      call zdscal(m,1d0/ruu,u,1)
      call zch1up(n,R,ldr,v,rw)
      do i = 1,n
        call zrot(m,Q(1,i),1,u,1,rw(i),conjg(v(i)))
      end do
      end subroutine
