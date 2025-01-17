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
      subroutine dgqvec(m,n,Q,ldq,u)
c purpose:      given an orthogonal m-by-n matrix Q, n < m, generates
c               a vector u such that Q'*u = 0 and norm(u) = 1.
c arguments:
c m (in)        number of rows of matrix Q.
c n (in)        number of columns of matrix Q.
c Q (in)        the orthogonal matrix Q.
c ldq (in)      leading dimension of Q.
c u (out)       the generated vector.
c
      integer m,n,ldq
      double precision Q(ldq,*),u(*)
      external ddot,daxpy,dnrm2,dscal
      double precision ddot,dnrm2,r
      integer info,i,j
c quick return if possible.
      if (m == 0) return
      if (n == 0) then
        u(1) = 1d0
        do i = 2,m
          u(i) = 0d0
        end do
        return
      end if
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (ldq < m) then
        info = 4
      end if
      if (info /= 0) then
        call xerbla('DGQVEC',info)
        return
      end if

      j = 1
 10   continue
c probe j-th canonical unit vector.
      do i = 1,m
        u(i) = 0d0
      end do
      u(j) = 1d0
c form u - Q*Q'*u
      do i = 1,n
        r = ddot(m,Q(1,i),1,u,1)
        call daxpy(m,-r,Q(1,i),1,u,1)
      end do
      r = dnrm2(m,u,1)
      if (r == 0d0) then
        j = j + 1
        if (j > n) then
c this is fatal, and in theory, it can't happen.
          stop 'fatal: impossible condition in DGQVEC'
        else
          j = j + 1
          goto 10
        end if
      end if
      call dscal(m,1d0/r,u,1)
      end subroutine
