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
      subroutine dqhqr(m,n,R,ldr,c,s)
c purpose:      given an m-by-n upper Hessenberg matrix R, this
c               subroutine updates R to upper trapezoidal form
c               using min(m-1,n) Givens rotations.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix R
c n (in)        number of columns of the matrix R
c R (io)        on entry, the upper Hessenberg matrix R
c               on exit, the updated upper trapezoidal matrix
c ldr (in)      leading dimension of R, >= m
c c(out)        rotation cosines, size at least min(m-1,n)
c s(out)        rotation sines, size at least min(m-1,n)
c
      integer m,n,ldr
      double precision R(ldr,*),c(*),s(*)
      external xerbla,dlartg
      double precision t
      integer info,i,ii,j
c quick return if possible.
      if (m == 0 .or. m == 1 .or. n == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (ldr < m) then
        info = 4
      end if
      if (info /= 0) then
        call xerbla('DQHQR',info)
        return
      end if
      do i = 1,n
c apply stored rotations, column-wise
        t = R(1,i)
        ii = min(m,i)
        do j = 1,ii-1
          R(j,i) = c(j)*t + s(j)*R(j+1,i)
          t = c(j)*R(j+1,i) - s(j)*t
        end do
        if (ii < m) then
c generate next rotation
          call dlartg(t,R(ii+1,i),c(i),s(i),R(ii,i))
          R(ii+1,i) = 0d0
        else
          R(ii,i) = t
        end if
      end do
      end subroutine
