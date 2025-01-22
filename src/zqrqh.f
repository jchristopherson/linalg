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
      subroutine zqrqh(m,n,R,ldr,c,s)
c purpose:      brings an upper trapezoidal matrix R into upper
c               Hessenberg form using min(m-1,n) Givens rotations.
c               (complex version)
c arguments:
c m (in)        number of rows of the matrix R
c n (in)        number of columns of the matrix R
c R (io)        on entry, the upper Hessenberg matrix R
c               on exit, the updated upper trapezoidal matrix
c ldr (in)      leading dimension of R, >= m
c c(in)         rotation cosines, size at least min(m-1,n)
c s(in)         rotation sines, size at least min(m-1,n)
c
      integer m,n,ldr
      double complex R(ldr,*),s(*)
      double precision c(*)
      external xerbla
      double complex t
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
        call xerbla('ZQRQH',info)
        return
      end if
      do i = 1,n
c apply stored rotations, column-wise
        ii = min(m-1,i)
        t = R(ii+1,i)
        do j = ii,1,-1
          R(j+1,i) = c(j)*t - conjg(s(j))*R(j,i)
          t = c(j)*R(j,i) + s(j)*t
        end do
        R(1,i) = t
      end do
      end subroutine
