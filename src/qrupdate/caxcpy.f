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
      subroutine caxcpy(n,a,x,incx,y,incy)
c purpose:      constant times a conjugated vector plus a vector.
c arguments:
c n (in)        vector length
c a (in)        complex factor
c x (in)        added vector
c incx (in)     x increments
c y (io)        accumulator vector
c incy (in)     y increments
c
      integer n,incx,incy
      complex a,x(*),y(*)
      integer i,ix,iy
c quick return if possible.
      if (n <= 0) return
      if (incx /= 1 .or. incy /= 1) then
c code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if (incx.lt.0) ix = (-n+1)*incx + 1
        if (incy.lt.0) iy = (-n+1)*incy + 1
        do i = 1,n
          y(iy) = y(iy) + a*conjg(x(ix))
          ix = ix + incx
          iy = iy + incy
        end do
      else
c code for both increments equal to 1
        do i = 1,n
           y(i) = y(i) + a*conjg(x(i))
        end do
      end if
      end subroutine
