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
      subroutine dqrtv1(n,u,w)
c purpose:      generates a sequence of n-1 Givens rotations that
c               eliminate all but the first element of a vector u.
c arguments:
c n (in)        the length of the vector u
c u (io)        on entry, the vector u.
c               on exit, u(2:n) contains the rotation sines, u(1)
c               contains the remaining element.
c w (o)         on exit, w contains the rotation cosines.
c
      integer n
      double precision u(*),w(*)
      external dlartg
      double precision rr,t
      integer i
c quick return if possible.
      if (n <= 0) return
      rr = u(n)
      do i = n-1,1,-1
        call dlartg(u(i),rr,w(i),u(i+1),t)
        rr = t
      end do
      u(1) = rr
      end subroutine
