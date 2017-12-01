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
      subroutine dqrot(dir,m,n,Q,ldq,c,s)
c purpose:      Apply a sequence of inv. rotations from right
c
c arguments:
c dir (in)      if 'B' or 'b', rotations are applied from backwards
c               if 'F' or 'f', from forwards.
c m (in)        number of rows of matrix Q
c n (in)        number of columns of the matrix Q
c Q (io)        on entry, the matrix Q
c               on exit, the updated matrix Q1
c ldq (in)      the leading dimension of Q
c c (in)        n-1 rotation cosines
c s (in)        n-1 rotation sines
c
      character dir
      integer m,n,ldq
      double precision Q(ldq,*),c(*),s(*)
      external drot,lsame
      logical lsame,fwd
      integer info,i
c quick return if possible
      if (m == 0 .or. n == 0 .or. n == 1) return
c check arguments.
      info = 0
      fwd = lsame(dir,'F')
      if (.not.(fwd .or. lsame(dir,'B'))) then
        info = 1
      else if (m < 0) then
        info = 2
      else if (n < 0) then
        info = 3
      else if (ldq < m) then
        info = 5
      end if
      if (info /= 0) then
        call xerbla('DQROT',info)
        return
      end if

      if (fwd) then
        do i = 1,n-1
          call drot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
      else
        do i = n-1,1,-1
          call drot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
      end if
      end subroutine
