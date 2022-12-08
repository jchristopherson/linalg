c Copyright (C) 2009  VZLU Prague, a.s., Czech Republic
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
      subroutine dlu1up(m,n,L,ldl,R,ldr,u,v)
c purpose:      updates an LU factorization after rank-1 modification
c               i.e., given an m-by-k lower-triangular matrix L with unit
c               diagonal and a k-by-n upper-trapezoidal matrix R,
c               where k = min(m,n),
c               this subroutine updates L -> L1 and R -> R1 so that
c               L is again lower unit triangular, R upper trapezoidal,
c               and L1*R1 = L*R + u*v.'.
c               (real version)
c arguments:
c m (in)        order of the matrix L.
c n (in)        number of columns of the matrix U.
c L (io)        on entry, the unit lower triangular matrix L.
c               on exit, the updated matrix L1.
c ldl (in)      the leading dimension of L. ldl >= m.
c R (io)        on entry, the upper trapezoidal m-by-n matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      the leading dimension of R. ldr >= min(m,n).
c u (io)        the left m-vector. On exit, if k < m, u is destroyed.
c v (io)        the right n-vector. On exit, v is destroyed.
c
c REMARK:       Algorithm is due to
c               J. Bennett: Triangular factors of modified matrices,
c                           Numerische Mathematik, 7 (1965)
c
      integer m,n,ldl,ldr
      double precision L(ldl,*),R(ldr,*),u(*),v(*)
      double precision ui,vi
      integer k,info,i,j
      external xerbla
c quick return if possible.
      k = min(m,n)
      if (k == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (ldl < m) then
        info = 4
      else if (ldr < k) then
        info = 6
      endif
      if (info /= 0) then
        call xerbla('DLU1UP',info)
        return
      end if

c The Bennett algorithm, modified for column-major access.
c The leading part.
      do i = 1,k
c prefetch        
        ui = u(i)
        vi = v(i)
c delayed R update        
        do j = 1,i-1
          R(j,i) = R(j,i) + u(j)*vi
          vi = vi - v(j)*R(j,i)
        end do
c diagonal update        
        R(i,i) = R(i,i) + ui*vi
        vi = vi/R(i,i)
c L update        
        do j = i+1,m
          u(j) = u(j) - ui*L(j,i)
          L(j,i) = L(j,i) + u(j)*vi
        end do
        u(i) = ui
        v(i) = vi
      end do

c Finish the trailing part of R if needed.
      do i = k+1,n
        vi = v(i)
        do j = 1,k
          R(j,i) = R(j,i) + u(j)*vi
          vi = vi - v(j)*R(j,i)
        end do
        v(i) = vi
      end do
      end subroutine
