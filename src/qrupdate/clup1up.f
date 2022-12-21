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
      subroutine clup1up(m,n,L,ldl,R,ldr,p,u,v,w)
c purpose:      updates a row-pivoted LU factorization after rank-1 modification
c               i.e., given an m-by-k lower-triangular matrix L with unit
c               diagonal, a k-by-n upper-trapezoidal matrix R, and a 
c               permutation matrix P, where k = min(m,n),
c               this subroutine updates L -> L1, R -> R1 and P -> P1 so that
c               L is again lower unit triangular, R upper trapezoidal,
c               P permutation and P1'*L1*R1 = P'*L*R + u*v.'.
c               (complex version)
c arguments:
c m (in)        order of the matrix L.
c n (in)        number of columns of the matrix U.
c L (io)        on entry, the unit lower triangular matrix L.
c               on exit, the updated matrix L1.
c ldl (in)      the leading dimension of L. ldl >= m.
c R (io)        on entry, the upper trapezoidal m-by-n matrix R.
c               on exit, the updated matrix R1.
c ldr (in)      the leading dimension of R. ldr >= min(m,n).
c p (in)        the permutation vector representing P
c u (in)        the left m-vector.
c v (in)        the right n-vector.
c w (work)      a workspace vector of size m.
c
c REMARK:       Algorithm is due to
c               A. Kielbasinski, H. Schwetlick, Numerische Lineare
c               Algebra, Verlag Harri Deutsch, 1988
c
      integer m,n,ldl,ldr,p(*)
      complex L(ldl,*),R(ldr,*),u(*),v(*),w(*)
      complex one,tmp
      real tau
      parameter (one = 1e0, tau = 1e-1)
      integer k,info,i,j,itmp
      external xerbla,ccopy,caxpy,ctrsv,cgeru,cgemv

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
        call xerbla('CLU1UP',info)
        return
      end if

c form L \ P*u.
      do i = 1,m
        w(i) = u(p(i))
      end do
      call ctrsv('L','N','U',k,L,ldl,w,1)
c if m > k = n, subtract the trailing part.
      if (m > k) then
        call cgemv('N',m-k,k,-one,L(k+1,1),ldl,w,1,one,w(k+1),1)
      end if
            
c work from bottom to top
      do j = k-1,1,-1
        if (abs(w(j)) < tau * abs(L(j+1,j)*w(j) + w(j+1))) then
c need pivoting. swap j and j+1
          tmp = w(j)
          w(j) = w(j+1)
          w(j+1) = tmp
c update p          
          itmp = p(j)
          p(j) = p(j+1)
          p(j+1) = itmp
c update L
          call cswap(m-j+1,L(j,j),1,L(j,j+1),1)
          call cswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
c update R
          call cswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
c make L lower triangular again
          tmp = -L(j,j+1)
          call caxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
c update R
          call caxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
c update w
          w(j) = w(j) - tmp*w(j+1)          
        end if
c eliminate w(j+1)
        tmp = w(j+1)/w(j)
        w(j+1) = 0
c update R.
        call caxpy(n-j+1,-tmp,R(j,j),ldr,R(j+1,j),ldr)
c update L.
        call caxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
      end do

c add a multiple of v to R
      call caxpy(n,w(1),v,1,R(1,1),ldr)

c forward sweep
      do j = 1,k-1
        if (abs(R(j,j)) < tau * abs(L(j+1,j)*R(j,j) + R(j+1,j))) then
c need pivoting. swap j and j+1
c update p          
          itmp = p(j)
          p(j) = p(j+1)
          p(j+1) = itmp
c update L
          call cswap(m-j+1,L(j,j),1,L(j,j+1),1)
          call cswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
c update R
          call cswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
c make L lower triangular again
          tmp = -L(j,j+1)
          call caxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
c update R
          call caxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
        end if
c eliminate R(j+1,j)
        tmp = R(j+1,j)/R(j,j)
c update R.
        R(j+1,j) = 0e0
        call caxpy(n-j,-tmp,R(j,j+1),ldr,R(j+1,j+1),ldr)
c update L.
        call caxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
      end do

c if m > k = n, complete the update by updating the lower part of L.
      if (m > k) then
        call ccopy(k,v,1,w,1)
        call ctrsv('U','T','N',k,R,ldr,w,1)
        call cgeru(m-k,k,one,w(k+1),1,w,1,L(k+1,1),ldl)
      endif
      end subroutine
      
            
