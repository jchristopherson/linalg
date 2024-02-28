      subroutine lstif3(ske, fe, xe, ye, xyk, funb, func, fung)
c---------------------------------------------------------------------------
c
c  This subroutine computes the local stiffness matrix for the 
c   Diffusion-Convection Equation with the 
c    variable cofficients, 'K(x,y), B(x,y), C(x,y) '
c
c  -Div( K(x,y) T(x,y)) + B(x,y) Tx + C(x,y) Ty = G
c
c   Here K(x,y) is a 2x2 Matrix, where each entry is a function of x and y.
c
c   K, B, C and G need to be supplied by user.
c    They need to be defined as externals in the calling routines.
c
c   PSI(i,x,y) : i-th shape fucntions on the standard triangle N, i=1, 2, 3
c      where N is the following.
c
c             (-1,1)
c                .
c                . .
c		 .   .
c		 .       .
c		 . . . . . . (1,-1)
c	       (-1,-1)
c
c   Local stiffness matrix is obtained by integral on the current 
c   element. To do so, change the current coordinates to N 
c    by Affine mapping, sending
c
c    (xe(1),ye(1))   ---> (-1,-1)
c    (xe(2),ye(2))   ---> (1,-1)
c    (xe(3),ye(3))   ---> (-1,1) .
c
c    Then we perform the integration on N
c     by Gaussian Quadrature with 9 points.
c
c---------------------------------------------------------------------------
c
c  on entry
c  ---------
c
c   xe       = x coordinates of the nodes in the current element.
c   ye       = y coordinates of the nodes in the current element.
c   xyk      = subroutine defining the function K(x,y). 
c   funb     = function defining the function b(x,y).
c   func     = function defining the function c(x,y).
c   fung     = function defining the function g(x,y).
c
c---------------------------------------------------------------------------
c
c  on return
c  ---------
c
c   ske : Local Stiffness Matrix.( 3x3 in this subroutine.)
c   fe  : Local Load Vector.
c
c---------------------------------------------------------------------------
        implicit real*8(a-h,o-z)
        dimension ske(3,3), fe(3), xe(3), ye(3), dn(3,2), 
     1            xyke(2,2), wei(9), gau1(9), gau2(9)
        external xyk, funb, func, fung

c  Gau1 and Gau2 are the Gaussian Quadrature Points for the Traingle N,
c   and Wei, are the corresponding weights.
c
c   They are derived from the 1-D case by Reiterated integrals.
c
        data gau1/-0.8,-0.1127016654,0.5745966692,-0.8872983346,-0.5, 
     1        -0.1127016654,  -0.9745966692, -0.8872983346, -0.8 /
        data gau2/3*-0.7745966692, 3*0., 3*0.7745966692 /
        data wei/0.2738575107,0.4381720172,0.2738551072,0.2469135803, 
     1     0.3950617284,0.2469135803,0.03478446464, 
     2     0.05565514341,0.03478446464 /

        npt = 9

c
c    Compute the Affine mappings from the current triangle to the
c     standard triangle N. Integration will be performed on that
c     triangle by Gaussian quadrature.
c
c    T = A X + B
c   
c    A11, A12, A21, A22, B1, B2 will denote the entries of
c     A & B.
c
        x1 = xe(1)
        x2 = xe(2)
        x3 = xe(3)
        y1 = ye(1)
        y2 = ye(2)
        y3 = ye(3)

        rj1 = (x3-x1)*(y2-y3) - (x2-x3)*(y3-y1)
        rj2 = (x3-x1)*(y1-y2) - (x1-x2)*(y3-y1)
        a11 = 2*(y1-y3)/rj1
        a12 = 2*(x3-x1)/rj1
        a21 = 2*(y1-y2)/rj2
        a22 = 2*(x2-x1)/rj2
        b1 = 1. - a11*x2 - a12*y2
        b2 = -1. - a21*x2 - a22*y2
c
c  Compute the first order partial derivatives of the shape functions.
c   dn(i,1) and dn(i,2) are the first order partial derivativ of i-th shape function
c    with respect to x and y, respectively.
c
        dn(1,1) = -0.5*(a11+a21)
        dn(1,2) = -0.5*(a12+a22)
        dn(2,1) = 0.5*a11
        dn(2,2) = 0.5*a12
        dn(3,1) = 0.5*a21
        dn(3,2) = 0.5*a22
c  Compute the Jacobian associated with T.
        Rja = a11*a22 - a12*a21
c
c  Find the inverse mapping of T
c

        u11 = a22/rja
        u12 = -a12/Rja
        u21 = -a21/rja
        u22 = a11/rja
        v1 = -u11*b1 - u12*b2
        v2 = -u21*b1 - u22*b2

        do 200 i = 1 , 3
          T4 = 0.
        do 220 j = 1 , 3
          T1 = 0.
          T2 = 0.
          T3 = 0.
          do 250 k = 1, npt
            r = gau1(k)
            s = gau2(k)
            w = wei(k)

            x = u11*r + u12*s + v1
            y = u21*r + u22*s + v2

            call xyk(xyke, x, y)

            derv2 = dn(i,1)*dn(j,1)*xyke(1,1) 
     1            + dn(i,2)*dn(j,2)*xyke(2,2)
     2            + dn(i,1)*dn(j,2)*xyke(1,2)
     3            + dn(i,2)*dn(j,1)*xyke(2,1)
            if(j .eq. 1) then
              t4 = t4 + w*fung(x,y)*psi(i,r,s)
            endif

            t1 = t1 + w*derv2
            t2 = t2 + w*funb(x,y)*psi(i,r,s)
            T3 = T3 + w*func(x,y)*psi(i,r,s)
250          continue
        
        ske(i,j) = (t1 + t2*dn(j,1) + t3*dn(j,2))/rja
220        continue        
        fe(i) = t4/rja
200        continue

        return
        end
c--- end of lstif3 -----------------------------------------------------
c-----------------------------------------------------------------------
       subroutine refall(nx, nelx,ijk,node,ndeg,x,y,
     *     ichild,iparnts,nodcode,nxmax,nelmax,ierr)
      implicit real*8  (a-h,o-z)
      integer nx, nelx, node, ndeg, nxmax, nelmax 
      integer ichild(ndeg,1),iparnts(2,nx),ijk(node,*), nodcode(nx)
      integer midnode(20),inod(20) 
      real*8  x(*),y(*)
c-------------------------------------------------------------
c refines a finite element grid using triangular elements.
c uses mid points to refine all the elements of the grid.
c
c nx	= number of nodes at input
c nelx	= number of elements at input
c ijk	= connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k. 
c node  = first dimension of array ijk [should be >=3] 
c ndeg	= first dimension of array ichild which is at least as large
c         as the max degree of each node
c x,y   = real*8 arrays containing the x(*) and y(*) coordinates 
c	  resp. of the nodes.
c ichild= list of the children of a node: ichild(1,k) stores 
c         the position in ichild(*,k)  of the last child so far.
c         (local use)
c iparnts= list of the 2 parents of each node.
c         (local use)
c nodcode= boundary information list for each node with the
c	   following meaning:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c corner elements are used only to generate the grid by refinement 
c since they do not  correspond to real elements. 
c nxmax  = maximum number of nodes allowed. If during the algorithm
c          the number of nodes being created exceeds nxmax then 
c	   refall  quits without modifying the (x,y) xoordinates
c	   and nx, nelx. ijk is modified. Also ierr is set to 1.
c nelmax = same as above for number of elements allowed. See ierr..
c ierr	 = error message: 
c	   0 --> normal return
c	   1 --> refall quit because nxmax  was exceeded.
c	   2 --> refall quit because nelmax was exceeded.
c--------------------------------------------------------------
c---------------------------------------------------------------  
c inilitialize lists of children and parents --
c data structure is as follows
c ichild(1,k) stores the position of last child of node k so far in list
c ichild(j,k) , j .ge. 2 = list of children of node k.
c iparnts(1,k) and iparnts(2,k) are the two parents of node k. 
c---------------------------------------------------------------  
c------ do a first check :
      if (nx .ge. nxmax) goto 800
      if (nelx .ge. nelmax) goto 900
c------ initialize
      do 1 k=1,nx
         do 2 j=2,ndeg
            ichild(j,k) = 0
 2       continue
         ichild(1,k) = 1
         iparnts(1,k)= 0
         iparnts(2,k)= 0
 1    continue
c------- initialize nelxnew and nxnew 
      nelxnew = nelx
      nxnew   = nx
      ierr    = 0
c--------------------------------------------------------------
c main loop: scan all elements
c--------------------------------------------------------------
c     do 100 nel = nelx,1,-1
      do 100 nel = 1, nelx
c note : interesting question which order is best for parallelism?
c alternative order: do 100 nel = nelx, 1, -1
c
c------ unpack nodes of element
         do 101 i=1,node
            inod(i) = ijk(i,nel)
c convention: node after last node = first node. 
            inod(node+i) = inod(i)
            midnode(i) = 0
 101     continue
c--------------------------------------------------------------
c for each new potential node determine if it has already been 
c numbered. a potential node is the middle of any two nodes ..
c-------------------------------------------------------------- 
         do 80 ii=1,node
                k1 = inod(ii)
                k2 = inod(ii+1)
c------- test for current pair :
                last = ichild(1,k1)
                do 21 k=2,last
                   jchild = ichild(k,k1) 
                   ipar1 = iparnts(1,jchild)
                   ipar2 = iparnts(2,jchild)
                   if( ((ipar1 .eq. k1 .and. ipar2 .eq. k2)) .or. 
     *                  ((ipar2 .eq. k1 .and. ipar1 .eq. k2))) then
c node has already been created and numbered ....
                      midnode(ii) = jchild
c... therefore it must be an internal node
                      nodcode(jchild) = 0
c... and no new node to create.
                      goto 80
                   endif
c-----------------------------------------------------            
 21         continue
c     
c else  create a new node
c
            nxnew = nxnew + 1
            if (nxnew .gt. nxmax) goto 800
c-------
            x(nxnew) = (x(k1) + x(k2))*0.5
            y(nxnew) = (y(k1) + y(k2))*0.5
            midnode(ii) = nxnew
c
c update nodcode information -- normally min0(nodcode(k1),nodcode(k2))
c 
            nodcode(nxnew) = min0(1,nodcode(k1),nodcode(k2))
c     
c update parents and children's lists
c 
            iparnts(1,nxnew) = k1
            iparnts(2,nxnew) = k2
c     
            last = last+1
            ichild(last,k1) = nxnew
            ichild(1,k1) = last
c     
            last = ichild(1,k2)+1
            ichild(last,k2) = nxnew
            ichild(1,k2) = last
c     
 80      continue                
c
c------- replace current element by new one
c
         do 81 i=1,node
            jnod = midnode(i)
            ijk(i,nel) = jnod
 81      continue
c-------create new elements
         do 82 ii=1, node
            nelxnew = nelxnew+1
            if (nelxnew .gt. nelmax) goto 900
            ijk(1,nelxnew) = inod(ii)
            k = ii
            do jj=2,node
               ijk(jj,nelxnew) = midnode(k)
               k = k+2
               if (k .gt. node) k =  k-node
           enddo
 82     continue
c------ done !
 100  continue
      nx = nxnew
      nelx = nelxnew
      return
 800  ierr = 1
      return
 900  ierr = 2
      return
      end
c
      subroutine checkref(nx,nelx,nodcode,nbound,nxnew,nelxnew) 
c-------------------------------------------------------------        
c returns the expected the new number of nodes and 
c elemnts of refall is applied to current grid once.
c
c nx	= number of nodes at input
c nelx	= number of elements at input
c [ijk	= connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k.] UNUSED here.
c nbound  = number of boundary points on entry - enter zero if
c	     unknown
c
c nodcode= boundary information list for each node with the
c	   following meaning:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c
c nxnew  = new number of nodes if refall were to be applied
c nelxnew = same for nelx.
c--------------------------------------------------------------
       integer nodcode(nx)
c
        nelxnew = nelx*4
c
c count the number of boundary nodes
c
        if (nbound .ne. 0) goto 2
        do 1 j=1, nx
        if (nodcode(j) .ge. 1) nbound = nbound+1
 1        continue
c number of edges=[3*(number of elmts) + number of bound nodes ]/ 2
 2        continue
        nxnew = nx + (3*nelx+nbound)/2
        nbound = 2*nbound
        return
        end                                
c----------------------------------------------------------------------- 
        subroutine unassbl (a,na,f,nx,nelx,ijk,node,x,y,ierr,xyk)
c----------------------------------------------------------------------- 
c a      = un-assembled matrix on output
c na	 = 1-st dimension of a.  a(na,node,node)
c
c f      = right hand side (global load vector) in un-assembled form
c nx     = number of nodes at input
c nelx	 = number of elements at input
c ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node	 = total number of nodal points in each element
c	   also second dimension of a.
c
c x,y   = real*8 arrays containing the $x$ and $y$ coordinates 
c	  resp. of the nodes.
c         K11, K22, and K12 at that element.
c ierr	= error message integer . 
c	  ierr = 0 --> normal return
c	  ierr = 1 --> negative area encountered (due to bad 
c	           numbering of nodes of an element)
c
c xyk	= subroutine defining the material properties at each 
c         element. Form: 
c 	call xyk(nel,xyke,x,y,ijk,node) 
c--------------------------------------------------------------
        implicit real*8 (a-h,o-z)
        dimension a(na,node,node),ijk(node,1),x(1),y(1),f(node,1),
     *            ske(3,3),fe(3),xe(3),ye(3),xyke(2,2)
        external xyk

c--------------------------------------------------------------
c   initialize
c--------------------------------------------------------------
        do 100 i=1, node 
        do 100 j=1, nx
               f(i,j) = 0.0d0
 100        continue
c---------------------------------------------------
c main loop
c--------------------------------------------------- 
        do 102 nel=1, nelx
c
c get coordinetes of nodal points
c 
        do 104 i=1, node
        j = ijk(i,nel)
        xe(i) = x(j)
        ye(i) = y(j)
 104        continue
c
c compute determinant
c
         det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
        if ( det .le. 0.) then
          print *, 'nel', nel, ' det = ' , det
          print *, xe(1), xe(2), xe(3)
          print *, ye(1), ye(2), ye(3)
        end if
c
c set material properties
c 
        call xyk(xyke,x,y)
c
c construct element stiffness matrix
c
        ierr = 0
        call estif3(nel,ske,det,xe,ye,xyke,ierr)
        if (ierr .ne. 0) then
          write (*,*) 'ERROR: estif3 gave an error',ierr
          return
        endif
c        write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
c assemble: add element stiffness matrix to global matrix
c 
        do 120 ka=1, node
            f(ka,nel) = fe(ka)
        do 108 kb = 1,node
            a(nel,ka,kb) = ske(ka,kb)
 108        continue
 120        continue
 102        continue
        return
        end
c-----------------------------------------------------------------------
        subroutine unassbl_lstif(a, na, f, nx, nelx, ijk, 
     *                     node, x, y, ierr, xyk, funb, func, fung)
c----------------------------------------------------------------------- 
c a      = un-assembled matrix on output
c
c na	 = 1-st dimension of a.  a(na,node,node)
c
c f      = right hand side (global load vector) in un-assembled form
c
c nx     = number of nodes at input
c
c nelx	 = number of elements at input
c
c ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c
c node	 = total number of nodal points in each element
c	   also second dimension of a.
c
c x,y   = real*8 arrays containing the $x$ and $y$ coordinates 
c	  resp. of the nodes.
c         K11, K22, and K12 at that element.
c
c ierr	= error message integer . 
c	  ierr = 0 --> normal return
c	  ierr = 1 --> negative area encountered (due to bad 
c	           numbering of nodes of an element)
c
c xyk	= subroutine defining the material properties at each 
c         element. Form:  	call xyk(xyke,x,y) 
c
c funb, = functions needed for the definition of lstif3 problem
c func,
c fung
c--------------------------------------------------------------
c moulitsa@cs.umn.edu : It uses lstif3 problem
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(na,node,node), ijk(node,*), x(*), y(*), f(node,*),
     &          ske(3,3), fe(3), xe(3), ye(3) 
      external xyk, funb, func, fung
c--------------------------------------------------------------
c   initialize
c--------------------------------------------------------------
      do i=1, node 
        do j=1, nx
          f(i,j) = 0.0d0
        end do
      end do

c---------------------------------------------------
c main loop
c--------------------------------------------------- 
      do nel=1, nelx
c
c get coordinetes of nodal points
c 
        do i=1, node
          j = ijk(i,nel)
          xe(i) = x(j)
          ye(i) = y(j)
        end do
c
c compute determinant
c
c	det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c       if ( det .le. 0.) then
c         print *, 'nel', nel, ' det = ' , det
c         print *, xe(1), xe(2), xe(3)
c         print *, ye(1), ye(2), ye(3)
c       end if
c
c construct element stiffness matrix
c
        ierr = 0
        call lstif3(ske, fe, xe, ye, xyk, funb, func, fung)
c        write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
c
c assemble: add element stiffness matrix to global matrix
c 
        do ka=1, node
          f(ka,nel) = fe(ka)
          do kb = 1,node
            a(nel,ka,kb) = ske(ka,kb)
          end do
        end do

      end do

      return
      end
c----------------------------------------------------------------------- 
        subroutine assmbo (nx, nelx, node, ijk, nodcode, x, y, a, ja, 
     *                  ia, f, iwk, jwk, ierr, xyk, funb, func, fung)
c----------------------------------------------------------------------- 
c nx     = number of nodes at input
c
c nelx	 = number of elements at input
c
c node	 = total number of nodal points in each element
c
c ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c
c nodcode= boundary information list for each node with the
c	   following meaning:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point (corner points
c
c x,y   = real arrays containing the $x$ and $y$ coordinates 
c	  resp. of the nodes.
c
c a,ja,ia= assembled matrix on output
c
c f      = right hand side (global load vector) 
c
c iwk,jwk = two integer work arrays.
c
c ierr	= error message integer . 
c	  ierr = 0 --> normal return
c	  ierr = 1 --> negative area encountered (due to bad 
c	           numbering of nodes of an element)
c
c xyk	= subroutine defining the material properties at each 
c         element. Form: 
c 	call xyk(nel,xyke,x,y,ijk,node) with on return
c         xyke =  material constant matrices. 
c         for each element nel, xyke(1,nel),xyke(2,nel) 
c         and xyke(3,nel) represent the constants
c         K11, K22, and K12 at that element.
c--------------------------------------------------------------
c moulitsa@cs.umn.edu : It has been modified so as to handle
c      more types of domains/meshes i.e.  |\ /|
c                                         | X |
c                                         |/ \|
c--------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension a(*),ijk(node,1),x(1),y(1),f(1),ske(3,3),fe(3),
     *      xe(3),ye(3),iwk(1),jwk(1)
      integer ia(1), ja(*), nodcode(1)
      external xyk, funb, func, fung

c--------------------------------------------------------------
c   initialize
c--------------------------------------------------------------
      do i=1,nx 
        f(i) = 0.0 
      end do
c initialize  pointer arrays.
      do k=1,nx+1
        ia(k) = 1
        jwk(k) = 0
      end do
      do k=1,nelx
        do j=1,node
          knod = ijk(j,k)
          ia(knod) = ia(knod) + 2
        end do
      end do
c---------------------------------------------------
      do k=1, nx
        if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
      end do
c
      ksav = ia(1)
      ia(1) = 1
      do j=2, nx+1
        ksavn = ia(j)
        ia(j) = ia(j-1) +  ksav
        iwk(j-1) = ia(j-1)-1
        ksav = ksavn
      end do

c-----------------
c main loop
c-----------------
      do nel=1, nelx
c
c get coordinates of nodal points
c
        do i=1, node
          j = ijk(i,nel)
          xe(i) = x(j)
          ye(i) = y(j)
        end do
c
c compute determinant
c
c         det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c
c set material properties
c 
c        call xyk(nel,xyke,x,y,ijk,node)
c
c construct element stiffness matrix
c
        ierr = 0
c
c        call evalg(nel, fe, xe, ye, fung, ierr)
c        call estif3(nel,ske,det,xe,ye,xyke,ierr)
        call lstif3(ske, fe, xe, ye, xyk, funb, func, fung)
        if (ierr .ne. 0) return
c
c assemble: add element stiffness matrix to global matrix
c 
        do ka=1, node
          ii = ijk(ka,nel)
          f(ii) = f(ii) + fe(ka)
c
c unpack row into jwk1
c
          irowst = ia(ii)
          ilast  = iwk(ii) 
          do k=irowst,ilast 
            jwk(ja(k)) = k
          end do
c
          do kb = 1,node
c
c column number = jj
c
            jj = ijk(kb,nel)
            k = jwk(jj)
            if (k .eq. 0) then
              ilast = ilast+1
              jwk(jj) = ilast
              ja(ilast) = jj
              a(ilast) = ske(ka,kb)
            else 
              a(k) = a(k) + ske(ka,kb)
            endif
          end do
c refresh jwk
          do k=irowst,ilast 
            jwk(ja(k)) = 0
          end do
          iwk(ii) = ilast 
        end do
c
      end do

c squeeze away the zero entries
c added so as to handle more type of domains/meshes
      do i=1, nx
        ista=ia(i)
        isto=ia(i+1)-1
        do j=ista, isto
          if (ja(j) .EQ. 0) then
            iwk(i)=j-ista
            go to 200
          end if
        end do
 200    continue
      end do

      do i=2, nx
        ksav=ia(i)
        ia(i)=ia(i-1)+iwk(i-1)
        ksavn=ia(i)
        do j=0, iwk(i)-1
          ja(ksavn+j)=ja(ksav+j)
          a(ksavn+j) = a(ksav+j)
        end do
      end do
      ia(nx+1)=ia(nx)+iwk(nx)

      return
      end
c-----------------------------------------------------------------------
        subroutine assmbo2 (nx, nelx, node, ijk, x, y, a, ja,
     *                   ia, f, iwk, jwk, ierr, xyk, funb, func, fung)
c----------------------------------------------------------------------- 
c nx     = number of nodes at input
c
c nelx	 = number of elements at input
c
c node	 = total number of nodal points in each element
c
c ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c
c x,y   = real arrays containing the $x$ and $y$ coordinates 
c	  resp. of the nodes.
c
c a,ja,ia= assembled matrix on output
c
c f      = right hand side (global load vector) 
c
c iwk,jwk = two integer work arrays.
c
c ierr	= error message integer . 
c	  ierr = 0 --> normal return
c	  ierr = 1 --> negative area encountered (due to bad 
c	           numbering of nodes of an element)
c
c xyk	= subroutine defining the material properties at each 
c         element. Form: 
c 	call xyk(nel,xyke,x,y,ijk,node) with on return
c         xyke =  material constant matrices. 
c         for each element nel, xyke(1,nel),xyke(2,nel) 
c         and xyke(3,nel) represent the constants
c         K11, K22, and K12 at that element.
c--------------------------------------------------------------
c
c moulitsa@cs.umn.edu : This routine yields the same results
c   as assmbo. It differs in that it constructs the ia array
c   by creating a list with the adjacent nodes for each node
c
c--------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension a(*),ijk(node,*),x(*),y(*),f(*),ske(3,3),fe(3),
     *      xe(3),ye(3),iwk(*),jwk(*), kwk(500)
      integer ia(*), ja(*) 
      external xyk, funb, func, fung
c--------------------------------------------------------------
c   initialize
c--------------------------------------------------------------
      do i=1,nx 
        f(i) = 0.0 
        iwk(i) = 0
        kwk(i) = 0
      end do

c iwk : how many elements a node belongs to
      do k=1,nelx
        do j=1,node
          knod = ijk(j,k)
          iwk(knod) = iwk(knod) + 1
        end do
      end do
c
c iwk : prepare for csr like format
      ksav=iwk(1)
      iwk(1)=1
      do j=2, nx+1
        ksavn = iwk(j)
        iwk(j) = iwk(j-1) +  ksav
        ksav = ksavn
      end do
c
c jwk : list of elements a node belongs to
      k=1
      do i=1,nelx
        do j=1,node
          knod = ijk(j,i)
          k=iwk(knod)
          jwk(k)=i
          iwk(knod)=iwk(knod)+1
        end do
      end do

c iwk : transform iwk back to what it was
      do i=nx+1,2,-1
        iwk(i)=iwk(i-1)
      end do
      iwk(1)=1

c kwk : mark edges that a node is associated with
      nedges=1
      ia(1)=1
      do i=1,nx
        kwk(i)=i
        do j=iwk(i), iwk(i+1)-1
          do k=1, node
            knod = ijk(k,jwk(j))
            if ( kwk(knod) .NE. i) then
              kwk(knod) = i
              nedges=nedges+1
            end if
          end do
        end do
        ia(i+1)=nedges
      end do
      do i=2,nx+1
        ia(i)=ia(i)+i-1
        iwk(i-1)=ia(i-1)-1
        jwk(i)=0
      end do
      jwk(1)=0
            
c-----------------
c main loop
c-----------------
      do nel=1, nelx
c
c get coordinates of nodal points
c
        do i=1, node
           j = ijk(i,nel)
          xe(i) = x(j)
          ye(i) = y(j)
        end do
c
c compute determinant
c
c         det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c
c set material properties
c 
c        call xyk(nel,xyke,x,y,ijk,node)
c
c construct element stiffness matrix
c
        ierr = 0
c
c        call evalg(nel, fe, xe, ye, fung, ierr)
c        call estif3(nel,ske,fe,det,xe,ye,xyke,ierr)
        call lstif3(ske, fe, xe, ye, xyk, funb, func, fung)
        if (ierr .ne. 0) return
c
c assemble: add element stiffness matrix to global matrix
c 
        do ka=1, node
          ii = ijk(ka,nel)
          f(ii) = f(ii) + fe(ka)
c
c unpack row into jwk1
c
          irowst = ia(ii)
          ilast  = iwk(ii) 
          do k=irowst,ilast 
            jwk(ja(k)) = k
          end do
c
          do kb = 1,node
c
c column number = jj
c
            jj = ijk(kb,nel)
            k = jwk(jj)
            if (k .eq. 0) then
              ilast = ilast+1
              jwk(jj) = ilast
              ja(ilast) = jj
              a(ilast) = ske(ka,kb)
            else 
              a(k) = a(k) + ske(ka,kb)
            endif
          end do
c refresh jwk
          do k=irowst,ilast 
            jwk(ja(k)) = 0
          end do
          iwk(ii) = ilast 
        end do
c
      end do

      return
      end
c----------------------------------------------------------------------- 
        subroutine chkelmt (x, y, nelx, ijk, node)
        implicit real*8 (a-h,o-z)
        dimension ijk(node,*),x(*),y(*)
c----------------------------------------------------------------------- 
c this subsourine checks the labeling within each elment and reorders
c the nodes in they ar not correctly ordered.
c----------------------------------------------------------------------- 
        do 1 nel =1, nelx
         det = x(ijk(2,nel))*(y(ijk(3,nel))-y(ijk(1,nel)))+
     *        x(ijk(3,nel))*(y(ijk(1,nel))-y(ijk(2,nel)))+
     *        x(ijk(1,nel))*(y(ijk(2,nel))-y(ijk(3,nel)))
c
c if determinant negative exchange last two nodes of elements.
c
        if (det .lt. 0.0d0) then
            j = ijk(2,nel)
            ijk(2,nel) = ijk(3,nel)
            ijk(3,nel) = j
        endif
 1        continue
c
        return                        
        end                         
c----------------------------------------------------------------------- 
      subroutine dlauny(x,y,nodes,elmnts,nemax,nelmnt)
      implicit real*8 (a-h,o-z)
c 
C code written by P.K. Sweby
c simple delauney triangulation routine (non optimal)
c
C     ******************************************************************
C     *                                                                *
C     * Performs a Delaunay triangularisation of a region given a set  *
C     * of mesh points.                                                *
C     *   X,Y    :- 1D arrays holding coordinates of mesh points.      *
C     *             dimensioned AT LEAST NODES+3.                      *
C     *   NODES  :- number of mesh points.                             *
C     *   ELMNTS :- INTEGER array, dimensioned NEMAX x 3, which on exit*
C     *             contains the index of global nodes associated with *
C     *             each element.                                      *
C     *   NELMNT :- on exit contains the number of elements in the     *
C     *             triangularisation.                                 *
C     *                                                                *
C     *                                   P.K.Sweby                    *
C     *                                                                *
C     ******************************************************************
C
      integer elmnts
      dimension x(nodes),y(nodes),elmnts(nemax,3)
      nn = 0 
c
c      pi=4.0*atan(1.0)
C
C     Calculate artificial nodes NODES+i i=1,2,3,4 and construct first
C     two (artificial) elements.
C
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      do 10 i=2,nodes
      xmin=min(xmin,x(i))
      xmax=max(xmax,x(i))
      ymin=min(ymin,y(i))
      ymax=max(ymax,y(i))
 10   continue
      dx=xmax-xmin
      dy=ymax-ymin
      xl=xmin-4.0*dx
      xr=xmax+4.0*dx
      yl=ymin-4.0*dy
      yr=ymax+4.0*dy
      x(nodes+1)=xl
      y(nodes+1)=yl
      x(nodes+2)=xl
      y(nodes+2)=yr
      x(nodes+3)=xr
      y(nodes+3)=yr
      x(nodes+4)=xr
      y(nodes+4)=yl
      elmnts(1,1)=nodes+1
      elmnts(1,2)=nodes+2
      elmnts(1,3)=nodes+3
      elmnts(2,1)=nodes+3
      elmnts(2,2)=nodes+4
      elmnts(2,3)=nodes+1
      nelmnt=2
      do 90 in=1,nodes
c
c     add one mesh point at a time and remesh locally if necessary
C
      ndel=0
      newel=0
      do 40 ie=1,nelmnt
c
c     is point in insided circumcircle of element IE ?
C
      i1=elmnts(ie,1)
      i2=elmnts(ie,2)
      i3=elmnts(ie,3)
      x2=x(i2)-x(i1)
      x3=x(i3)-x(i1)
      y2=y(i2)-y(i1)
      y3=y(i3)-y(i1)
      z=(x2*(x2-x3)+y2*(y2-y3))/(y2*x3-y3*x2)
      cx=0.5*(x3-z*y3)
      cy=0.5*(y3+z*x3)
      r2=cx**2+cy**2
      rn2=((x(in)-x(i1)-cx)**2+(y(in)-y(i1)-cy)**2)
      if(rn2.gt.r2)goto 40
c
c     yes it is inside,create new elements and mark old for deletion.
C
      do 30 j=1,3
      do 20 k=1,3
      elmnts(nelmnt+newel+j,k)=elmnts(ie,k)
 20   continue
      elmnts(nelmnt+newel+j,j)=in
 30   continue
      newel=newel+3
      elmnts(ie,1)=0
      ndel=ndel+1
c
 40   continue
C
C     If IN was inside circumcircle of more than 1 element then will
C     have created 2 identical new elements: delete them both.
C
      if(ndel.gt.1)then
          do 60 ie=nelmnt+1,nelmnt+newel-1
          do 60 je=ie+1,nelmnt+newel
          match=0
          do 50 k=1,3
          do 50 l=1,3
          if(elmnts(ie,k).eq.elmnts(je,l))match=match+1
 50       continue
          if(match.eq.3)then
              elmnts(ie,1)=0
              elmnts(je,1)=0
              ndel=ndel+2
          endif
 60       continue
      endif
c
c     delete any elements
C
      nn=nelmnt+newel
      ie=1
 70   continue
      if(elmnts(ie,1).eq.0)then
          do 80 j=ie,nn-1
          do 80 k=1,3
          elmnts(j,k)=elmnts(j+1,k)
 80       continue
          nn=nn-1
          ie=ie-1
      endif
      ie=ie+1
      if(ie.le.nn)goto 70
      nelmnt=nn
 90   continue
c
c     finally remove elements containing artificial nodes
C
      ie=1
 100  continue
      nart=0
      do 110 l=1,3
      if(elmnts(ie,l).gt.nodes)nart=nart+1
 110  continue
      if(nart.gt.0)then
          do 120 j=ie,nn-1
          do 120 k=1,3
          elmnts(j,k)=elmnts(j+1,k)
 120      continue
          nelmnt=nelmnt-1
          ie=ie-1
      endif
      ie=ie+1
      if(ie.le.nelmnt)goto 100
      return
      end
c-----------------------------------------------------------------------
      subroutine estif3(nel,ske,det,xe,ye,xyke,ierr)
c----------------------------------------------------------------------- 
c this subroutine constructs the element stiffness matrix for heat
c condution problem
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary
c
c using 3-node triangular elements arguments:
c nel        = element number
c ske        = element stiffness matrix
c [fe        = element load vector] unused here.
c det        = 2*area of the triangle
c xy, ye= coordinates of the three nodal points in an element.
c xyke  = material constants (kxx, kxy, kyx, kyy)
c
c------------------------------------------------------------------------
        implicit real*8 (a-h,o-z)
        dimension ske(3,3), xe(3), ye(3), dn(3,2),xyke(2,2)
c
c initialize
c
        area = 0.5*det
c
        do 200 i=1,3
        do 200 j=1,3
        ske(i,j) = 0.0d0
 200        continue
c
c get first gradient of shape function
c
        call gradi3(nel,xe,ye,dn,det,ierr)
        if (ierr .ne. 0) return
c
        do 100 i=1,3
        do 100 j=1,3
        t = 0.0d0
        do 102 k=1,2
        do 102 l=1,2
 102        t = t+xyke(k,l)*dn(i,k)*dn(j,l)
 100        ske(i,j) = t*area
c
        return
        end
c-------------------------------------------------------
        subroutine gradi3(nel, xe, ye, dn, det,ierr)
c-------------------------------------------------------
c constructs the first derivative of the shape functions.
c arguments:
c nel       = element nuumber
c xy, ye    = coordinates of the three nodal points in an element.
c dn        = gradients (1-st derivatives) of the shape functions.
c area      = area of the triangle
c
c------------------------------------------------------- 
        implicit real*8 (a-h,o-z)
        dimension xe(3), ye(3), dn(3,2)
        data eps/1.d-17/
c compute area
        ierr = 0 
        if (det .le. eps) goto 100
c
        dn(1,1) = (ye(2)-ye(3))/det
        dn(2,1) = (ye(3)-ye(1))/det
        dn(3,1) = (ye(1)-ye(2))/det
        dn(1,2) = (xe(3)-xe(2))/det
        dn(2,2) = (xe(1)-xe(3))/det
        dn(3,2) = (xe(2)-xe(1))/det
c
        return
c
 100        continue
        ierr = 3
        print *, 'ERROR:negative area encountered at elmt: ',nel
c        write(iout,*) det,(xe(i),ye(i),i=1,3)
        return
        end
c----------------------------------------------------------------------- 
        subroutine hsourc (indic,nelx,node,x,y,ijk,fs,f) 
        implicit real*8 (a-h,o-z) 
        real*8 x(*),y(*),fs(*),f(*),xe(3),ye(3),det,areao3
        integer ijk(node,*)
c
c generates the load vector f in assembled/unassembled form from the
c the element contributions fs. 
c indic = indicates if f is to be assembled (1) or not (zero) 
c note: f(*) not initilazed. because might use values from boundary 
c conditions.
c 
        jnod = 0
        do 130 nel = 1,nelx
c
c get coordinates of nodal points
c        
        do 104 i=1, node
        j = ijk(i,nel)
        xe(i) = x(j)
        ye(i) = y(j)
 104        continue
c
c compute determinant
c
        det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c area3 = area/3 
        areao3 = det/6.0
c 
c contributions to nodes in the element
c 
        if (indic .eq. 0) then
           do 115 ka=1,node
           jnod = jnod+1
           f(jnod) = fs(nel)*areao3
 115           continue
        else
            do 120 ka=1, node
            ii = ijk(ka,nel)
            f(ii) = f(ii) + fs(nel)*areao3
 120            continue
        endif 
c
 130        continue
        return
        end
c----- end of hsourc --------------------------------------------------- 
c-----------------------------------------------------------------------
      subroutine bound (nx,nelx,ijk,nodcode,node,nint,iperm,
     *             x,y,wk,iwk)
c-----------------------------------------------------------------------
c this routine counts the number of boundary points and 
c reorders the points in such a way that the boundary nodes
c are last.
c 
c nx, nelx, ijk, nodcode, node: see other subroutines
c iperm = permutation array from old orderin to new ordering,
c iwk   = reverse permutation array or return.
c wk	= real work array
c On return
c x, y, nodecode, are permuted
c ijk  is updated according to new oerdering.
c nint = number of interior points.
c 
c-----------------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension ijk(node,1),x(1),y(1),wk(1),iwk(1),iperm(1),
     *     nodcode(1)

c     put all boundary points at the end, backwards
      nint = 1
      nbound = nx
      do 1 j=1, nx
         if (nodcode(j) .eq. 0) then
            iperm(nint) = j
            nint = nint+1
         else
            iperm(nbound) = j
            nbound = nbound-1
         endif
 1    continue
c-------------------------------------------------------------------
      nint = nint-1
c     
c permute x's
c   
      do 2 k=1, nx
         wk(k) = x(k)
 2    continue
      do 3 k=1,nx
         x(k) = wk(iperm(k))
 3    continue
c 
c     permute the y's
c
      do 4 k=1, nx
         wk(k) = y(k)
 4    continue
      do 5 k=1, nx
         y(k) = wk(iperm(k))
 5    continue
c 
c     permute the boundary information
c
      do 6 k=1, nx
         iwk(k) = nodcode(k)
 6    continue 
      do 7 k=1,nx
         nodcode(k) = iwk(iperm(k))
 7    continue
c
c     get reverse permutation
c
      do 8 k=1, nx
         iwk(iperm(k)) = k
 8    continue
c
c     update the elements connectivity matrix
c
      do 10 nel = 1, nelx
         do 9 j=1, node
            knod = ijk(j,nel)
            ijk(j,nel) = iwk(knod) 
 9       continue
 10   continue
      return
      end
c-----------------------------------------------------------------------
      subroutine symbound (nx,nelx,ijk,nodcode,node,nint,
     *     iperm,iwk)
c-----------------------------------------------------------------------
c     this routine is a symbolic version of routine bound.
c     
c   nx, nelx, ijk, nodcode, node: see other subroutines
c   iperm = permutation array from old orderin to new ordering,
c   iwk   = reverse permutation array or return.
c   [wk	  = real work array] unused
c   On return
c   ijk   = is updated according to new oerdering.
c   nint  = number of interior points.
c 
c-----------------------------------------------------------------------
        implicit real*8  (a-h,o-z)
        dimension ijk(node,*),iwk(*),iperm(*),nodcode(*)

c put all boundary points at the end, backwards
        nint = 1
        nbound = nx
        do 1 j=1, nx
        if (nodcode(j) .eq. 0) then
          iperm(nint) = j
          nint = nint+1
        else
        iperm(nbound) = j
        nbound = nbound-1
        endif
 1        continue
c-------------------------------------------------------------------
        nint = nint-1
c 
c permute the boundary information
c
        do 6 k=1, nx
           iwk(k) = nodcode(k)
 6      continue 
        do 7 k=1,nx
           nodcode(k) = iwk(iperm(k))
 7        continue
c
c get reverse permutation
c
        do 8 k=1, nx
           iwk(iperm(k)) = k
 8        continue
c
c update the elements connectivity matrix
c
        do 10 nel = 1, nelx
           do 9 j=1, node
              knod = ijk(j,nel)
              ijk(j,nel) = iwk(knod) 
 9           continue
 10        continue
        return
        end    
c----------------------------------------------------------------------- 
        subroutine diric (nint,a,ja,ia)
c--------------------------------------------------------------
c this routine takes into account the boundary conditions
c and removes the unnecessary boundary points.
c--------------------------------------------------------------
        implicit real*8  (a-h,o-z)
        dimension a(*),ia(*),ja(*)
c call extract from UNARY
        call submat (1,1,nint,1,nint,a,ja,ia,nr,nc,a,ja,ia)
        write (*,*) 'nr=',nr,'nc=',nc
        return
c----------- end of diric ------------------------------------- 
        end
c-----------------------------------------------------------------------
        subroutine symdiric (nint,a,ja,ia)
c--------------------------------------------------------------
c this routine takes into account the boundary conditions
c and removes the unnecessary boundary points.
c--------------------------------------------------------------
        implicit real*8  (a-h,o-z)
        dimension a(*),ia(*),ja(*)
c     call submat from UNARY, with job = 0,
c     meaning no movement of real values.
        call submat (0,1,nint,1,nint,a,ja,ia,nr,nc,a,ja,ia)
        return
c----------- end of symdiric ------------------------------------- 
        end
c-----------------------------------------------------------------------
      subroutine cleannods (nx,x,y,nelx,ijk,node,nodcode,iperm) 
c      implicit none 
      integer nx,nelx,node,ijk(node,nelx),nodcode(*),iperm(nx) 
      real*8 x(nx),y(nx)
c-----------------------------------------------------------------------
c     this routine removes the nodes that do not belong to any element
c     (spurious points) and relabels the ijk array accordingly.
c-----------------------------------------------------------------------
      integer nel,i,k,j,indx
c
      do j=1, nx
         iperm(j) = 0
      enddo
c     
      do nel = 1, nelx
         do i=1,node
            k = ijk(i,nel) 
            iperm(k) = nel 
         enddo
      enddo
c
      indx = 0 
      do j =1, nx
         if (iperm(j) .ne. 0) then
            indx = indx+1
            iperm(indx) = j
            x(indx) = x(j)
            y(indx) = y(j) 
            nodcode(indx) = nodcode(j) 
         endif
      enddo
c     
c     update nx
c     
      nx = indx
c     
c     old number to new numbers
c     
      do j =1, nx 
         iperm(nx+iperm(j)) = j 
      enddo 
c     
c     
c     change all node numbers in ijk
c     
      do nel = 1, nelx
         do i=1,node
            k = ijk(i,nel) 
            k = iperm(nx+k) 
            ijk(i,nel) = k
         enddo
      enddo
      return
c-----------------------------------------------------------------------
c-----end-of-cleannod---------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine cleanel (nelx,ijk,node,nodcode,nodexc) 
c      implicit none
      integer nelx,node,nodexc,ijk(node,nelx),nodcode(*)
c-----------------------------------------------------------------------
c     this routine remove certain types of elements from the mesh
c     An element whose nodes are all labelled by the same label 
c     nodexc are removed. nelx is changed accordingly on return.
c-----------------------------------------------------------------------
      logical exclude
      integer nel, i,k 
      nel = 1
 1    continue
      exclude = .true. 
      do i=1,node
         k = ijk(i,nel) 
         exclude = (exclude .and. nodcode(k).eq. nodexc) 
      enddo
c     
      if (exclude) then
         do i=1,node
            ijk(i,nel) = ijk(i,nelx) 
         enddo
         nelx = nelx - 1
      else
         nel = nel+1
      endif
      if (nel .le. nelx) goto 1
      return
c-----------------------------------------------------------------------
c-----end-of-cleanel---------------------------------------------------- 
      end
c---------------------------------------------------------------------------
C  Piecewise linear fucntions on triangle.
      function psi(i,r,s)
      implicit real*8(a-h,o-z)
      
      goto (100,200,300) ,i
 100  psi = -(r+s)/2.
      return
 200  psi = (r+1.)/2.
      return
 300  psi = (s+1.)/2.
      return
      end
