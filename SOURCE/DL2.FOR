c**** new **************************************************************
      subroutine iclear(ia,m)
c-ZACE 2005.08
c
c.... program to clear an integer array
c
      dimension ia(*)
c
      do 100 i=1,m
      ia(i) = 0
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine igen(ia,m)
c
c.... program to read and generate integer nodal data
c
c        ia = input array
c         m = number of rows in ia
c         n = node number
c        ne = end node in generation sequence
c        ng = generation increment
c
      dimension ia(m,1),ib(13)
      common /iounit/ iin,iout,irsin,irsout
c
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)
      if (n.eq.0) return
      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
c
      do 200 i=n,ne,ng
      call imove(ia(1,i),ib,m)
  200 continue
c
      go to 100
c
 1000 format(16i5)
      end
c**** new **************************************************************
      subroutine imove(ia,ib,n)
c
c.... program to move an integer array
c
      dimension ia(*),ib(*)
c
      do 100 i=1,n
      ia(i)=ib(i)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine input(f,ndof,numnp,j,nlvect,iprtin,time)
c
c.... program to read, generate and write nodal input data
c
c        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
c                             = initial displacements (j=1)
c                             = initial velocities(j=2)
c                             = initial accelerations (j=3)
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical lzero
      dimension f(ndof,numnp,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      call clear(f,nlvect*numnp*ndof)
c
      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof)
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
c
      if (iprtin.eq.0) then
c
         if (lzero) then
            if (j.eq.0) write(iout,1000) nlv
            if (j.eq.1) write(iout,2000)
            if (j.eq.2) write(iout,3000)
            if (j.eq.3) write(iout,4000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
c
            if (j.eq.1)
     &      call printd(' i n i t i a l   d i s p l a c e m e n t s  ',
     &                  f,ndof,numnp,0,time)
c
            if (j.eq.2)
     &      call printd(' i n i t i a l   v e l o c i t i e s        ',
     &                  f,ndof,numnp,0,time)
c
            if (j.eq.3)
     &      call printd(' i n i t i a l   a c c e l e r a t i o n s  ',
     &                  f,ndof,numnp,0,time)
c
         endif
      endif
c
  100 continue
c
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',
     &    'kinematic boundary conditions for load vector number ',i5)
 2000 format('1'//,' there are no nonzero initial displacements')
 3000 format('1'//,' there are no nonzero initial velocities   ')
 4000 format('1'//,' there are no nonzero initial accelerations')
      end
c**** new **************************************************************
      subroutine interp(x,y,xx,yy,n)
c
c.... program to perform linear interpolation
c
c        x(i) = abscissas
c        y(i) = ordinates
c          xx = input abscissa
c          yy = output ordinate
c           n = total number of data points (1.le.i.le.n)
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension x(1),y(1)
c
      if (xx.le.x(1)) then
         yy = y(1)
c
      else if (xx.ge.x(n)) then
         yy = y(n)
c
      else
         do 100 i=1,n
         if (x(i).ge.xx) then
            yy = y(i-1) + (xx - x(i-1))*(y(i) - y(i-1))/(x(i) - x(i-1))
            return
         endif
  100    continue
c
      endif
c
      return
      end
c**** new **************************************************************
      subroutine iterup(id,d,dpred,vpred,a,brhs,ndof,numnp,ldyn)
c
c.... program to perform intermediate update of displacements,
c        velocities and accelerations during iterative loop in
c        predictor/corrector algorithm
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      dimension id(ndof,1),d(ndof,1),dpred(ndof,1),vpred(ndof,1),
     &          a(ndof,1),brhs(1)
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
c
      if (ldyn) then
c
         do 200 i=1,ndof
c
         do 100 j=1,numnp
         k = id(i,j)
         if (k.gt.0) then
            dpred(i,j) = dpred(i,j) + coeff5*brhs(k)
            vpred(i,j) = vpred(i,j) + coeff4*brhs(k)
            a(i,j) = a(i,j) + brhs(k)
         endif
  100    continue
c
  200    continue
c
      else
c
         do 400 i=1,ndof
c
         do 300 j=1,numnp
         k = id(i,j)
         if (k.gt.0) d(i,j) = brhs(k)
  300    continue
c
  400    continue
c
      endif
c
      return
      end
c**** new **************************************************************
      subroutine lfac(g,t,g1,nltftn,nptslf)
c
c.... program to compute load factors at time t
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension g(nptslf,2,1),g1(1)
c
      do 100 nlf=1,nltftn
      call interp(g(1,1,nlf),g(1,2,nlf),t,g1(nlf),nptslf)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine load(id,f,brhs,g1,ndof,numnp,nlvect)
c
c.... program to accumulate nodal forces and transfer into
c        right-hand-side vector
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension id(ndof,1),f(ndof,numnp,1),brhs(1),g1(1)
c
      do 300 i=1,ndof
c
      do 200 j=1,numnp
      k = id(i,j)
      if (k.gt.0) then
c
         do 100 nlv=1,nlvect
         brhs(k) = brhs(k) + f(i,j,nlv)*g1(nlv)
  100    continue
c
      endif
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine local(ien,x,xl,nen,nrowx,nrowxl)
c
c.... program to localize a global array
c
c        note: it is assumed nrowxl.le.nrowx
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension ien(1),x(nrowx,1),xl(nrowxl,1)
c
      do 200 j=1,nen
      node = ien(j)
c
      do 100 i=1,nrowxl
      xl(i,j)= x(i,node)
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      function lout(i,j)
c
c.... program to determine logical switch
c
      logical lout
c
      lout = .false.
      if (j.eq.0) return
      if (mod(i,j).eq.0) lout = .true.
c
      return
      end
c**** new **************************************************************
      subroutine ltimef(g,nptslf,nltftn,iprtin)
c
c.... program to read, write and store load-time functions
c
c        g(i,1,l) = time i for load-time function l
c        g(i,2,l) = load factor at time i for load-time function l
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension g(nptslf,2,nltftn)
      common /iounit/ iin,iout,irsin,irsout
c
      do 200 l=1,nltftn
c
      do 100 i=1,nptslf
      read(iin,1000) g(i,1,l),g(i,2,l)
  100 continue
c
  200 continue
c
      if (iprtin.eq.1) return
c
      write(iout,2000) nltftn
      do 400 l=1,nltftn
c
      do 300 i=1,nptslf
      if (mod(i,50).eq.1) write(iout,3000) l
      write(iout,4000) g(i,1,l),g(i,2,l)
  300 continue
c
  400 continue
c
      return
c
 1000 format(2f10.0)
 2000 format('1',' l o a d - t i m e   f u n c t i o n   d a t a ',//5x,
     &' number of load-time funtions . . . . . (nltftn  ) = ',i5       )
 3000 format(///5x,' function number ',i5,//
     &       16x,'time',13x,'load factor'/)
 4000 format(5x,2(1pe20.8))
c
      end
c**** new **************************************************************
      subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
c
c.... program to add rectangular matrices
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(ma,1),b(mb,1),c(mc,1)
c
      go to (1000,2000,3000),iopt
c
c.... iopt = 1, add entire matrices
c
 1000 do 1200 j=1,n
c
      do 1100 i=1,m
      c(i,j) = a(i,j) + b(i,j)
 1100 continue
c
 1200 continue
      return
c
c.... iopt = 2, add lower triangular and diagonal elements
c
 2000 do 2200 j=1,n
c
      do 2100 i=j,m
      c(i,j) = a(i,j) + b(i,j)
 2100 continue
c
 2200 continue
      return
c
c.... iopt = 3, add upper triangular and diagonal elements
c
 3000 do 3200 j=1,n
c
      do 3100 i=1,j
      c(i,j) = a(i,j) + b(i,j)
 3100 continue
c
 3200 continue
      return
c
      end
c**** new **************************************************************
      subroutine meansh(shgbar,w,det,r,shg,nen,nint,iopt,nesd,nrowsh)
c
c.... program to calculate mean values of shape function
c        global derivatives for b-bar method
c
c        note: if iopt.eq.2, det(l) = det(l)*r(l) upon entry
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shgbar(3,1),w(1),det(1),r(1),shg(nrowsh,nen,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      call clear(shgbar,3*nen)
c
      volinv = one/coldot(w,det,nint)
c
      do 300 l=1,nint
      temp1 = w(l)*det(l)*volinv
      if (iopt.eq.2) temp2 = temp1/r(l)
c
      do 200 j=1,nen
c
      do 100 i=1,nesd
      shgbar(i,j) = shgbar(i,j) + temp1*shg(i,j,l)
  100 continue
c
      if (iopt.eq.2) shgbar(3,j) = shgbar(3,j) + temp2*shg(3,j,l)
  200 continue
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine minmax(x,xmax,xmin,l,m,n)
c
c.... program to compute the min and max in the row of a matrix
c
c        x = matrix
c        l = number of rows in x
c        m = number of columns in x
c        n = row number
c
      dimension x(l,1)
c
      xmax = x(n,1)
      xmin = x(n,1)
c
      do 100 i = 2,m
        if (x(n,i).gt.xmax) xmax = x(n,i)
        if (x(n,i).lt.xmin) xmin = x(n,i)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine move(a,b,n)
c
c.... program to move a floating-point array
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1),b(1)
c
      do 100 i=1,n
      a(i) = b(i)
  100 continue
c
      return
      end
c**** new **************************************************************
      function mpoint(name,ndim1,ndim2,ndim3,ipr)
c
c.... program to calculate storage pointer
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      character*4 name
      dimension name(2)
      common /bpoint/ mfirst,mlast,mtot,iprec
c
      mpoint = mfirst
      if ( iprec.eq.2 .and. mod(mpoint,2).eq.0 ) mpoint = mpoint + 1
      call dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast)
      mfirst = mpoint + ndim1*max0(1,ndim2)*max0(1,ndim3)*ipr
      if (mfirst.ge.mlast) call serror(name,mfirst-mlast)
c
      return
      end
c**** new **************************************************************
      subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
c
c.... program to multiply two matrices
c
c        l = range of dot-product index
c        m = number of active rows in c
c        n = number of active columns in c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(ma,1),b(mb,1),c(mc,1)
c
      go to (1000,2000,3000,4000),iopt
c
c.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
c
 1000 do 1200 i=1,m
c
      do 1100 j=1,n
      c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
 1100 continue
c
 1200 continue
      return
c                                            t
c.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
c
 2000 do 2200 i=1,m
c
      do 2100 j=1,n
      c(i,j) = coldot(a(1,i),b(1,j),l)
 2100 continue
c
 2200 continue
      return
c                                                t
c.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
c
 3000 do 3200 i=1,m
c
      do 3100 j=1,n
      c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
 3100 continue
c
 3200 continue
      return
c                                            t    t
c.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
c
 4000 do 4200 i=1,m
c
      do 4100 j=1,n
      c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
 4100 continue
c
 4200 continue
c
      return
      end
c**** new **************************************************************
      subroutine pivots(a,idiag,neq,nsq,*)
c
c.... program to determine the number of zero and negative terms in
c        array d of factorization a = u(transpose) * d * u
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1),idiag(1)
      common /iounit/ iin,iout,irsin,irsout
c
      iz = 0
      in = 0
c
      do 100 n=1,neq
      i = idiag(n)
      if (a(i).eq.0.) iz = iz + 1
      if (a(i).lt.0.) in = in + 1
  100 continue
c
      write(iout,1000) nsq,iz,in
c
      return 1
c
 1000 format(' ',
     &' zero and/or negative pivots encountered                ', ///5x,
     &' time sequence number   . . . . . . . . . . . (nsq  ) = ',i5//5x,
     &' number of zeroes . . . . . . . . . . . . . . . . . . = ',i5//5x,
     &' number of negatives  . . . . . . . . . . . . . . . . = ',i5//5x)
c
      end
c**** new **************************************************************
      subroutine predct(d,v,a,dpred,vpred,ndof,numnp)
c
c.... program to calculate predictor for displacements,velocities
c        and accelerations
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension d(ndof,1),v(ndof,1),a(ndof,1),
     &          dpred(ndof,1),vpred(ndof,1)
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 200 i=1,ndof
c
      do 100 j=1,numnp
      dpred(i,j) = d(i,j) + coeff6*v(i,j) + coeff7*a(i,j)
      vpred(i,j) = v(i,j) + coeff8*a(i,j)
      a(i,j) = zero
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine princ(n,s,p)
c
c.... program to compute principal values of symmetric 2nd-rank tensor
c
c        s = symmetric second-rank tensor stored as a vector
c        n = number of dimensions (2 or 3)
c        p = vector of principal values
c
c.... the components of s must be stored in the following orders
c
c        2-d problems: s11,s22,s12
c        3-d problems: s11,s22,s33,s12,s23,s31
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension s(1),p(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      data rt2/1.41421356237309/,pi23/2.09439510239321/,
     &     two2/22.50/,four5/45.0/
c
      if (n.eq.2) then
c
c....... 2-d problem
c
         a = two2/atan(one)
         x = pt5*(s(1) + s(2))
         y = pt5*(s(1) - s(2))
         r = sqrt(y*y + s(3)*s(3))
         p(1) = x + r
         p(2) = x - r
         p(3) = r
         p(4) = four5
         if (y.ne.zero .or. s(3).ne.zero) p(4) = a*atan2(s(3),y)
      endif
c
      if (n.eq.3) then
c
c....... 3-d problem
c
  100    r = zero
         x = (s(1) + s(2) + s(3))/three
         y = s(1)*(s(2) + s(3)) + s(2)*s(3)
     &       - s(4)*s(4) - s(6)*s(6) - s(5)*s(5)
         z = s(1)*s(2)*s(3) - two*s(4)*s(6)*s(5) - s(1)*s(5)*s(5)
     &       - s(2)*s(6)*s(6) - s(3)*s(4)*s(4)
         t = three*x*x - y
         u = zero
         if (t.ne.zero) then
            u = sqrt(two*t/three)
            ucubed = u*u*u
            if (ucubed.ne.zero) then
                a = (z + (t - x*x)*x)*rt2/ucubed
                r = sqrt(abs(one - a*a))
                if (r.ne.zero.or.a.ne.zero) then
                    r = atan2(r,a)/three
                else
                    r = zero
                endif
            else
                u = zero
            endif
         endif
         p(1) = x + u*rt2*cos(r)
         p(2) = x + u*rt2*cos(r - pi23)
         p(3) = x + u*rt2*cos(r + pi23)
      endif
c
      return
      end
c**** new **************************************************************

      subroutine printd(name,dva,ndof,numnp,ntstep,time)
c
c.... program to print kinematic data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical lzero,lskip
c ZACE correction 21.05.1997
	character*(*) name
      dimension dva(ndof,1) ! ,name(11)
      common /iounit/ iin,iout,irsin,irsout
c
      nn = 0
      lskip = .true.
c
      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1)
     &      write(iout,1000) name,ntstep,time,(i,i=1,ndof)
         write(iout,2000) n,(dva(i,n),i=1,ndof)
         lskip = .false.
      endif
  100 continue
c
      if (lskip) then
         write(iout,1000) name,ntstep,time,(i,i=1,ndof)
         write(iout,3000)
      endif
c
      return
c
 1000 format('1',a44//5x, !ZACE correction 21.05.1997   11a4//5x, 
     &' step number = ',i10//5x,
     &' time        = ',1pe10.3///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i5,10x,6(1pe15.8,2x))
 3000 format(' ',//,' there are no nonzero components')
      end
c**** new **************************************************************
      subroutine printf(f,ndof,numnp,nlv)
c
c.... program to print prescribed force and boundary condition data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical lzero
      dimension f(ndof,numnp,1)
      common /iounit/ iin,iout,irsin,irsout
c
      nn = 0
c
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1)
     &      write(iout,1000) nlv,(i,i=1,ndof)
         write(iout,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
c
      return
c
 1000 format('1',
     &' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',
     &'  b o u n d a r y   c o n d i t i o n s'//5x,
     &' load vector number = ',i5///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i5,10x,6(1pe15.8,2x))
      end
c**** new **************************************************************
      subroutine printp(a,idiag,neq,nsq,*)
c
c.... program to print array d after Crout factorization
c        a = u(transpose) * d * u
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension a(1),idiag(1)
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,neq
      if (mod(n,50).eq.1) write(iout,1000) nsq
      i = idiag(n)
      write(iout,2000) n,a(i)
  100 continue
c
      return 1
c
 1000 format('1',' array d of factorization',
     &' a = u(transpose) * d * u ',                                //5x,
     &' time sequence number   . . . . . . . . . . . . (nsq) = ',i5//5x)
 2000 format(1x,i5,4x,1pe20.8)
      end
c**** new **************************************************************
      subroutine prntel(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(1),ien(nen,1)
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iout,1000) (i,i=1,nen)
      write(iout,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format('1',
     &' e l e m e n t   d a t a',//5x,
     &' element   material',9('  node ',i1,:,2x),/5x,
     &'  number    number'//)
 2000 format(6x,i5,9(5x,i5))
      end
c**** new **************************************************************
      subroutine prop2d(rho,rdampm,rdampk,th,c,numat,iopt,nrowb)
c
c.... program to read, write and store properties for two-dimensional
c        continuum elements
c
c        note: this routine is presently restricted to the
c              isotropic linearly-elastic case
c
c              iopt = 0; plane stress
c                   = 1; plane strain
c                   = 2; torsionless axisymmetric
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension rho(1),rdampm(1),rdampk(1),th(1),c(nrowb,nrowb,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,numat
      if (mod(n,50).eq.1) write(iout,1000) numat
      read(iin,2000) m,e,pois,rho(m),rdampm(m),rdampk(m),th(m)
      if (th(m).eq.zero) th(m) = one
      write(iout,3000) m,e,pois,rho(m),rdampm(m),rdampk(m),th(m)
c
c.... set material constants for out-of-plane components
c
      amu2 = e/(one + pois)
      alam = amu2*pois/(one - two*pois)
c
      c(1,4,m) = alam
      c(2,4,m) = alam
      c(3,4,m) = zero
      c(4,4,m) = alam + amu2
c
      c(4,1,m) = c(1,4,m)
      c(4,2,m) = c(2,4,m)
      c(4,3,m) = c(3,4,m)
c
c.... set material constants for in-plane components
c
      if (iopt.eq.0) alam = alam*amu2/(alam + amu2)
c
      c(1,1,m) = alam + amu2
      c(1,2,m) = alam
      c(2,2,m) = c(1,1,m)
      c(1,3,m) = zero
      c(2,3,m) = zero
      c(3,3,m) = pt5*amu2
c
      c(2,1,m) = c(1,2,m)
      c(3,1,m) = c(1,3,m)
      c(3,2,m) = c(2,3,m)
c
  100 continue
c
      return
c
 1000 format('1',
     &' m a t e r i a l   s e t   d a t a                      '   //5x,
     &' number of material sets . . . . . . . . . . (numat ) = ',i5///,
     & 7x,'set',5x,'young''s',4x,'poisson''s',5x,'mass',8x,'mass',
     & 6x,'stiffness',3x,'thickness',/6x,'number',3x,'modulus',
     & 6x,'ratio',6x,'density',5x,'damping',5x,'damping',/)
 2000 format(i5,5x,7f10.0)
 3000 format(4x,i5,3x,6(2x,1pe10.4))
      end
c**** new **************************************************************
      subroutine prtdc
c
c.... program to print memory-pointer dictionary
c
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /iounit/ iin,iout,irsin,irsout
c-ZACE-2005.08
      include 'memory_size.inc'
      common ia(max_size)
c
      n = (mtot-mlast)/7
      j = mtot + 1
c
      do 100 i=1,n
      if (mod(i,50).eq.1) write(iout,1000)
      j = j - 7
      call prtdc1(i,ia(j),ia(j+2),ia(j+3),ia(j+4),ia(j+5),ia(j+6))
  100 continue
c
      return
c
 1000 format('1',
     &' d y n a m i c   s t o r a g e    a l l o c a t i o n',
     &'   i n f o r m a t i o n '//
     &  12x,'array no.',5x,'array',8x,'address',6x,'dim1',6x,'dim2',
     &  6x, 'dim3',6x,'prec.'/)
c
      end
c**** new **************************************************************
      subroutine prtdc1(i,iname,iadd,ndim1,ndim2,ndim3,ipr)
c
c.... program to print memory-pointer information for an array
c
      character*4 name(2),nelpar,lefths
      dimension iname(2)
      common /iounit/ iin,iout,irsin,irsout
cZACE MODIF (RD  6 jan 1987) -----ordre save/data + character name
      save neg
      data nelpar,lefths/'npar','alhs'/
c
      if (i.eq.1) neg = 1
      call intchr (iname(1),name(1))
      call intchr (iname(2),name(2))
      if (name(1).eq.nelpar) then
        write (iout,1000) neg
        neg = neg + 1
      endif
      if (name(1).eq.lefths) write (iout,2000)
      write(iout,3000) i,name(1),name(2),iadd,ndim1,ndim2,ndim3,ipr
c
      return
c
 1000 format(/14x,'*****',7x,'begin element group number',i5/' ')
 2000 format(/14x,'*****',7x,'end element group data',/' ')
 3000 format(14x,i5,7x,2a4,1x,6i10)
      end
c**** new **************************************************************
      subroutine prts2d(xint,stress,pstrs,strain,pstrn,
     &                  nn,nntot,neg,nel,lint)
c
c.... program to print stress, strain, and principal values
c        for two-dimensional continuum elements
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension xint(2),stress(4),pstrs(4),strain(4),pstrn(4)
      common /iounit/ iin,iout,irsin,irsout
c
      nn = nn+1
      if (mod(nn,nntot).eq.1) write(iout,1000) neg
      write(iout,2000) nel,lint,xint,stress,pstrs,strain,pstrn
c
      return
c
 1000 format('1',
     &' e l e m e n t  s t r e s s e s   a n d   s t r a i n s ',  //5x,
     &' element group number  . . . . . . . . . . . (neg   ) = ',i5///
     &'  element  int. pt.         x1        x2    ',5x,
     &'  stress    stress    stress    stress  ',
     &'  princ.    princ.    shear     stress  ',  /,
     &'   number   number                          ',5x,
     &'    11        22        12        33    ',
     &' stress 1  stress 2   stress     angle  ',//,49x,
     &'  strain    strain    strain    strain  ',
     &'  princ.    princ.    shear     strain  ', /,49x,
     &'    11        22        12        33    ',
     &' strain 1  strain 2   strain     angle  ')
 2000 format(/2x,i5,6x,i2,8x,2(1pe10.2),5x,8(1pe10.2)/48x,8(1pe10.2))
      end
c**** new **************************************************************
      function rcdot(a,b,ma,n)
c
c.... program to compute the dot product of a vector stored row-wise
c        with a vector stored column-wise
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(ma,1),b(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      rcdot = zero
c
      do 100 i=1,n
      rcdot = rcdot + a(1,i)*b(i)
  100 continue
c
      return
      end
c**** new **************************************************************
      function rowdot(a,b,ma,mb,n)
c
c.... program to compute the dot product of vectors stored row-wise
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension a(ma,1),b(mb,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      rowdot = zero
c
      do 100 i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine rsin(d,v,a,ndof,numnp,ntstep,time)
c
c.... program to read restart file
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension d(ndof,1),v(ndof,1),a(ndof,1)
      common /iounit/ iin,iout,irsin,irsout
c
      read(irsin,1000) ntstep,time,ijunk
c
      do 100 j=1,numnp
      read(irsin,2000) (d(i,j),v(i,j),a(i,j),i=1,ndof)
 100  continue
c
      return
c
 1000 format(//,15x,i5,/,15x,e12.5/,i1)
 2000 format(3(6e16.8/))
      end
c**** new **************************************************************
      subroutine rsout(d,v,a,ndof,numnp,ntstep,time)
c
c.... program to write restart file
c
      double precision a,d,v,time
c
c.... deactivate above card(s)s for single precision operation
c
      character*4 title
      dimension d(ndof,1),v(ndof,1),a(ndof,1)
      common /iounit/ iin,iout,irsin,irsout
      common /titlec/ title(20)
c
      write(irsout,1000) title,ntstep,time
c
      do 100 j=1,numnp
      write(irsout,2000) (d(i,j),v(i,j),a(i,j),i=1,ndof)
 100  continue
c
      return
c
 1000 format(' ',20a4//,' step number = ',i5/,
     &                  ' time        = ',1pe12.5/' ')
 2000 format(3(6e16.8/))
      end
c**** new **************************************************************
      subroutine serror(name,i)
c
c.... program to print error message if available storage is exceeded
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      character*4 name
      dimension name(2)
      common /iounit/ iin,iout,irsin,irsout
c
      call prtdc
      write(iout,1000) i,name(1),name(2)
      stop
c
 1000 format(1x,5('*'),'storage exceeded by ',i10,
     &' words in attempting to store array ',2a4)
      end
c***** end *************************************************************
      subroutine setupd(c,dmat,const,nstr,nrowb)
c
c.... program to calculate the d matrix
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension c(nrowb,1),dmat(nrowb,1)
c
      do 200 j=1,nstr
c
      do 100 i=1,j
      dmat(i,j) = const*c(i,j)
      dmat(j,i) = dmat(i,j)
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine smult(a,b,c,mb,mc,m,n,iopt)
c
c.... program to perform scalar multiplication of a matrix
c
c        c(i,j) = a*b(i,j)
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension b(mb,1),c(mc,1)
c
      go to (1000,2000,3000),iopt
c
c.... iopt = 1, multiply entire matrix
c
 1000 do 1200 j=1,n
c
      do 1100 i=1,m
      c(i,j) = a*b(i,j)
 1100 continue
c
 1200 continue
      return
c
c.... iopt = 2, multiply lower triangular and diagonal elements
c
 2000 do 2200 j=1,n
c
      do 2100 i=j,m
      c(i,j) = a*b(i,j)
 2100 continue
c
 2200 continue
      return
c
c.... iopt = 3, multiply upper triangular and diagonal elements
c
 3000 do 3200 j=1,n
c
      do 3100 i=1,j
      c(i,j) = a*b(i,j)
 3100 continue
c
 3200 continue
      return
c
      end
c**** new **************************************************************
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      subroutine stored(idhist,d,v,a,dout,ndof,ndout,ndout1)
c
c.... program to store nodal time histories as single-precision data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      real dout(ndout1,1)
      dimension idhist(3,1),d(ndof,1),v(ndof,1),a(ndof,1)
      common /hplotc/ nplpts,locplt,time
c
      dout(1,locplt) = real(time)
c
      do 100 i=1,ndout
      node = idhist(1,i)
      idof = idhist(2,i)
      idva = idhist(3,i)
      if (idva .eq. 1) dout(i+1,locplt) = real(d(idof,node))
      if (idva .eq. 2) dout(i+1,locplt) = real(v(idof,node))
      if (idva .eq. 3) dout(i+1,locplt) = real(a(idof,node))
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine timcon(nsq,nstep ,ndprt ,nsprt ,nhplt ,niter ,
     &                      nstep1,ndprt1,nsprt1,nhplt1,niter1,
     &                      alpha ,beta  ,gamma ,dt    )
c
c.... program to compute current time sequence parameters
c        and time-integration coefficients
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension nstep(1),ndprt(1),nsprt(1),nhplt(1),niter(1),
     &          alpha(1),beta(1) ,gamma(1) ,dt(1)
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      nstep1 = nstep(nsq)
      ndprt1 = ndprt(nsq)
      nsprt1 = nsprt(nsq)
      nhplt1 = nhplt(nsq)
      niter1 = niter(nsq)
      alpha1 = alpha(nsq)
      beta1  = beta(nsq)
      gamma1 = gamma(nsq)
      dt1    = dt(nsq)
c
      coeff1 = one + alpha1
      coeff2 = gamma1*dt1
      coeff3 = beta1*dt1*dt1
      coeff4 = coeff1*coeff2
      coeff5 = coeff1*coeff3
      coeff6 = coeff1*dt1
      coeff7 = pt5*coeff1*(one - two*beta1)*dt1*dt1
      coeff8 = coeff1*(one - gamma1)*dt1
c
      return
      end
c**** new **************************************************************
      subroutine timing(time)
c
c.... program to determine elapsed cpu time
c.... **** this is a system-dependent routine ****
c....     note: can only access clock time on vax/vms
c      time = secnds(0.0)

      integer*4 hr,min,sec,hund
      ! microsoft Fortran PowerStation v.4.0
      call gettim(hr,min,sec,hund)
      time = 3600.*hr + 60.*min + sec + hund * 0.01

      ! Fortran Lahey 77/90
c       call timer ( iticks )
c       time = 0.01 * iticks

      return
      end
c**** new **************************************************************
      subroutine timlog
c
c.... program to print log of execution times
c
      character*4 title
      common /etimec/ etime(7)
      common /iounit/ iin,iout,irsin,irsout
      common /titlec/ title(20)
c
      subtot = 0.0
      do 100 i=3,7
      subtot = subtot + etime(i)
  100 continue
c
      write(iout,1000) title,etime,subtot
c
      return
c
 1000 format('1',20a4///5x,
     &' e x e c u t i o n   t i m i n g   i n f o r m a t i o n'  ///5x,
     &' i n i t i a l i z a t i o n   p h a s e        = ',1pe10.3///5x,
     &' s o l u t i o n   p h a s e                    = ',1pe10.3///5x,
     &'     formation of left-hand-side matrices       = ',1pe10.3 //5x,
     &'     factorizations                             = ',1pe10.3 //5x,
     &'     formation of right-hand-side vectors       = ',1pe10.3 //5x,
     &'     forward reductions/back substitutions      = ',1pe10.3 //5x,
     &'     calculation of element output              = ',1pe10.3  /5x,
     &51x,'_________',//5x,
     &'     subtotal                                   = ',1pe10.3     )
c
      end
c**** new **************************************************************
      subroutine tseqin(nstep,ndprt,nsprt,nhplt,niter,alpha,beta,
     &                  gamma,dt,numseq,nplpts,ldyn)
c
c.... program to read, write and store time sequence data
c
c        note: "nplpts" is passed to subroutine hplot by way of
c               common /hplotc/
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical ldyn
      dimension nstep(1),ndprt(1),nsprt(1),nhplt(1),niter(1)
     &         ,alpha(1),beta(1),gamma(1),dt(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      nplpts = 1
c
      do 100 i=1,numseq
      read(iin,1000) n,nstep(n),ndprt(n),nsprt(n),nhplt(n),niter(n),
     &                 alpha(n),beta(n) ,gamma(n),dt(n)
      if (nhplt(n).gt.0) nplpts = nplpts + nstep(n)/nhplt(n)
  100 continue
c
c.... set default sequence parameters for static analysis
c
      if (.not.ldyn) then
         do 200 i=1,numseq
         nstep(i) = max0(1,nstep(i))
         ndprt(i) = 1
         nsprt(i) = 1
         nhplt(i) = 0
         niter(i) = 1
         alpha(i) = zero
         beta(i)  = one
         gamma(i) = zero
         dt(i)    = one
  200    continue
      endif
c
      do 300 n=1,numseq
      if (mod(n,2).eq.1) write(iout,2000) numseq
      write(iout,3000) n,nstep(n),ndprt(n),nsprt(n),nhplt(n),niter(n),
     &                   alpha(n),beta(n) ,gamma(n),dt(n)
  300 continue
c
      return
c
 1000 format(6i5,4f10.0)
 2000 format('1',' t i m e   s e q u e n c e   d a t a      ',     //5x,
     &' number of time sequences . . . . . . (numseq  ) = ',    i5///  )
 3000 format(5x,
     &' time sequence number . . . . . . . . (n       ) = ',     i5//5x,
     &' number of time steps . . . . . . . . (nstep(n)) = ',     i5//5x,
     &' kinematic print increment  . . . . . (ndprt(n)) = ',     i5//5x,
     &' stress/strain print increment  . . . (nsprt(n)) = ',     i5//5x,
     &' time history plot increment  . . . . (nhplt(n)) = ',     i5//5x,
     &' number of iterations . . . . . . . . (niter(n)) = ',     i5//5x,
     &' first integration parameter  . . . . (alpha(n)) = ',1pe12.5//5x,
     &' second integration parameter . . . . (beta(n) ) = ',1pe12.5//5x,
     &' third integration parameter  . . . . (gamma(n)) = ',1pe12.5//5x,
     &' time step  . . . . . . . . . . . . . (dt(n)   ) = ',1pe12.5////)
c
      end
c**** new **************************************************************
      subroutine ztest(a,n,lzero)
c
c.... program to determine if an array contains only zero entries
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension a(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      logical lzero
c
      lzero = .true.
c
      do 100 i=1,n
      if (a(i).ne.zero) then
         lzero = .false.
         return
      endif
  100 continue
c
      return
      end
