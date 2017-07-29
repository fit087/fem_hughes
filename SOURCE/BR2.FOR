c**** new **************************************************************
      subroutine bkcb(shg,shgbar,b,nrowsh,nrowb,nen,ibbar)
c
c.... program to set up the strain-displacement matrix "b" for
c        three-dimensional continuum elements
c
c        ibbar = 0, standard b-matrix
c
c        ibbar = 1, mean-dilatational b-matrix
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shg(nrowsh,1),shgbar(nrowsh,1),b(nrowb,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 100 j=1,nen
c
      j3   = 3*j
      j3m1 = j3 - 1
      j3m2 = j3 - 2
c
      b(1,j3m2) = shg(1,j)
      b(1,j3m1) = zero
      b(1,j3  ) = zero
      b(2,j3m2) = zero
      b(2,j3m1) = shg(2,j)
      b(2,j3  ) = zero
      b(3,j3m2) = zero
      b(3,j3m1) = zero
      b(3,j3  ) = shg(3,j)
      b(4,j3m2) = shg(2,j)
      b(4,j3m1) = shg(1,j)
      b(4,j3  ) = zero
      b(5,j3m2) = zero
      b(5,j3m1) = shg(3,j)
      b(5,j3  ) = shg(2,j)
      b(6,j3m2) = shg(3,j)
      b(6,j3m1) = zero
      b(6,j3  ) = shg(1,j)
c
  100 continue
c
      if (ibbar.eq.0) return
c
c.... add contributions to form b-bar
c
      constb = one/three
c
      do 200 j=1,nen
c
      j3   = 3*j
      j3m1 = j3 - 1
      j3m2 = j3 - 2
c
      temp1 = constb*(shgbar(1,j) - shg(1,j))
      temp2 = constb*(shgbar(2,j) - shg(2,j))
      temp3 = constb*(shgbar(3,j) - shg(3,j))
c
      b(1,j3m2) = b(1,j3m2) + temp1
      b(1,j3m1) = b(1,j3m1) + temp2
      b(1,j3  ) = b(1,j3  ) + temp3
      b(2,j3m2) = b(2,j3m2) + temp1
      b(2,j3m1) = b(2,j3m1) + temp2
      b(2,j3  ) = b(2,j3  ) + temp3
      b(3,j3m2) = b(3,j3m2) + temp1
      b(3,j3m1) = b(3,j3m1) + temp2
      b(3,j3  ) = b(3,j3  ) + temp3
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine bkck(shgbar,w,det,shg,b,c,dmat,db,elstif,constk,
     &                ibbar,nen,nint,nesd,nrowsh,nrowb,nstr,nee)
c
c.... program to form stiffness matrix for a continuum element
c        with "nen" nodes
c
c        note: the b-bar option is restricted to the mean-dilatation
c              formulation. to generalize to other formulations,
c              redimension array "shgbar", and replace routines
c              "meansh" and "bkcb".
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shgbar(nrowsh,nen,1),w(1),det(1),shg(nrowsh,nen,1),
     &          b(nrowb,1),c(nrowb,1),dmat(nrowb,1),db(nrowb,1),
     &          elstif(nee,1)
c
c.... calculate mean values of shape function global derivatives
c        for mean-dilatational b-bar formulation
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,0,nesd,nrowsh)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call bkcb(shg(1,1,l),shgbar,b,nrowsh,nrowb,nen,ibbar)
c
c.... set up the constitutive matrix
c
      call setupd(c,dmat,temp,nstr,nrowb)
c
c.... multiply d*b
c
      call multab(dmat,b,db,nrowb,nrowb,nrowb,nstr,nstr,nee,1)
c
c.... multiply b(transpose) * db, taking account of symmetry,
c        and accumulate in elstif
c
      call btdb(elstif,b,db,nee,nrowb,nstr)
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine bkckd(shgbar,w,det,shg,b,dl,strain,c,stress,work,
     &                 elresf,constk,ibbar,nen,nint,nrowsh,
     &                 nesd,nrowb,nee,nstr)
c
c.... program to form internal force ("-k*d") for a continuum element
c        with "nen" nodes
c
c        note: the b-bar option is restricted to the mean-dilatation
c              formulation. to generalize to other formulations,
c              redimension array "shgbar", and replace routines
c              "meansh" and "bkcb".
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shgbar(nrowsh,nen,1),w(1),det(1),shg(nrowsh,nen,1),
     &          b(nrowb,1),dl(1),strain(1),c(nrowb,1),stress(1),
     &          work(1),elresf(1)
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,0,nesd,nrowsh)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call bkcb(shg(1,1,l),shgbar,b,nrowsh,nrowb,nen,ibbar)
c
c.... calculate strains
c
      call multab(b,dl,strain,nrowb,nee,nstr,nee,nstr,1,1)
c
c.... calculate stresses
c
      call multab(c,strain,stress,nrowb,nstr,nstr,nstr,nstr,1,1)
c
c.... calculate element internal force
c
      call smult(temp,stress,stress,nstr,nstr,nstr,1,1)
      call multab(b,stress,work,nrowb,nstr,nee,nstr,nee,1,2)
      call matadd(elresf,work,elresf,nee,nee,nee,nee,1,1)
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine bkcrsf(ielno,iside,press,nsurf)
c
c.... program to read, write and store surface force data for the
c        eight-node brick, elastic continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension ielno(1),iside(1),press(4,1)
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,nsurf
      if (mod(n,50).eq.1) write(iout,1000) nsurf
      read(iin,2000) ielno(n),iside(n),press(1,n),press(2,n),
     &               press(3,n),press(4,n)
      write(iout,3000) ielno(n),iside(n),press(1,n),press(2,n),
     &                press(3,n),press(4,n)
  100 continue
c
      return
c
 1000 format('1',
     &' e l e m e n t   s u r f a c e   f o r c e   d a t a    ',  //5x,
     &' number of surface force cards . . . . . . . (nsurf ) = ',i5///
     &5x,' element    side    ',4('   pressure   '),/,
     &5x,2('  number  '),'    node i        node j    ',
     &                   '    node k        node l    ',/)
 2000 format(2i5,4f10.0)
 3000 format(6x,i5,7x,i2,3x,4(2x,e12.4))
      end
c**** new **************************************************************
      subroutine bkcshg(xl,det,shl,shg,nint,nel,neg)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a eight-node brick element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local ("zeta") derivative of shape function
c        shl(4,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = z-derivative of shape function
c        shg(4,i,l) = shl(4,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 8
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension xl(3,1),det(1),shl(4,8,1),shg(4,8,1),xs(3,3),xsj(3,3)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      call move(shg,shl,32*nint)
c
      do 700 l=1,nint
c
      do 300 j=1,3
      do 200 i=1,3
      xs(i,j) = rowdot(shg(i,1,l),xl(j,1),4,3,8)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*(xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))
     &        +xs(1,2)*(xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3))
     &        +xs(1,3)*(xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))
      if (det(l).le.zero) then
         write(iout,1000) nel,neg
         stop
      endif
c
c     inversion
c
      xsj(1,1) = (xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))/det(l)
      xsj(1,2) = (xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3))/det(l)
      xsj(1,3) = (xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))/det(l)
      xsj(2,1) = (xs(1,3)*xs(3,2)-xs(1,2)*xs(3,3))/det(l)
      xsj(2,2) = (xs(1,1)*xs(3,3)-xs(1,3)*xs(3,1))/det(l)
      xsj(2,3) = (xs(1,2)*xs(3,1)-xs(1,1)*xs(3,2))/det(l)
      xsj(3,1) = (xs(1,2)*xs(2,3)-xs(1,3)*xs(2,2))/det(l)
      xsj(3,2) = (xs(1,3)*xs(2,1)-xs(1,1)*xs(2,3))/det(l)
      xsj(3,3) = (xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1))/det(l)
c
      do 600 i=1,8
      temp1=xsj(1,1)*shg(1,i,l)+xsj(2,1)*shg(2,i,l)+xsj(3,1)*shg(3,i,l)
      temp2=xsj(1,2)*shg(1,i,l)+xsj(2,2)*shg(2,i,l)+xsj(3,2)*shg(3,i,l)
      shg(3,i,l) =
     &      xsj(1,3)*shg(1,i,l)+xsj(2,3)*shg(2,i,l)+xsj(3,3)*shg(3,i,l)
      shg(1,i,l) = temp1
      shg(2,i,l) = temp2
  600 continue
c
  700 continue
c
      return
c
 1000 format('1','non-positive determinant in element number  ',i5,
     &          ' in element group  ',i5)
      end
c**** new **************************************************************
      subroutine bkcshl(shl,w,nint)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a eight-node brick element
c
c             r,s,t = local element coord ("xi", "eta", "zeta" resp.)
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local ("zeta") derivative of shape function
c        shl(4,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 8
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shl(4,8,1),w(1),ra(8),sa(8),ta(8)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      data ra/-0.50, 0.50, 0.50,-0.50,-0.50, 0.50, 0.50,-0.50/
      data sa/-0.50,-0.50, 0.50, 0.50,-0.50,-0.50, 0.50, 0.50/
      data ta/-0.50,-0.50,-0.50,-0.50, 0.50, 0.50, 0.50, 0.50/
c
      g = zero
      w(1) = four + four
      if (nint.eq.8) then
         g = two/sqrt(three)
         w(1) = one
         w(2) = one
         w(3) = one
         w(4) = one
         w(5) = one
         w(6) = one
         w(7) = one
         w(8) = one
      endif
c
      do 200 l=1,nint
      r = g*ra(l)
      s = g*sa(l)
      t = g*ta(l)
c
      do 100 i=1,8
      tempr = pt5 + ra(i)*r
      temps = pt5 + sa(i)*s
      tempt = pt5 + ta(i)*t
      shl(1,i,l) = ra(i)*temps*tempt
      shl(2,i,l) = tempr*sa(i)*tempt
      shl(3,i,l) = tempr*temps*ta(i)
      shl(4,i,l) = tempr*temps*tempt
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine bkcstr(shg,shgbar,b,dl,strain,c,stress,pstrn,pstrs,
     &                  nrowsh,nesd,nrowb,ibbar,nen,ned,nee,nstr)
c
c.... program to calculate stress, strain and principal values at an
c        integration point for a three-dimensional continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shg(nrowsh,1),shgbar(nrowsh,1),b(nrowb,1),dl(ned,1),
     &          strain(1),c(nrowb,1),stress(1),pstrn(1),pstrs(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
c.... set up strain-displacement matrix
c
      call bkcb(shg,shgbar,b,nrowsh,nrowb,nen,ibbar)
c
c.... calculate strains
c
      call multab(b,dl,strain,nrowb,nee,nstr,nee,nstr,1,1)
c
c.... calculate stresses
c
      call multab(c,strain,stress,nrowb,nstr,nstr,nstr,nstr,1,1)
c
c.... calculate principal strains; account for engineering shear strain
c
      strain(4) = pt5*strain(4)
      strain(5) = pt5*strain(5)
      strain(6) = pt5*strain(6)
      call princ(nesd,strain,pstrn)
      strain(4) = two*strain(4)
      strain(5) = two*strain(5)
      strain(6) = two*strain(6)
c
c.... calculate principal stress
c
      call princ(nesd,stress,pstrs)
c
      return
      end
c**** new **************************************************************
      subroutine bkcsuf(ielno,ien,x,xl,iside,press,elresf,
     &                  brhs,lm,fac,nsurf,nen,nsd,nesd,ned,nee)
c
c.... program to compute consistent surface loads for the
c        eight-node brick, elastic continuum element
c
c        note: two by two points gaussian quadrature is employed
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension work(2),ielno(1),ien(nen,1),x(nsd,1),xl(nesd,1),
     &          iside(1),mat(1),press(4,1),shl(3,4,4),w(4),det(4),
     &          elresf(ned,1),brhs(1),lm(ned,nen,1),ienl(4,6),sigp(6),
     &          p(3),xlf(3,4),xs(3,3)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      data nint/4/
      data ienl /1,2,3,4 ,
     &           5,6,7,8 ,
     &           1,2,6,5 ,
     &           2,3,7,6 ,
     &           3,4,8,7 ,
     &           4,1,5,8 /
      data sigp /1.,-1.,-1.,-1.,-1.,-1./
c
      do 300 k=1,nsurf
      nel = ielno(k)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call clear(elresf,nee)
c
      i = iside(k)
      call local(ienl(1,i),xl,xlf,4,nsd,nesd)
c
      call qdcshl (shl,w,nint)
c
      do 200 l=1,nint
c
      call sufdet (xlf,shl(1,1,l),nel,neg,i
     &            ,xs,det(l))
      pn = rcdot(shl(3,1,l),press(1,k),3,4)
      pn = pn*fac*sigp(i)*w(l)*det(l)
      p(1)=xs(3,1)*pn
      p(2)=xs(3,2)*pn
      p(3)=xs(3,3)*pn
c
      do 100 n=1,3
      do 100 j=1,4
      jj=ienl(j,i)
      elresf(n,jj) = elresf(n,jj) + shl(3,j,l)*p(n)
  100 continue
c
  200 continue
      call addrhs(brhs,elresf,lm(1,1,nel),nee)
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine prop3d(rho,rdampm,rdampk,c,numat,nrowb)
c
c.... program to read, write and store properties for three-dimensional
c        continuum elements
c
c        note: this routine is presently restricted to the
c              isotropic linearly-elastic case
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension rho(1),rdampm(1),rdampk(1),c(nrowb,nrowb,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,numat
      if (mod(n,50).eq.1) write(iout,1000) numat
      read(iin,2000) m,e,pois,rho(m),rdampm(m),rdampk(m)
      write(iout,3000) m,e,pois,rho(m),rdampm(m),rdampk(m)
c
c.... set material constants for out-of-plane components
c
      amu2 = e/(one + pois)
      alam = e*pois/((one+pois)*(one-two*pois))
c
      call clear (c(1,1,m),nrowb*nrowb)
      c(1,1,m) = alam + amu2
      c(2,2,m) = alam + amu2
      c(3,3,m) = alam + amu2
      c(2,1,m) = alam
      c(3,1,m) = alam
      c(3,2,m) = alam
      c(1,2,m) = alam
      c(1,3,m) = alam
      c(2,3,m) = alam
      c(4,4,m) = pt5*amu2
      c(5,5,m) = pt5*amu2
      c(6,6,m) = pt5*amu2
c
  100 continue
c
      return
c
 1000 format('1',
     &' m a t e r i a l   s e t   d a t a                      '   //5x,
     &' number of material sets . . . . . . . . . . (numat ) = ',i5///,
     & 7x,'set',5x,'young''s',4x,'poisson''s',5x,'mass',8x,'mass',
     & 6x,'stiffness',3x,/6x,'number',3x,'modulus',
     & 6x,'ratio',6x,'density',5x,'damping',5x,'damping',/)
 2000 format(i5,5x,7f10.0)
 3000 format(4x,i5,3x,6(2x,1pe10.4))
      end
c**** new **************************************************************
      subroutine prts3d(xint,stress,pstrs,strain,pstrn,
     &                  nn,nntot,neg,nel,lint)
c
c.... program to print stress, strain, and principal values
c        for three-dimensional continuum elements
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension xint(3),stress(6),pstrs(3),strain(6),pstrn(3)
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
     &'  element  int. pt.         x1        x2        x3    ',5x,
     &'  stress    stress    stress    stress  ',
     &'  stress    stress  ',/,
     &'   number   number                          ',15x,
     &'    11        22        33        12        23        31',/,59x,
     &' stress 1  stress 2   stress 3  ',//,59x,
     &'  strain    strain    strain    strain  ',
     &'  strain    strain  ',/,59x,
     &'    11        22        33        12        23        31',/,59x,
     &' strain 1  strain 2   strain 3  ')
 2000 format(/2x,i5,6x,i2,8x,3(1pe10.2),5x,6(1pe10.2)/58x,3(1pe10.2)
     &,/58x,6(1pe10.2),/58x,3(1pe10.2))
      end
c**** new **************************************************************
      subroutine sufdet(xl,shl,nel,neg,if
     &                 ,xs,det)
c
c.... program to calculate jacobian and his determinant
c     for a four-node quadrilateral surface
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension xl(3,1),shl(3,1),xs(3,3)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
c
      do 300 j=1,3
      do 200 i=1,2
      xs(i,j) = rowdot(shl(i,1),xl(j,1),3,3,4)
  200 continue
  300 continue
c
c     third direction
c
      xs(3,1)= xs(1,2)*xs(2,3)-xs(2,2)*xs(1,3)
      xs(3,2)= xs(2,1)*xs(1,3)-xs(1,1)*xs(2,3)
      xs(3,3)= xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
      zn     = sqrt(xs(3,1)**2+xs(3,2)**2+xs(3,3)**2)
      do 400 i=1,3
      xs(3,i)= xs(3,i)/zn
  400 continue
c
      det    = xs(1,1)*(xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))
     &        +xs(1,2)*(xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3))
     &        +xs(1,3)*(xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))
      if (det.le.zero) then
         write(iout,1000) if,nel,neg
         stop
      endif
c
      return
c
 1000 format('1','non-positive determinant in side  ',i5,
     &           ' in element number  ',i5,
     &           ' in element group  ',i5)
      end
