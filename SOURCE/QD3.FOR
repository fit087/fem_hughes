c**** new **************************************************************
      subroutine qdct1(shl   ,w     ,rho   ,rdampm,rdampk,th    ,
     &                 c     ,grav  ,ien   ,mat   ,id    ,lm    ,
     &                 idiag ,ielno ,iside ,press ,shear ,ishist,
     &                 ntype ,numel ,numat ,nsurf ,nsout ,iopt  ,
     &                 istprt,lfsurf,lfbody,nicode,nint  ,ibbar ,
     &                 imass ,impexp,nrowsh,nrowb ,nesd  ,nen   ,
     &                 ndof  ,ned   ,iprtin,ldyn  )
c
c.... program to read, generate and write data for the
c        four-node quadrilateral, elastic continuum element
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      dimension shl(nrowsh,nen,1),w(1),rho(1),rdampm(1),rdampk(1),
     &          th(1),c(nrowb,nrowb,1),grav(nesd),ien(nen,1),mat(1),
     &          id(ndof,1),lm(ned,nen,1),idiag(1),ielno(1),iside(1),
     &          press(2,1),shear(2,1),ishist(3,1)
      common /iounit/ iin,iout,irsin,irsout
c
      write(iout,1000) ntype,numel,numat,nsurf,nsout,iopt,istprt,
     &                 lfsurf,lfbody
      write(iout,2000) nicode,ibbar
      if (ldyn) write(iout,3000) imass,impexp
c
      call qdcshl(shl,w,nint)
c
      call prop2d(rho,rdampm,rdampk,th,c,numat,iopt,nrowb)
c
      read (iin,4000) grav
      write (iout,5000) grav
c
      call genel(ien,mat,nen)
c
      if (iprtin.eq.0) call prntel(mat,ien,nen,numel)
c
      call formlm(id,ien,lm,ndof,ned,nen,numel)
c
      if ( (.not.ldyn) .or. (impexp.eq.0) .or. (imass.eq.0) )
     &   call colht(idiag,lm,ned,nen,numel)
c
      if (nsurf.gt.0) call qdcrsf(ielno,iside,press,shear,nsurf)
c
      if (nsout.gt.0) call shist(ishist,nsout,ntype)
c
      return
c
 1000 format(' ',' f o u r - n o d e  ',
     &' q u a d r i l a t e r a l   e l e m e n t s',              //5x,
     &' element type number . . . . . . . . . . . . (ntype ) = ',i5//5x,
     &' number of elements  . . . . . . . . . . . . (numel ) = ',i5//5x,
     &' number of element material sets . . . . . . (numat ) = ',i5//5x,
     &' number of surface force cards . . . . . . . (nsurf ) = ',i5//5x,
     &' number of stress/strain time histories  . . (nsout ) = ',i5//5x,
     &' analysis option . . . . . . . . . . . . . . (iopt  ) = ',i5//5x,
     &'    eq.0, plane stress                                  ',   /5x,
     &'    eq.1, plane strain                                  ',   /5x,
     &'    eq.2, axisymmetric                                  ',  //5x,
     &' stress output print code  . . . . . . . . . (istprt) = ',i5//5x,
     &'    eq.0, stress output printed                         ',   /5x,
     &'    eq.1, stress output not printed                     ',  //5x,
     &' surface force load-time function number . . (lfsurf) = ',i5//5x,
     &' body force load-time function number  . . . (lfbody) = ',i5 /5x)
 2000 format(5x,
     &' numerical integration code  . . . . . . . . (nicode) = ',i5//5x,
     &'    eq.0, 2 x 2 gaussian quadrature                     ',   /5x,
     &'    eq.1, 1-point gaussian quadrature                   ',  //5x,
     &' strain-displacement option  . . . . . . . . (ibbar ) = ',i5//5x,
     &'    eq.0, standard formulation                          ',   /5x,
     &'    eq.1, b-bar formulation                             ',   /5x)
 3000 format(5x,
     &' mass type code  . . . . . . . . . . . . . . (imass ) = ',i5//5x,
     &'    eq.0, consistent mass matrix                        ',   /5x,
     &'    eq.1, lumped mass matrix                            ',   /5x,
     &'    eq.2, no mass matrix                                ',  //5x,
     &' implicit/explicit code  . . . . . . . . . . (impexp) = ',i5//5x,
     &'    eq.0, implicit element group                        ',   /5x,
     &'    eq.1, explicit element group                        ',  //5x)
 4000 format(8f10.0)
 5000 format(////' ',
     &' g r a v i t y   v e c t o r   c o m p o n e n t s      ',  //5x,
     &' x-1 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x,
     &' x-2 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x)
c
      end
c**** new **************************************************************
      subroutine qdct2(eleffm,ien   ,x     ,xl    ,mat   ,det   ,
     &                 shl   ,shg   ,r     ,rdampm,rdampk,th    ,
     &                 rho   ,w     ,work  ,shgbar,b     ,c     ,
     &                 dmat  ,db    ,alhs  ,idiag ,lm    ,
     &                 impexp,imass ,numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,ldyn  ,ned   ,
     &                 iopt  ,ibbar ,nrowb ,nstr  ,nee   )
c
c.... program to calculate effective mass matrix for the
c        four-node quadrilateral, elastic continuum element and
c        assemble into the global left-hand-side matrix
c
c         impexp = 0, implicit time integration
c                = 1, explicit time integration
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,ldiag,lquad
      dimension eleffm(nee,1),ien(nen,1),x(nsd,1),xl(nesd,1),mat(1),
     &          det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),r(1),
     &          rdampm(1),rdampk(1),th(1),rho(1),w(1),work(1),
     &          shgbar(3,1),b(nrowb,1),c(nrowb,nrowb,1),dmat(nrowb,1),
     &          db(nrowb,1),alhs(1),idiag(1),lm(ned,nen,1)
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      ldiag = .false.
      if ( (impexp.eq.1) .and. (imass.eq.1) ) ldiag = .true.
c
      do 200 nel=1,numel
c
      call clear(eleffm,neesq)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      m = mat(nel)
      lquad = .true.
      if (ien(3,nel).eq.ien(4,nel)) lquad = .false.
      call qdcshg(xl,det,shl,shg,nint,nel,neg,lquad)
c
      if (iopt.eq.2) then
c
         do 100 l=1,nint
         r(l) = rowdot(shg(nrowsh,1,l),xl,nrowsh,nesd,nen)
         det(l) = det(l)*r(l)
  100    continue
c
      endif
c
      if ( ldyn .and. (imass.ne.2) ) then
c
c....... form mass matrix
c
         constm = (one + rdampm(m)*coeff4)*th(m)*rho(m)
         if (constm.ne.zero) call contm(shg,xl,w,det,eleffm,work,
     &             constm,imass,nint,nrowsh,nesd,nen,ned,nee,.false.)
c
      endif
c
      if ( (.not.ldyn) .or. (impexp.eq.0) ) then
c
c....... form stiffness matrix
c
         constk = (coeff4*rdampk(m) + coeff5)*th(m)
         call qdck(shgbar,w,det,r,shg,b,c(1,1,m),dmat,db,eleffm,constk,
     &            ibbar,nen,nint,iopt,nesd,nrowsh,nrowb,nstr,nee)
c
      endif
c
c.... assemble element effective mass matrix into global
c        left-hand-side matrix
c
      call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,ldiag)
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine qdct3(mat   ,ien   ,dpred ,dl    ,vpred ,vl    ,
     &                 a     ,al    ,rdampk,rdampm,rho   ,grav  ,
     &                 elresf,x     ,xl    ,det   ,shl   ,shg   ,
     &                 r     ,g1    ,work  ,th    ,w     ,eleffm,
     &                 shgbar,b     ,strain,c     ,stress,brhs  ,
     &                 lm    ,ielno ,iside ,press ,shear ,
     &                 numel ,ned   ,nen   ,ndof  ,ldyn  ,nee   ,
     &                 imass ,nesd  ,lfbody,nsd   ,nint  ,nrowsh,
     &                 neg   ,iopt  ,nrowb ,nstr  ,ibbar , nsurf,
     &                 lfsurf)
c
c.... program to calculate residual-force vector for the
c        four-node quadrilateral, elastic continuum element and
c        assemble into the global right-hand-side vector
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,formma,formkd,zeroal,zerodl,zerog,lquad
      dimension mat(1),ien(nen,1),dpred(ndof,1),dl(ned,1),vpred(ndof,1),
     &          vl(ned,1),a(ndof,1),al(ned,1),rdampk(1),rdampm(1),
     &          rho(1),grav(1),elresf(1),x(nsd,1),xl(nesd,1),det(1),
     &          shl(nrowsh,nen,1),shg(nrowsh,nen,1),r(1),g1(1),work(1),
     &          th(1),w(1),eleffm(nee,1),shgbar(3,1),b(nrowb,1),
     &          strain(1),c(nrowb,nrowb,1),stress(1),brhs(1),
     &          lm(ned,nen,1),ielno(1),iside(1),press(2,1),shear(2,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 600 nel=1,numel
c
      formma = .false.
      formkd = .false.
      m = mat(nel)
c
c.... note: for static analysis mpdprd = mpd, hence reference to
c           array "dpred" will access the contents of array "d".
c
      call local(ien(1,nel),dpred,dl,nen,ndof,ned)
      if (ldyn) then
c
         call local(ien(1,nel),vpred,vl,nen,ndof,ned)
         call local(ien(1,nel),a,al,nen,ndof,ned)
c
         do 200 j=1,nen
c
         do 100 i=1,ned
         dl(i,j) = dl(i,j) + rdampk(m)*vl(i,j)
         al(i,j) = al(i,j) + rdampm(m)*vl(i,j)
  100    continue
c
  200    continue
c
         call ztest(al,nee,zeroal)
         if ( (.not.zeroal) .and. (imass.ne.2) .and. (rho(m).ne.zero) )
     &      formma = .true.
c
      else
c
         call clear(al,nee)
c
      endif
c
      call ztest(dl,nee,zerodl)
      if (.not.zerodl) formkd = .true.
      call ztest(grav,nesd,zerog)
c
      if ((.not.zerog) .and. (lfbody.ne.0) .and. (rho(m).ne.zero)
     &   .and. (imass.ne.2)) then
         formma = .true.
         do 400 i=1,ned
         temp = grav(i)*g1(lfbody)
c
         do 300 j=1,nen
         al(i,j) = al(i,j) - temp
  300    continue
c
  400    continue
c
      endif
c
      if (formma.or.formkd) then
c
         call clear(elresf,nee)
         call local(ien(1,nel),x,xl,nen,nsd,nesd)
         lquad = .true.
         if (ien(3,nel).eq.ien(4,nel)) lquad = .false.
         call qdcshg(xl,det,shl,shg,nint,nel,neg,lquad)
c
         if (iopt.eq.2) then
            do 500 l=1,nint
            r(l) = rowdot(shg(nrowsh,1,l),xl,nrowsh,nesd,nen)
            det(l) = det(l)*r(l)
  500       continue
         endif
c
         if (formma) then
c
c.......... form inertial and/or body force
c
            constm = - th(m)*rho(m)
            call contma(shg,xl,w,det,al,eleffm,work,elresf,constm,imass,
     &                  nint,nrowsh,nesd,nen,ned,nee)
         endif
c
         if (formkd) then
c
c.......... form internal force
c
            constk = - th(m)
            call qdckd(shgbar,w,det,r,shg,b,dl,strain,c(1,1,m),stress,
     &                 work,elresf,constk,ibbar,nen,nint,iopt,nrowsh,
     &                 nesd,nrowb,nee,nstr)
         endif
c
         call addrhs(brhs,elresf,lm(1,1,nel),nee)
c
      endif
c
  600 continue
c
c.... form surface force
c
c        note: assembly of surface loads is performed inside qdcsuf
c
      if ( (nsurf.gt.0) .and. (lfsurf.gt.0) )
     &   call qdcsuf(ielno,ien,x,xl,iside,mat,th,press,shear,elresf,
     &               brhs,lm,g1(lfsurf),nsurf,nen,nsd,nesd,ned,nee,iopt)
c
      return
      end
c**** new **************************************************************
      subroutine qdct4(mat   ,ien   ,d     ,dl    ,x     ,xl    ,
     &                 det   ,shl   ,shg   ,xint  ,r     ,shgbar,
     &                 w     ,b     ,strain,c     ,stress,pstrn ,
     &                 pstrs ,
     &                 nint  ,numel ,nen   ,ndof  ,ned   ,nsd   ,
     &                 nesd  ,nrowsh,neg   ,iopt  ,ibbar ,nrowb ,
     &                 nee   ,nstr  )
c
c.... program to calculate and print stress, strain and
c        principal values for the four-node quadrilateral,
c        elastic continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical lquad
      dimension mat(1),ien(nen,1),d(ndof,1),dl(ned,1),x(nsd,1),
     &          xl(nesd,1),det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),
     &          xint(nesd,1),r(1),shgbar(3,1),w(1),b(nrowb,1),strain(1),
     &          c(nrowb,nrowb,1),stress(1),pstrn(1),pstrs(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      nntot = 16
      nn = 0
c
      do 300 nel=1,numel
c
      m = mat(nel)
      call local(ien(1,nel),d,dl,nen,ndof,ned)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      lquad = .true.
      if (ien(3,nel).eq.ien(4,nel)) lquad = .false.
      call qdcshg(xl,det,shl,shg,nint,nel,neg,lquad)
c
c.... calculate coordinates of integration points
c
      do 100 l=1,nint
      xint(1,l) = rowdot(shg(nrowsh,1,l),xl(1,1),nrowsh,nesd,nen)
      xint(2,l) = rowdot(shg(nrowsh,1,l),xl(2,1),nrowsh,nesd,nen)
      if (iopt.eq.2) then
         r(l) = xint(1,l)
         det(l) = det(l)*r(l)
      endif
  100 continue
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,iopt,nesd,nrowsh)
c
c.... loop over integration points
c
      do 200 l=1,nint
c
c.... calculate stress, strain and principal values
c
      call qdcstr(shg(1,1,l),shgbar,b,r(l),dl,strain,c(1,1,m),stress,
     &      pstrn,pstrs,nrowsh,nesd,nrowb,ibbar,nen,ned,nee,nstr,iopt)
c
c.... print stress, strain and principal values
c
      call prts2d(xint(1,l),stress,pstrs,strain,pstrn,
     &            nn,nntot,neg,nel,l)
  200 continue
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine qdct5(ishist,sout  ,mat   ,ien   ,d     ,dl    ,
     &                 x     ,xl    ,det   ,shl   ,shg   ,r     ,
     &                 shgbar,w     ,b     ,strain,c     ,stress,
     &                 pstrn ,pstrs ,work  ,
     &                 nsout ,nen   ,ndof  ,ned   ,nsd   ,nesd  ,
     &                 nrowsh,nint  ,neg   ,iopt  ,ibbar ,nrowb ,
     &                 nee   ,nstr  ,nsout1)
c
c.... program to calculate and store element time-histories for the
c        four-node quadrilateral, elastic continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      real sout
      logical lquad
      dimension ishist(3,1),sout(nsout1,1),mat(1),ien(nen,1),d(ndof,1),
     &          dl(ned,1),x(nsd,1),xl(nesd,1),det(1),shl(nrowsh,nen,1),
     &          shg(nrowsh,nen,1),r(1),shgbar(3,nen,1),w(1),b(nrowb,1),
     &          strain(1),c(nrowb,nrowb,1),stress(1),pstrn(1),pstrs(1),
     &          work(1)
      common /hplotc/ nplpts,locplt,time
c
      sout(1,locplt) = real(time)
c
      do 300 i=1,nsout
c
      nel   = ishist(1,i)
      intpt = ishist(2,i)
      ncomp = ishist(3,i)
c
      m = mat(nel)
      call local(ien(1,nel),d,dl,nen,ndof,ned)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      lquad = .true.
      if (ien(3,nel).eq.ien(4,nel)) lquad = .false.
      call qdcshg(xl,det,shl,shg,nint,nel,neg,lquad)
c
      if (iopt.eq.2) then
c
         do 100 l=1,nint
         r(l) = rowdot(shg(nrowsh,1,l),xl,nrowsh,nesd,nen)
         det(l) = det(l)*r(l)
  100    continue
c
      endif
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,iopt,nesd,nrowsh)
c
c.... calculate stress, strain and principal values
c
      call qdcstr(shg(1,1,intpt),shgbar,b,r(intpt),dl,strain,c(1,1,m),
     &            stress,pstrn,pstrs,nrowsh,nesd,nrowb,ibbar,nen,ned,
     &            nee,nstr,iopt)
c
      do 200 j=1,4
      work(j     ) = stress(j)
      work(j +  4) = pstrs(j)
      work(j +  8) = strain(j)
      work(j + 12) = pstrn(j)
  200 continue
c
      sout(i+1,locplt) = real(work(ncomp))
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine qdcb(shg,shgbar,b,r,iopt,nrowsh,nrowb,nen,ibbar)
c
c.... program to set up the strain-displacement matrix "b" for
c        two-dimensional continuum elements
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
      dimension shg(nrowsh,1),shgbar(3,1),b(nrowb,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 100 j=1,nen
c
      j2   = 2*j
      j2m1 = j2 - 1
c
      b(1,j2m1) = shg(1,j)
      b(1,j2  ) = zero
      b(2,j2m1) = zero
      b(2,j2  ) = shg(2,j)
      b(3,j2m1) = shg(2,j)
      b(3,j2  ) = shg(1,j)
c
      if (iopt.eq.2) then
         b(4,j2m1) = shg(3,j)/r
         b(4,j2  ) = zero
      endif
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
      j2   = 2*j
      j2m1 = j2 - 1
c
      if (iopt.eq.2) then
         temp3 = constb*(shgbar(3,j) - shg(3,j)/r)
         b(1,j2m1) = b(1,j2m1) + temp3
         b(2,j2m1) = b(2,j2m1) + temp3
         b(4,j2m1) = b(4,j2m1) + temp3
      else
         b(4,j2m1) = zero
         b(4,j2  ) = zero
      endif
c
      temp1 = constb*(shgbar(1,j) - shg(1,j))
      temp2 = constb*(shgbar(2,j) - shg(2,j))
c
      b(1,j2m1) = b(1,j2m1) + temp1
      b(1,j2  ) = b(1,j2  ) + temp2
      b(2,j2m1) = b(2,j2m1) + temp1
      b(2,j2  ) = b(2,j2  ) + temp2
      b(4,j2m1) = b(4,j2m1) + temp1
      b(4,j2  ) = b(4,j2  ) + temp2
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine qdck(shgbar,w,det,r,shg,b,c,dmat,db,elstif,constk,
     &                ibbar,nen,nint,iopt,nesd,nrowsh,nrowb,nstr,nee)
c
c.... program to form stiffness matrix for a continuum element
c        with "nen" nodes
c
c        note: the b-bar option is restricted to the mean-dilatation
c              formulation. to generalize to other formulations,
c              redimension array "shgbar", and replace routines
c              "meansh" and "qdcb".
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shgbar(3,nen,1),w(1),det(1),r(1),shg(nrowsh,nen,1),
     &          b(nrowb,1),c(nrowb,1),dmat(nrowb,1),db(nrowb,1),
     &          elstif(nee,1)
c
c.... calculate mean values of shape function global derivatives
c        for mean-dilatational b-bar formulation
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,iopt,nesd,nrowsh)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call qdcb(shg(1,1,l),shgbar,b,r(l),iopt,nrowsh,nrowb,nen,ibbar)
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
      subroutine qdckd(shgbar,w,det,r,shg,b,dl,strain,c,stress,work,
     &                 elresf,constk,ibbar,nen,nint,iopt,nrowsh,
     &                 nesd,nrowb,nee,nstr)
c
c.... program to form internal force ("-k*d") for a continuum element
c        with "nen" nodes
c
c        note: the b-bar option is restricted to the mean-dilatation
c              formulation. to generalize to other formulations,
c              redimension array "shgbar", and replace routines
c              "meansh" and "qdcb".
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shgbar(3,nen,1),w(1),det(1),r(1),shg(nrowsh,nen,1),
     &          b(nrowb,1),dl(1),strain(1),c(nrowb,1),stress(1),
     &          work(1),elresf(1)
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,iopt,nesd,nrowsh)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call qdcb(shg(1,1,l),shgbar,b,r(l),iopt,nrowsh,nrowb,nen,ibbar)
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
      subroutine qdcrsf(ielno,iside,press,shear,nsurf)
c
c.... program to read, write and store surface force data for the
c        four-node quadrilateral, elastic continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension ielno(1),iside(1),press(2,1),shear(2,1)
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,nsurf
      if (mod(n,50).eq.1) write(iout,1000) nsurf
      read(iin,2000) ielno(n),iside(n),press(1,n),press(2,n),
     &               shear(1,n),shear(2,n)
      write(iout,3000) ielno(n),iside(n),press(1,n),press(2,n),
     &                shear(1,n),shear(2,n)
  100 continue
c
      return
c
 1000 format('1',
     &' e l e m e n t   s u r f a c e   f o r c e   d a t a    ',  //5x,
     &' number of surface force cards . . . . . . . (nsurf ) = ',i5///
     &5x,' element    side    ',2('   pressure   '),
     &       2('    shear     '),/
     &5x,2('  number  '),2('    node i        node j    '),/)
 2000 format(2i5,4f10.0)
 3000 format(6x,i5,7x,i2,3x,4(2x,e12.4))
      end
c**** new **************************************************************
      subroutine qdcshg(xl,det,shl,shg,nint,nel,neg,lquad)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a four-node quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
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
      logical lquad
      dimension xl(2,1),det(1),shl(3,4,1),shg(3,4,1),xs(2,2)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      call move(shg,shl,12*nint)
c
      do 700 l=1,nint
c
      if (.not.lquad) then
         do 100 i=1,3
         shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shg(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shg(i,1,l),xl(j,1),3,2,4)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det(l).le.zero) then
         write(iout,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det(l)
  400 continue
  500 continue
c
      do 600 i=1,4
        temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
        shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
        shg(1,i,l) = temp
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
      subroutine qdcshl(shl,w,nint)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shl(3,4,1),w(1),ra(4),sa(4)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      data ra/-0.50,0.50,0.50,-0.50/,sa/-0.50,-0.50,0.50,0.50/
c
      g = zero
      w(1) = four
      if (nint.eq.4) then
         g = two/sqrt(three)
         w(1) = one
         w(2) = one
         w(3) = one
         w(4) = one
      endif
c
      do 200 l=1,nint
      r = g*ra(l)
      s = g*sa(l)
c
      do 100 i=1,4
      tempr = pt5 + ra(i)*r
      temps = pt5 + sa(i)*s
      shl(1,i,l) = ra(i)*temps
      shl(2,i,l) = tempr*sa(i)
      shl(3,i,l) = tempr*temps
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine qdcstr(shg,shgbar,b,r,dl,strain,c,stress,pstrn,pstrs,
     &                  nrowsh,nesd,nrowb,ibbar,nen,ned,nee,nstr,iopt)
c
c.... program to calculate stress, strain and principal values at an
c        integration point for a two-dimensional continuum element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shg(nrowsh,1),shgbar(3,1),b(nrowb,1),dl(ned,1),
     &          strain(1),c(nrowb,1),stress(1),pstrn(1),pstrs(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
c.... set up strain-displacement matrix
c
      call qdcb(shg,shgbar,b,r,iopt,nrowsh,nrowb,nen,ibbar)
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
      strain(3) = pt5*strain(3)
      call princ(nesd,strain,pstrn)
      strain(3) = two*strain(3)
      pstrn(3)  = two*pstrn(3)
c
c.... calculate principal stress
c
      call princ(nesd,stress,pstrs)
c
      if (iopt.eq.0) then
         stress(4) = zero
         strain(4) = - ( c(4,1)*strain(1) + c(4,2)*strain(2)
     &                 + c(4,3)*strain(3) )/c(4,4)
      endif
c
      if ( (iopt.eq.1) .and. (ibbar.eq.0) ) then
         strain(4) = zero
         stress(4) = c(4,1)*strain(1) + c(4,2)*strain(2)
     &             + c(4,3)*strain(3)
      endif
c
      return
      end
c**** new **************************************************************
      subroutine qdcsuf(ielno,ien,x,xl,iside,mat,th,press,shear,elresf,
     &                  brhs,lm,fac,nsurf,nen,nsd,nesd,ned,nee,iopt)
c
c.... program to compute consistent surface loads for the
c        four-node quadrilateral, elastic continuum element
c
c        note: two-point gaussian quadrature is employed
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension z(2),work(2),ielno(1),ien(nen,1),x(nsd,1),xl(nesd,1),
     &          iside(1),mat(1),th(1),press(2,1),shear(2,1),
     &          elresf(ned,1),brhs(1),lm(ned,nen,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      z(2) = pt5/sqrt(three)
      z(1) = - z(2)
c
      do 300 k=1,nsurf
      nel = ielno(k)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call clear(elresf,nee)
      i = iside(k)
      j = i + 1
      if (j.eq.5) j = 1
      dx = xl(1,j) - xl(1,i)
      dy = xl(2,j) - xl(2,i)
      m = mat(nel)
      temp = pt5*fac*th(m)
c
      do 200 l=1,2
      shi = pt5 - z(l)
      shj = pt5 + z(l)
      p = shi*press(1,k) + shj*press(2,k)
      s = shi*shear(1,k) + shj*shear(2,k)
c
      if (iopt.eq.2) then
         r = shi*xl(1,i) + shj*xl(1,j)
         p = p*r
         s = s*r
      endif
c
      work(1) = temp*( - p*dy + s*dx)
      work(2) = temp*(   p*dx + s*dy)
c
      do 100 n=1,2
      elresf(n,i) = elresf(n,i) + shi*work(n)
      elresf(n,j) = elresf(n,j) + shj*work(n)
  100 continue
c
  200 continue
c
      call addrhs(brhs,elresf,lm(1,1,nel),nee)
c
  300 continue
c
      return
      end
