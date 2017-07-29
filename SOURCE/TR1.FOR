c**** new **************************************************************
      subroutine trust1(shl   ,w     ,rho   ,rdampm,rdampk,area  ,
     &                  c     ,grav  ,ien   ,mat   ,id    ,lm    ,
     &                  idiag ,ishist,
     &                  ntype ,numel ,numat ,nen   ,nsout ,istprt,
     &                  lfbody,nint  ,imass ,impexp,nrowsh,nrowb ,
     &                  nesd  ,ndof  ,ned   ,iprtin,ldyn  )
c
c.... program to read, generate and write element data for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical ldyn
      dimension shl(nrowsh,nen,1),w(1),rho(1),rdampm(1),rdampk(1),
     &          area(1),c(nrowb,nrowb,1),grav(nesd),ien(nen,1),mat(1),
     &          id(ndof,1),lm(ned,nen,1),idiag(1),ishist(3,1)
      common /iounit/ iin,iout,irsin,irsout
c
      write(iout,1000) ntype,numel,numat,nen,nsout,istprt,lfbody,nint
      if (ldyn) write(iout,2000) imass,impexp
c
      call trushl(shl,w,nint,nen)
c
      call truspr(rho,rdampm,rdampk,area,c,numat)
c
      read(iin,3000) grav
      write(iout,4000) grav
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
      if (nsout.gt.0) call shist(ishist,nsout,ntype)
c
      return
c
 1000 format(' ',
     &' t w o / t h r e e - n o d e   t r u s s   e l e m e n t s',//5x,
     &' element type number . . . . . . . . . . . . (ntype ) = ',i5//5x,
     &' number of elements  . . . . . . . . . . . . (numel ) = ',i5//5x,
     &' number of element material sets . . . . . . (numat ) = ',i5//5x,
     &' number of element nodes . . . . . . . . . . (nen   ) = ',i5//5x,
     &' number of stress/strain time histories  . . (nsout ) = ',i5//5x,
     &' stress output print code  . . . . . . . . . (istprt) = ',i5//5x,
     &'    eq.0, stress output printed                         ',   /5x,
     &'    eq.1, stress output not printed                     ',  //5x,
     &' body force load-time function number  . . . (lfbody) = ',i5//5x,
     &' integration code  . . . . . . . . . . . . . (nint  ) = ',i5//5x,
     &'    eq.1, 1-point gaussian quadrature                   ',   /5x,
     &'    eq.2, 2-point gaussian quadrature                   ',   /5x,
     &'    eq.3, 3-point gaussian quadrature                   ',   /5x)
 2000 format(' ',/5x,
     &' mass type code  . . . . . . . . . . . . . . (imass ) = ',i5//5x,
     &'    eq.0, consistent mass matrix                        ',   /5x,
     &'    eq.1, lumped mass matrix                            ',   /5x,
     &'    eq.2, no mass matrix                                ',  //5x,
     &' implicit/explicit code  . . . . . . . . . . (impexp) = ',i5//5x,
     &'    eq.0, implicit element group                        ',   /5x,
     &'    eq.1, explicit element group                        ',  //5x)
 3000 format(8f10.0)
 4000 format(////' ',
     &' g r a v i t y   v e c t o r   c o m p o n e n t s      ',  //5x,
     &' x-1 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x,
     &' x-2 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x,
     &' x-3 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x)
c
      end
c**** new **************************************************************
      subroutine trust2(eleffm,ien   ,x     ,xl    ,mat   ,det   ,
     &                  shl   ,shg   ,rdampm,rdampk,area  ,rho   ,
     &                  w     ,work  ,b     ,c     ,dmat  ,db    ,
     &                  alhs  ,idiag ,lm    ,xs    ,
     &                  impexp,imass ,numel ,neesq ,nen   ,nsd   ,
     &                  nesd  ,nint  ,neg   ,nrowsh,ldyn  ,ned   ,
     &                  nrowb ,nstr  ,nee   )
c
c.... program to calculate effective mass matrix for the
c        three-dimensional, elastic truss element and
c        assemble into the global left-hand-side matrix
c
c         impexp = 0, implicit time integration
c                = 1, explicit time integration
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,ldiag,lnode3
      dimension eleffm(nee,1),ien(nen,1),x(nsd,1),xl(nesd,1),mat(1),
     &          det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),rdampm(1),
     &          rdampk(1),area(1),rho(1),w(1),work(1),b(nrowb,1),
     &          c(nrowb,nrowb,1),dmat(nrowb,1),db(nrowb,1),alhs(1),
     &          idiag(1),lm(ned,nen,1),xs(nesd,1)
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      ldiag = .false.
      if ( (impexp.eq.1) .and. (imass.eq.1) ) ldiag = .true.
c
      do 100 nel=1,numel
c
      call clear(eleffm,neesq)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      m = mat(nel)
      lnode3 = .true.
      if ( nen .eq. 3 .and. ien(2,nel).eq.ien(3,nel) ) lnode3 = .false.
      call trushg(xl,det,shl,shg,xs,nen,nint,nel,neg,lnode3)
c
      if ( ldyn .and. (imass.ne.2) ) then
c
c....... form mass matrix
c
         constm = (one + rdampm(m)*coeff4)*area(m)*rho(m)
         if (constm.ne.zero) call contm(shg,xl,w,det,eleffm,work,
     &             constm,imass,nint,nrowsh,nesd,nen,ned,nee,.false.)
c
      endif
c
      if ( (.not.ldyn) .or. (impexp.eq.0) ) then
c
c....... form stiffness matrix
c
         constk = (coeff4*rdampk(m) + coeff5)*area(m)
         call trusk(w,det,shg,xs,xl,b,c(1,1,m),dmat,db,eleffm,
     &            constk,nen,nint,nesd,nrowsh,nrowb,nstr,nee)
c
      endif
c
c.... assemble element effective mass matrix into global
c        left-hand-side matrix
c
      call addlhs(alhs,eleffm,idiag,lm(1,1,nel),nee,ldiag)
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine trust3(mat   ,ien   ,dpred ,dl    ,vpred ,vl    ,
     &                  a     ,al    ,rdampk,rdampm,rho   ,grav  ,
     &                  elresf,x     ,xl    ,det   ,shl   ,shg   ,
     &                  g1    ,work  ,area  ,w     ,eleffm,b     ,
     &                  strain,c     ,dmat  ,stress,brhs  ,lm    ,
     &                  xs    ,
     &                  numel ,ned   ,nen   ,ndof  ,ldyn  ,nee   ,
     &                  imass ,nesd  ,lfbody,nsd   ,nint  ,nrowsh,
     &                  neg   ,nrowb )
c
c.... program to calculate residual-force vector for the
c        three-dimensional, elastic truss element and
c        assemble into the global right-hand-side vector
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,formma,formkd,zeroal,zerodl,zerog,lnode3
      dimension mat(1),ien(nen,1),dpred(ndof,1),dl(ned,1),vpred(ndof,1),
     &          vl(ned,1),a(ndof,1),al(ned,1),rdampk(1),rdampm(1),
     &          rho(1),grav(1),elresf(1),x(nsd,1),xl(nesd,1),det(1),
     &          shl(nrowsh,nen,1),shg(nrowsh,nen,1),g1(1),work(1),
     &          area(1),w(1),eleffm(nee,1),b(nrowb,1),strain(1),
     &          c(nrowb,nrowb,1),dmat(nrowb,1),stress(1),brhs(1),
     &          lm(ned,nen,1),xs(nesd,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 500 nel=1,numel
c
      formma = .false.
      formkd = .false.
      m = mat(nel)
c
c.... note: for static analysis mpdprd = mpd, hence reference
c           to array "dpred" will access the contents of array "d".
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
         work(i) = grav(i)*g1(lfbody)
c
         do 300 j=1,nen
         al(i,j) = al(i,j) - work(i)
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
         lnode3 = .true.
         if (nen.eq.3 .and. ien(2,nel).eq.ien(3,nel)) lnode3 = .false.
         call trushg(xl,det,shl,shg,xs,nen,nint,nel,neg,lnode3)
c
         if (formma) then
c
c.......... form inertial and/or body force
c
            constm = - area(m)*rho(m)
            call contma(shg,xl,w,det,al,eleffm,work,elresf,constm,imass,
     &                  nint,nrowsh,nesd,nen,ned,nee)
         endif
c
         if (formkd) then
c
c.......... form internal force
c
            constk = - area(m)
            call truskd(w,det,shg,xs,xl,b,dl,strain,c(1,1,m),dmat,
     &                  stress,work,elresf,constk,nen,nint,nrowsh,
     &                  nesd,nrowb,nee)
         endif
c
         call addrhs(brhs,elresf,lm(1,1,nel),nee)
c
      endif
c
  500 continue
c
      return
      end
c**** new **************************************************************
      subroutine trust4(mat   ,ien   ,d     ,dl    ,x     ,xl    ,
     &                  det   ,shl   ,shg   ,xs    ,xint  ,b     ,
     &                  strain,c     ,stress,force ,area  ,
     &                  nint  ,numel ,nen   ,ndof  ,ned   ,nsd   ,
     &                  nesd  ,nrowsh,neg   ,nrowb ,nee   )
c
c.... program to calculate and print stress, strain and force for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical lnode3
      dimension mat(1),ien(nen,1),d(ndof,1),dl(ned,1),x(nsd,1),
     &          xl(nesd,1),det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),
     &          xs(nesd,1),xint(nesd,1),b(nrowb,1),strain(1),
     &          c(nrowb,nrowb,1),stress(1),force(1),area(1)
c
      nntot = 24
      nn = 0
c
      do 300 nel=1,numel
c
      m = mat(nel)
      call local(ien(1,nel),d,dl,nen,ndof,ned)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      lnode3 = .true.
      if ( nen .eq. 3 .and. ien(2,nel).eq.ien(3,nel) ) lnode3 = .false.
      call trushg(xl,det,shl,shg,xs,nen,nint,nel,neg,lnode3)
c
c.... loop over integration points
c
      do 200 l=1,nint
c
c.... calculate coordinates of integration points
c
      do 100 i=1,nesd
      xint(i,l) = rowdot(shg(nrowsh,1,l),xl(i,1),nrowsh,nesd,nen)
  100 continue
c
c.... calculate stress, strain and force
c
      call trustr(shg(1,1,l),xs(1,l),b,dl,strain,c(1,1,m),stress,
     &            force,area(m),nrowsh,nesd,nrowb,nen,nee)
c
c.... print stress, strain and force
c
      call truspt(xint(1,l),stress,force,strain,nn,nntot,neg,nel,l)
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine trust5(ishist,sout  ,mat   ,ien   ,d     ,dl    ,
     &                  x     ,xl    ,det   ,shl   ,shg   ,xs    ,
     &                  b     ,strain,c     ,stress,force ,area  ,
     &                  work  ,
     &                  nsout ,nen   ,ndof  ,ned   ,nsd   ,nesd  ,
     &                  nrowsh,nint  ,neg   ,nrowb ,nee   ,nsout1)
c
c.... program to calculate and store element time-histories for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      real sout
      logical lnode3
      dimension ishist(3,1),sout(nsout1,1),mat(1),ien(nen,1),d(ndof,1),
     &          dl(ned,1),x(nsd,1),xl(nesd,1),det(1),shl(nrowsh,nen,1),
     &          shg(nrowsh,nen,1),xs(nesd,1),b(nrowb,1),strain(1),
     &          c(nrowb,nrowb,1),stress(1),force(1),area(1),work(1)
      common /hplotc/ nplpts,locplt,time
c
      sout(1,locplt) = real(time)
c
      do 100 i=1,nsout
c
      nel   = ishist(1,i)
      intpt = ishist(2,i)
      ncomp = ishist(3,i)
c
      m = mat(nel)
      call local(ien(1,nel),d,dl,nen,ndof,ned)
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      lnode3 = .true.
      if ( nen .eq. 3 .and. ien(2,nel).eq.ien(3,nel) ) lnode3 = .false.
      call trushg(xl,det,shl,shg,xs,nen,nint,nel,neg,lnode3)
c
      call trustr(shg(1,1,intpt),xs(1,intpt),b,dl,strain,c(1,1,m),
     &            stress,force,area(m),nrowsh,nesd,nrowb,nen,nee)
c
      work(1) = stress(1)
      work(2) = force(1)
      work(3) = strain(1)
      sout(i+1,locplt) = real(work(ncomp))
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine trusb(b,shg,xs,nen,nesd,nrowb,nrowsh)
c
c.... program to set up the strain-displacement matrix for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension b(nrowb,1),shg(nrowsh,1),xs(1)
c
      do 200 j=1,nen
      k = (j - 1)*nesd
c
      do 100 i=1,nesd
      b(1,k+i) = shg(1,j)*xs(i)
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine trushg(xl,det,shl,shg,xs,nen,nint,nel,neg,lnode3)
c
c.... program to calculate global derivatives of shape functions
c        and jacobian determinants for the three-dimensional,
c        elastic truss element
c
c           xl(j,l) = global coordinates of nodal points
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c        shg(1,i,l) = global ("arc-length") derivative of shape ftn
c        shg(2,i,l) = shl(2,i,l)
c           xs(j,l) = jth component of the local derivative
c                        of the position vector; then scaled to
c                        direction cosine
c            det(l) = euclidean length of xs
c                 i = local node number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical lnode3
      dimension xl(3,1),det(1),shl(2,nen,1),shg(2,nen,1),xs(3,1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      call move(shg,shl,2*nen*nint)
c
      do 400 l=1,nint
      if (.not.lnode3) then
         temp = pt5*shg(1,3,l)
         shg(1,1,l) = shg(1,1,l) + temp
         shg(1,2,l) = shg(1,2,l) + temp
         temp = pt5*shg(2,3,l)
         shg(2,1,l) = shg(2,1,l) + temp
         shg(2,2,l) = shg(2,2,l) + temp
      endif
      det(l) = zero
c
      do 100 j=1,3
      xs(j,l) = rowdot(shl(1,1,l),xl(j,1),2,3,nen)
      det(l) = det(l) + xs(j,l)**2
  100 continue
c
      det(l) = sqrt(det(l))
c
      if (det(l).le.zero) then
         write(iout,1000) nel,neg
         stop
      endif
c
      do 200 j=1,3
      xs(j,l) = xs(j,l)/det(l)
  200 continue
c
      do 300 i=1,nen
      shg(1,i,l) = shl(1,i,l)/det(l)
  300 continue
c
  400 continue
c
      return
c
 1000 format('1','non-positive determinant in element number  ',i5,
     &          ' in element group  ',i5)
      end
c**** new **************************************************************
      subroutine trushl(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a two or three node,
c        one-dimensional element
c
c                 r = local element coordinate ("xi")
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration-point number
c              nint = number of integration points, eq. 1, 2 or 3
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shl(2,nen,1),w(1),ra(3)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      data ra/-1.00,1.00,0.00/,
     &     five9/0.5555555555555555/,eight9/0.8888888888888888/
c
      if (nint.eq.1) then
         w(1) = two
         g = zero
      endif
c
      if (nint.eq.2) then
         w(1) = one
         w(2) = one
         g = one/sqrt(three)
      endif
c
      if (nint.eq.3) then
         w(1) = five9
         w(2) = five9
         w(3) = eight9
         g = sqrt(three/five)
      endif
c
      do 100 l=1,nint
      r = g*ra(l)
c
      shl(1,1,l) = - pt5
      shl(1,2,l) =   pt5
      shl(2,1,l) =   pt5*(one - r)
      shl(2,2,l) =   pt5*(one + r)
c
      if (nen.eq.3) then
         shl(1,3,l) = - two*r
         shl(2,3,l) = one - r**2
c
         temp = - pt5*shl(2,3,l)
         shl(1,1,l) = shl(1,1,l) + r
         shl(1,2,l) = shl(1,2,l) + r
         shl(2,1,l) = shl(2,1,l) + temp
         shl(2,2,l) = shl(2,2,l) + temp
c
      endif
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine trusk(w,det,shg,xs,xl,b,c,dmat,db,elstif,constk,
     &                 nen,nint,nesd,nrowsh,nrowb,nstr,nee)
c
c.... program to form stiffness matrix for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension w(1),det(1),shg(nrowsh,nen,1),xs(nesd,1),xl(nesd,1),
     &          b(nrowb,1),db(nrowb,1),elstif(nee,1)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp1 = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call trusb(b,shg(1,1,l),xs(1,l),nen,nesd,nrowb,nrowsh)
c
c.... set up the constitutive "matrix"
c
      dmat = c*temp1
c
c.... multiply dmat * b
c
      call smult(dmat,b,db,nrowb,nrowb,nstr,nee,1)
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
      subroutine truskd(w,det,shg,xs,xl,b,dl,strain,c,dmat,stress,work,
     &                  elresf,constk,nen,nint,nrowsh,nesd,nrowb,nee)
c
c.... program to form internal force ("-k*d") for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension w(1),det(1),shg(nrowsh,nen,1),xs(nesd,1),xl(nesd,1),
     &          b(nrowb,1),dl(1),work(1),elresf(1)
c
c.... loop on integration points
c
      do 100 l=1,nint
      temp = constk*w(l)*det(l)
c
c.... set up the strain-displacement matrix
c
      call trusb(b,shg(1,1,l),xs,nen,nesd,nrowb,nrowsh)
c
c.... calculate strain
c
      strain = rcdot(b,dl,nrowb,nee)
c
c.... calculate stress
c
      stress = c*strain
c
c.... calculate element internal force
c
      stress = temp*stress
      call smult(stress,b,work,nrowb,1,1,nee,1)
      call matadd(elresf,work,elresf,nee,nee,nee,nee,1,1)
c
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine truspr(rho,rdampm,rdampk,area,c,numat)
c
c.... program to read, write and store properties for
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension rho(1),rdampm(1),rdampk(1),area(1),c(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /iounit/ iin,iout,irsin,irsout
c
      do 100 n=1,numat
      if (mod(n,50).eq.1) write(iout,1000) numat
      read(iin,2000) m,e,rho(m),rdampm(m),rdampk(m),area(m)
      write(iout,3000) m,e,rho(m),rdampm(m),rdampk(m),area(m)
      c(m) = e
  100 continue
c
      return
c
 1000 format('1',
     &' m a t e r i a l   s e t   d a t a                      '   //5x,
     &' number of material sets . . . . . . . . . . (numat ) = ',i5///,
     &  7x,'set',5x,'young''s',6x,'mass',8x,'mass',
     &  6x,'stiffness',6x,'area',/6x,'number',3x,'modulus',
     &  5x,'density',5x,'damping',5x,'damping',/)
 2000 format(i5,5x,7f10.0)
 3000 format(4x,i5,3x,5(2x,1pe10.4))
      end
c**** new **************************************************************
      subroutine truspt(xint,stress,force,strain,nn,nntot,neg,nel,lint)
c
c.... program to print stress, strain and force for the
c        three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension xint(3)
      common /iounit/ iin,iout,irsin,irsout
c
      nn = nn + 1
c
      if (mod(nn,nntot).eq.1) then
         write(iout,1000) neg
         nn = 1
      endif
c
      write(iout,2000) nel,lint,xint,stress,force,strain
c
      return
c
 1000 format('1',
     &' e l e m e n t   s t r e s s e s   a n d   s t r a i n s',  //5x,
     &' element group number  . . . . . . . . . . . (neg   ) = ',i5///,
     &'  element  int. pt.          x1        x2        x3    ',5x,
     &'  stress    force     strain  ',/
     &'   number   number         coord.    coord.    coord.  ')
 2000 format(/2x,i5,7x,i2,8x,3(1pe10.2),5x,3(1pe10.2))
      end
c**** new **************************************************************
      subroutine trustr(shg,xs,b,dl,strain,c,stress,force,area,
     &                  nrowsh,nesd,nrowb,nen,nee)
c
c.... program to calculate stress, strain and force at an integration
c        point for the three-dimensional, elastic truss element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension shg(nrowsh,1),xs(nesd,1),b(nrowb,1),dl(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
c.... set up strain-displacement matrix
c
      call trusb(b,shg,xs,nen,nesd,nrowb,nrowsh)
c
c.... calculate strain
c
      strain = rcdot(b,dl,nrowb,nee)
c
c.... calculate stress
c
      stress = c*strain
c
c.... calculate forces
c
      force = area*stress
c
      return
      end
