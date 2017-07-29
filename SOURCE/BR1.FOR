c**** new **************************************************************
      subroutine brickc(itask,npar,mp,neg)
c
c.... program to set storage and call tasks for the
c        eight-node brick, elastic element
c
      double precision time
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      dimension npar(1),mp(1)
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /hplotc/ nplpts,locplt,time
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
      common /spoint/ mpd,mpx,mpid,mpf,mpg,mpg1,mpdiag,mpngrp,
     &                mpalhs,mpbrhs
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
c
      mw     = 1
      mdet   = 2
      mshl   = 3
      mshg   = 4
      mshgbr = 5
      mrho   = 6
      mrdpm  = 7
      mrdpk  = 8
      mc     = 9
      mgrav  = 10
      mien   = 11
      mmat   = 12
      mlm    = 13
      mielno = 14
      miside = 15
      mpress = 16
      mishst = 17
      msout  = 18
      melefm = 19
      mxl    = 20
      mwork  = 21
      mb     = 22
      mdmat  = 23
      mdb    = 24
      mvl    = 25
      mal    = 26
      melres = 27
      mdl    = 28
      mstrn  = 29
      mstrs  = 30
      mpstrn = 31
      mpstrs = 32
c
      ntype  = npar( 1)
      numel  = npar( 2)
      numat  = npar( 3)
      nsurf  = npar( 4)
      nsout  = npar( 5)
      iopt   = npar( 6)
      istprt = npar( 7)
      lfsurf = npar( 8)
      lfbody = npar( 9)
      nicode = npar(10)
      ibbar  = npar(11)
      imass  = npar(12)
      impexp = npar(13)
c
c.... set element parameters
c
      nen    = 8
      ned    = 3
      nee    = nen*ned
      nesd   = 3
      nrowsh = 4
      neesq  = nee*nee
      nrowb  = 6
      nstr   = 6
      if (ibbar.eq.1) nstr = 4
      nint   = 1
      if (nicode.eq.0) nint = 8
      nrint  = 1
      if (itask.eq.1) then
c
c....... set memory pointers
c
c
c        note:  the mp array is stored directly after the npar array,
c               beginning at location mpnpar + 16 of blank common.
c               the variable "junk" is not used subsequently.
c
         junk       = mpoint('mp      ',32     ,0     ,0     ,1)
c
         mp(mw    ) = mpoint('w       ',nint   ,0     ,0     ,iprec)
         mp(mdet  ) = mpoint('det     ',nint   ,0     ,0     ,iprec)
         mp(mshl  ) = mpoint('shl     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshg  ) = mpoint('shg     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshgbr) = mpoint('shgbar  ',nrowsh ,nen   ,nrint ,iprec)
         mp(mrho  ) = mpoint('rho     ',numat  ,0     ,0     ,iprec)
         mp(mrdpm ) = mpoint('rdampm  ',numat  ,0     ,0     ,iprec)
         mp(mrdpk ) = mpoint('rdampk  ',numat  ,0     ,0     ,iprec)
         mp(mc    ) = mpoint('c       ',nrowb  ,nrowb ,numat ,iprec)
         mp(mgrav ) = mpoint('grav    ',nesd   ,0     ,0     ,iprec)
         mp(mien  ) = mpoint('ien     ',nen    ,numel ,0     ,1)
         mp(mmat  ) = mpoint('mat     ',numel  ,0     ,0     ,1)
         mp(mlm   ) = mpoint('lm      ',ned    ,nen   ,numel ,1)
         mp(mielno) = mpoint('ielno   ',nsurf  ,0     ,0     ,1)
         mp(miside) = mpoint('iside   ',nsurf  ,0     ,0     ,1)
         mp(mpress) = mpoint('press   ',4      ,nsurf ,0     ,iprec)
c
      if (nsout.eq.0) then
         mp(mishst) = junk
         mp(msout ) = junk
      else
         mp(mishst) = mpoint('ishist  ',3      ,nsout ,0     ,1)
         mp(msout ) = mpoint('sout    ',nsout+1,nplpts,0     ,1)
      endif
c
         mp(melefm) = mpoint('eleffm  ',nee    ,nee   ,0     ,iprec)
         mp(mxl   ) = mpoint('xl      ',nesd   ,nen   ,0     ,iprec)
         mp(mwork ) = mpoint('work    ',24     ,0     ,0     ,iprec)
         mp(mb    ) = mpoint('b       ',nrowb  ,nee   ,0     ,iprec)
         mp(mdmat ) = mpoint('dmat    ',nrowb  ,nrowb ,0     ,iprec)
         mp(mdb   ) = mpoint('db      ',nrowb  ,nee   ,0     ,iprec)
         mp(mvl   ) = mpoint('vl      ',ned    ,nen   ,0     ,iprec)
         mp(mal   ) = mpoint('al      ',ned    ,nen   ,0     ,iprec)
         mp(melres) = mpoint('elresf  ',nee    ,0     ,0     ,iprec)
         mp(mdl   ) = mpoint('dl      ',ned    ,nen   ,0     ,iprec)
         mp(mstrn ) = mpoint('strain  ',nrowb  ,0     ,0     ,iprec)
         mp(mstrs ) = mpoint('stress  ',nrowb  ,0     ,0     ,iprec)
         mp(mpstrn) = mpoint('pstrn   ',nrowb  ,0     ,0     ,iprec)
         mp(mpstrs) = mpoint('pstrs   ',nrowb  ,0     ,0     ,iprec)
      endif
c
c.... task calls
c
      if (itask.gt.6) return
      go to (100,200,300,400,500,600),itask
c
  100 continue
c
c.... input element data ('input___')
c
      call bkct1(a(mp(mshl  )),a(mp(mw    )),a(mp(mrho  )),
     &           a(mp(mrdpm )),a(mp(mrdpk )),
     &           a(mp(mc    )),a(mp(mgrav )),a(mp(mien  )),
     &           a(mp(mmat  )),a(mpid      ),a(mp(mlm   )),
     &           a(mpdiag    ),a(mp(mielno)),a(mp(miside)),
     &           a(mp(mpress)),              a(mp(mishst)),
     &           ntype ,numel ,numat ,nsurf ,nsout ,
     &           istprt,lfsurf,lfbody,nicode,nint  ,ibbar ,
     &           imass ,impexp,nrowsh,nrowb ,nesd  ,nen   ,
     &           ndof  ,ned   ,iprtin,ldyn  )
c
      return
c
  200 continue
c
c.... form element effective mass and assemble into global
c        left-hand-side matrix  ('form_lhs')
c
      call bkct2(a(mp(melefm)),a(mp(mien  )),a(mpx       ),
     &           a(mp(mxl   )),a(mp(mmat  )),a(mp(mdet  )),
     &           a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mrdpm )),a(mp(mrdpk )),
     &           a(mp(mrho  )),a(mp(mw    )),a(mp(mwork )),
     &           a(mp(mshgbr)),a(mp(mb    )),a(mp(mc    )),
     &           a(mp(mdmat )),a(mp(mdb   )),a(mpalhs    ),
     &           a(mpdiag    ),a(mp(mlm   )),
     &           impexp,imass ,numel ,neesq ,nen   ,nsd   ,
     &           nesd  ,nint  ,neg   ,nrowsh,ldyn  ,ned   ,
     &           ibbar ,nrowb ,nstr  ,nee   )
c
      return
c
  300 continue
c
c.... form element residual-force vector and assemble into global
c        right-hand-side vector ('form_rhs')
c
      call bkct3(a(mp(mmat  )),a(mp(mien  )),a(mpdprd    ),
     &           a(mp(mdl   )),a(mpvprd    ),a(mp(mvl   )),
     &           a(mpa       ),a(mp(mal   )),a(mp(mrdpk )),
     &           a(mp(mrdpm )),a(mp(mrho  )),a(mp(mgrav )),
     &           a(mp(melres)),a(mpx       ),a(mp(mxl   )),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &                         a(mpg1      ),a(mp(mwork )),
     &                         a(mp(mw    )),a(mp(melefm)),
     &           a(mp(mshgbr)),a(mp(mb    )),a(mp(mstrn )),
     &           a(mp(mc    )),a(mp(mstrs )),a(mpbrhs    ),
     &           a(mp(mlm   )),a(mp(mielno)),a(mp(miside)),
     &           a(mp(mpress)),
     &           numel ,ned   ,nen   ,ndof  ,ldyn  ,nee   ,
     &           imass ,nesd  ,lfbody,nsd   ,nint  ,nrowsh,
     &           neg   ,nrowb ,nstr  ,ibbar , nsurf,
     &           lfsurf)
c
      return
c
  400 continue
c
c.... calculate and print element stress/strain output ('str_prnt')
c
      if (istprt.eq.0)
     &   call bkct4(a(mp(mmat  )),a(mp(mien  )),a(mpd       ),
     &              a(mp(mdl   )),a(mpx       ),a(mp(mxl   )),
     &              a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &              a(mp(mwork ))              ,a(mp(mshgbr)),
     &              a(mp(mw    )),a(mp(mb    )),a(mp(mstrn )),
     &              a(mp(mc    )),a(mp(mstrs )),a(mp(mpstrn)),
     &              a(mp(mpstrs)),
     &              nint  ,numel ,nen   ,ndof  ,ned   ,nsd   ,
     &              nesd  ,nrowsh,neg   ,ibbar ,nrowb ,
     &              nee   ,nstr  )
c
      return
c
  500 continue
c
c.... calculate and store element time-histories ('str_stor')
c
      if (nsout.gt.0)
     &   call bkct5(a(mp(mishst)),a(mp(msout )),a(mp(mmat  )),
     &              a(mp(mien  )),a(mpd       ),a(mp(mdl   )),
     &              a(mpx       ),a(mp(mxl   )),a(mp(mdet  )),
     &              a(mp(mshl  )),a(mp(mshg  )),
     &              a(mp(mshgbr)),a(mp(mw    )),a(mp(mb    )),
     &              a(mp(mstrn )),a(mp(mc    )),a(mp(mstrs )),
     &              a(mp(mpstrn)),a(mp(mpstrs)),a(mp(mwork )),
     &              nsout ,nen   ,ndof  ,ned   ,nsd   ,nesd  ,
     &              nrowsh,nint  ,neg   ,ibbar ,nrowb ,
     &              nee   ,nstr  ,nsout+1)
c
      return
c
  600 continue
c
c.... plot element time-histories ('str_plot')
c
      if (nsout.gt.0)
     &   call hplot(a(mp(mishst)),a(mp(msout )),nsout ,3,ntype )
      return
c
      end
c**** new **************************************************************
      subroutine bkct1(shl   ,w     ,rho   ,rdampm,rdampk,
     &                 c     ,grav  ,ien   ,mat   ,id    ,lm    ,
     &                 idiag ,ielno ,iside ,press ,ishist,
     &                 ntype ,numel ,numat ,nsurf ,nsout ,
     &                 istprt,lfsurf,lfbody,nicode,nint  ,ibbar ,
     &                 imass ,impexp,nrowsh,nrowb ,nesd  ,nen   ,
     &                 ndof  ,ned   ,iprtin,ldyn  )
c
c.... program to read, generate and write data for the
c        eight-node brick, elastic element
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      dimension shl(nrowsh,nen,1),w(1),rho(1),rdampm(1),rdampk(1),
     &          c(nrowb,nrowb,1),grav(nesd),ien(nen,1),mat(1),
     &          id(ndof,1),lm(ned,nen,1),idiag(1),ielno(1),iside(1),
     &          press(4,1),ishist(3,1)
      common /iounit/ iin,iout,irsin,irsout
c
      write(iout,1000) ntype,numel,numat,nsurf,nsout,iopt,istprt,
     &                 lfsurf,lfbody
      write(iout,2000) nicode,ibbar
      if (ldyn) write(iout,3000) imass,impexp
c
      call bkcshl(shl,w,nint)
c
      call prop3d(rho,rdampm,rdampk,c,numat,nrowb)
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
      if (nsurf.gt.0) call bkcrsf(ielno,iside,press,nsurf)
c
      if (nsout.gt.0) call shist(ishist,nsout,ntype)
c
      return
c
 1000 format(' ',' e i g h t - n o d e  ',
     &' b r i c k   e l e m e n t s',                              //5x,
     &' element type number . . . . . . . . . . . . (ntype ) = ',i5//5x,
     &' number of elements  . . . . . . . . . . . . (numel ) = ',i5//5x,
     &' number of element material sets . . . . . . (numat ) = ',i5//5x,
     &' number of surface force cards . . . . . . . (nsurf ) = ',i5//5x,
     &' number of stress/strain time histories  . . (nsout ) = ',i5//5x,
     &' free option . . . . . . . . . . . . . . . . (iopt  ) = ',i5//5x,
     &' stress output print code  . . . . . . . . . (istprt) = ',i5//5x,
     &'    eq.0, stress output printed                         ',   /5x,
     &'    eq.1, stress output not printed                     ',  //5x,
     &' surface force load-time function number . . (lfsurf) = ',i5//5x,
     &' body force load-time function number  . . . (lfbody) = ',i5 /5x)
 2000 format(5x,
     &' numerical integration code  . . . . . . . . (nicode) = ',i5//5x,
     &'    eq.0, 2 x 2 x 2 gaussian quadrature                 ',   /5x,
     &'    eq.1, 1-point   gaussian quadrature                 ',  //5x,
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
     &' x-2 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x,
     &' x-3 direction . . . . . . . . . . . . . . = ',      1pe15.8//5x)
c
      end
c**** new **************************************************************
      subroutine bkct2(eleffm,ien   ,x     ,xl    ,mat   ,det   ,
     &                 shl   ,shg   ,rdampm,rdampk,
     &                 rho   ,w     ,work  ,shgbar,b     ,c     ,
     &                 dmat  ,db    ,alhs  ,idiag ,lm    ,
     &                 impexp,imass ,numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,ldyn  ,ned   ,
     &                 ibbar ,nrowb ,nstr  ,nee   )
c
c.... program to calculate effective mass matrix for the
c        eight-node brick, elastic element and
c        assemble into the global left-hand-side matrix
c
c         impexp = 0, implicit time integration
c                = 1, explicit time integration
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,ldiag
      dimension eleffm(nee,1),ien(nen,1),x(nsd,1),xl(nesd,1),mat(1),
     &          det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),
     &          rdampm(1),rdampk(1),rho(1),w(1),work(1),
     &          shgbar(nrowsh,1),b(nrowb,1),c(nrowb,nrowb,1),
     &          dmat(nrowb,1),
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
      call bkcshg(xl,det,shl,shg,nint,nel,neg)
c
      if ( ldyn .and. (imass.ne.2) ) then
c
c....... form mass matrix
c
         constm = (one + rdampm(m)*coeff4)*rho(m)
         if (constm.ne.zero) call contm(shg,xl,w,det,eleffm,work,
     &             constm,imass,nint,nrowsh,nesd,nen,ned,nee,.false.)
c
      endif
c
      if ( (.not.ldyn) .or. (impexp.eq.0) ) then
c
c....... form stiffness matrix
c
         constk = (coeff4*rdampk(m) + coeff5)
         call bkck(shgbar,w,det,shg,b,c(1,1,m),dmat,db,eleffm,constk,
     &            ibbar,nen,nint,nesd,nrowsh,nrowb,nstr,nee)
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
      subroutine bkct3(mat   ,ien   ,dpred ,dl    ,vpred ,vl    ,
     &                 a     ,al    ,rdampk,rdampm,rho   ,grav  ,
     &                 elresf,x     ,xl    ,det   ,shl   ,shg   ,
     &                 g1    ,work  ,w     ,eleffm,
     &                 shgbar,b     ,strain,c     ,stress,brhs  ,
     &                 lm    ,ielno ,iside ,press ,
     &                 numel ,ned   ,nen   ,ndof  ,ldyn  ,nee   ,
     &                 imass ,nesd  ,lfbody,nsd   ,nint  ,nrowsh,
     &                 neg   ,nrowb ,nstr  ,ibbar , nsurf,
     &                 lfsurf)
c
c.... program to calculate residual-force vector for the
c        eight-node brick, elastic element and
c        assemble into the global right-hand-side vector
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,formma,formkd,zeroal,zerodl,zerog
      dimension mat(1),ien(nen,1),dpred(ndof,1),dl(ned,1),vpred(ndof,1),
     &          vl(ned,1),a(ndof,1),al(ned,1),rdampk(1),rdampm(1),
     &          rho(1),grav(1),elresf(1),x(nsd,1),xl(nesd,1),det(1),
     &          shl(nrowsh,nen,1),shg(nrowsh,nen,1),g1(1),work(1),
     &          w(1),eleffm(nee,1),shgbar(nrowsh,1),b(nrowb,1),
     &          strain(1),c(nrowb,nrowb,1),stress(1),brhs(1),
     &          lm(ned,nen,1),ielno(1),iside(1),press(4,1)
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
         call bkcshg(xl,det,shl,shg,nint,nel,neg)
c
         if (formma) then
c
c.......... form inertial and/or body force
c
            constm = - rho(m)
            call contma(shg,xl,w,det,al,eleffm,work,elresf,constm,imass,
     &                  nint,nrowsh,nesd,nen,ned,nee)
         endif
c
         if (formkd) then
c
c.......... form internal force
c
            constk = - one
            call bkckd(shgbar,w,det,shg,b,dl,strain,c(1,1,m),stress,
     &                 work,elresf,constk,ibbar,nen,nint,nrowsh,
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
c        note: assembly of surface loads is performed inside bkcsuf
c
      if ( (nsurf.gt.0) .and. (lfsurf.gt.0) )
     &   call bkcsuf(ielno,ien,x,xl,iside,press,elresf,
     &               brhs,lm,g1(lfsurf),nsurf,nen,nsd,nesd,ned,nee)
c
      return
      end
c**** new **************************************************************
      subroutine bkct4(mat   ,ien   ,d     ,dl    ,x     ,xl    ,
     &                 det   ,shl   ,shg   ,xint  ,shgbar,
     &                 w     ,b     ,strain,c     ,stress,pstrn ,
     &                 pstrs ,
     &                 nint  ,numel ,nen   ,ndof  ,ned   ,nsd   ,
     &                 nesd  ,nrowsh,neg   ,ibbar ,nrowb ,
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
      dimension mat(1),ien(nen,1),d(ndof,1),dl(ned,1),x(nsd,1),
     &          xl(nesd,1),det(1),shl(nrowsh,nen,1),shg(nrowsh,nen,1),
     &          xint(nesd,1),shgbar(nrowsh,1),w(1),b(nrowb,1),
     &          strain(1),
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
      call bkcshg(xl,det,shl,shg,nint,nel,neg)
c
c.... calculate coordinates of integration points
c
      do 100 l=1,nint
      xint(1,l) = rowdot(shg(nrowsh,1,l),xl(1,1),nrowsh,nesd,nen)
      xint(2,l) = rowdot(shg(nrowsh,1,l),xl(2,1),nrowsh,nesd,nen)
  100 continue
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,0,nesd,nrowsh)
c
c.... loop over integration points
c
      do 200 l=1,nint
c
c.... calculate stress, strain and principal values
c
      call bkcstr(shg(1,1,l),shgbar,b,dl,strain,c(1,1,m),stress,
     &      pstrn,pstrs,nrowsh,nesd,nrowb,ibbar,nen,ned,nee,nstr)
c
c.... print stress, strain and principal values
c
      call prts3d(xint(1,l),stress,pstrs,strain,pstrn,
     &            nn,nntot,neg,nel,l)
  200 continue
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine bkct5(ishist,sout  ,mat   ,ien   ,d     ,dl    ,
     &                 x     ,xl    ,det   ,shl   ,shg   ,
     &                 shgbar,w     ,b     ,strain,c     ,stress,
     &                 pstrn ,pstrs ,work  ,
     &                 nsout ,nen   ,ndof  ,ned   ,nsd   ,nesd  ,
     &                 nrowsh,nint  ,neg   ,ibbar ,nrowb ,
     &                 nee   ,nstr  ,nsout1)
c
c.... program to calculate and store element time-histories for the
c        eight-node brick, elastic element
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      real sout
      dimension ishist(3,1),sout(nsout1,1),mat(1),ien(nen,1),d(ndof,1),
     &          dl(ned,1),x(nsd,1),xl(nesd,1),det(1),shl(nrowsh,nen,1),
     &          shg(nrowsh,nen,1),shgbar(nrowsh,nen,1),w(1),b(nrowb,1),
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
      call bkcshg(xl,det,shl,shg,nint,nel,neg)
c
      if (ibbar.eq.1)
     &   call meansh(shgbar,w,det,r,shg,nen,nint,0,nesd,nrowsh)
c
c.... calculate stress, strain and principal values
c
      call bkcstr(shg(1,1,intpt),shgbar,b,dl,strain,c(1,1,m),
     &            stress,pstrn,pstrs,nrowsh,nesd,nrowb,ibbar,nen,ned,
     &            nee,nstr)
c
      do 200 j=1,6
      work(j     ) = stress(j)
      work(j +  6) = pstrs(j)
      work(j + 12) = strain(j)
      work(j + 18) = pstrn(j)
  200 continue
c
      sout(i+1,locplt) = real(work(ncomp))
c
  300 continue
c
      return
      end

