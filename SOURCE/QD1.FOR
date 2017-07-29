c**** new **************************************************************
      subroutine quadc(itask,npar,mp,neg)
c
c.... program to set storage and call tasks for the
c        four-node quadrilateral, elastic continuum element
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
      mr     = 3
      mshl   = 4
      mshg   = 5
      mshgbr = 6
      mrho   = 7
      mrdpm  = 8
      mrdpk  = 9
      mth    = 10
      mc     = 11
      mgrav  = 12
      mien   = 13
      mmat   = 14
      mlm    = 15
      mielno = 16
      miside = 17
      mpress = 18
      mshear = 19
      mishst = 20
      msout  = 21
      melefm = 22
      mxl    = 23
      mwork  = 24
      mb     = 25
      mdmat  = 26
      mdb    = 27
      mvl    = 28
      mal    = 29
      melres = 30
      mdl    = 31
      mstrn  = 32
      mstrs  = 33
      mpstrn = 34
      mpstrs = 35
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
      nen    = 4
      ned    = 2
      nee    = nen*ned
      nesd   = 2
      nrowsh = 3
      neesq  = nee*nee
      nrowb  = 4
      nstr   = 3
      if ( (iopt.eq.2) .or. (ibbar.eq.1) ) nstr = 4
      nint   = 1
      if (nicode.eq.0) nint = 4
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
         junk       = mpoint('mp      ',35     ,0     ,0     ,1)
c
         mp(mw    ) = mpoint('w       ',nint   ,0     ,0     ,iprec)
         mp(mdet  ) = mpoint('det     ',nint   ,0     ,0     ,iprec)
         mp(mr    ) = mpoint('r       ',nint   ,0     ,0     ,iprec)
         mp(mshl  ) = mpoint('shl     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshg  ) = mpoint('shg     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshgbr) = mpoint('shgbar  ',nrowsh ,nen   ,nrint ,iprec)
         mp(mrho  ) = mpoint('rho     ',numat  ,0     ,0     ,iprec)
         mp(mrdpm ) = mpoint('rdampm  ',numat  ,0     ,0     ,iprec)
         mp(mrdpk ) = mpoint('rdampk  ',numat  ,0     ,0     ,iprec)
         mp(mth   ) = mpoint('th      ',numat  ,0     ,0     ,iprec)
         mp(mc    ) = mpoint('c       ',nrowb  ,nrowb ,numat ,iprec)
         mp(mgrav ) = mpoint('grav    ',nesd   ,0     ,0     ,iprec)
         mp(mien  ) = mpoint('ien     ',nen    ,numel ,0     ,1)
         mp(mmat  ) = mpoint('mat     ',numel  ,0     ,0     ,1)
         mp(mlm   ) = mpoint('lm      ',ned    ,nen   ,numel ,1)
         mp(mielno) = mpoint('ielno   ',nsurf  ,0     ,0     ,1)
         mp(miside) = mpoint('iside   ',nsurf  ,0     ,0     ,1)
         mp(mpress) = mpoint('press   ',2      ,nsurf ,0     ,iprec)
         mp(mshear) = mpoint('shear   ',2      ,nsurf ,0     ,iprec)
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
         mp(mwork ) = mpoint('work    ',16     ,0     ,0     ,iprec)
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
      call qdct1(a(mp(mshl  )),a(mp(mw    )),a(mp(mrho  )),
     &           a(mp(mrdpm )),a(mp(mrdpk )),a(mp(mth   )),
     &           a(mp(mc    )),a(mp(mgrav )),a(mp(mien  )),
     &           a(mp(mmat  )),a(mpid      ),a(mp(mlm   )),
     &           a(mpdiag    ),a(mp(mielno)),a(mp(miside)),
     &           a(mp(mpress)),a(mp(mshear)),a(mp(mishst)),
     &           ntype ,numel ,numat ,nsurf ,nsout ,iopt  ,
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
      call qdct2(a(mp(melefm)),a(mp(mien  )),a(mpx       ),
     &           a(mp(mxl   )),a(mp(mmat  )),a(mp(mdet  )),
     &           a(mp(mshl  )),a(mp(mshg  )),a(mp(mr    )),
     &           a(mp(mrdpm )),a(mp(mrdpk )),a(mp(mth   )),
     &           a(mp(mrho  )),a(mp(mw    )),a(mp(mwork )),
     &           a(mp(mshgbr)),a(mp(mb    )),a(mp(mc    )),
     &           a(mp(mdmat )),a(mp(mdb   )),a(mpalhs    ),
     &           a(mpdiag    ),a(mp(mlm   )),
     &           impexp,imass ,numel ,neesq ,nen   ,nsd   ,
     &           nesd  ,nint  ,neg   ,nrowsh,ldyn  ,ned   ,
     &           iopt  ,ibbar ,nrowb ,nstr  ,nee   )
c
      return
c
  300 continue
c
c.... form element residual-force vector and assemble into global
c        right-hand-side vector ('form_rhs')
c
      call qtsk3(itask,npar,mp,neg)
c
      return
c
  400 continue
c
c.... calculate and print element stress/strain output ('str_prnt')
c
      if (istprt.eq.0)
     &   call qdct4(a(mp(mmat  )),a(mp(mien  )),a(mpd       ),
     &              a(mp(mdl   )),a(mpx       ),a(mp(mxl   )),
     &              a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &              a(mp(mwork )),a(mp(mr    )),a(mp(mshgbr)),
     &              a(mp(mw    )),a(mp(mb    )),a(mp(mstrn )),
     &              a(mp(mc    )),a(mp(mstrs )),a(mp(mpstrn)),
     &              a(mp(mpstrs)),
     &              nint  ,numel ,nen   ,ndof  ,ned   ,nsd   ,
     &              nesd  ,nrowsh,neg   ,iopt  ,ibbar ,nrowb ,
     &              nee   ,nstr  )
c
      return
c
  500 continue
c
c.... calculate and store element time-histories ('str_stor')
c
      if (nsout.gt.0)
     &   call qdct5(a(mp(mishst)),a(mp(msout )),a(mp(mmat  )),
     &              a(mp(mien  )),a(mpd       ),a(mp(mdl   )),
     &              a(mpx       ),a(mp(mxl   )),a(mp(mdet  )),
     &              a(mp(mshl  )),a(mp(mshg  )),a(mp(mr    )),
     &              a(mp(mshgbr)),a(mp(mw    )),a(mp(mb    )),
     &              a(mp(mstrn )),a(mp(mc    )),a(mp(mstrs )),
     &              a(mp(mpstrn)),a(mp(mpstrs)),a(mp(mwork )),
     &              nsout ,nen   ,ndof  ,ned   ,nsd   ,nesd  ,
     &              nrowsh,nint  ,neg   ,iopt  ,ibbar ,nrowb ,
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
