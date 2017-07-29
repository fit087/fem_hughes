c**** new **************************************************************
      subroutine qtsk3(itask,npar,mp,neg)
c
c.... program to set storage and call task 3
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
c
c.... task 3 call
c
  300 continue
c
c.... form element residual-force vector and assemble into global
c        right-hand-side vector ('form_rhs')
c
      call qdct3(a(mp(mmat  )),a(mp(mien  )),a(mpdprd    ),
     &           a(mp(mdl   )),a(mpvprd    ),a(mp(mvl   )),
     &           a(mpa       ),a(mp(mal   )),a(mp(mrdpk )),
     &           a(mp(mrdpm )),a(mp(mrho  )),a(mp(mgrav )),
     &           a(mp(melres)),a(mpx       ),a(mp(mxl   )),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mr    )),a(mpg1      ),a(mp(mwork )),
     &           a(mp(mth   )),a(mp(mw    )),a(mp(melefm)),
     &           a(mp(mshgbr)),a(mp(mb    )),a(mp(mstrn )),
     &           a(mp(mc    )),a(mp(mstrs )),a(mpbrhs    ),
     &           a(mp(mlm   )),a(mp(mielno)),a(mp(miside)),
     &           a(mp(mpress)),a(mp(mshear)),
     &           numel ,ned   ,nen   ,ndof  ,ldyn  ,nee   ,
     &           imass ,nesd  ,lfbody,nsd   ,nint  ,nrowsh,
     &           neg   ,iopt  ,nrowb ,nstr  ,ibbar , nsurf,
     &           lfsurf)
c
      return
c
      end
