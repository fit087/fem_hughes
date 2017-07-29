c**** new **************************************************************
      subroutine tseq
c
c.... program to set memory pointers for time sequence and
c        nodal time history data arrays
c
      double precision time
      logical ldyn
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /hplotc/ nplpts,locplt,time
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
c
      mpstep = mpoint('nstep   ',numseq,0,0,1)
      mpdprt = mpoint('ndprt   ',numseq,0,0,1)
      mpsprt = mpoint('nsprt   ',numseq,0,0,1)
      mphplt = mpoint('nhplt   ',numseq,0,0,1)
      mpiter = mpoint('niter   ',numseq,0,0,1)
      mpalph = mpoint('alpha   ',numseq,0,0,iprec)
      mpbeta = mpoint('beta    ',numseq,0,0,iprec)
      mpgamm = mpoint('gamma   ',numseq,0,0,iprec)
      mpdt   = mpoint('dt      ',numseq,0,0,iprec)
c
      call tseqin(a(mpstep),a(mpdprt),a(mpsprt),a(mphplt),
     &            a(mpiter),a(mpalph),a(mpbeta),a(mpgamm),
     &            a(mpdt  ),numseq,nplpts,ldyn)
c
      if (ndout.eq.0) then
         mpidhs = 1
         mpdout = 1
      else
         mpidhs = mpoint('idhist  ',3      ,ndout ,0,1)	
         mpdout = mpoint('dout    ',ndout+1,nplpts,0,1)
      endif
c
      return
      end
c**** new **************************************************************
      subroutine statin(neq)
c
c.... program to set memory pointers for static analysis data arrays,
c        and call associated input routines
c
      double precision zero,pt1667,pt25,pt5,one,two,three,four,five
c
c.... deactivate above card(s) for single precision operation
c
      logical ldyn
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
      common /spoint/ mpd,mpx,mpid,mpf,mpg,mpg1,mpdiag,mpngrp,
     &                mpalhs,mpbrhs
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
c
      mpd    = mpoint('d       ',ndof  ,numnp ,0    ,iprec)
      if (.not.ldyn) mpdprd = mpd
      mpx    = mpoint('x       ',nsd   ,numnp ,0     ,iprec)
      mpid   = mpoint('id      ',ndof  ,numnp ,0     ,1)
c
      if (nlvect.eq.0) then
         mpf = 1
      else
         mpf = mpoint('f       ',ndof  ,numnp ,nlvect,iprec)
      endif
c
      if (nltftn.eq.0) then
         mpg  = 1
         mpg1 = 1
      else
         mpg  = mpoint('g       ',nptslf,2     ,nltftn,iprec)
         mpg1 = mpoint('g1      ',nltftn,0     ,0     ,iprec)
      endif
c
c.... input coordinate data
c
      call coord(a(mpx),nsd,numnp,iprtin)
c
c.... input boundary condition data and establish equation numbers
c
      call bc(a(mpid),ndof,numnp,neq,iprtin)
c
c.... input nodal force and prescribed kinematic boundary-value data
c
      if (nlvect.gt.0) call input(a(mpf),ndof,numnp,0,nlvect,
     &                            iprtin,zero)
c
c.... input load-time functions
c
      if (nltftn.gt.0) call ltimef(a(mpg),nptslf,nltftn,iprtin)
c
c.... allocate memory for idiag array and clear
c
      mpdiag = mpoint('idiag   ',neq   ,0     ,0     ,1)
      call iclear(a(mpdiag),neq)
c
      mpngrp = mpoint('ngrp    ',numeg ,0     ,0     ,1)
c
      return
      end

