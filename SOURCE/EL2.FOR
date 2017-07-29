c**** new **************************************************************
      block data
c
c.... program to define output labels and numerical constants
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
c
      character*4 labeld,label1,label2,label3
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /labels/ labeld(3),label1(16),label2(3),label3(24)
c
c        labeld(3)  = displacement, velocity and acceleration labels
c        label1(16) = output labels for element-type 1
c        label2(3)  = output labels for element-type 2
c        label2(24) = output labels for element-type 3
c
c.... note: add label arrays for any additional elements
c
      data   zero,pt1667,pt25,pt5
     &      /0.00,0.1666666666666667,0.25,0.50/,
     &       one,two,three,four,five
     &      /1.00,2.00,3.00,4.00,5.00/
c
      data labeld/'disp','vel ','acc '/
c
      data label1/'s 11','s 22','s 12','s 33','ps 1','ps 2',
     &            'tau ','sang','e 11','e 22','g 12','e 33',
     &            'pe 1','pe 2','gam ','eang'/
c
      data label2/'strs','forc','strn'/
c
      data label3/'s 11','s 22','s 33','s 12','s 23','s 31',
     &            'ps 1','ps 2','ps 3','----','----','----',
     &            'e 11','e 22','e 33','g 12','g 23','g 31',
     &            'pe 1','pe 2','pe 3','----','----','----'/
c
      end
c**** new **************************************************************
      subroutine elmlib(ntype,mpnpar,itask,neg)
c
c.... program to call element routines
c
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
c
      go to (100,200,300),ntype
c
  100 continue
      call quadc(itask,a(mpnpar),a(mpnpar+16),neg)
      return
c
  200 continue
      call truss(itask,a(mpnpar),a(mpnpar+16),neg)
      return
c
  300 continue
      call brickc(itask,a(mpnpar),a(mpnpar+16),neg)
      return
c
c.... add additional elements for fun and valuable prizes
c
      end
c**** new **************************************************************
      subroutine dhist(idhist,ndout)
c
c.... program to read, write and store nodal time-history input data
c
      dimension idhist(3,1)
      common /iounit/ iin,iout,irsin,irsout
      common /labels/ labeld(3),label1(16),label2(3),label3(24)
c
      do 100 n=1,ndout
      read(iin,2000) node,idof,idva
      if (mod(n,50).eq.1) write(iout,1000) ndout
      write(iout,3000) node,idof,labeld(idva)
      idhist(1,n) = node
      idhist(2,n) = idof
      idhist(3,n) = idva
  100 continue
c
      return
c
 1000 format('1',' n o d a l   t i m e - h i s t o r y  ',
     &' i n f o r m a t i o n'//5x,
     &' number of nodal time histories  . . . . . . (ndout ) = ',i5///
     &5x,'    node      dof     kinematic ',/
     &5x,'   number    number     type    ',/)
 2000 format(3i5)
 3000 format(7x,i5,5x,i5,7x,a4)
c
      end
c**** new **************************************************************
      subroutine dlearn
c
c.... DLEARN - a linear static and dynamic finite element
c                analysis program: global driver
c
      double precision
     &       time,zero,pt1667,pt25,pt5,one,two,three,four,five,tempf,
     &       coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &       coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      character*4 title,ciao
c
c.... catalog of common statements
c
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /colhtc/ neq
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /etimec/ etime(7)
      common /genelc/ n,nel(3),incel(3),inc(3)
      common /genflc/ tempf(6,20),nf,numgpf,nincf(3),incf(3)
      common /hplotc/ nplpts,locplt,time
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
      common /iounit/ iin,iout,irsin,irsout
      common /labels/ labeld(3),label1(16),label2(3),label3(24)
      common /spoint/ mpd,mpx,mpid,mpf,mpg,mpg1,mpdiag,mpngrp,
     &                mpalhs,mpbrhs
      common /titlec/ title(20)
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
      dimension ia(MAX_SIZE)
      equivalence(a(1),ia(1))
      data ciao/'*end'/
c
c.... input phase
c
      call echo
  100 continue
      do 200 i=1,7
  200 etime(i) = 0.0
      call timing(t1)
      read(iin,1000) title
      if (title(1).eq.ciao) return
      read(iin,2000) iexec,iacode,ireadr,iwritr,iprtin,irank,numseq,
     &               ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,numeg
      ldyn = .true.
      if (iacode.eq.1) ldyn = .false.
      write(iout,3000) title , iexec,iacode,ireadr,iwritr,iprtin
      write(iout,4000) irank ,numseq, ndout,   nsd, numnp,  ndof,
     &                 nlvect,nltftn,nptslf, numeg
c
c.... initialization phase
c
      call tseq
      call dynpts 
      if (ndout.gt.0) call dhist(ia(mpidhs),ndout)
      call statin(neq)
c
      ntstep = 0
      time   = zero
      if (ldyn) then
c
         if (ireadr.eq.1) then
c
c.......... read initial conditions from restart file
c
            call rsin(a(mpd),a(mpv),a(mpa),ndof,numnp,ntstep,time)
            if (iprtin.eq.0) then
            call printd(' r e s t a r t   d i s p l a c e m e n t s  ',
     &                  a(mpd),ndof,numnp,ntstep,time)
            call printd(' r e s t a r t   v e l o c i t i e s        ',
     &                  a(mpv),ndof,numnp,ntstep,time)
            call printd(' r e s t a r t   a c c e l e r a t i o n s  ',
     &                  a(mpa),ndof,numnp,ntstep,time)
            endif
c
         else
c
c.......... read initial conditions from input file
c
            call input(a(mpd),ndof,numnp,1,1,iprtin,time)
            call input(a(mpv),ndof,numnp,2,1,iprtin,time)
            call input(a(mpa),ndof,numnp,3,1,iprtin,time)
c
         endif
c
      endif
c
c.... store initial kinematic data for time histories
c
      if ( ldyn .and. ndout.gt.0 ) then
         locplt = 1
         call stored(a(mpidhs),a(mpd),a(mpv),a(mpa),a(mpdout),
     &               ndof,ndout,ndout+1)
      endif
c
c.... input element data
c
      call elemnt('input___',a(mpngrp))
c
c.... store initial stress/strain data for element time histories
c
      if (ldyn) then
         locplt = 1
         call elemnt('str_stor',a(mpngrp))
      endif
c
c.... allocate memory for global equation system
c
      call eqset(neq,nalhs)
      call timing(t2)
      etime(1) = t2 - t1
c
c.... solution phase
c
      if (iexec.eq.1) call driver(ntstep,neq,nalhs)
c
c.... write restart file
c
      if ( ldyn .and. (iwritr.eq.1) )
     &   call rsout(a(mpd),a(mpv),a(mpa),ndof,numnp,ntstep,time)
c
c.... print memory-pointer dictionary
c
      call prtdc
c
      call timing(t1)
      etime(2) = t1 - t2
c
c.... print elapsed time summary
c
      call timlog
      go to 100
c
 1000 format(20a4)
 2000 format(16i5)
 3000 format('1',20a4///
     &' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,
     &' execution code  . . . . . . . . . . . . . . (iexec ) = ',i5//5x,
     &'    eq. 0, data check                                   ',   /5x,
     &'    eq. 1, execution                                    ',  //5x,
     &' analysis code . . . . . . . . . . . . . . . (iacode) = ',i5//5x,
     &'    eq. 0, dynamic analysis                             ',   /5x,
     &'    eq. 1, static analysis                              ',  //5x,
     &' read restart file code  . . . . . . . . . . (ireadr) = ',i5//5x,
     &'    eq. 0, do not read restart file                     ',   /5x,
     &'    eq. 1, read restart file                            ',  //5x,
     &' write restart file code . . . . . . . . . . (iwritr) = ',i5//5x,
     &'    eq. 0, do not write restart file                    ',   /5x,
     &'    eq. 1, write restart file                           ',  //5x,
     &' input data print code . . . . . . . . . . . (iprtin) = ',i5//5x,
     &'    eq. 0, print nodal and element input data           ',   /5x,
     &'    eq. 1, do not print nodal and element input data    ',   /5x)
 4000 format(5x,
     &' rank check code . . . . . . . . . . . . . . (irank ) = ',i5//5x,
     &'    eq. 0, do not perform rank check                    ',   /5x,
     &'    eq. 1, print numbers of zero and negative pivots    ',   /5x,
     &'    eq. 2, print all pivots                             ',  //5x,
     &' number of time sequences  . . . . . . . . . (numseq) = ',i5//5x,
     &' number of nodal output time-histories . . . (ndout ) = ',i5//5x,
     &' number of space dimensions  . . . . . . . . (nsd   ) = ',i5//5x,
     &' number of nodal points  . . . . . . . . . . (numnp ) = ',i5//5x,
     &' number of nodal degrees-of-freedom  . . . . (ndof  ) = ',i5//5x,
     &' number of load vectors  . . . . . . . . . . (nlvect) = ',i5//5x,
     &' number of load-time functions . . . . . . . (nltftn) = ',i5//5x,
     &' number of points on load-time functions . . (nptslf) = ',i5//5x,
     &' number of element groups  . . . . . . . . . (numeg ) = ',i5//5x)
c
      end
c**** new **************************************************************

      subroutine hplot(ih,xt,nplots,nrows,io)
c
c.... program to plot output histories
c
c           ih(nrows,nplots) = dof/component information
c        xt(nplots+1,nplpts) = output history data
c        xt(       1,nplpts) = time record
c                     nplots = number of histories to be plotted
c                     nplpts = number of time points at which
c                              data is to be plotted
c                      nrows = number of rows in ih array
c                         io = output code
c                                 eq.0, nodal output histories
c                                 eq.n.gt.0, element output histories
c                                    (n = ntype in calling routine)
c
      double precision time
c
c.... deactivate above card(s) for single-precision operation
c
      character*1 iblank,istar,line(53)
      character*4 title,labeld,label1,label2,label3
      dimension ih(nrows,1),xt(nplots+1,1)
      common /hplotc/ nplpts,locplt,time
      common /iounit/ iin,iout,irsin,irsout
      common /labels/ labeld(3),label1(16),label2(3),label3(24)
      common /titlec/ title(20)
c
      data iblank,istar/' ','*'/,nchar/53/
c
      do 300 i=1,nplots
c
      i1 = ih(1,i)
      i2 = ih(2,i)
      i3 = ih(3,i)
c
      if (io.eq.0) write(iout,1000) title,i1,i2,labeld(i3)
      if (io.eq.1) write(iout,2000) title,i1,i2,label1(i3)
      if (io.eq.2) write(iout,2000) title,i1,i2,label2(i3)
      if (io.eq.3) write(iout,2000) title,i1,i2,label3(i3)
c
c.... add if/write statements as above for additional element types
c
      call minmax(xt,xmax,xmin,nplots+1,nplpts,i+1)
      if (xmax.eq.xmin) then
c
         write(iout,3000) xmax
c
      else
c
         scale = xmax - xmin
         write(iout,4000) xmin,xmax
c
         do 200 j=1,nplpts
         t = xt(1,j)
c
         do 100 k = 1,nchar
         line(k) = iblank
  100    continue
c
         xk = ((xt(i+1,j) - xmin)/scale)*nchar
         k  = xk + 1
         if (k.gt.nchar) k = nchar
         line(k) = istar
         write(iout,5000) t,xt(i+1,j),(line(k),k=1,nchar)
  200    continue
      endif
c
  300 continue
c
      return
c
 1000 format('1',20a4///
     &' node number = ',i5//
     &' dof number  = ',i5,5x,'output: ',a4//5x)
 2000 format('1',20a4///
     &' element number            = ',i5//
     &' integration point number  = ',i5,5x,'output: ',a4//5x)
 3000 format(' ',
     &' value is constant ( = ',1pe11.4,' ), plot omitted')
 4000 format(' ',4x,'time',8x,'value',6x,1pe11.4,31x,1pe11.4/
     &2x,10('-'),3x,10('-'),3x,53('-'))
 5000 format(' ',1pe11.4,2x,1pe11.4,3x,53a1)
      end


c**** new **************************************************************
      subroutine shist(ishist,nsout,ntype)
c
c.... program to read, write and store element time-history input data
c
      character*4 labeld,label1,label2,label3
      dimension ishist(3,1)
      common /iounit/ iin,iout,irsin,irsout
      common /labels/ labeld(3),label1(16),label2(3),label3(24)
c
      do 100 n=1,nsout
      if (mod(n,50).eq.1) write(iout,1000) nsout
      read(iin,2000) nel,intpt,ncomp
      if (intpt.eq.0) intpt = 1
c
      if (ntype.eq.1) write(iout,3000) nel,intpt,label1(ncomp)
      if (ntype.eq.2) write(iout,3000) nel,intpt,label2(ncomp)
      if (ntype.eq.3) write(iout,3000) nel,intpt,label3(ncomp)
c
c.... add if/write statements as above for additional element types
c
      ishist(1,n) = nel
      ishist(2,n) = intpt
      ishist(3,n) = ncomp
  100 continue
c
      return
c
 1000 format('1',
     &' e l e m e n t   t i m e   h i s t o r y  ',
     &' i n f o r m a t i o n '//5x,
     &' number of stress/strain time histories  . . (nsout ) = ',i5///
     &5x,'   element   int pt   component',/,
     &5x,'   number    number            ',/)
 2000 format(3i5)
 3000 format(7x,i5,5x,i5,7x,1a4)
      end
