c**** new **************************************************************
      subroutine driver(ntstep,neq,nalhs)
c
c.... solution driver program
c
      double precision
     &       time,coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &            coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn,lout
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /etimec/ etime(7)
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
c.... time sequence loop
c
      do 300 nsq=1,numseq
c
c.... set current time sequence parameters
c
      call timcon(nsq,a(mpstep),a(mpdprt),a(mpsprt),a(mphplt),a(mpiter),
     &                   nstep1,   ndprt1,   nsprt1,   nhplt1,   niter1,
     &                a(mpalph),a(mpbeta),a(mpgamm),a(mpdt  ))
c
c.... form effective mass matrix
c
      call clear(a(mpalhs),nalhs)
      call timing(t1)
      call elemnt('form_lhs',a(mpngrp))
      call timing(t2)
      etime(3) = etime(3) + t2 - t1
c
c.... perform factorization of effective mass matrix
c
      call factor(a(mpalhs),a(mpdiag),neq)
      call timing(t1)
      etime(4) = etime(4) + t1 - t2
c
c.... rank check (note: subroutines "pivots" and "printp"
c        return to statement 300)
c
      if (irank.eq.1) call pivots(a(mpalhs),a(mpdiag),neq,nsq,*300)
      if (irank.eq.2) call printp(a(mpalhs),a(mpdiag),neq,nsq,*300)
c
c.... time step loop
c
      do 200 n=1,nstep1
c
      time = time + dt1
      ntstep = ntstep + 1
c
      if (ldyn) then
c
c....... predictor update of all degrees-of-freedom
c
         call predct(a(mpd),a(mpv),a(mpa),a(mpdprd),a(mpvprd),
     &               ndof,numnp)
      else
         call clear(a(mpd),ndof*numnp)
      endif
c
c
c.... evaluate load-time functions at time n+1
c
      if (nlvect.gt.0)
     &   call lfac(a(mpg),time,a(mpg1),nltftn,nptslf)
c
c
c....... overwrite predictors to account for kinematic
c           boundary conditions
c
      if (dt1.ne.0)
     &   call compbc(a(mpid),a(mpd),a(mpv),a(mpa),a(mpdprd),a(mpvprd),
     &               a(mpf),a(mpg1),ndof,numnp,nlvect,ldyn)
c
c.... multi-corrector iteration loop
c
      do 100 i=1,niter1
c
      call clear(a(mpbrhs),neq)
c
c
c....... evaluate load-time functions at time n+1+alpha
c
      if (nltftn.gt.0)
     &   call lfac(a(mpg),time+alpha1*dt1,a(mpg1),nltftn,nptslf)
c
c
c....... form nodal contribution to residual force vector
c
      if (nlvect.gt.0)
     &   call load(a(mpid),a(mpf),a(mpbrhs),a(mpg1),ndof,numnp,nlvect)
c
c.... form element contribution to residual force vector
c
      call timing(t1)
      call elemnt('form_rhs',a(mpngrp))
      call timing(t2)
      etime(5) = etime(5) + t2 - t1
c
c.... solve equation system
c
      call back(a(mpalhs),a(mpbrhs),a(mpdiag),neq)
      call timing(t1)
      etime(6) = etime(6) + t1 - t2
c
c.... perform intermediate update of active degrees-of-freedom
c
      call iterup(a(mpid),a(mpd),a(mpdprd),a(mpvprd),a(mpa),a(mpbrhs),
     &            ndof,numnp,ldyn)
  100 continue
c
c
c.... perform corrector update of all degrees-of-freedom
c
      if (ldyn)
     &   call corrct(a(mpd),a(mpv),a(mpdprd),a(mpvprd),ndof,numnp)
c
      if (lout(n,ndprt1)) then
c
c....... write kinematic output
c
         call printd(' d i s p l a c e m e n t s                  ',
     &               a(mpd),ndof,numnp,ntstep,time)
         if (ldyn) then
            call printd(' v e l o c i t i e s                        ',
     &                  a(mpv),ndof,numnp,ntstep,time)
            call printd(' a c c e l e r a t i o n s                  ',
     &                  a(mpa),ndof,numnp,ntstep,time)
         endif
      endif
c
      call timing(t1)
c
c.... calculate and write element output
c
      if (lout(n,nsprt1)) call elemnt('str_prnt',a(mpngrp))
c
      if (ldyn.and.lout(n,nhplt1)) then
         locplt = locplt + 1
c
c....... note: variables "locplt" and "time" are passed into sub-
c              routine "stored" and element routines performing
c              task5 ('str_stor') by way of common /hplotc/.
c
c....... store kinematic time history data
c
cZACE MODIF (RD 15 may 1987) -------(T.J.R. Hughes tel. call)-----------
         IF (NDOUT.GT.0)
     &    call stored(a(mpidhs),a(mpd),a(mpv),a(mpa),a(mpdout),
     &               ndof,ndout,ndout+1)
c
c....... calculate and store element time history data
c
         call elemnt('str_stor',a(mpngrp))
      endif
      call timing(t2)
      etime(7) = etime(7) + t2 - t1
c
  200 continue
c
  300 continue
      if (ldyn) then
c
c....... plot nodal time-histories
c
         if (ndout.gt.0) 
     X   call hplot (a(mpidhs),a(mpdout),ndout,3,0)
c
c....... plot element time-histories
c
         call elemnt('str_plot',a(mpngrp))
      endif
c
      return
      end
c**** new **************************************************************
      subroutine addlhs(alhs,eleffm,idiag,lm,nee,ldiag)
c
c.... program to add element left-hand-side matrix to
c        global left-hand-side matrix
c
c        ldiag = .true.,  add diagonal element matrix
c
c        ldiag = .false., add upper triangle of full element matrix
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldiag
      dimension alhs(1),eleffm(nee,1),idiag(1),lm(1)
c
      if (ldiag) then
c
         do 100 j=1,nee
         k = lm(j)
         if (k.gt.0) then
            l = idiag(k)
            alhs(l) = alhs(l) + eleffm(j,j)
         endif
  100    continue
c
      else
c
         do 300 j=1,nee
         k = lm(j)
         if (k.gt.0) then
c
            do 200 i=1,j
            m = lm(i)
            if (m.gt.0) then
               if (k.ge.m) then
                  l = idiag(k) - k + m
               else
                  l = idiag(m) - m + k
               endif
               alhs(l) = alhs(l) + eleffm(i,j)
            endif
  200       continue
c
         endif
  300    continue
c
      endif
c
      return
      end
c**** new **************************************************************
      subroutine addrhs(brhs,elresf,lm,nee)
c
c.... program to add element residual-force vector to
c        global right-hand-side vector
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension brhs(1),elresf(1),lm(1)
c
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) brhs(k) = brhs(k) + elresf(j)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine back(a,b,idiag,neq)
c
c.... program to perform forward reduction and back substitution
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1),b(1),idiag(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
c.... forward reduction
c
      jj = 0
c
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1)
     &   b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neq
      ajj = a(idiag(j))
c
c.... warning: diagonal scaling is not performed if ajj equals zero
c
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neq.eq.1) return
      jjnext = idiag(neq)
c
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return
      end
c**** new **************************************************************
      subroutine bc(id,ndof,numnp,neq,iprtin)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
      dimension id(ndof,1)
c
      common /iounit/ iin,iout,irsin,irsout
      logical pflag
c
      call iclear(id,ndof*numnp)
      call igen(id,ndof)
c
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
c
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
c
         if (pflag) then
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iout,1000) (i,i=1,ndof)
            write(iout,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
c
c.... establish equation numbers
c
      neq = 0
c
      do 400 n=1,numnp
c
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
c
  300 continue
c
  400 continue
c
      return
c
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n   c o
     & d e s'///
     & 5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i5,5x,6(5x,i5))
c
      end
c**** new **************************************************************
      subroutine btdb(elstif,b,db,nee,nrowb,nstr)
c
c.... program to multiply b(transpose) * db taking account of symmetry
c        and accumulate into element stiffness matrix
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension elstif(nee,1),b(nrowb,1),db(nrowb,1)
c
      do 200 j=1,nee
c
      do 100 i=1,j
      elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine clear(a,m)
c
c.... program to clear a floating-point array
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 100 i=1,m
      a(i) = zero
  100 continue
c
      return
      end
c**** new **************************************************************
      function coldot(a,b,n)
c
c.... program to compute the dot product of vectors stored column-wise
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1),b(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      coldot = zero
c
      do 100 i=1,n
      coldot = coldot + a(i)*b(i)
  100 continue
c
      return
      end
c**** new **************************************************************
      subroutine colht(idiag,lm,ned,nen,numel)
c
c.... program to compute column heights in global left-hand-side matrix
c
      dimension idiag(1),lm(ned,nen,1)
      common /colhtc/ neq
c
      do 500 k=1,numel
      min = neq
c
      do 200 j=1,nen
c
      do 100 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) min = min0(min,num)
  100 continue
c
  200 continue
c
      do 400 j=1,nen
c
      do 300 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) then
         m = num - min
         if (m.gt.idiag(num)) idiag(num) = m
      endif
c
  300 continue
c
  400 continue
c
  500 continue
c
      return
      end
c**** new **************************************************************
      subroutine compbc(id,d,v,a,dpred,vpred,f,g1,
     &                  ndof,numnp,nlvect,ldyn)
c
c.... program to compute displacement, velocity and
c        acceleration boundary conditions
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical ldyn
      dimension id(ndof,1),d(ndof,1),v(ndof,1),a(ndof,1),
     &          dpred(ndof,1),vpred(ndof,1),f(ndof,numnp,1),g1(1)
c
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt1
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      do 700 i=1,ndof
c
      do 600 j=1,numnp
c
      k = id(i,j)
      if (k.gt.0) go to 500
      val = zero
      do 100 lv=1,nlvect
      val = val + f(i,j,lv)*g1(lv)
  100 continue
c
      m = 1 - k
      go to (200,300,400),m
c
  200 continue
      if (ldyn) then
         temp = coeff1*val - alpha1*d(i,j)
         a(i,j) = (temp - dpred(i,j))/coeff5
         dpred(i,j) = temp
         vpred(i,j) = vpred(i,j) + coeff4*a(i,j)
      else
         d(i,j) = val
      endif
      go to 500
c
  300 temp = coeff1*val - alpha1*v(i,j)
      a(i,j) = (temp - vpred(i,j))/coeff4
      vpred(i,j) = temp
      dpred(i,j) = dpred(i,j) + coeff5*a(i,j)
      go to 500
c
  400 dpred(i,j) = dpred(i,j) + coeff5*val
      vpred(i,j) = vpred(i,j) + coeff4*val
      a(i,j) = val
c
  500 continue
c
  600 continue
  700 continue
      return
      end
c**** new **************************************************************
      subroutine contm(shg,xl,w,det,elmass,work,constm,imass,nint,
     &                 nrowsh,nesd,nen,ned,nee,column)
c
c.... program to form mass matrix for a continuum element
c        with "nen" nodes
c
c        imass = mass code, eq. 0, consistent mass
c                           eq. 1, lumped mass
c                           otherwise return
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      logical column
      dimension shg(nrowsh,nen,1),xl(nesd,1),w(1),det(1),
     &          elmass(nee,1),work(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      if (imass.eq.0) then
c
c....... consistent mass
c
         do 400 l=1,nint
         temp1 = constm*w(l)*det(l)
c
         do 300 j=1,nen
         n = (j - 1)*ned
c
         do 200 i=1,j
         m = (i - 1)*ned
         temp2 = temp1*shg(nrowsh,i,l)*shg(nrowsh,j,l)
c
         do 100 k=1,ned
         elmass(m + k,n + k) = elmass(m + k,n + k) + temp2
  100    continue
c
  200    continue
c
  300    continue
c
  400    continue
c
      endif
c
      if (imass.eq.1) then
c
c....... lumped mass
c
         dsum   = zero
         totmas = zero
         call clear(work,nen)
c
         do 600 l=1,nint
         temp1 = constm*w(l)*det(l)
         totmas = totmas + temp1
c
         do 500 j=1,nen
         temp2 = temp1*shg(nrowsh,j,l)**2
         dsum = dsum + temp2
         work(j) = work(j) + temp2
  500    continue
c
  600    continue
c
c....... scale diagonal to conserve total mass
c
         temp1 = totmas/dsum
c
         if (column) then
c
c.......... store terms in first column of matrix
c
            do 800 j=1,nen
            temp2 = temp1*work(j)
            n = (j - 1)*ned
c
            do 700 k=1,ned
            elmass(n + k,1) = elmass(n + k,1) + temp2
  700       continue
c
  800       continue
c
         else
c
c.......... store terms along diagonal of matrix
c
            do 1000 j=1,nen
            temp2 = temp1*work(j)
            n = (j - 1)*ned
c
            do 900 k=1,ned
            elmass(n + k,n + k) = elmass(n + k,n + k) + temp2
  900       continue
c
 1000       continue
c
         endif
c
      endif
c
      return
      end
c**** new **************************************************************
      subroutine contma(shg,xl,w,det,al,elmass,work,elresf,constm,imass,
     &                  nint,nrowsh,nesd,nen,ned,nee)
c
c.... program to calculate inertial and gravity/body force ("-m*(a-g)")
c        for a continuum element with "nen" nodes
c
c        imass = mass code, eq. 0, consistent mass
c                           eq. 1, lumped mass
c                           otherwise return
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension shg(nrowsh,nen,1),xl(nesd,1),w(1),det(1),
     &          al(ned,1),elmass(nee,1),work(1),elresf(ned,1)
c
      if (imass.eq.0) then
c
c....... consistent mass
c
         do 300 l=1,nint
         temp = constm*w(l)*det(l)
c
         do 200 i=1,ned
         acc = rowdot(shg(nrowsh,1,l),al(i,1),nrowsh,ned,nen)
c
         do 100 j=1,nen
         elresf(i,j) = elresf(i,j) + temp*acc*shg(nrowsh,j,l)
  100    continue
c
  200    continue
c
  300    continue
c
      endif
c
      if (imass.eq.1) then
c
c....... lumped mass
c
         call clear(elmass,nee)
         call contm(shg,xl,w,det,elmass,work,constm,imass,nint,
     &              nrowsh,nesd,nen,ned,nee,.true.)
c
         do 500 j=1,nen
         k = (j - 1)*ned
c
         do 400 i=1,ned
         elresf(i,j) = elresf(i,j) + al(i,j)*elmass(k + i,1)
  400    continue
c
  500    continue
c
      endif
c
      return
      end
c**** new **************************************************************
      subroutine coord(x,nsd,numnp,iprtin)
c
c.... program to read, generate and write coordinate data
c
c        x(nsd,numnp) = coordinate array
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension x(nsd,1)
      common /iounit/ iin,iout,irsin,irsout
c
      call genfl(x,nsd)
c
      if (iprtin.eq.1) return
c
      do 100 n=1,numnp
      if (mod(n,50).eq.1) write(iout,1000) (i,i=1,nsd)
      write(iout,2000) n,(x(i,n),i=1,nsd)
  100 continue
c
      return
c
 1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,
     &' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i5,10x,3(1pe15.8,2x))
      end
c**** new **************************************************************
      subroutine corrct(d,v,dpred,vpred,ndof,numnp)
c
c.... program to perform corrector update of displacements
c        and velocities
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension d(ndof,1),v(ndof,1),dpred(ndof,1),vpred(ndof,1)
c
      common /coeffs/ coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,
     &                coeff7,coeff8,alpha1,beta1 ,gamma1,dt
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      temp = one/coeff1
c
      do 200 i=1,ndof
c
      do 100 j=1,numnp
      dn = d(i,j)
      vn = v(i,j)
      d(i,j) = (dpred(i,j) - dn)*temp + dn
      v(i,j) = (vpred(i,j) - vn)*temp + vn
  100 continue
c
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast)
c
c.... program to store pointer information in dictionary
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      character*4 name
      dimension name(2)
c-ZACE-2005.08
      include 'memory_size.inc'
      common ia(max_size)
c
      mlast = mlast - 7
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
c      ia(mlast+1) = name(1)
c      ia(mlast+2) = name(2)
      call chrint(ia(mlast+1),name(1))
      call chrint(ia(mlast+2),name(2))
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      ia(mlast+3) = mpoint
      ia(mlast+4) = ndim1
      ia(mlast+5) = ndim2
      ia(mlast+6) = ndim3
      ia(mlast+7) = ipr
c
      return
      end
c**** new **************************************************************
      subroutine diag(idiag,neq,n)
c
c.... program to compute diagonal addresses of left-hand-side matrix
c
      dimension idiag(1)
c
      n = 1
      idiag(1) = 1
      if (neq.eq.1) return
c
      do 100 i=2,neq
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neq)
c
      return
      end
c**** new **************************************************************
      subroutine dynpts 
c
c.... program to set memory pointers for dynamic analysis data arrays
c
      logical ldyn
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /dpoint/ mpstep,mpdprt,mpsprt,mphplt,mpiter,mpalph,mpbeta,
     &                mpgamm,mpdt  ,mpidhs,mpdout,mpvprd,mpdprd,mpa,mpv
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
c-ZACE-2005.08
      if (ldyn) then
        mpvprd = mpoint('vpred   ',ndof   ,numnp ,0,iprec)
        mpdprd = mpoint('dpred   ',ndof   ,numnp ,0,iprec)
        mpa    = mpoint('a       ',ndof   ,numnp ,0,iprec)
        mpv    = mpoint('v       ',ndof   ,numnp ,0,iprec)
      else
        mpvprd = mpoint('vpred   ',0   ,0 ,0,iprec)
        mpdprd = mpoint('dpred   ',0   ,0 ,0,iprec)
        mpa    = mpoint('a       ',0   ,0 ,0,iprec)
        mpv    = mpoint('v       ',0   ,0 ,0,iprec)
      end if
c
      return
      end
c**** new **************************************************************
      subroutine echo
c
c.... program to echo input data
c
      dimension ia(20)
      common /iounit/ iin,iout,irsin,irsout
c
      read(iin,1000) iecho
      if (iecho.eq.0) return
c
      write(iout,2000) iecho
      backspace iin
c
      do 100 i=1,100000
      read(iin,3000,end=200) ia
      if (mod(i,50).eq.1.and.i.gt.1) write(iout,4000)
      write(iout,5000) ia
  100 continue
c
  200 continue
      rewind iin
      read(iin,1000) iecho
c
      return
c
 1000 format(16i5)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,
     &' echo print code . . . . . . . . . . . . . . (iecho ) = ',i5//5x,
     &'    eq. 0, no echo of input data                        ',   /5x,
     &'    eq. 1, echo input data                              ',   ///,
     & 8('123456789*'),//)
 3000 format(20a4)
 4000 format('1',8('123456789*'),//)
 5000 format(' ',20a4)
      end
c**** new **************************************************************
      subroutine elemnt(task,ngrp)
c
c.... program to calculate element task number
c
      character*8 task,eltask(6)
      dimension ngrp(1)
      common /info  / iexec,iacode,ldyn,ireadr,iwritr,iprtin,irank,
     &                numseq,ndout,nsd,numnp,ndof,nlvect,nltftn,nptslf,
     &                numeg
      include 'memory_size.inc'
c-ZACE-2005.08
      common ia(max_size)
      data ntask,    eltask
     &    /    6,'input___',
     &           'form_lhs',
     &           'form_rhs',
     &           'str_prnt',
     &           'str_stor',
     &           'str_plot'/
c
      do 100 i=1,ntask
      if (task.eq.eltask(i)) itask = i
  100 continue
c
      do 200 neg=1,numeg
c
      if (itask.eq.1) then
         mpnpar = mpoint('npar    ',16   ,0,0,1)
         ngrp(neg) = mpnpar
         call elcard(ia(mpnpar),neg)
      else
         mpnpar = ngrp(neg)
      endif
c
      ntype  = ia(mpnpar)
      call elmlib(ntype,mpnpar,itask,neg)
  200 continue
c
      return
      end
c**** new **************************************************************
      subroutine elcard(npar,neg)
c
c.... program to read element group control card
c
      dimension npar(1)
      common /iounit/ iin,iout,irsin,irsout
c
      read(iin,1000) (npar(i),i=1,16)
      write(iout,2000) neg
c
      return
c
 1000 format(16i5)
 2000 format('1',' e l e m e n t   g r o u p   d a t a         ',  //5x,
     &' element group number  . . . . . . . . . . . (neg   ) = ',i5/// )
c
      end
c**** new **************************************************************
      subroutine eqset(neq,nalhs)
c
c.... program to allocate storage for global equation system
c
      character*4 title
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /iounit/ iin,iout,irsin,irsout
      common /spoint/ mpd,mpx,mpid,mpf,mpg,mpg1,mpdiag,mpngrp,
     &                mpalhs,mpbrhs
      common /titlec/ title(20)
c-ZACE-2005.08
      include 'memory_size.inc'
      common a(MAX_SIZE)
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      dimension iadum(MAX_SIZE)
      equivalence(a(1),iadum(1))
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
c
c.... determine addresses of diagonals in left-hand-side matrix
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
      call diag(iadum(mpdiag),neq,nalhs)
      mpalhs = mpoint('alhs    ',nalhs,0,0,iprec)
      mpbrhs = mpoint('brhs    ',neq  ,0,0,iprec)
      meanbw = nalhs/neq
      nwords = mtot - mlast + mfirst - 1
c
c.... write equation system data
c
      write(iout,1000) title,neq,nalhs,meanbw,nwords
c
      return
 1000 format('1',20a4///
     &' e q u a t i o n    s y s t e m    d a t a              ',  //5x,
     &' number of equations . . . . . . . . . . . . (neq   ) = ',i8//5x,
     &' number of terms in left-hand-side matrix  . (nalhs ) = ',i8//5x,
     &' mean half bandwidth . . . . . . . . . . . . (meanbw) = ',i8//5x,
     &' total length of blank common required . . . (nwords) = ',i8    )
c
      end
c**** new **************************************************************
      subroutine factor(a,idiag,neq)
c
c.... program to perform Crout factorization: a = u(transpose) * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(1),idiag(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      jj = 0
c
      do 300 j=1,neq
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         jlngth = i - istart + 1
         length = min0(icolht-1,jlngth)
         if (length.gt.0)
     &      a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)
         ij = ij + 1
  100    continue
c
      endif
c
      if (jcolht.ge.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c
         do 200 ij=jjlast+1,jj-1
c
         ii = idiag(jtemp + ij)
c
c....... warning: the following calculations are skipped
c                 if a(ii) equals zero
c
         if (a(ii).ne.zero) then
            temp  = a(ij)
            a(ij) = temp/a(ii)
            a(jj) = a(jj) - temp*a(ij)
         endif
  200    continue
c
      endif
c
  300 continue
c
      return
      end
c**** new **************************************************************
      subroutine formlm (id,ien,lm,ndof,ned,nen,numel)
c
c.... program to form lm array
c
      dimension id(ndof,1),ien(nen,1),lm(ned,nen,1)
c
      do 300 k=1,numel
c
      do 200 j=1,nen
      node=ien(j,k)
c
      do 100 i=1,ned
      lm(i,j,k) = id(i,node)
  100 continue
c
  200 continue
c
  300 continue
c
      return
      end
c***********************************************************************
      subroutine genel(ien,mat,nen)
c
c.... program to read and generate element node and material numbers
c
c         ien(nen,numel) = element node numbers
c         mat(numel)     = element material numbers
c         nen            = number of element nodes (le.27)
c         n              = element number
c         ng             = generation parameter
c         nel(i)         = number of elements in direction i
c         incel(i)       = element number increment for direction i
c         inc(i)         = node number increment for direction i
c
      dimension ien(nen,1),mat(1),itemp(27)
      common /iounit/ iin,iout,irsin,irsout
      common /genelc/ n,nel(3),incel(3),inc(3)
c
  100 continue
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng
      if (n.eq.0) return
      call imove(ien(1,n),itemp,nen)
      mat(n)=m
      if (ng.ne.0) then
c
c....... generate data
c
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)
         call genel1(ien,mat,nen)
      endif
      go to 100
c
 1000 format(16i5,10x,14i5)
c
      end
c***********************************************************************
      subroutine genel1(ien,mat,nen)
c
c.... program to generate element node and material numbers
c
      dimension ien(nen,1),mat(1)
      common /genelc/ n,nel(3),incel(3),inc(3)
c
c.... set defaults
c
      call geneld
c
c.... generation algorithm
c
      ie = n
      je = n
      ke = n
c
      ii = nel(1)
      jj = nel(2)
      kk = nel(3)
c
      do 300 k=1,kk
c
      do 200 j=1,jj
c
      do 100 i=1,ii
c
      if (i.ne.ii) then
         le = ie
         ie = le + incel(1)
         call geneli(ien(1,ie),ien(1,le),inc(1),nen)
         mat(ie) = mat(le)
      endif
  100 continue
c
      if (j.ne.jj) then
         le = je
         je = le + incel(2)
         call geneli(ien(1,je),ien(1,le),inc(2),nen)
         mat(je) = mat(le)
         ie = je
      endif
  200 continue
c
      if (k.ne.kk) then
         le = ke
         ke = le + incel(3)
         call geneli(ien(1,ke),ien(1,le),inc(3),nen)
         mat(ke) = mat(le)
         ie = ke
         je = ke
      endif
  300 continue
c
      return
      end
c***********************************************************************
      subroutine geneld
c
c.... program to set defaults for element node
c        and material number generation
c
      common /genelc/ n,nel(3),incel(3),inc(3)
c
      if (nel(1).eq.0) nel(1) = 1
      if (nel(2).eq.0) nel(2) = 1
      if (nel(3).eq.0) nel(3) = 1
c
      if (incel(1).eq.0) incel(1) = 1
      if (incel(2).eq.0) incel(2) = nel(1)
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)
c
      if (inc(1).eq.0) inc(1) = 1
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)
c
      return
      end
c***********************************************************************
      subroutine geneli(ien2,ien1,inc,nen)
c
c.... program to increment element node numbers
c
      dimension ien1(1),ien2(1)
c
      do 100 i=1,nen
      if (ien1(i).eq.0) then
         ien2(i) = 0
      else
         ien2(i) = ien1(i) + inc
      endif
  100 continue
c
      return
      end
c***********************************************************************
      subroutine genfl(a,nra)
c
c.... program to read and generate floating-point nodal data
c
c         a       = input array
c         nra     = number of rows in a (le.6)
c         n       = node number
c         numgp   = number of generation points
c         ninc(i) = number of increments for direction i
c         inc(i)  = increment for direction i
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(nra,1)
      common /iounit/ iin,iout,irsin,irsout
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)
c
  100 continue
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra)
      if (n.eq.0) return
      call move(a(1,n),temp,nra)
      if (numgp.ne.0) then
         do 200 j=2,numgp
c
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)
         if (mgen.ne.0) call move(temp(1,j),a(1,m),nra)
c
  200    continue
         read(iin,2000) (ninc(i),inc(i),i=1,3)
         call genfl1(a,nra)
      endif
      go to 100
c
 1000 format(2i5,6f10.0)
 2000 format(16i5)
c
      end
c***********************************************************************
      subroutine genfl1(a,nra)
c
c.... program to generate floating-point nodal data
c        via isoparametric interpolation
c
c         iopt = 1, generation along a line
c              = 2, generation over a surface
c              = 3, generation within a volume
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(nra,1),sh(20)
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      iopt = 3
      if (ninc(3).eq.0) iopt = 2
      if (ninc(2).eq.0) iopt = 1
c
      dr = zero
      ds = zero
      dt = zero
c
      if (ninc(1).ne.0) dr = two/ninc(1)
      if (ninc(2).ne.0) ds = two/ninc(2)
      if (ninc(3).ne.0) dt = two/ninc(3)
c
      ii = ninc(1)+1
      jj = ninc(2)+1
      kk = ninc(3)+1
c
      ni = n
      nj = n
      nk = n
c
      t = -one
      do 300 k=1,kk
c
      s = -one
      do 200 j=1,jj
c
      r = -one
      do 100 i=1,ii
c
      call gensh(r,s,t,sh,numgp,iopt)
      call multab(temp,sh,a(1,ni),6,20,nra,numgp,nra,1,1)
      ni = ni + inc(1)
      r = r + dr
  100 continue
c
      nj = nj + inc(2)
      ni = nj
      s = s + ds
  200 continue
c
      nk = nk + inc(3)
      ni = nk
      t = t + dt
  300 continue
c
      return
      end
c***********************************************************************
      subroutine gensh(r,s,t,sh,numgp,iopt)
c
c.... program to call shape function routines
c        for isoparametric generation
c
      implicit double precision (a-h,o-z)
c
c.... modify above card for single-precision operation
c
      dimension sh(1)
c
      go to (100,200,300),iopt
c
  100 call gensh1(r,sh,numgp)
      return
c
  200 call gensh2(r,s,sh,numgp)
      return
c
  300 call gensh3(r,s,t,sh,numgp)
      return
c
      end
c***********************************************************************
      subroutine gensh1(r,sh,n)
c
c.... program to compute 1d shape functions
c        for isoparametric generation
c
      implicit double precision (a-h,o-z)
c
c.... modify above card(s) for single-precision operation
c
      dimension sh(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      sh(2) = pt5*r
      sh(1) = pt5 - sh(2)
      sh(2) = pt5 + sh(2)
      if (n.eq.3) then
         sh(3) = one - r*r
         sh(1) = sh(1) - pt5*sh(3)
         sh(2) = sh(2) - pt5*sh(3)
      endif
c
      return
      end
c***********************************************************************
      subroutine gensh2(r,s,sh,n)
c
c.... program to compute 2d shape functions
c        for isoparametric generation
c
      implicit double precision (a-h,o-z)
c
c.... modify above card for single-precision operation
c
      dimension sh(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      r2 = pt5*r
      r1 = pt5 - r2
      r2 = pt5 + r2
      s2 = pt5*s
      s1 = pt5 - s2
      s2 = pt5 + s2
      sh(1) = r1*s1
      sh(2) = r2*s1
      sh(3) = r2*s2
      sh(4) = r1*s2
      if (n.eq.4) return
c
      r3 = one - r*r
      s3 = one - s*s
      sh(5) = r3*s1
      sh(6) = s3*r2
      sh(7) = r3*s2
      sh(8) = s3*r1
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
c
      return
      end
c***********************************************************************
      subroutine gensh3(r,s,t,sh,n)
c
c.... program to compute 3d shape functions
c        for isoparametric generation
c
      implicit double precision (a-h,o-z)
c
c.... modify above card for single-precision operation
c
      dimension sh(1)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five
c
      r2 = pt5*r
      r1 = pt5 - r2
      r2 = pt5 + r2
      s2 = pt5*s
      s1 = pt5 - s2
      s2 = pt5 + s2
      t2 = pt5*t
      t1 = pt5 - t2
      t2 = pt5 + t2
c
      rs1 = r1*s1
      rs2 = r2*s1
      rs3 = r2*s2
      rs4 = r1*s2
      sh(1) = rs1*t1
      sh(2) = rs2*t1
      sh(3) = rs3*t1
      sh(4) = rs4*t1
      sh(5) = rs1*t2
      sh(6) = rs2*t2
      sh(7) = rs3*t2
      sh(8) = rs4*t2
      if (n.eq.8) return
c
      r3 = one - r*r
      s3 = one - s*s
      t3 = one - t*t
      sh(17) = t3*rs1
      sh(18) = t3*rs2
      sh(19) = t3*rs3
      sh(20) = t3*rs4
      rs1 = r3*s1
      rs2 = s3*r2
      rs3 = r3*s2
      rs4 = s3*r1
      sh( 9) = rs1*t1
      sh(10) = rs2*t1
      sh(11) = rs3*t1
      sh(12) = rs4*t1
      sh(13) = rs1*t2
      sh(14) = rs2*t2
      sh(15) = rs3*t2
      sh(16) = rs4*t2
c
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
c
      return
      end
