      program main
c
c
c     ************************************************************
c     *                                                          *
c     *               * * *   D L E A R N   * * *                *
c     *                                                          *
c     *                                                          *
c     *               A LINEAR STATIC AND DYNAMIC                *
c     *                                                          *
c     *             FINITE ELEMENT ANALYSIS PROGRAM              *
c     *                                                          *
c     *                                                          *
c     *               T. J. R.   H U G H E S                     *
c     *                                                          *
c     *                  R. M.   F E R E N C Z                   *
c     *                                                          *
c     *                  A. M.   R A E F S K Y                   *
c     *                                                          *
c     *                                                          *
c     *                 1 5   M a y   1 9 8 6                    *
c     *                                                          *
c     ************************************************************
c     *                                                          *
c     *   Reference:  Thomas J. R. Hughes, "The Finite Element   *
c     *   Method -- Linear Static and Dynamic Finite Element     *
c     *   Analysis," Prentice-Hall, Englewood Cliffs, N.J. 1987  *
c     *                                                          *
c     ************************************************************
c     *                                                          *
c     *         NOTICE TO USERS OF THIS PROGRAM PACKAGE          *
c     *                                                          *
c     *   No responsibility or legal liability is assumed by     *
c     *   the authors for the accuracy, completeness, or use-    *
c     *   fulness of any information or process contained in     *
c     *   this computer program package, or for any errors,      *
c     *   mistakes, or misrepresentations that may occur from    *
c     *   the use of this computer program package.  All soft-   *
c     *   ware provided are in "as is, with all defects" con-    *
c     *   dition.  No warranties of any kind, whether statu-     *
c     *   tory, written, oral, expressed, or implied (includ-    *
c     *   ing warranties of fitness and merchantability) shall   *
c     *   apply.                                                 *
c     *                                                          *
c     *   Distribution of this program is made possible by the   *
c     *   authors with the stipulation that the program neither  *
c     *   be sold in whole or in part for direct profit nor      *
c     *   royalties or development charges made for its use.     *
c     *   By acceptance of delivery of this program package,     *
c     *   the purchaser understands the restrictions on the use  *
c     *   and distribution of the program.  The fee paid repre-  *
c     *   sent a charge for materials, duplication, packaging,   *
c     *   mailing, communications, and assorted labor costs.     *
c     *   The legal ownership of the program remains with the    *
c     *   developers.                                            *
c     *                                                          *
c     *   Copyright: 1987 by T.J.R. Hughes, all rights reserved. *
c     *                                                          *
c     ************************************************************
c
c
c.... program to set storage capacity, precision and input/output units
c
      common /bpoint/ mfirst,mlast,mtot,iprec
      common /iounit/ iin,iout,irsin,irsout
      common a(50000)

c        mfirst = address of first available word in blank common
c        mlast  = address of last available word in blank common
c        mtot   = total storage allocated to blank common
c        iprec  = precision flag; eq.1, single precision
c                                 eq.2, double precision
c
      mfirst = 1
      mlast  = 50000
      mtot   = 50000
      iprec  = 2
c
c        iin    = input unit number
c        iout   = output unit number
c        irsin  = restart input unit number
c        irsout = restart output unit number
c
      iin    = 5
      iout   = 6
      irsin  = 7
      irsout = 8
c
c.... system-dependent unit/file specifications
c
c        the following lines are appropriate for the VMS operating
c        system -- change as necessary for other operating systems
c
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
C      call assign(   iin,   'input.dat', 9)
C      call assign(  iout,  'output.dat',10)
C      call assign( irsin, 'rsinput.dat',11)
C      call assign(irsout,'rsoutput.dat',12)
      open(unit=   iin,file='input.dat',status='old')
      open(unit=  iout,file='output.dat')
cZACE MODIF (RD 30 mar 1987) -------------------------------------------
c
      call dlearn
c
c.... system-dependent unit/file specifications
c
c        the following lines are appropriate for the VMS operating
c        system -- change as necessary for other operating systems
c
      call close(iin)
      call close(iout)
c
      stop
      end
