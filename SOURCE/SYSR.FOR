c--------------------------------------------------------------------------
      subroutine chrint (ia,a1)
c... converts char*4 to integer*4
      character*4 a1
      ia=0
      do 10 i=1,4
 10   ia=ia*256 + ichar(a1(i:i))
      return
      end
      subroutine intchr (ia,a1)
c... converts integer*4 to char*4
      character*4 a1
      ival=ia
      do 10 i=4,1,-1
      ich=mod(ival,256)
      a1(i:i)=char(ich)
      ival=int( (ival-ich)/256 )
 10   continue
      return
      end
      SUBROUTINE CLOSE (IUNIT)
C... INTERFACE FOR DLEARN PROGRAM (T.J.R. HUGHES) 
C              TO CLOSE AND KEEP A FILE IF EXIST
C
      CLOSE (UNIT=IUNIT,STATUS='KEEP')
      RETURN
      END
      FUNCTION SECNDS (T)
C... INTERFACE FOR DLEARN PROGRAM (T.J.R. HUGHES) TO GET CPU TIME
C
      CALL SYSCPU (T0)
      SECNDS = T0-T
      RETURN
      END
      SUBROUTINE SYSCPU (T)
C... GET CPU TIME ON IBM PC MACHINES USING MS-FORTRAN
C**************** NOT ACTIVE ************************
      T=0.0
      RETURN
      END

