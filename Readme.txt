
______________________________________________________________________________
DLEARN : A  LINEAR  STATIC AND DYNAMIC FINITE ELEMENT ANALYSIS PROGRAM
______________________________________________________________________________

Input Data Instructions and examples are formally introduced in the
textbook " The Finite Element Method : linear static and dynamic finite
element analysis ", Prentice Hall, by   T.J.R HUGHES

Last program UPDATE : 10.06.97
Program has been compiled with Fortran Power Station v.4.0

_____________________________________________________________________________

*****************************************************************************
*             THIS IS A VERSION FOR IBM AT AND COMPATIBLES.                 *
*                                                                           *
*  THE EXECUTABLE FILES WERE CREATED USING THE COMPAQ VISUAL FORTRAN 6.6C COMPILER.  *
*****************************************************************************


THIS HIGH DENSITY MINI-FLOPPY DISK CONTAINS TWO SUBDIRECTORIES :
       
	   1.  \DLEARNPC        2.  \SOURCE
	    
Contents :
__________

README                            / GENERAL INSTRUCTIONS
DLEARNPC     <DIR>                / Executable + Data Files of Program DLEARN
SOURCE       <DIR>                / Source   + Object Files of Program DLEARN 
      
-------------------------------------------------------------------------------


****************************
  SUBDIRECTORY:  \PRJ
****************************
Contains project files dlearn.dsw and dlearn.dsp required by COMPAQ VISUAL FORTRAN 6.6C compiler

****************************
  SUBDIRECTORY:  \DLEARNPC
****************************

Contains the Executable files, Run  files and examples of program DLEARN, identified as DLEARN for the original version and DLEARNB
for the version that includes the Brick element.The source files of these
versions can be found in sub-directory \SOURCE .

The following batch files are included to run the executable files.

            RUN.BAT     for  DLEARN.EXE 


RUN procedure : 
    a) To run DLEARN    :        
	          Type   RUN  and the Name of the DATAFILE
			  
Contents :
__________
README
RUN      BAT   / Batch file to run the program  
DLEARN   EXE   / Executable file
TH1      DAT   / Data file for an 8-node brick element (patch test)
TH2      DAT   / Data file for a  brick cantilever beam
WE1      DAT   / Data file for a brick element problem
TU1      DAT   / Data file for truss tower problem
TU2      DAT   / Data file for a cantilever beam (Static case )
TU3      DAT   / Data file for a cantilever beam (Dynamic case)
TU4      DAT   / Data file for an implicit-explicit dynamic case
TU5      DAT   / Data file for 4-node quad element
TU6      DAT   / Data file for 4-node quad element (rank check)

--------------------------------------------------------------------------------
**************************
  SUBDIRECTORY:  \SOURCE            DLEARN SOURCE LIBRARY 
**************************


The following batch files are included to link the object files
			
contents :
___________
READMES
DL       FOR   /MAIN PROGRAM
DL0      FOR   /general routines 
DL1      FOR   /general routines 
DL2      FOR   /general routines 
EL2      FOR   /general routines (modified for BRICK)
SYSR     FOR   /SYSTEM DEPENDENT ROUTINES
TR0      FOR   /truss element
TR1      FOR   /truss element
QD1      FOR   /quad element
QD2      FOR   /quad element
QD3      FOR   /quad element
BR1      FOR   /brick element
BR2      FOR   /brick element

ALL OBJECT FILES OF THE ABOVE SOURCE FILES
--------------------------------------------------------------------------------

THE ABOVE SOURCE FILES CONTAIN THE FOLLOWING SUBROUTINES :
__________________________________________________________________________

  DL.FOR      :      MAIN PROGRAM DLEARN
  
  DL0.FOR     :      TSEQ
                     STATIN
					
  DL1.FOR     :      DRIVER        ADDLHS        ADDRHS        BC 					
                     BTDB          COLHT         COMPBC        CONTM
					 CONTMA        COORD         CORRCT        DCTNRY
					 DIAG          DYNPTS        ECHO          ELEMNT
					 ELCARD        EQSET         FACTOR        FORMLM
					 GENEL         GENEL1        GENELD        GENELI
					 GENFL         GENFL1        GENSH         GENSH1
					 GENSH2        GENSH3
 
  DL2.FOR     :      ICLEAR        IGEN          IMOVE         INPUT
                     INTERP        ITERUP        LFAC          LOAD
					 LOCAL         LOUT          LTIMEF        MATADD
					 MEANSH        MINMAX        MOVE          MPOINT
					 MULTAB        PIVOTS        PREDCT        PRINC
					 PRINTD        PRINTF        PRINTP        PRNTEL
					 PROP2D        PRTDC         PRTDC1        PRTS2D 
					 RCDOT         RSIN          RSOUT         SERROR
					 SETUPD        SMULT         STORED        TIMCOM
					 TIMING        TIMLOG        TSEQUIN       ZTEST
					 
 (EL2.FOR)           HPLOT         SHIST 					  


  SYSR.FOR    :      CHRINT        INTCHR        CLOSE         SECNDS
                     SYSCPU
					 
  TR0         :      TRUSS
  
  TR1         :      TRUST1        TRUST2        TRUST3        TRUST4
                     TRUST5        TRUSB         TRUSHG        TRUSHL
					 TRUSK         TRUSKD        TRUSPR        TRUSPT
					 TRUSTR
					 
  QD1         :      QUADC
  
  QD2         :      QTSK3
  
  QD3         :      QDCT1         QDCT2         QDCT3         QDCT4
                     QDCT5         QDCB          QDCK          QDCKD
					 QDRSF         QDCSHG        QDCSHL        QDCSTR
					 QDCSUF
					 
  BR1         :      BRICKC        BKCT1         BKCT2         BKCT3
                     BKCT4         BKCT5
					 
  BR2         :      BKCB          BKCK          BKCKD         BKCRSF
                     BKCSHG        BKCSHL        BKCSTR        BKCSUF
					 PROP3D        PRTS3D        SUFDET
					 
					 					 				 					 					 


