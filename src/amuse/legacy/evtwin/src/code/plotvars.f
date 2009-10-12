      MODULE PLOTVARIABLES
C     This module can be approached by FUNCS1, the subroutine were all
C     physics is calculated, and by the various subroutines in PRINTB,
C     were a selected set of physcial variables are written to output
C     files such as file.plt1 etcetera
C     
C     If you (temporarily) need any additional output in the plt or mdl
C     files, which is available in FUNCS1, this is the module you're
C     looking for. 
C
C     To do
C     - make this module a substitute for the commonblock /PLOT/
C
C     Hist:
C     SdM 30jan09
      IMPLICIT NONE

C Variables previously stored in the common block /PLOT/ as created by MvdS
      DOUBLE PRECISION ::WML, MTR, DHDT, DHSP(2), DHMB(2), DHSO(2), DHML,
     &     DHGW, DHMT(2), QCNV, TEFF 


      END MODULE    

