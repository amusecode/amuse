      SUBROUTINE output(time)
C=========================================================================
C===
C===  task: write output files
C===  -----
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER    i,k
      DOUBLE PRECISION time
      CHARACTER*79     bofile

C--------------------------------------------------------------------
C---  output of body data
C--------------------------------------------------------------------

C---  get filename for current step
      IF (lbod) THEN
         WRITE(bofile,500) diskstep
 500     FORMAT(I5.5,'.dat')

         OPEN(UNIT=ubo,FILE=bofile)

         WRITE(ubo,'(I10)')      diskstep
         WRITE(ubo,'(I10)')      N
         WRITE(ubo,'(1P,E15.6)') time

         DO i=1,N
            WRITE(ubo,510) ind(i)-1,m(i),(x(k,i),k=1,3),(v(k,i),k=1,3)
 510        FORMAT(I10.10,1P,E15.6,2X,3E14.6,2X,3E14.6)
         ENDDO

         CLOSE(ubo)
      ENDIF

C---------------------------------------------------------------------
C---  output to tot file
C---------------------------------------------------------------------
      IF (ltot) THEN              
         WRITE(utot,520) time,Timesteps,E_tot,E_pot,E_kin,xcm,vcm,
     $                   mom
 520     FORMAT(1P,E15.6,I8,12E15.6)
         CALL flush(utot)
      ENDIF


      RETURN
      END
