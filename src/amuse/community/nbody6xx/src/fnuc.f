      SUBROUTINE FNUC(XI,XIDOT,FM,FD)
*
*
*       Galaxy point-mass force.
*       ------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
*
*
*       Evaluate force and first derivative from global variables.
      R2 = XI(1)**2 + XI(2)**2 + XI(3)**2
      RRD = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3))/R2
      H3 = GMG/(R2*SQRT(R2))
*
      DO 10 K = 1,3
          FM(K) = -H3*XI(K)
          FD(K) = -H3*(XIDOT(K) - RRD*XI(K))
   10 CONTINUE
*
      RETURN
*
      END
