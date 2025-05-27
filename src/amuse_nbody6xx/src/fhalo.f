      SUBROUTINE FHALO(XI,XIDOT,FM,FD)
*
*
*       Logarithmic halo force.
*       -----------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
*
*
*       Obtain force and first derivative for logarithmic potential.
      R2 = XI(1)**2 + XI(2)**2 + XI(3)**2 + RL2
      RRD = 2.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3))/R2
      H2 = V02/R2
*
*       Note softening in the logarithmic potential (cf. Binney & Tremaine).
      DO 10 K = 1,3
          FM(K) = -H2*XI(K)
          FD(K) = -H2*(XIDOT(K) - RRD*XI(K))
   10 CONTINUE
*
      RETURN
*
      END
