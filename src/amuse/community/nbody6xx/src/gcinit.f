      SUBROUTINE GCINIT
*
*
*       Initialization of 3D cluster orbit.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
*
* 
*     Obtain initial force and first derivative from point-mass galaxy.
      IF (KZ(14).GT.2) THEN
         RIN2 = 1.0/(RG(1)**2 + RG(2)**2 + RG(3)**2)
         RIN3 = RIN2*SQRT(RIN2)
         RGVG = 3.0*(RG(1)*VG(1) + RG(2)*VG(2) + RG(3)*VG(3))*RIN2
*     
         DO 1 K = 1,3
            FG(K) = -GMG*RG(K)*RIN3
            FGD(K) = -GMG*(VG(K) - RGVG*RG(K))*RIN3
    1    CONTINUE
*
      END IF
      TG = TIME + TOFF
*
      RETURN
*
      END
