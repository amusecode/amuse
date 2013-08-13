      SUBROUTINE FDISK(XI,XIDOT,FM,FD)
*
*
*       Miyamoto disk force.
*       --------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
*
*
*       Obtain force & derivative for global variables XI & XIDOT.
      R2 = XI(1)**2 + XI(2)**2
      BZ = SQRT(B**2 + XI(3)**2)
      AZ = SQRT(R2 + (A + BZ)**2)
      AZ3 = DISK/AZ**3
*       Note missing square root sign in (b^2 + z^2) of Book eq. (8.52).
      AZDOT = XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                         (A + BZ)*XI(3)*XIDOT(3)/BZ
      FM(1) = -AZ3*XI(1)
      FM(2) = -AZ3*XI(2)
      FM(3) = -AZ3*XI(3)*(A + BZ)/BZ
*     RDOT = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2))/SQRT(R2)
      FD(1) = -AZ3*(XIDOT(1) - 3.0*AZDOT*XI(1)/AZ**2)
      FD(2) = -AZ3*(XIDOT(2) - 3.0*AZDOT*XI(2)/AZ**2)
*       Note wrong sign in first term of eq. (8.52) (see printer errata).
      Y1 = 3.0*(A + BZ)*XI(3)*AZDOT/(AZ**2*BZ)
*       Note printer mistake of dZ/dt in denominator (see errata).
      Y2 = (A*B**2 + BZ**3)*XIDOT(3)/BZ**3
      FD(3) = AZ3*(Y1 - Y2)
*       Refer to Miyamoto & Nagai (PASJ 27, 533) and Book eq. (8.52).
*
      RETURN
*
      END
