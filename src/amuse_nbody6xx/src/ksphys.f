      SUBROUTINE KSPHYS(Q,P,XR,PR)
*
*
*       Transformation from KS to physical variables.
*       ---------------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Q(4),P(4),XR(3),PR(3)
*
*
*       Obtain physical coordinates & momenta from KS.
      XR(1)=Q(1)**2-Q(2)**2-Q(3)**2+Q(4)**2
      XR(2)=2.D0*(Q(1)*Q(2)-Q(3)*Q(4))
      XR(3)=2.D0*(Q(1)*Q(3)+Q(2)*Q(4))
*
      R=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2
      A=0.5D0/R
      PR(1)=(Q(1)*P(1)-Q(2)*P(2)-Q(3)*P(3)+Q(4)*P(4))*A
      PR(2)=(Q(2)*P(1)+Q(1)*P(2)-Q(4)*P(3)-Q(3)*P(4))*A
      PR(3)=(Q(3)*P(1)+Q(4)*P(2)+Q(1)*P(3)+Q(2)*P(4))*A
*
      RETURN
      END
