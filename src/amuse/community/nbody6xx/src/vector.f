      SUBROUTINE VECTOR(Q,P,W)
*
*
*       KS matrix of Q times P.
*       -----------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Q(4),P(4),W(4)
*
*
*       Evaluate W = L(Q)P; the KS matrix of Q times P.
      W(1)=(+Q(1)*P(1)+Q(2)*P(2)+Q(3)*P(3)+Q(4)*P(4))
      W(2)=(-Q(2)*P(1)+Q(1)*P(2)+Q(4)*P(3)-Q(3)*P(4))
      W(3)=(-Q(3)*P(1)-Q(4)*P(2)+Q(1)*P(3)+Q(2)*P(4))
      W(4)=(+Q(4)*P(1)-Q(3)*P(2)+Q(2)*P(3)-Q(1)*P(4))
*
      RETURN
      END
