      REAL*8 FUNCTION PFAC(A,Z)

*       Precession factor for hierarchy.
*       --------------------------------
*
*       Maths by Douglas Heggie; function by Rosemary Mardling (11/96).
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      REAL*8  m,m1,A,Z,C,C1,C2,C3,D1,D2,D3,Z1

*The function is singular for m1=0, but the following holds:
*For 0 < m1 <=1 (m1 never > 1),  \infty > PFAC(m1) > 1.6, with P(0.01)=3.6.
*So taking it as 1 isn't too bad for an order of magnitude estimate. 

 
      C1 = 0.5*(Z + A**2)
      Z1 = 5.0 - Z + 4*A**2
      C2 = (Z1 + SQRT(Z1))/6.0
      C3 = (Z1 - SQRT(Z1))/6.0
      IF (C1.GE.C2) THEN
          D1 = C1
          D2 = C2
          D3 = C3
      ELSE IF (C2.GE.C1.AND.C1.GE.C3) THEN
          D1 = C2
          D2 = C1
          D3 = C3
      ELSE
          D1 = C2
          D2 = C3
          D3 = C1
      END IF
*
      C = 1.0/SQRT(D1 - D3)
      m = (D2 - D3)/(D1 - D3)
*
*       Evaluate elliptic integral by approximate expression.
      a0=1.3862944
      a1=0.1119723
      a2=0.0725296
      b0=0.5
      b1=0.1213478
      b2=0.0288729
      m1 = 1.0 - m
      m1 = MAX(m1,1.0D-10)

      PFAC = (a0+a1*m1+a2*m1**2)+(b0+b1*m1+b2*m1**2)*log(1.0/m1)
*
*       Include the square root of (D1 - D3).
      PFAC = C*PFAC

      END
