      SUBROUTINE STATUS(X,I1,I2,I3,I4)
*
*
*       Current configuration.
*       ----------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  X(12),RR(6)
      INTEGER  II(4,4),IN(6),JN(6),IJ(6)
      LOGICAL  SWITCH
      COMMON/IND6/  IND(6)
      COMMON/CONFIG/  R2(4,4),J1,J2,J3,J4
      DATA  II/0,3,2,2,4,0,1,1,4,4,0,1,3,3,2,0/
*
*
*       Form mutual square distances and determine closest pair indices.
      R2MIN = 1.0D+10
      L = 0
      DO 10 I = 1,3
          DO 5 J = I+1,4
              IL = 3*(I - 1)
              JL = 3*(J - 1)
              R2(I,J) = (X(IL+1) - X(JL+1))**2 + (X(IL+2) - X(JL+2))**2
     &                                         + (X(IL+3) - X(JL+3))**2
              R2(J,I) = R2(I,J)
              L = L + 1
              RR(L) = R2(I,J)
              IN(L) = I
              JN(L) = J
              IF (R2(I,J).LT.R2MIN) THEN
                  R2MIN = R2(I,J)
*       Set closest pair indices in J1 & J2.
                  J1 = I
                  J2 = J
              END IF
    5     CONTINUE
   10 CONTINUE
*
*       Set indices of the two distant particles.
      J3 = II(J1,J2)
      J4 = II(J2,J1)
*       Ensure that body #J3 is closest to #J1 & J2.  
      IF (R2(J1,J3) + R2(J2,J3).GT.R2(J1,J4) + R2(J2,J4)) THEN
          JDUM = J3
          J3 = J4
          J4 = JDUM
      END IF
*
*       Sort square distances.
      CALL RSORT(RR,SWITCH,IND)
      DO 20 K = 3,5
          IJ(1) = IN(IND(1))
          IJ(2) = JN(IND(1))
          IJ(3) = IN(IND(2))
          IJ(4) = JN(IND(2))
          IJ(5) = IN(IND(K))
          IJ(6) = JN(IND(K))
*       Form particle chain for regularization.
          CALL ICHAIN(IJ,KO,I1,I2,I3,I4)
          IF (KO.EQ.1) GO TO 30
   20 CONTINUE
*
   30 RETURN
*
      END
