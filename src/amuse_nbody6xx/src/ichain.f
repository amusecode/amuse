      SUBROUTINE ICHAIN(IJ,KO,I1,I2,I3,I4)
*
*
*       Determination of chain vector.
*       ------------------------------
*
      INTEGER  IJ(6)
*
*
      ICON12 = 0
      IF   (IJ(1).EQ.IJ(3).OR.IJ(1).EQ.IJ(4)
     &  .OR.IJ(2).EQ.IJ(3).OR.IJ(2).EQ.IJ(4)) ICON12 = 1
      ICON23 = 0
      IF   (IJ(3).EQ.IJ(5).OR.IJ(3).EQ.IJ(6)
     &  .OR.IJ(4).EQ.IJ(5).OR.IJ(4).EQ.IJ(6)) ICON23 = 1
      ICON14 = 0
      IF   (IJ(1).EQ.IJ(5).OR.IJ(1).EQ.IJ(6)
     &  .OR.IJ(2).EQ.IJ(5).OR.IJ(2).EQ.IJ(6)) ICON14 = 1
*
      IF (ICON12 + ICON23 + ICON14.NE.2) THEN
          KO = 0
          RETURN
      END IF
*
      KO = 1
      IF (ICON12.EQ.0) THEN
          IF (IJ(3).EQ.IJ(5).OR.IJ(3).EQ.IJ(6)) THEN
              I2 = IJ(3)
              I1 = IJ(4)
          ELSE
              I1 = IJ(3)
              I2 = IJ(4)
          END IF
          IF (IJ(1).EQ.IJ(5).OR.IJ(1).EQ.IJ(6)) THEN
              I3 = IJ(1)
              I4 = IJ(2)
          ELSE
              I3 = IJ(2)
              I4 = IJ(1)
          END IF
          RETURN
*
      ELSE IF (ICON23.EQ.0) THEN
          IF (IJ(5).EQ.IJ(1).OR.IJ(5).EQ.IJ(2)) THEN
              I2 = IJ(5)
              I1 = IJ(6)
          ELSE
              I1 = IJ(5)
              I2 = IJ(6)
          END IF
          IF (IJ(3).EQ.IJ(1).OR.IJ(3).EQ.IJ(2)) THEN
              I3 = IJ(3)
              I4 = IJ(4)
          ELSE
              I4 = IJ(3)
              I3 = IJ(4)
          END IF
          RETURN
*
      ELSE IF (ICON14.EQ.0) THEN
          IF (IJ(1).EQ.IJ(3).OR.IJ(1).EQ.IJ(4)) THEN
              I2 = IJ(1)
              I1 = IJ(2)
          ELSE
              I1 = IJ(1)
              I2 = IJ(2)
          END IF
          IF (IJ(5).EQ.IJ(3).OR.IJ(5).EQ.IJ(4)) THEN
              I3 = IJ(5)
              I4 = IJ(6)
          ELSE
              I4 = IJ(5)
              I3 = IJ(6)
          END IF
          RETURN
      END IF
*
      END
