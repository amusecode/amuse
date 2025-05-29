      SUBROUTINE CLINT
*
*
*       Cloud integration.
*       ------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
      REAL*8  FC(3)
*
*
*       Check the time for advancing clouds.
      IF (TCL + STEPCL.GT.TIME + TOFF) GO TO 30 
*
*       Set time-step and update the cloud reference time.
      DT = TIME + TOFF - TCL
      DT1 = 0.5D0*DT
      TCL = TIME + TOFF
*
*       Integrate cloud orbits in rotating coordinates with tidal effects.
      DO 20 ICL = 1,NCL
          FC(1) = XCL(1,ICL)*TIDAL(1) + XDOTCL(2,ICL)*TIDAL(4)
          FC(2) = -XDOTCL(1,ICL)*TIDAL(4)
          FC(3) = XCL(3,ICL)*TIDAL(3)
          XCL2 = 0.0
*
          DO 10 K = 1,3
              XCL(K,ICL) = (FC(K)*DT1 + XDOTCL(K,ICL))*DT + XCL(K,ICL)
              XDOTCL(K,ICL) = FC(K)*DT + XDOTCL(K,ICL)
              XCL2 = XCL2 + (XCL(K,ICL) - RDENS(K))**2
   10     CONTINUE
*
*       Determine the minimum impact parameter (squared).
          PCL2 = MIN(XCL2,PCL2)
*
*       Check for modification of the cloud mass (increase or reduction).
          IF (XCL2.LT.RB2) THEN
              IF (BODYCL(ICL).LT.CLM(ICL)) THEN
                  BODYCL(ICL) = BODYCL(ICL) + CLMDOT(ICL)*DT
              END IF
          ELSE
*       Reduce cloud mass gradually by 'sun-set' procedure.
              BODYCL(ICL) = BODYCL(ICL) - CLMDOT(ICL)*DT
          END IF
*
*       Initialize a new cloud when current mass becomes negative.
          IF (BODYCL(ICL).LE.0.0) CALL CLOUD(ICL)
   20 CONTINUE
*
   30 RETURN
*
      END
