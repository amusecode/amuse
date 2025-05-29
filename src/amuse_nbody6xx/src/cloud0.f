      SUBROUTINE CLOUD0
*
*
*       Cloud initialization.
*       ---------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
*
*
*       Initialize cloud variables.
      NCL = 0
      TCL = 0.0D0
      NEWCL = 0
      PCL2 = 0.0
*
*       Read the cloud parameters (pc, km/sec, Msun & pc).
      if (amusein.eq.0) READ (5,*)  NCL, RB2, VCL, SIGMA,
     &            (CLM(J),J=1,NCL), (RCL2(J),J=1,NCL)
      if(rank.eq.0)then
      WRITE (6,1)  NCL, RB2, VCL, SIGMA
    1 FORMAT (/,12X,'CLOUDS:    NCL =',I4,'  RB =',F5.1,
     &                      '   MEAN VELOC =',F5.1,'  DISP =',F5.1)
      WRITE (6,2)  (CLM(J),J=1,NCL)
    2 FORMAT (/,12X,'CLOUD MASSES:   ',1P,10E9.1)
      WRITE (6,3)  (RCL2(J),J=1,NCL)
    3 FORMAT (/,12X,'PLUMMER RADII:  ',1P,10E9.1)
      end if
      RBAR1 = RBAR
      IF (RBAR.EQ.0.0) RBAR1 = 1.0
*       Set cloud parameters in scaled units.
      RB2 = RB2/RBAR1
      A1 = 0.047*SQRT(ZMASS*ZMBAR/RBAR1)
*       Rms velocity of cluster members in km/sec.
      A2 = A1/SQRT(0.5D0*ZMASS)
*       Velocity unit.
      VCL = VCL/A2
*       Cloud velocity in scaled units.
      SIGMA = SIGMA/A2
*       Specify conservative cloud integration step using crossing time.
      STEPCL = 0.002*TCR*RB2/VCL
*
*       Adopt a quantized value.
      CALL STEPK(STEPCL,DTN)
      STEPCL = DTN
*
*       Scale radii & masses to model units.
      DO 10 J = 1,NCL
          RCL2(J) = RCL2(J)/RBAR1
          CLM(J) = CLM(J)/ZMBAR
   10 CONTINUE
*
      if(rank.eq.0)then
      WRITE (6,15)  RB2, VCL, SIGMA, STEPCL
   15 FORMAT (/,12X,'SCALING:     RB =',F6.1,'  VCL =',F5.1,
     &                        '  SIGMA =',F5.1,'  STEP =',1P,E10.2)
      WRITE (6,20)  (CLM(J),J=1,NCL)
   20 FORMAT (/,12X,'SCALED MASSES:  ',10F7.2)
      WRITE (6,25)  (RCL2(J),J=1,NCL)
   25 FORMAT (/,12X,'SCALED RADII:   ',10F7.2)
      end if
      CLDOT = 0.1*RB2/VCL
*       Time scale for 'sun-rise' is 0.05 of the cloud crossing time.
      CLDOT = 1.0/CLDOT
*
*       Define the square of cloud half-mass radii & growth times.
      DO 30 J = 1,NCL
          RCL2(J) = RCL2(J)**2
          CLMDOT(J) = CLM(J)*CLDOT
   30 CONTINUE
*
*       Set square boundary radius & impact parameter.
      RB2 = RB2**2
      PCL2 = RB2
*       Define density centre for routine CLOUD.
      DO 40 K = 1,3
          RDENS(K) = 0.0
   40 CONTINUE
*
*       Initialize new clouds on the boundary.
      DO 50 ICL = 1,NCL
          CALL CLOUD(ICL)
   50 CONTINUE
*
      RETURN
*
      END
