      SUBROUTINE ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
*
*
*       Assessment of hierarchical stability.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  XX(3,3),VV(3,3)
      SAVE ITIME
      DATA ITIME /0/
*
*
*       Define indices for components and c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
*
*       Determine inclination between inner relative motion and outer orbit.
      RV = 0.0
      DO 10 K = 1,3
          XX(K,1) = XREL(K,IM)
          XX(K,2) = 0.0
          XX(K,3) = X(K,I2)
          VV(K,1) = VREL(K,IM)
          VV(K,2) = 0.0
          VV(K,3) = XDOT(K,I2)
          RV = RV + XREL(K,IM)*VREL(K,IM)
   10 CONTINUE
*
      CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ANGLE)
*
*       Form inner eccentricity (neglect radial velocity for minimum).
      RIN = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
      SEMI0 = -0.5*BODY(I1)/HM(IM)
      ECC2 = (1.0 - RIN/SEMI0)**2 + RV**2/(BODY(I1)*SEMI0)
      ECC0 = SQRT(ECC2)
      PMIN = SEMI*(1.0 - ECC)
*
*       Evaluate the general stability function.
      IF (ECC.LT.1.0) THEN
          EOUT = ECC
*       Modify outer eccentricity for consistency with acceptance.
          IF (EOUT.GT.0.80) THEN
              DE = 0.5*(1.0 - EOUT)
              DE = MIN(DE,0.01D0)
              EOUT = ECC - DE
              PMIN = SEMI*(1.0 - EOUT)
          END IF
          NST = NSTAB(SEMI0,SEMI,ECC0,EOUT,ANGLE,CM(1,IM),
     &                                    CM(2,IM),BODY(I2))
          IF (NST.EQ.0) THEN
              PCRIT = 0.98*PMIN
              PCR = stability(CM(1,IM),CM(2,IM),BODY(I2),
     &                                          ECC0,ECC,ANGLE)*SEMI0
*        Reduce termination distance if old criterion < PMIN/2.
              IF (PCR.LT.0.5*PMIN) THEN
                  PCRIT = 0.75*PMIN
              END IF
              ITIME = ITIME + 1
              IF (ITIME.GT.2000000000) ITIME = 0
              IF (MOD(ITIME,1000).EQ.0) THEN
                  ALPH = 360.0*ANGLE/TWOPI
                  WRITE (6,20)  ECC0, ECC, ALPH, SEMI, PCRIT, PCR,
     &                          R0(IPAIR)
   20             FORMAT (' ASSESS    E0 E1 INC A1 PCR PC0 R0 ',
     &                                2F7.3,F7.1,1P,4E9.1)
              END IF
          ELSE
              PCRIT = 1.01*PMIN
          END IF
      ELSE
          PCRIT = stability(CM(1,IM),CM(2,IM),BODY(I2),
     &                                        ECC0,ECC,ANGLE)*SEMI0
      END IF
*
*       Set new stability distance or define termination.
      IF (PCRIT.LT.PMIN) THEN
          ITERM = 0
          R0(IPAIR) = PCRIT
      ELSE
          ITERM = 1
      END IF
*
      RETURN
*
      END
