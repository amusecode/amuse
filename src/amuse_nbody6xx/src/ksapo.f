      SUBROUTINE KSAPO(IPAIR)
*
*
*       Apocentre/pericentre/random KS variables.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RAN2
*
*
*       Specify half regularized period (PI/2) or random phase (unperturbed).
      IF (IPAIR.GT.0) THEN
          THETA = 0.25D0*TWOPI
          IKICK = 0
      ELSE
          THETA = 0.5*TWOPI*RAN2(IDUM1)
          IPAIR = -IPAIR
          IKICK = 1
          IF (LIST(1,2*IPAIR-1).GT.0) THETA = 0.0D0
*       Initialize time because HDOT & TDOT2 not updated for RESOLV.
          T0(2*IPAIR-1) = TIME
*       Skip hyperbolic orbit (i.e. kick for second binary component).
          IF (H(IPAIR).GT.0.0) GO TO 30
      END IF
*
*       Form transformation coefficients (Stiefel & Scheifele p. 85).
      XC = COS(THETA)
      YS = SIN(THETA)
      FF = SQRT(0.5D0*ABS(H(IPAIR)))
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
*
*       Generate analytical solutions for U & UDOT using old U0 & UDOT.
      DO 10 K = 1,4
          U(K,IPAIR) = U0(K,IPAIR)*XC + UDOT(K,IPAIR)*YS/FF
          UDOT(K,IPAIR) = UDOT(K,IPAIR)*XC - U0(K,IPAIR)*YS*FF
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
   10 CONTINUE
*
*       Impose R' < 0 for apocentre procedures (IKICK = 0).
      SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
      IF (TDOT2(IPAIR).GT.0.0D0.AND.R(IPAIR).GT.SEMI) THEN
          IF (IKICK.EQ.0) THEN
              TDOT2(IPAIR) = -1.0E-20
          END IF
      END IF
*
*       Include diagnostic check that correct apocentre has been set.
*     SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
*     if(rank.eq.0)
*    &WRITE (6,20)  SEMI, R(IPAIR), H(IPAIR), GAMMA(IPAIR)
*  20 FORMAT (' APOCENTRE:    A RA H G ',1P,2E10.2,2E10.1)
*
*       Save KS parameters for WD or neutron star kick (routine FCORR).
   30 IF (IKICK.GT.0) THEN
          CALL KICK(IPAIR,0,0)
      END IF
*
      RETURN
*
      END
