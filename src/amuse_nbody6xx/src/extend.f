      SUBROUTINE EXTEND(ISUB)
*
*
*       Size of (un)perturbed subsystem.
*       --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Set global index for chain c.m. body.
      NCM = NAMES(1,ISUB)
      IF (ISYS(ISUB).GE.3) THEN
          ICM = ICH
          IF (NAME(ICH).EQ.0) GO TO 20
          GO TO 12
      END IF
*
*       Find global index for c.m. body (unperturbed triple or quad).
      DO 10 I = IFIRST,N
          IF (NCM.EQ.NAME(I)) THEN
              ICM = I
              GO TO 20
          END IF
   10 CONTINUE
*
*       Include safety procedure in case c.m. body not identified.
   12 if(rank.eq.0) WRITE (6,15)  ISUB, NCM
   15 FORMAT (5X,'WARNING!   SUBSYSTEM TERMINATION   ISUB =',I3,
     &                                                     '  NCM =',I5)
      STEPS(ISUB) = 0.0D0
      if(rank.eq.0) WRITE (6,60)  ISUB, NCH, ICH, NAME(ICH)
   60 FORMAT (' EXTEND:   ISUB NCH ICH NAME(ICH) ',4I5)
      GO TO 40
*
*       Determine the largest perturbing force (M/R**3).
   20 PMAX = 0.0
      NNB = LIST(1,ICM)
      DO 30 L = 2,NNB+1
          J = LIST(L,ICM)
          RIJ2 = (X(1,J) - X(1,ICM))**2 + (X(2,J) - X(2,ICM))**2 +
     &                                    (X(3,J) - X(3,ICM))**2
          PIJ = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PIJ.GT.PMAX) THEN
              PMAX = PIJ
          END IF
   30 CONTINUE
*
*       Choose maximum of dominant perturber and all neighbours at boundary.
      PMAX = MAX(PMAX,FLOAT(NNB)*BODYM/RS(ICM)**3)
*
*       Specify maximum size of (un)perturbed system (limit of 100*GSTAR).
      IF (ISYS(ISUB).LE.2) THEN
          GSTAR = GMIN
      ELSE
          GSTAR = 0.01*GMAX
      END IF
      RMAXS(ISUB) = (100.0*GSTAR*BODY(ICM)/(2.0*PMAX))**0.3333
*
*       Update time limit interval unless termination has been signalled.
      IF (STEPS(ISUB).GT.0.0D0) THEN
          STEPS(ISUB) = STEP(ICM)
      ELSE
*       Set phase indicator for termination (step reduction in STEPS).
          IF (ISYS(ISUB).LE.2) IPHASE = 4
      END IF
*
*     if(rank.eq.0)
*    &WRITE (6,35)  ICM,NAME(ICM),RMAXS(ISUB),STEPS(ISUB),STEP(ICM)
*  35 FORMAT (' EXTEND:  ICM NM RMAX STEPS STEP ',2I5,1P,3E10.2)
   40 RETURN
*
      END
