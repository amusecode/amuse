      SUBROUTINE PERMIT(PERIM,IGO)
*
*
*       Check on existing multiple regularization.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
*
*
*       Search any existing subsystem.
      ISUB = 0
      ICHSUB = 1
      DO 10 L = 1,NSUB
*       Identify chain pointer for possible reduction of STEPS.
          IF (ISYS(L).EQ.3) ICHSUB = ISYS(L)
*       Distinguish between triple & quad case (denoted ISUB = 1 or 2).
          IF (JCOMP.LE.N.AND.NAMES(4,L).EQ.0) THEN
              ISUB = 1
      ELSE IF (JCOMP.GT.N.AND.NAMES(4,L).GT.0) THEN
              ISUB = 2
          END IF
   10 CONTINUE
*
*       Do not allow a second regularization of the same type.
      IF (ISUB.GT.0) THEN
*       See whether the case ISUB = 1 or 2 is used already.
          DO 20 L = 1,NSUB
              IF (ISUB.EQ.ISYS(L)) IGO = 1
   20     CONTINUE
*       Enforce chain termination at next extension if new system < RSUM/2.
          IF (PERIM.LT.0.5*RSUM.AND.IGO.GT.0) THEN
              STEPS(ICHSUB) = 0.0
          END IF
      END IF
*
      RETURN
*
      END
