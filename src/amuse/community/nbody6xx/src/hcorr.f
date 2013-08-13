      SUBROUTINE HCORR(I,DM,RNEW)
*
*
*       Mass loss correction of KS orbit.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
*
*
*       Copy pair index and determine secondary component (note I1 = I here).
      IPAIR = KSPAIR
      I1 = I
      I2 = 2*IPAIR - 1
      IF (I2.EQ.I1) THEN
          I2 = I2 + 1
      END IF
*
*       Define c.m. index and old semi-major axis.
      J = N + IPAIR
      SEMI0 = -0.5D0*BODY(J)/H(IPAIR)
*
*       Distinguish different cases (KS, triple or outer ghost of quadruple).
      IF (BODY(I).GT.0.0D0) THEN
          ZMU0 = BODY(I1)*BODY(I2)/BODY(J)
      ELSE
          IM = 0
          JJ = I
          IF (I.LT.IFIRST) JJ = N + KVEC(I)
*       Determine merger index for inner and outer binary masses.
          DO 1 K = 1,NMERGE
              IF (NAMEG(K).EQ.NAME(JJ)) IM = K
    1     CONTINUE
          IF (IM.EQ.0) THEN
              if(rank.eq.0)
     &        WRITE (6,5) NAME(I), NAME(JJ), (NAMEG(K),K=1,NMERGE)
    5         FORMAT (' DANGER HCORR!    NMI NMJJ NMG ',12I6)
              STOP
          END IF
*       Form old reduced mass from inner and outer binary.
          ZM1 = CM(1,IM) + CM(2,IM)
          ZM2 = CM(3,IM) + CM(4,IM)
          ZMU0 = ZM1*ZM2/BODY(J)
      END IF
*
*       Obtain energy change from MK*A = const and H = -M/(2*A) (DCH 8/96).
      DH = DM/SEMI0*(1.0 - 0.5*DM/BODY(J))
*
*       Reduce mass of c.m. body and subtract energy correction terms.
      BODY(J) = BODY(J) - DM
      IF (BODY(I).GT.0.0D0) THEN
          ZMU1 = (BODY(I1) - DM)*BODY(I2)/BODY(J)
      ELSE
          ZMU1 = ZM1*(ZM2 - DM)/BODY(J)
      END IF
      EMDOT = EMDOT - ZMU1*DH - (ZMU1 - ZMU0)*H(IPAIR)
*       Retain binding energy gained in past interactions.
      EGRAV = EGRAV - ZMU1*DH - (ZMU1 - ZMU0)*H(IPAIR)
*
*       Include diagnostic warning for large relative energy change.
      IF (rank.eq.0.and.DH.GT.0.2*ABS(H(IPAIR))) THEN
          WRITE (6,10)  NAME(I), DH, H(IPAIR), R(IPAIR)/SEMI0, DM*SMU
   10     FORMAT (' WARNING!    LARGE CORRECTION    NM DH H R/A DMS ',
     &                                              I6,2F8.2,2F7.2)
      END IF
*
*       Update the binding energy due to mass loss DM and set new radius.
      H(IPAIR) = H(IPAIR) + DH
      RADIUS(I) = RNEW
*
*       Modify KS variables at constant eccentricity.
      CALL EXPAND(IPAIR,SEMI0)
*
      RETURN
*
      END
