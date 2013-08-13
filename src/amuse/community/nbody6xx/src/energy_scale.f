      SUBROUTINE ENERGY_SCALE
*
*
*       Total energy for special scaling with initial binaries.
*       -------------------------------------------------------
*       Assume all NBIN0 binaries will be regularised afterwards.
*       For KyotoII it is assumed that the binaries are at the beginning.
*
      INCLUDE 'common6.h'
      REAL*8 GRAV,PARSEC,SUNM
      DATA GRAV/6.6732D-8/
      DATA PARSEC/3.0857D18/
      DATA SUNM/1.989D33/
*
      GRAV = GRAV*SUNM/PARSEC/1.D10*ZMBAR*FLOAT(N)
*       Sum the total energy of pairs.
      NPAIRS = NBIN0
      IFIRST = 2*NPAIRS + 1
      NTOT = N + NPAIRS
*
      if(rank.eq.0)
     *   PRINT*,' Start Energy Scale with N,NTOT,NPAIRS,IFIRST=',
     *   N,NTOT,NPAIRS,IFIRST,GRAV
*
      XMBTOT = 0.D0
      IF (NBIN0.GT.0) THEN
          EBIN = 0.0D0
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
              RIJ2 = 0.0D0
              RDOT = 0.0D0
              VIJ2 = 0.0D0
              J = N + IPAIR
              XMB = BODY(ICOMP) + BODY(JCOMP)
              BODY(J) = XMB
              XMBTOT = XMBTOT + XMB
              DO 75 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
                  RDOT = RDOT + (X(K,ICOMP) - X(K,JCOMP))*
     &           (XDOT(K,ICOMP) - XDOT(K,JCOMP))
                  VIJ2 = VIJ2 + (XDOT(K,ICOMP) - XDOT(K,JCOMP))**2
                  X(K,J) = (BODY(ICOMP)*X(K,ICOMP) + 
     &                     BODY(JCOMP)*X(K,JCOMP))/XMB 
                  XDOT(K,J) = (BODY(ICOMP)*XDOT(K,ICOMP) + 
     &                        BODY(JCOMP)*XDOT(K,JCOMP))/XMB
   75         CONTINUE
              RIJ = DSQRT(RIJ2)
*
              SEMI = 1.D0/(2.0D0/RIJ - VIJ2/(GRAV*XMB))
              ECC2 = (1.0D0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*GRAV*XMB)
              ECC = DSQRT(ECC2)
              EBIN = EBIN + 0.5D0*GRAV*BODY(ICOMP)*BODY(JCOMP)/SEMI
*      if(rank.eq.0)WRITE(67,111)IPAIR,SEMI,ECC,XMB,(X(K,J),K=1,3),
*    &    (XDOT(K,J),K=1,3)
*       IF(RIJ.LT.RMIN.and.rank.eq.0)
*    &     WRITE(68,112)IPAIR,RIJ,RDOT,DSQRT(VIJ2),
*    &     GRAV*XMB/SEMI,(1.0D0 - RIJ/SEMI)**2,
*    &     RDOT**2/(SEMI*GRAV*XMB)
   50     CONTINUE
*111  FORMAT(1X,I5,1P,12D15.5)
*112  FORMAT(1X,I5,1P,12D15.5)
      END IF
*
*       Calculate the potential energy.
      ZKIN = 0.D00
      POT = 0.0
      GRAV = 1.D0
*
      DO 20 I = 1,N-1
      JMIN = I + 1
      IF (I.LE.2*NPAIRS) THEN
*       Binding energy of regularized pairs is included explicitly above.
          IPAIR = KVEC(I)
          JMIN = 2*IPAIR + 1
      END IF
*
      IPAIR = 0
      POTJ = 0.D00
*
      DO 30 J = JMIN,N
      IF (J.EQ.I .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
     *    BODY(J).EQ.0.0D0 .OR. BODY(I).EQ.0.0D0) GO TO 30
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          A4 = GRAV*BODY(J)/DSQRT (A1*A1 + A2*A2 + A3*A3)
      POTJ = POTJ + A4
*
   30 CONTINUE
      POT = POT + BODY(I)*POTJ
   20 CONTINUE
*
*       Sum the kinetic energy (include c.m. bodies but not components).
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
*
      ZKIN = 0.5D0*ZKIN
*
      if(rank.eq.0)then
      PRINT*,' Energies Absolute Scale ',ZKIN,POT,EBIN
      PRINT*,' Energies per mass ',ZKIN/ZMASS,POT/ZMASS,EBIN/ZMASS
      PRINT*,' Total Mass in Binaries ',XMBTOT
      end if
*       Obtain the tidal potential if external field is present.
      ETIDE = 0.0D0
      IF (KZ(14).GT.0) THEN
          CALL XTRNLV(1,N)
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Total energy = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL.
*
      EBIN = 0
      NPAIRS = 0
      IFIRST = 1
      NTOT = N
*
      RETURN
*
      END
