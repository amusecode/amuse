      SUBROUTINE HIPOP
*
*
*       Initial hierarchical population.
*       --------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3),
     &        BS(MMAX),XS(3,MMAX),VS(3,MMAX)
      REAL*8  RAN2
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
* recover random number sequence from COMMON variable IDUM1.
      KDUM = IDUM1
*
*       Read input parameters (same usage as routine BINPOP).
      if((rank.eq.0) .AND. (amusein.eq.0)) READ (5,*)  SEMI0, ECC0,
     &      RATIO, RANGE, ICIRC
#if MPIINIT
      CALL MPI_BCAST(SEMI0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ECC0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RATIO,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RANGE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ICIRC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      NHI = NHI0
      WRITE (6,1)  NHI, SEMI0, ECC0, RATIO, RANGE, ICIRC
    1 FORMAT (/,12X,'HIERARCHIES:    NHI =',I4,'  A =',F9.6,
     &              '  E =',F6.2,'  RATIO =',F4.1,'  RANGE =',F6.1,
     &              '  ICIRC =',I2,/)
*
      IF (NHI.GT.NBIN0) THEN
          WRITE (6,2)  NBIN0, NHI
    2     FORMAT (5X,'FATAL ERROR!   NBIN0 =',I4,'  NHI =',I4)
          STOP
      END IF
*
*       Introduce binary components by splitting the primary.
      DO 40 I = 1,NHI
*
*       Randomize perihelion, node & inclination.
          PI = TWOPI*RAN2(KDUM)
          OMEGA = TWOPI*RAN2(KDUM)
          ZI = 0.25*TWOPI*RAN2(KDUM)
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PI)*COS(OMEGA) - SIN(PI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PI)*COS(OMEGA) - COS(PI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PI)*SIN(OMEGA) + SIN(PI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PI)*SIN(OMEGA) + COS(PI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PI)*SIN(ZI)
          QX(3) = COS(PI)*SIN(ZI) 
*
*       Determine two-body elements for original binary.
          I1 = 2*I - 1
          I2 = 2*I
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 5 K = 1,3
              XREL(K) = X(K,I1) - X(K,I2)
              VREL(K) = XDOT(K,I1) - XDOT(K,I2)
              RIJ2 = RIJ2 + XREL(K)**2
              VIJ2 = VIJ2 + VREL(K)**2
              RDOT = RDOT + XREL(K)*VREL(K)
    5     CONTINUE
          RIJ = SQRT(RIJ2)
          ZMB1 = BODY(I1) + BODY(I2)
          A1 = 2.0/RIJ - VIJ2/ZMB1
          SEMI1 = 1.0/A1
          ECC2 = (1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*ZMB1)
          ECC1 = SQRT(ECC2)
          PMIN = SEMI1*(1.0 - ECC1)
*
*       Specify component masses (primary fraction range 0.5 - 0.9).
          Q0 = 0.5 + 0.4*RAN2(KDUM)
          IF (RATIO.EQ.1.0) Q0 = 0.5
          BS(I) = BODY(I2)
          BODY(I1) = Q0*BODY(I1)
          BODY(I2) = BODY(I1)*(1.0 - Q0)/Q0
          ZMB = BODY(I1) + BODY(I2)
*
*       Choose random (thermalized) or fixed eccentricity.
          IF (ECC0.LT.0.0) THEN
              ECC2 = RAN2(KDUM)
              ECC = SQRT(ECC2)
          ELSE
              ECC = ECC0
          END IF
*
*       Select semi-major axis from uniform distribution in log(A) or SEMI0.
          ITER = 0
   10     IF (RANGE.GT.0.0) THEN
              IF (ITER.LE.5) THEN
                  EXP = RAN2(KDUM)*LOG10(RANGE)
                  SEMI = SEMI0/10.0**EXP
              ELSE
*       Shrink by factor 2 if no success after 5 iterations.
                  SEMI = 0.5*SEMI
              END IF
          ELSE
               SEMI = SEMI0
          END IF
*
*       Check stability criterion (maximum 12 tries with inclination effect).
          PCRIT = stability(BODY(I1),BODY(I2),BS(I),ECC,ECC1,ZI)
          PCRIT = PCRIT*SEMI
          ITER = ITER + 1
          IF (PMIN.LT.PCRIT.AND.ITER.LE.12.AND.SEMI.LT.SEMI0) GO TO 10
*
          P0 = DAYS*SEMI*SQRT(SEMI/ZMB)
          P1 = DAYS*SEMI1*SQRT(SEMI1/ZMB1)
      if(rank.eq.0)
     &    WRITE (6,20)  ITER, ECC, ECC1, PMIN, PCRIT, P0, P1
   20     FORMAT (' HIERARCHY:    IT E E1 PMIN PCRIT P0 P1 ',
     &                            I4,2F7.3,1P,4E9.1)
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMI*(1.0 + ECC)
          XORB(2) = 0.0
          VORB(1) = 0.0
          VORB(2) = SQRT(ZMB*(1.0D0 - ECC)/(SEMI*(1.0D0 + ECC)))
          EBIN0 = EBIN0 - 0.5*BODY(I1)*BODY(I2)/SEMI
*
*       Transform to relative variables.
          DO 25 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   25     CONTINUE
*
*       Save old secondary and set global variables for each component.
          DO 30 K = 1,3
              XS(K,I) = X(K,I2)
              VS(K,I) = XDOT(K,I2)
              X(K,I1) = X(K,I1) + BODY(I2)*XREL(K)/ZMB
              X(K,I2) = X(K,I1) - XREL(K)
              XDOT(K,I1) = XDOT(K,I1) + BODY(I2)*VREL(K)/ZMB
              XDOT(K,I2) = XDOT(K,I1) - VREL(K)
   30     CONTINUE
   40 CONTINUE
*
*       Move single particles down by NHI to make room for outer components.
      L = N
      DO 50 J = 2*NBIN0+1,N
          I = NHI + L
          BODY(I) = BODY(L)
          DO 45 K = 1,3
              X(K,I) = X(K,L)
              XDOT(K,I) = XDOT(K,L)
   45     CONTINUE
          L = L - 1
   50 CONTINUE
*
*       Place hierarchical components immediately after the binaries.
      DO 60 L = 1,NHI
          I = 2*NBIN0 + L
          BODY(I) = BS(L)
          DO 55 K = 1,3
              X(K,I) = XS(K,L)
              XDOT(K,I) = VS(K,L)
   55     CONTINUE
   60 CONTINUE
*
*       Update the particle number and reset NHI (might be used elsewhere).
      N = N + NHI
      NZERO = N
      NTOT = N
      NHI = 0
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
      RETURN
*
      END

