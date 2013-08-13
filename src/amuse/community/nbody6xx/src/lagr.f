      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common6.h'
      INTEGER NPARTC(11), NCORE
      REAL*8  R2(NMAX),RSNGL(NMAX),RBIN(NMAX)
      REAL*8  C(3), FLAGR(11), RLAGR(11), AVMASS(11), AVMRC,VROTC
      REAL*8  RSLAGR(11),RBLAGR(11),SIGR2(11),SIGT2(11),VROT(11)
      INTEGER ISLIST(NMAX),IBLIST(NMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*       Lagrangian radii at 1,2,5,7.5,10,15,20,50,70,90,100 % of total mass.
      DATA FLAGR/0.01,0.02,0.05,0.10,0.20,0.30,0.40,0.50,0.70,0.90,1.00/
*  correspond to:  1%, 2% , 5% , 10%, 20%, 30%, 40%, 50%, 70%, 90%, 100%
*
*       All quantities refer to initial mass/particle distribution
      IF(FIRST) ZMASS0 = ZMASS
*
*       Set square radii of all single particles & c.m. bodies.
      NP = 0
      DO 10 I = 1,N
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)
*
*       Set square radii of singles and binaries for primordial binaries
      IF (KZ(8).GT.0 .OR. NBIN0.GT.0) THEN
*
          IF (FIRST) THEN
          ZMB0 = 0.0D0
          ZMS0 = 0.0D0
              DO 101 I = IFIRST,N
 101          ZMS0 = ZMS0 + BODY(I)
              DO 102 I = N+1,NTOT
 102          ZMB0 = ZMB0 + BODY(I)
          END IF
*
          NSNGL = 0
          DO 103 I = IFIRST,N
              NSNGL = NSNGL + 1
              RSNGL(NSNGL) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                            (X(3,I) - C(3))**2
              ISLIST(NSNGL) = I
 103      CONTINUE
*       Sort square distances with respect to the centre C.
      CALL SORT1(NSNGL,RSNGL,ISLIST)
*
          NBIN  = 0
          DO 104 I = N+1, NTOT
              NBIN = NBIN + 1
              RBIN(NBIN) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                          (X(3,I) - C(3))**2
              IBLIST(NBIN) = I
 104      CONTINUE
*       Sort square distances with respect to the centre C.
      CALL SORT1(NBIN,RBIN,IBLIST)
*
      END IF
*
*  Determine the Lagrangian radii for specified mass fractions.
*     RLAGR = Lagrangian radius
*     AVMASS = average mass of a spherical shell with radius R2(I)
*     AVMRC = average mass inside of core radius RC
*     NPARTC = particle counter within a shell
*     NCORE = particle counter for the core
*     RC = Core radius (calculated in core.f)
      SIGR2C = 0.0D0
      SIGT2C = 0.0D0
      VROTC = 0.0D0
      AVMRC = 0.D0
      RC2 = RC*RC
      NCORE = 0
      ZM = 0.0D0
      ZMH = 0.5D0*ZMASS0
      I = 0
*
      DO 15 J = 1,11
         SIGR2(J) = 0.0D0
         SIGT2(J) = 0.0D0
         VROT(J) = 0.0D0
         AVMASS(J) = 0.0D0
         NPARTC(J) = 0
*
 20      I = I + 1
         IM = JLIST(I)
         ZM = ZM + BODY(IM)
*
         VRI = 0.D0
         DO 25 K = 1,3
 25      VRI = VRI + XDOT(K,IM)*(X(K,IM)-C(K))/DSQRT(R2(I))
         VR2I = 0.D0
         VT2I = 0.D0
         DO 26 K = 1,3
         VRADK = VRI*(X(K,IM)-C(K))/DSQRT(R2(I))
         VR2I = VR2I + VRADK**2
 26      VT2I = VT2I + (XDOT(K,IM) - VRADK)**2
*
         SIGR2(J) = SIGR2(J) + BODY(IM)*VR2I
         SIGT2(J) = SIGT2(J) + BODY(IM)*VT2I/2.D0
*
         RR12 = X(1,IM)**2 + X(2,IM)**2
         XR12 = 0.D0
         DO 27 K = 1,2
 27      XR12 = XR12 + XDOT(K,IM)*(X(K,IM)-C(K))
*
         VROT1 = XDOT(1,IM) - XR12/RR12*X(1,IM)
         VROT2 = XDOT(2,IM) - XR12/RR12*X(2,IM)
*
         XSIGN = VROT1*X(2,IM)/DSQRT(RR12) - VROT2*X(1,IM)/DSQRT(RR12)
         VROTM = DSQRT(VROT1**2+VROT2**2)
*
         IF(XSIGN.GT.0.D0) THEN
         VROT(J) = VROT(J) + BODY(IM)*VROTM
         ELSE
         VROT(J) = VROT(J) - BODY(IM)*VROTM
         END IF
*
         AVMASS(J) = AVMASS(J) + BODY(IM)
         NPARTC(J) = NPARTC(J) + 1
*
*        Independent determination of mass in core radius.
         IF (R2(I) .LT. RC2) THEN
            SIGR2C = SIGR2C + BODY(IM)*VR2I
            SIGT2C = SIGT2C + BODY(IM)*VT2I/2.D0
            VROTC = VROTC + BODY(IM)*VROTM
            AVMRC = AVMRC + BODY(IM)
            NCORE = NCORE + 1
         END IF
*
*        Check whether mass within Langrangian radius is complete.
         IF (I.LT.N.AND.ZM.LT.FLAGR(J)*ZMASS0) GOTO 20
*
         RLAGR(J) = SQRT(R2(I))
         SIGR2(J) = SIGR2(J)/AVMASS(J)
         SIGT2(J) = SIGT2(J)/AVMASS(J)
         VROT(J) = VROT(J)/AVMASS(J)
         AVMASS(J) = AVMASS(J)/NPARTC(J)
 15   CONTINUE
      SIGR2C = SIGR2C/AVMRC
      SIGT2C = SIGT2C/AVMRC
      VROTC = VROTC/AVMRC
      AVMRC = AVMRC/NCORE
*
*  Determine half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5D0*ZMASS0
      I = 0
 30   I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 30
*  Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*
      IF (KZ(8).GT.0 .OR. NBIN0.GT.0) THEN
*  Determine the Lagrangian radii for singles and binaries separately.
      ZMS = 0.0D0
      I = 0
*
      DO 251 J = 1, 11
 201     I = I + 1
         IM = ISLIST(I)
         ZMS = ZMS + BODY(IM)
         IF (I.LT.N-2*NPAIRS.AND.ZMS.LT.FLAGR(J)*ZMS0) GOTO 201
         RSLAGR(J) = SQRT(RSNGL(I))
 251  CONTINUE
*
      ZMB = 0.0D0
      I = 0
*
      DO 351 J = 1, 11
 301     I = I + 1
         IM = IBLIST(I)
         ZMB = ZMB + BODY(IM)
         IF (I.LT.NPAIRS.AND.ZMB.LT.FLAGR(J)*ZMB0) GOTO 301
         RBLAGR(J) = SQRT(RBIN(I))
 351  CONTINUE
*
      END IF
*
      if(rank.eq.0)then
*        Write on diagnostics (only until r_h).
         IF (KZ(7).EQ.2 .OR. KZ(7).EQ.4) THEN
          WRITE (6,40) ' TIME ',(FLAGR(K),K=1,11),' <RC'
   40     FORMAT (/,3X,A6,' M/MT:  ',1P,11(1X,D9.2),A4)
          WRITE (6,41) TTOT, (RLAGR(K),K=1,11),RC
   41     FORMAT (3X,F6.1,' RLAGR: ',1P,12(1X,D9.2))
*
         IF (KZ(8).GT.0 .OR. NBIN0.GT.0) THEN
          WRITE (6,401) TTOT, (RSLAGR(K),K=1,11)
  401     FORMAT (3X,F6.1,' RSLAGR: ',1P,11(1X,D9.2))
          WRITE (6,402) TTOT, (RBLAGR(K),K=1,11)
  402     FORMAT (3X,F6.1,' RBLAGR: ',1P,11(1X,D9.2))
         END IF
*
          WRITE (6,42) TTOT, (AVMASS(K),K=1,11),AVMRC
   42     FORMAT (3X,F6.1,' AVMASS:',1P,12(1X,D9.2))
          WRITE (6,43) TTOT, (NPARTC(K),K=1,11),NCORE
   43     FORMAT (3X,F6.1,' NPARTC:',12I10)
          WRITE (6,45) TTOT, (SIGR2(K),K=1,11),SIGR2C
   45     FORMAT (3X,F6.1,' SIGR2: ',1P,12(1X,D9.2))
          WRITE (6,46) TTOT, (SIGT2(K),K=1,11),SIGT2C
   46     FORMAT (3X,F6.1,' SIGT2: ',1P,12(1X,D9.2))
          WRITE (6,47) TTOT, (VROT(K),K=1,11),VROTC
   47     FORMAT (3X,F6.1,' VROT:  ',1P,12(1X,D9.2))


         END IF
      end if
*        Write all data in binary format on unit 7.
*
      if(rank.eq.0)then
         IF (KZ(7).GE.3) THEN
            WRITE (7) TTOT, (RLAGR(K),K=1,11),
     &         (AVMASS(K),K=1,11),(NPARTC(K),K=1,11),RC,AVMRC,NCORE
         CALL FLUSH(7)
         END IF
      end if
*
      RETURN
*
      END


