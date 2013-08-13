      SUBROUTINE MTRACE(J,DM)
*
*
*       Orbit diagnostics for mass loss.
*       --------------------------------
*
      INCLUDE 'common6.h'
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Open unit #14 the first time.
      IF (rank.eq.0.and.FIRST) THEN
          OPEN (UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE='MDOT')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF (NRG.EQ.0) THEN
              WRITE (14,1)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                      NBIN0, NZERO
    1         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I4,'  N0 =',I6,//)
          END IF
*
          WRITE (14,2)
    2     FORMAT ('     TIME  NAME   K*   M    DM     r/Rc   VR',
     &            '     EI        N')
      END IF
*
*       Evaluate potential energy on GRAPE (#I: single star or binary c.m.).
      IF (J.GE.IFIRST) THEN
          I = J
          CALL POTI(I,POTJ)
      ELSE
          KS = KVEC(J)
          I = N + KS
          CALL POTI(I,POTJ)
      END IF
*
*       Obtain central distance (scaled by core radius) and velocities.
      RI2 = 0.0
      VI2 = 0.0
      VR = 0.0
      DO 10 K = 1,3
          RI2 = RI2 + (X(K,I) - RDENS(K))**2
          VI2 = VI2 + XDOT(K,I)**2
          VR = VR + (X(K,I) - RDENS(K))*XDOT(K,I)
   10 CONTINUE
      RI = SQRT(RI2)/RC
      VR = VR/SQRT(RI2)
*
*       Form binding energy per unit mass (note: POTJ < 0).
      EI = 0.5*VI2 + POTJ
*
*       Include optional external tidal field (note: HT includes mass).
      IF (KZ(14).GT.0) THEN
          CALL XTRNLV(I,I)
          EI = EI + HT/(BODY(I) + 1.0E-20)
      END IF
*
      if(rank.eq.0)
     &WRITE (14,20)  TTOT, NAME(J), KSTAR(J), BODY(J)*ZMBAR, DM*ZMBAR,
     &               RI, VR, EI, N
   20 FORMAT (1X,F8.1,I6,I4,F7.1,F6.2,F7.2,F7.2,F8.3,I6)
      CALL FLUSH(14)
*
      RETURN
*
      END
