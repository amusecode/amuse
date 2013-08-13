      SUBROUTINE DELAY(KCH,KS)
*
*
*       Delay of multiple regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SAVEIT/  IPH0,KS0,KCH0,KS20,JC0,JCL0
*
*
*       Check for saving KS pointers from IMPACT or termination at TBLOCK.
      IF (KS.GE.0) THEN
          IPH0 = IPHASE
          KS0 = KSPAIR
          KCH0 = KCH
          KS20 = KS
          JC0 = JCOMP
          JCL0 = JCLOSE
      ELSE
*       Copy relevant indices for delayed KS termination.
          IPHASE = IPH0
          KSPAIR = KS0
          KCHAIN = KCH0
          KS2 = KS20
          JCOMP = JC0
          JCLOSE = JCL0
*      Preserve contents of KSAVE during chain regularization.
          IF (NCH.EQ.0) THEN
              KSAVE(1) = 0
              KSAVE(2) = 0
          END IF
          KSKIP = 0
*
*      Exit in case of new merger or merger termination.
          IF (IPHASE.EQ.6.OR.IPHASE.EQ.7) THEN
              GO TO 10
          END IF
*
*       Include the case of two interacting KS pairs.
          IF (JCOMP.GT.N) THEN
              IF (KCHAIN.GT.0.AND.KSTAR(N+KSPAIR).NE.0) THEN
                  KSAVE(1) = KSTAR(N+KSPAIR)
                  KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
                  KSKIP = 1
              END IF
*       Terminate smallest pair first and copy second pair index.
              CALL KSTERM
              KSPAIR = KS2
*       Specify JCOMP < 0 to prevent spurious prediction second KSTERM call.
              JCOMP = -1
          END IF
*
*       Save KSTAR (> 0) and sum of component names (for chain termination).
      IF (KCHAIN.GT.0.AND.KSTAR(N+KSPAIR).NE.0.AND.KSKIP.EQ.0) THEN
          KSAVE(1) = KSTAR(N+KSPAIR)
          KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
      END IF
*
*       Terminate binary in triple or widest binary-binary collision pair.
          CALL KSTERM
*
*       See whether chain regularization indicator should be switched on.
          IF (KCHAIN.GT.0) THEN
              IPHASE = 8
          END IF
      END IF
*
   10 RETURN
*
      END
