      SUBROUTINE KICK2(I)
*
*
*       Velocity kick for Roche NS & BH stars.
*       --------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Randomize orbital phase and set disruption velocity in routine KICK.
      KSPAIR = KVEC(I)
      JPAIR = -KSPAIR
      KSTAR(I) = -KSTAR(I)
      CALL KSAPO(JPAIR)
*
*       Update the physical variables for termination.
      CALL RESOLV(KSPAIR,2)
*
*       Specify new single star location and terminate regularization.
      I = I + 2*(NPAIRS - KSPAIR)
      IPHASE = 2
      JCOMP = 0
      CALL KSTERM
*
*       Implement kick velocity for single component (mass loss in ROCHE).
      KW = KSTAR(I)
      CALL KICK(I,1,KW)
*
*       Re-initialize the KS regularization.
      ICOMP = IFIRST
      JCOMP = IFIRST + 1
      IPHASE = 1
      CALL KSREG
*
*       Check neighbour step reduction to compensate for velocity increase.
      NNB = ILIST(1)
      DO 10 L = 2,NNB+1
          J = ILIST(L)
          IF (T0(J) + 0.5*STEP(J).GT.TIME) THEN
              STEP(J) = 0.5*STEP(J)
              TIMENW(J) = T0(J) + STEP(J)
          END IF
   10 CONTINUE
*
*       Specify negative index for new particle sequence.
      IPHASE = -1
*
      RETURN
*
      END
