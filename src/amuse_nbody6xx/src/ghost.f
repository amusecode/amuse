      SUBROUTINE GHOST(J)
*
*
*       Initialization of ghost particle.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*
*
*       Form ghost and initialize integration variables for body #J.
      BODY(J) = 0.0D0
      T0(J) = 1.0E+06
      DO 10 K = 1,3
          X0DOT(K,J) = 0.0D0
          XDOT(K,J) = 0.0D0
          F(K,J) = 0.0D0
          FDOT(K,J) = 0.0D0
          D0(K,J) = 0.0D0
          D1(K,J) = 0.0D0
          D2(K,J) = 0.0D0
          D3(K,J) = 0.0D0
          D0R(K,J) = 0.0D0
          D1R(K,J) = 0.0D0
          D2R(K,J) = 0.0D0
          D3R(K,J) = 0.0D0
   10 CONTINUE
*
*       Set large X0 & X to avoid perturber selection (no escape).
      X0(1,J) = 1.0E+06
      X(1,J) = 1.0E+06
*
*       Specify zero neighbour membership.
      LIST(1,J) = 0
*
*       Copy neighbour list for routine NBREM.
      NNB = LIST(1,J)
      DO 20 L = 2,NNB+1
          JPERT(L-1) = LIST(L,J)
   20 CONTINUE
*
*       Remove ghost from lists of neighbours containing body #ICH.
      JLIST(1) = J
      CALL NBREM(ICH,1,NNB)
*
      RETURN
*
      END
