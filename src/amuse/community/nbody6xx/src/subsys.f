      SUBROUTINE SUBSYS(NSYS,CM)
*
*
*       Initialization of subsystem.
*       ----------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  CM(7)
*
*
*       Increase subsystem counter and set current index.
      NSUB = NSUB + 1
      ISUB = NSUB
*
*       Set zero name & mass in last component to distinguish triple case.
      IF (NSYS.EQ.3) THEN
          NAMES(4,ISUB) = 0
          BODYS(4,ISUB) = 0.0D0
      END IF
*
*       Save global masses & names of new subsystem components.
      BX = 0.0
      DO 10 L = 1,NSYS
          J = JLIST(L)
          BODYS(L,ISUB) = BODY(J)
          NAMES(L,ISUB) = NAME(J)
          IF (BODY(J).GT.BX) THEN
              BX = BODY(J)
              LX = L
          END IF
   10 CONTINUE
*
*       Form ghosts and initialize integration variables (NSYS = 3 or 4).
      DO 20 L = 1,NSYS
          J = JLIST(L)
          BODY(J) = 0.0D0
          T0(J) = 1.0E+06
          LIST(1,J) = 0
          DO 15 K = 1,3
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
   15     CONTINUE
*       Set large X0 & X to avoid perturber selection (no escape).
          X0(1,J) = 1.0E+06
          X(1,J) = 1.0E+06
   20 CONTINUE
*
*       Place c.m. of subsystem in first location (ICOMP may switch!).
      ICOMP = JLIST(1)
*       Ensure heaviest body is selected in case of ARchain.
      IF (ISYS(ISUB).EQ.4) THEN
          ICOMP = JLIST(LX)  
      END IF
*       Define zero name for identification (only chain or ARchain c.m.).
      IF (ISYS(ISUB).GE.3) NAME(ICOMP) = 0
*
*       Initialize c.m in ICOMP (redefined above).
      T0(ICOMP) = TIME
      BODY(ICOMP) = CM(7)
      DO 30 K = 1,3
          X(K,ICOMP) = CM(K)
          X0(K,ICOMP) = CM(K)
          XDOT(K,ICOMP) = CM(K+3)
          X0DOT(K,ICOMP) = CM(K+3)
   30 CONTINUE
*
*       Predict coordinates & velocities for all neighbours (order FDOT).
      NNB = LIST(1,ICOMP)
      CALL XVPRED(ICOMP,NNB)
*
*       Obtain new neighbour list (to ensure membership > 0).
      RS0 = RS(ICOMP)
      CALL NBLIST(ICOMP,RS0)
*
*       Construct force polynomials for c.m. motion.
      CALL FPOLY1(ICOMP,ICOMP,0)
      CALL FPOLY2(ICOMP,ICOMP,0)
*
*       Initialize decision-making variables for multiple regularization.
      T0S(ISUB) = TIME
      TS(ISUB) = TIME
      STEPS(ISUB) = STEP(ICOMP)
*
*       Obtain maximum unperturbed separation based on dominant neighbour.
      CALL EXTEND(ISUB)
*
      RETURN
*
      END
