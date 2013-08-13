      SUBROUTINE REINIT(ISUB)
*
*
*       Re-initialization of chain system.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,ANG(3),FIRR(3),FD(3)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/ECHAIN/  ECH
*
*
*       Set phase indicator for step reduction (routine STEPS).
      IPHASE = 8
*
*       Obtain current coordinates & velocities of neighbours (for FPOLY).
      NNB = LIST(1,ICH)
      CALL XVPRED(ICH,NNB)
*
*       Initialize force polynomials (include differential force correction).
      CALL FPOLY1(ICH,ICH,0)
      DO 5 K = 1,3
          FIRR(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE
      CALL CHFIRR(ICH,0,X(1,ICH),XDOT(1,ICH),FIRR,FD)
      DO 10 K = 1,3
          F(K,ICH) = F(K,ICH) + FIRR(K)
          FDOT(K,ICH) = FDOT(K,ICH) + FD(K)
   10 CONTINUE
      CALL FPOLY2(ICH,ICH,0)
*
*       Obtain maximum unperturbed separation based on dominant neighbour.
      CALL EXTEND(ISUB)
*
*       Ensure new global coordinates are available in chpert.f.
      CALL XCPRED(0)
*
*       Update decision-making variables for chain regularization.
      TS(ISUB) = TIME
      STEPS(ISUB) = 0.01*STEP(ICH)
*
*       Initialize perturber list for chain (use modified value of RSUM).
      CALL CHLIST(ICH)
*
*       Re-calculate new energy of chain system (just in case).
      CALL CONST(XCH,VCH,M,NN,ECH,ANG,ALAG)
*
*       Set phase indicator < 0 to ensure new time-step list in INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
