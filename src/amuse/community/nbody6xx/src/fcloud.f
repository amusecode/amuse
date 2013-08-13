      SUBROUTINE FCLOUD(I,FREG,FDR,KCASE)
*
*
*       Force due to interstellar clouds.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
      REAL*8  FREG(3),FDR(3),A(6)
*
*
*       Distinguish between standard case and initialization.
      IF (KCASE.EQ.1) THEN
          RB3 = RB2*SQRT(RB2)
*
*       Sum over all clouds.
          DO 10 ICL = 1,NCL
              DO 5 K = 1,3
                  A(K) = XCL(K,ICL) - X(K,I)
                  A(K+3) = XDOTCL(K,ICL) - XDOT(K,I)
    5         CONTINUE
*
              RIJ2 = A(1)**2 + A(2)**2 + A(3)**2 + RCL2(ICL)
              A7 = BODYCL(ICL)/RB3
              A8 = BODYCL(ICL)/(RIJ2*SQRT(RIJ2))
              A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))/RIJ2
*
*       Add cloud force & derivative but subtract the average field effect.
              DO 8 K = 1,3
                  FREG(K) = FREG(K) + A(K)*A8 + (X(K,I) - RDENS(K))*A7
                  FDR(K) = FDR(K) + (A(K+3) - A(K)*A9)*A8 + XDOT(K,I)*A7
    8         CONTINUE
   10     CONTINUE
      ELSE
*
*       Obtain total cloud force & first derivative (for routine FPOLY1).
          DO 20 ICL = 1,NCL
              DO 15 K = 1,3
                  A(K) = XCL(K,ICL) - X(K,I)
                  A(K+3) = XDOTCL(K,ICL) - XDOT(K,I)
   15         CONTINUE
*
              RIJ2 = A(1)**2 + A(2)**2 + A(3)**2 + RCL2(ICL)
              A8 = BODYCL(ICL)/(RIJ2*SQRT(RIJ2))
              A9 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))/RIJ2
*
*       Include direct cloud contributions but ignore the average field.
              DO 18 K = 1,3
                  FR(K,I) = FR(K,I) + A(K)*A8
                  D1R(K,I) = D1R(K,I) + (A(K+3) - A(K)*A9)*A8
   18         CONTINUE
   20     CONTINUE
      END IF
*
      RETURN
*
      END
