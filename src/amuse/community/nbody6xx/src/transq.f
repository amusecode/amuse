      SUBROUTINE TRANSQ
*
*
*       Transformations from chain variables to KS.
*       -------------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8  KSCH
      LOGICAL KSLOW,KCOLL
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
*
*
*        Centre of mass.
      DO K=1,3
      CMX(K)=0.0
      CMV(K)=0.0
      END DO
      MASS=0.0
      DO I=1,N
      L=3*(I-1)
      MC(I)=M(INAME(I))
      MASS=MASS+MC(I)
      DO K=1,3
      CMX(K)=CMX(K)+M(I)*X(L+K)
      CMV(K)=CMV(K)+M(I)*V(L+K)
      END DO
      END DO
*
      DO K=1,3
      CMX(K)=CMX(K)/MASS
      CMV(K)=CMV(K)/MASS
      END DO
*        Auxiliary quantities.
      DO I=1,N-1
      TKK(I)=.5D0*(1./MC(I)+1./MC(I+1))
      TK1(I)=-1./MC(I)
      MKK(I)=MC(I)*MC(I+1)
      DO J=I+1,N
      MIJ(I,J)=MC(I)*MC(J)
      MIJ(J,I)=MIJ(I,J)
      END DO
      END DO
*
*       Initialize the whole slow-down array TK2.
      DO I=0,N
      TK2(I)=0.0
      END DO
*
*       Pysical momenta.
      DO I=1,N
      L=3*(I-1)
      LF=3*INAME(I)-3
      DO K=1,3
      XI(L+K)=X(LF+K)
      PI(L+K)=M(INAME(I))*(V(LF+K)-CMV(K))
      END DO
      END DO
*
*        Chain momenta.
      L=3*(N-2)
      DO K=1,3
      WC(K)=-PI(K)
      WC(L+K)=PI(L+K+3)
      END DO
      DO I=2,N-2
      L=3*(I-1)
      DO K=1,3
      WC(L+K)=WC(L+K-3)-PI(L+K)
      END DO
      END DO
*
*       Include artificial STOP to get round compiler bug.
      IF (XI(4)-XI(1).EQ.0.0D0) STOP
*
*       Chain coordinates.
      DO I=1,N-1
      L=3*(I-1)
      DO K=1,3
      XC(L+K)=XI(L+K+3)-XI(L+K)
      END DO
      END DO
*
*        KS-transformations.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL PHYSKS(XC(L1),WC(L1),Q(KS1),P(KS1))
      END DO
*
      RETURN
      END
