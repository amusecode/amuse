      SUBROUTINE TRANSX
*
*
*       Transformations from KS to chain variables.
*     ---------------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8 Q0(3)
*
*
*	Transform to chain coordinates.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL KSPHYS(Q(KS1),P(KS1),XC(L1),WC(L1))
      END DO
*
*	Obtain physical variables from chain quantities.
      L=3*(N-2)
      DO K=1,3
      PI(K)=-WC(K)
      PI(L+K+3)=WC(L+K)
      END DO
      DO I=2,N-1
      L=3*(I-1)
      DO K=1,3
      PI(L+K)=WC(L+K-3)-WC(L+K)
      END DO
      END DO
*
      DO K=1,3
      XI(K)=0.0
      Q0(K)=0.0
      END DO
*
      DO I=1,N-1
      L=3*(I-1)
      DO K=1,3
      XI(L+3+K)=XI(L+K)+XC(L+K)
      END DO
      END DO
*
      DO I=1,N
      L=3*(I-1)
      DO K=1,3
      Q0(K)=Q0(K)+XI(L+K)*MC(I)/MASS
      END DO
      END DO
*
*	Rearrange according to INAME(I) and add CM (if required).
      DO I=1,N
      L=3*(I-1)
      LF=3*(INAME(I)-1)
      DO K=1,3
      X(LF+K)=XI(L+K)-Q0(K)
*    &       +CMX(K)
      V(LF+K)=PI(L+K)/MC(I)
*    &       +CMV(K)
      END DO
      END DO
*
      RETURN
      END
