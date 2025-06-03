      SUBROUTINE XTF(FW,FV,XCM,CHTIME)
*
*
*       External chain force.
*       ---------------------
*
      INCLUDE 'commonc.h'
      REAL*8 FW(NMX3),FV(3),XCM(3)
      REAL*8 ACC(NMX3),DDP(NMX3),Q0(3)
*
*
      DO K=1,3
      Q0(K)=0.0
      END DO
*
      DO I=1,N
      L=3*(I-1)
      DO K=1,3
      Q0(K)=Q0(K)+XI(L+K)*MC(I)/MASS
      END DO
      END DO
*
*       Re-arrange according to INAME(I) and add c.m. (if needed).
      DO I=1,N
      L=3*(I-1)
      LF=3*(INAME(I)-1)
      DO K=1,3
      X(LF+K)=XI(L+K)-Q0(K)
*    &       +XCM(K)
      END DO
      END DO
*
*       Obtain the perturbing accelerations.
      CALL XTPERT(ACC,CHTIME)
*
*	Centre of mass force.
      DO K=1,3
      FV(K)=0.0
      END DO
      DO I=1,N
      L=3*(I-1)
      DO K=1,3
      FV(K)=FV(K)+M(I)/MASS*ACC(L+K)
      END DO
      END DO
*
*       Physical chain forces.
      DO I=1,N
      L=3*(I-1)
      LF=3*(INAME(I)-1)
      DO K=1,3
      DDP(L+K)=MC(I)*(ACC(LF+K)-FV(K))
      END DO
      END DO
*
      L=3*(N-2)
      DO K=1,3
      FW(K)=-DDP(K)
      FW(L+K)=DDP(L+K+3)
      END DO
      DO I=2,N-2
      L=3*(I-1)
      DO K=1,3
      FW(L+K)=FW(L+K-3)-DDP(L+K)
      END DO
      END DO
*
      RETURN
      END
