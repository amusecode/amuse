      SUBROUTINE SWITCH(Y)
*
*
*       Switching of chain.
*       -------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      LOGICAL  KSLOW,KCOLL,Itest
      REAL*8  Y(NMX8),XCNEW(NMX3),KSCH,ksnew(nmx)
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      INTEGER  IOLD(NMX)
*
*
*       Copy Y-array to COMMON.
      CALL YSAVE(Y)
*
*	First transform to chain coordinates.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL KSPHYS(Q(KS1),P(KS1),XC(L1),WC(L1))
      END DO
*
      L2=3*(INAME(1)-1)
      DO K=1,3
      X(L2+K)=0.0
      END DO
*
*       Set X for determining new chain indices.
      DO I=1,N-1
      L=3*(I-1)
      L1=L2
      L2=3*(INAME(I+1)-1)
      DO K=1,3
      X(L2+K)=X(L1+K)+XC(L+K)
      END DO
      END DO
*
*	Save the old chain indices.
      DO I=1,N
      IOLD(I)=INAME(I)
      END DO
*
*       Select new indices.
      CALL SELECT
*
*       EXit if chain vectors are unchanged.
      ISW=0
      DO I=1,N
      IF(INAME(I).NE.IOLD(I))ISW=ISW+1
      END DO
      IF(ISW.EQ.0)RETURN
*
*	Transform chain momenta.
      L1=3*(IOLD(1)-1)
      LN=3*(IOLD(N)-1)
      L=3*(N-2)
      DO K=1,3
      PI(L1+K)=-WC(K)
      PI(LN+K)=WC(L+K)
      END DO
      DO I=2,N-1
      L=3*(I-1)
      LI=3*(IOLD(I)-1)
      DO K=1,3
      PI(LI+K)=WC(L+K-3)-WC(L+K)
      END DO
      END DO
      L1=3*(INAME(1)-1)
      LN=3*(INAME(N)-1)
      L=3*(N-2)
      DO K=1,3
      WC(K)=-PI(L1+K)
      WC(L+K)=PI(LN+K)
      END DO
      DO I=2,N-2
      L=3*(I-1)
      LI=3*(INAME(I)-1)
      DO K=1,3
      WC(L+K)=WC(L+K-3)-PI(LI+K)
      END DO
      END DO
*
*       Construct new chain coordinates.
      DO I=1,3*(N-1)
      XCNEW(I)=0.0
      END DO
*       Transformation matrix (old to new) has only coefficients -1, 0 or +1.
      DO ICNEW=1,N-1
*       Find K0 & K1 for IOLD(K0) = INAME(ICNEW) & IOLD(K1) = INAME(ICNEW+1).
      LNEW=3*(ICNEW-1)
      DO I=1,N
      IF(IOLD(I).EQ.INAME(ICNEW))K0=I
      IF(IOLD(I).EQ.INAME(ICNEW+1))K1=I
      END DO
      DO ICOLD=1,N-1
      LOLD=3*(ICOLD-1)
      IF((K1.GT.ICOLD).AND.(K0.LE.ICOLD))THEN
*       Add.
      DO K=1,3
      XCNEW(LNEW+K)=XCNEW(LNEW+K)+XC(LOLD+K)
      END DO
      ELSE IF((K1.LE.ICOLD).AND.(K0.GT.ICOLD))THEN
*	Subtract.
      DO K=1,3
      XCNEW(LNEW+K)=XCNEW(LNEW+K)-XC(LOLD+K)
      END DO
      END IF
      END DO
      END DO
*
*	Perform KS-transformations and update chain coordinates & RINV.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS1=4*(I-1)+1
      CALL PHYSKS(XCNEW(L1),WC(L1),Q(KS1),P(KS1))
      DO K=1,3
      XC(L1+K-1)=XCNEW(L1+K-1)
      END DO
      RINV(I)=1.0/(Q(KS1)**2+Q(KS1+1)**2+Q(KS1+2)**2+Q(KS1+3)**2)
      END DO
*
*       Swop Ksch acording to Iold and Iname.
      do i=1,n-1
          ksnew(i)=1.0
      end do
      do i=1,n-1
          if(ksch(i).ne.1.0d0)then
              i1=Iold(i)
              i2=Iold(i+1)
              do j=1,n-1
              Itest=((i1.eq.iname(j)).and.(i2.eq.iname(j+1))).or.
     &              ((i2.eq.iname(j)).and.(i1.eq.iname(j+1)))
              if(Itest)then
                  ksnew(j)=ksch(i)
              end if
              end do
          end if
      end do
      do i=1,n-1
          ksch(i)=ksnew(i)
      end do
*       Define auxiliary quantities.
      MASS=0.0
      DO I=1,N
      L=3*(I-1)
      MC(I)=M(INAME(I))
      MASS=MASS+MC(I)
      END DO
*
      do i=1,n
          tk1(i)=-1./MC(I)
      end do
      do i=1,n-1
          if(ksch(i).ne.1.0d0)then
              tk1(i)=tk1(i)/Ksch(i)
              tk1(i+1)=tk1(i+1)/Ksch(i)
          end if
      end do
      DO I=1,N-1
          TKK(I)=.5D0*(-tk1(i)-tk1(i+1))
          MKK(I)=MC(I)*MC(I+1)/Ksch(i)
          DO J=I+1,N
              MIJ(I,J)=MC(I)*MC(J)
              MIJ(J,I)=MIJ(I,J)
           END DO
      END DO
      do i=1,n-1
          m12=mc(i)+mc(i+1)
          dt12=0.5d0*(1.d0-1.d0/Ksch(i))/m12
          if(i.gt.1) TKK(i-1)=tkk(i-1)+dt12
          if(i.lt.n-1) TKK(i+1)=tkk(i+1)+dt12
          if(i.gt.1.and.i.lt.n-1) TK2(i)=-2.0d0*dt12
      end do
*       Note: TK2(0) & TK2(N) should be zero but never used.
      TK2(0) = 0.0
      TK2(N) = 0.0
*
*       Copy Y-array from COMMON.
      CALL YCOPY(Y)
*
      RETURN
      END
