      SUBROUTINE DIFSY3(N,EPS,H,X,Y)
*
*
*       Bulirsch-Stoer integrator.
*       --------------------------
*
*       Works if Gamma = (H - E)/L. For other time transformations 'eps'
*       must be scaled appropriately such that the test is esentially
*       of the form (H - E)/L < EPS. This is achieved by using TFAC.
*       Convergence test: ABS(Q'*DP) < EPS*TFAC & ABS(P'*DQ) < EPS*TFAC.
*       Reference: Mikkola (1987). In 'The Few Body Problem' p. 311.
*       This works only if eqs. are canonical and we have P:s & Q:s.
*       One additional eq. is allowed (e.g. for t'=...??) but not checked.
*
*
      PARAMETER  (NMX=17,NMX7=7*NMX)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      REAL*8  Y(N),YA(NMX),YL(NMX),YM(NMX),DY(NMX),DZ(NMX),DT(NMX,7),
     &        D(7),X,XN,H,G,B,B1,U,V,C,TA,W,DABS
*
      LOGICAL  KONV,BO,KL,GR,FYBAD
      DATA DT /NMX7*0.0D0/
*     DIMENSION  EP(4)
*     DATA  EP/0.4D-1,0.16D-2,0.64D-4,0.256D-5/
*
*
      SAVE
      
      NHALF2=(N/2)*2
      JTI=0
      FY=1.
      DO I=1,N
      YA(I)=Y(I)
      END DO
      CALL DERQP3(Y(1),Y(9),DZ(1),DZ(9),DZ(17))
      IF (JC.GT.0)H=DSC
   10 XN=X+H
      BO=.FALSE.
      M=1
      JR=2
      JS=3
      DO J=1,10
      IF(BO)THEN
      D(2)=16D0/9.d0
      D(4)=64.d0/9.d0
      D(6)=256.D0/9.d0
      ELSE
      D(2)=2.25D0
      D(4)=9.D0
      D(6)=3.6D1
      END IF
      IF(J.GT.7)THEN
      L=7
      D(7)=6.4D1
      ELSE
      L=J
      D(L)=M*M
      END IF
      KONV=L.GT.3
      M=M+M
      G=H/M
      B=G+G
      M=M-1
      DO I=1,N
      YL(I)=YA(I)
      YM(I)=YA(I)+G*DZ(I)
      END DO
      DO K=1,M
      CALL DERQP3(YM(1),YM(9),DY(1),DY(9),DY(17))
      DO I=1,N
      U=YL(I)+B*DY(I)
      YL(I)=YM(I)
      YM(I)=U
      END DO
      END DO
*       Switch on tolerance scaling indicator.
      ITFAC = 1
      CALL DERQP3(YM(1),YM(9),DY(1),DY(9),DY(17))
      KL=L.LT.2
      GR=L.GT.5
      FS=0.
      DO I=1,N
      V=DT(I,1)
      C=(YM(I)+YL(I)+G*DY(I))*0.5D0
      DT(I,1)=C
      TA=C
      IF(.NOT.KL)THEN
      DO K=2,L
      B1=D(K)*V
      B=B1-C
      W=C-V
      U=V
      IF(B.NE.0.D0)THEN
      B=W/B
      U=C*B
      C=B1*B
      END IF
      V=DT(I,K)
      DT(I,K)=U
      TA=U+TA
      END DO
      IS=I+N/2
      IF(IS.GT.NHALF2)IS=I-(N/2)
*       Use modified EPS for convergence test.
      DYIS=DABS(DY(IS))/TFAC
      IF(I.GT.NHALF2)DYIS=0.0
      IF(KONV)THEN
      TEST=DABS( (Y(I)-TA)*DYIS )
      IF(TEST.GT.EPS) KONV=.FALSE.
      END IF
      IF(.NOT.GR)THEN
      FV=DABS(W)*DYIS
      IF(FS.LT.FV) FS=FV
      END IF
      END IF
      Y(I)=TA
      END DO
      IF(FS.NE.0.D0)THEN
      FA=FY
      K=L-1
      FY=(EP(K)/FS)**(1./(L+K))
      FA7=0.7*FA
      IF(L.EQ.2)FA7=0.0
      FYBAD=.NOT.((FA7.GT.FY).OR.(FY.GT.0.7))
      IF(FYBAD)THEN
      H=H*FY
      JTI=JTI+1
      IF(JTI.GT.5)THEN
      H=0.0
      DO I=1,N
      Y(I)=YA(I)
      END DO
      RETURN
      END IF
      GOTO 10
      END IF
      END IF
      IF(KONV)THEN
      X=XN
      H=H*FY
      RETURN
      END IF
      D(3)=4.D0
      D(5)=1.6D1
      BO=.NOT.BO
      M=JR
      JR=JS
      JS=M+M
      END DO
      H=0.5*H
      GOTO 10
      END
