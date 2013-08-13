      SUBROUTINE SELECT
*
*
*       Select chain indices.
*       ---------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8 RIJ2(NMXM)
      INTEGER IC(NMX2),IJ(NMXM,2),IND(NMXM)
      LOGICAL USED(NMXM),SUC,LOOP
      SAVE
      
*
*
      L=0
      DO I=1,N-1
      DO J=I+1,N
      L=L+1
      RIJ2(L) = (X(3*I-2)-X(3*J-2))**2+(X(3*I-1)-X(3*J-1))**2+
     &          (X(3*I)-X(3*J))**2
*       Include masses for defining dominant interaction (09/06).
      RIJ2(L) = RIJ2(L)/(M(I) + M(J))
      IJ(L,1)=I
      IJ(L,2)=J
      USED(L)=.FALSE.
      END DO
      END DO
      CALL HPSORT(L,RIJ2,IND)
      LMIN=1+NMX
      LMAX=2+NMX
      IC(LMIN)=IJ(IND(1),1)
      IC(LMAX)=IJ(IND(1),2)
      USED(IND(1))=.TRUE.
1     DO I=2,L
      LI=IND(I)
      IF(.NOT.USED(LI))THEN
      CALL REDRAW(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
      IF(SUC)THEN
      USED(LI)=.TRUE.
      GOTO 2
      ELSE
      USED(LI)=LOOP
      END IF
      END IF
      END DO
2     IF(LMAX-LMIN+1.LT.N)GO TO 1
      L=0
      DO I=LMIN,LMAX
      L=L+1
      INAME(L)=IC(I)
      END DO
*
      RETURN
      END
