      SUBROUTINE YCOPY(Y)
*
*
*       Copy Y-array from COMMON.
*       -------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8 Y(NMX8)
*
*
      NC=N-1
      L=0
      DO I=1,4*NC
      L=L+1
      Y(L)=Q(I)
      END DO
      DO I=1,3
      L=L+1
      Y(L)=CMX(I)
      END DO
      L=L+1
      Y(L)=ENERGY
      DO I=1,4*NC
      L=L+1
      Y(L)=P(I)
      END DO
      DO I=1,3
      L=L+1
      Y(L)=CMV(I)
      END DO
      L=L+1
      Y(L)=CHTIME
*
      RETURN
      END
