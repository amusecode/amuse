      SUBROUTINE HPSORT(N,Array,Indx)
*
*
*       Heapsort algorithm (from Press).
*       --------------------------------
*
      implicit real*8 (a-h,o-z)
*     dimension Array(1),Indx(1)
      REAL*8 Array(N)
      INTEGER Indx(N)
*
*
      do 11 j=1,N
      Indx(j)=J
11    continue
      if(N.lt.2)RETURN
      l=N/2+1
      ir=N
10    CONTINUE
      IF(l.gt.1)THEN
      l=l-1
      Indxt=Indx(l)
      q=Array(Indxt)
      ELSE
      Indxt=Indx(ir)
      q=Array(Indxt)
      Indx(ir)=Indx(1)
      ir=ir-1
      IF(ir.eq.1)THEN
      Indx(1)=Indxt
      RETURN
      END IF
      END IF
      i=l
      j=l+l
20    IF(j.le.ir)THEN
          IF(j.lt.ir)THEN
             IF(Array(Indx(j)).lt.Array(Indx(j+1)))j=j+1
          END IF
          IF(q.lt.Array(Indx(j)))THEN
             Indx(i)=Indx(j)
             i=j
             j=j+j
          ELSE
             j=ir+1
          END IF
       GOTO 20
       END IF
       Indx(i)=Indxt
       GO TO 10
       END
