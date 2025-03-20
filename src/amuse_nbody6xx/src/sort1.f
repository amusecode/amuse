      SUBROUTINE SORT1(N,RA,RB)
*     
*     
*     Heapsort method (Press p. 231).
*     -------------------------------
*     
      INTEGER  RB(N),RRB
      REAL*8  RA(N),RRA
*     
*     
      IF (N.LE.1) RETURN
*     !     bug fix Nov 2007.
      L = N/2+1
      IR=N
 10   CONTINUE
      IF (L.GT.1) THEN
         L=L-1
         RRA=RA(L)
         RRB=RB(L)
      ELSE
         RRA=RA(IR)
         RRB=RB(IR)
         RA(IR)=RA(1)
         RB(IR)=RB(1)
         IR=IR-1
         IF (IR.EQ.1) THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
         END IF
      END IF
      I=L
      J=L+L
 20   IF (J.LE.IR) THEN
         IF (J.LT.IR) THEN
            IF (RA(J).LT.RA(J+1)) J=J+1
         END IF
         IF (RRA.LT.RA(J)) THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         END IF
         GO TO 20
      END IF
      RA(I)=RRA
      RB(I)=RRB
      GO TO 10
*     
      END
