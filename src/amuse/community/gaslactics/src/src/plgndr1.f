
      FUNCTION PLGNDR1(L,X)
      PMM=1.
      IF(L.EQ.0) THEN
        PLGNDR1=PMM
      ELSE
        PMMP1=X*PMM
        IF(L.EQ.1) THEN
          PLGNDR1=PMMP1
        ELSE
          DO 12 LL=2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL-1)*PMM)/(LL)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR1=PLL
        ENDIF
      ENDIF
      RETURN
      END
