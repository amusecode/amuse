      BLOCK DATA BLOCK
*
*
*       Run time initializations.
*       -------------------------
*
      REAL*8  EP,TFAC
      COMMON/BSSAVE/  EP(4),TFAC,ITFAC,JC,NHALF2
*
*
*       Initialize Bulirsch-Stoer variables.
      DATA  EP  /0.04D0,0.0016D0,0.64D-4,0.256D-5/
      DATA  ITFAC,JC,NHALF2  /0,-1,16/
*
      END
