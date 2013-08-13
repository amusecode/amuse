      BLOCK DATA BLOCK
*
*
*       Run time initializations.
*       -------------------------
*
      INCLUDE 'params.h'
      REAL*8  EP,DSC,FACM,TFAC,RANGE
      COMMON/RAND2/  IY,IDUM2,IV(32),IXYZ(65)
      COMMON/ICPU0/  ICPU
      COMMON/IND6/  IND(6)
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/COUNTS/  NCOUNT(168)
*
*
*       Initialize COMMON indices & B-S data array.
      DATA  ICPU,IY,IDUM2,IV,IXYZ  /0,0,123456789,32*0,65*0/
      DATA  IND  /1,2,3,4,5,6/
      DATA  EP  /0.04D0,0.0016D0,0.64D-4,0.256D-5/
      DATA  RANGE  /50.0D0/
      DATA  ISLOW  /1,2,4,8,16,32,64,128,256,512/
      DATA  NCOUNT  /168*0/
*
      END
