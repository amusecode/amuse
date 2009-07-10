*
* const_bse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      INTEGER bhflag
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 sigma,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp
*
