*
* const_mobse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32),counter_ran3
      COMMON /RAND3/ idum2,iy,ir,counter_ran3
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag,piflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag,piflag
      INTEGER bhflag
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 sigma1,sigma2,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma1,sigma2,bhflag
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,33)
      COMMON /BINARY/ bcm,bpp
*
