! simple code for random number sequences
! from numerical recipes..
function ranfb(ix)
  real :: y,ranfb
  integer :: ix,ix1,ix2,ia1,ia2,ic,j1,j2,j3

  ia1=273
  ia2=25301
  ic=226908345
  
  if (ix.lt.0) ix=-ix
  if (ix.gt.2**30) ix=MOD(ix,2**30)

  ix1=ix/(2**15)
  ix2=ix-ix1*(2**15)
  j1=MOD(ix1*ia2,2**15)*2**15
  j2=MOD(ix2*ia1,2**15)*2**15
  j3=ix2*ia2
  ix=MOD(j1+j2,2**30)
  ix=MOD(ix+j3,2**30)
  ix=MOD(ix+ic,2**30)
  y=FLOAT(ix)

  ranfb=y/2**30

end function

function ran2(ix)
 real ran2,ranfb
 integer ix
 ran2=ranfb(ix)
end function 

FUNCTION gasdev(iseed)
 REAL :: gasdev,x(2),r2,fac,ran2
 INTEGER :: iseed
 INTEGER, SAVE :: iset=0
 REAL, SAVE :: saved
 IF(iset.EQ.0) THEN
   r2=2.
   DO WHILE(r2.GE.1.OR.r2.EQ.0)
     x(1)=ran2(iseed)
     x(2)=ran2(iseed)
     x=2*x-1
     r2=sum(x**2)
   ENDDO
   fac=sqrt(-2.*log(r2)/r2)
   saved=x(1)*fac
   gasdev=x(2)*fac
   iset=1
 ELSE
   iset=0
   gasdev=saved
 ENDIF
END FUNCTION

!
   SUBROUTINE makefractal(nstar,rstar,vstar,fdim,idum)
!
! makes fractal of nstar particles of dimension fdim, using ndiv 
! subunits forcing the number of cells if force=.true.
! Simon Goodwin & Ant Whitworth (2004, A&A, 413, 929)
!
! fdim (real) = fractal dimension
! ndiv (int) = number of divisions of box (>=2)
! nstar (int) = number of stars in *final* distribution
! rstar (3,n) = positions to return
! vstar (3,n) = velocities to return
! force (logical) = make box dimension such that integer no. children
! vtype (logical) = .true. - inherit velocities
! vir (real) = virial ratio (virialised = 0.5)
! idum = random number seed
! outfile = output file name 
   IMPLICIT NONE
   INTEGER :: ndiv,nstar,idum,isnext,iforce
   INTEGER :: is,nsmax,islast,flagabort,k
   INTEGER :: k1,k2,k3,i,isfirst,nswap,ndo,nsnow
   REAL :: fdim,dum,delta,dr(3),prob,rad,cov(3)
   REAL :: how,how2,ran2,com(3)
   REAL :: gasdev
   LOGICAL :: force,vtype
   INTEGER :: n,isurv,level
   REAL :: rstar(1:3,nstar),vstar(1:3,nstar)
   INTEGER, DIMENSION(:), ALLOCATABLE :: forced
   REAL, DIMENSION(:,:), ALLOCATABLE :: r,v
!
! set initial parameters
   ndiv=2
   force=.TRUE.
   vtype=.TRUE.
   ALLOCATE(r(1:3,1:nstar*50))
   ALLOCATE(v(1:3,1:nstar*50))
   ALLOCATE(forced(1:ndiv**3))
!
! ===================================================================
! probability of a child surviving
   prob=REAL(ndiv)**(fdim-3.)
!
! list characteristics
   WRITE(6,*) 'fractal dimension', fdim
   WRITE(6,*) 'ndiv', ndiv
   WRITE(6,*) 'number of stars', nstar
   WRITE(6,*) 'survival chance', prob
   WRITE(6,*) 'forcing?', force
!
! warning
   how=REAL(ndiv)**fdim
   WRITE(6,*) '# of children expected', how
   how2=how - INT(how)
   IF (how2>0.2 .AND. how2<0.8) THEN
     WRITE(6,*) 'WARNING:'
     IF (how2<0.5) THEN
       dum=LOG(REAL(INT(how)))/LOG(REAL(ndiv))
       WRITE(6,*) 'could look like fdim=', dum
     ELSE
       dum=LOG(REAL(INT(how)+1))/LOG(REAL(ndiv))
       WRITE(6,*) 'could look like fdim=', dum
     END IF
   END IF
!
   IF (force) THEN
     iforce=INT(how + 0.5)
     WRITE(6,*) 'forcing to survivors=', iforce
     dum=LOG(REAL(iforce))/LOG(REAL(ndiv))
     WRITE(6,*) 'which means fdim=', dum
   END IF
!
! build fractal cube
   r=0.
   v=0.
   nsmax=1
   islast=1
   flagabort=0
   delta=1./REAL(ndiv)
   level=0
   DO
     level=level + 1
     isfirst=nsmax + 1
     isnext=isfirst
! loop over all parents
     DO is=1,nsmax
!
! if forced decide on survivors now
       IF (force) THEN
         forced=0
         DO n=1,iforce
4          isurv=INT(ran2(idum)*REAL(ndiv**3)) + 1
           IF (forced(isurv)==1) GOTO 4
           forced(isurv)=1
         END DO
       END IF
! loop over all ndiv**3 possible children
       n=0
       DO k1=1,ndiv
         dr(1)=REAL(2*k1-1-ndiv)*delta
         DO k2=1,ndiv
           dr(2)=REAL(2*k2-1-ndiv)*delta
           DO k3=1,ndiv
             n=n+1
             dr(3)=REAL(2*k3-1-ndiv)*delta
             r(1:3,isnext)=r(1:3,is) + dr(1:3)
             IF (vtype) THEN
! keep parent's velocity
               v(1:3,isnext)=v(1:3,is)
             END IF
! do children survive?
! forced or not forced here
             IF (force) THEN
               IF (forced(n)==1) THEN
                 islast=isnext
                 isnext=isnext + 1
               END IF
             ELSE
               IF (ran2(idum)<prob) THEN
                 islast=isnext
                 isnext=isnext + 1
               END IF
             END IF
           END DO
         END DO
       END DO
     END DO
! are there children (we hope so)
     IF (islast/=nsmax) THEN
! delete parents
       DO is=isfirst,islast
         r(1:3,is-nsmax)=r(1:3,is)
         v(1:3,is-nsmax)=v(1:3,is)
       END DO
! new total number
       nsmax=islast - nsmax
       DO is=nsmax+1,islast
         r(1:3,is)=0.
         v(1:3,is)=0.
       END DO
! resize delta
       delta=delta/REAL(ndiv)
     ELSE
       flagabort=1
       EXIT
     END IF
     IF (flagabort==1) STOP 'no children'
! add noise to positions and velocity
     DO is=1,nsmax
       DO k=1,3
         r(k,is)=r(k,is) + (2.*ran2(idum) - 1.)*delta
!         v(k,is)=v(k,is) + 2.**(fdim-1.)*gasdev(idum)
         v(k,is)=v(k,is) + delta*gasdev(idum)
       END DO
     END DO
! exit if enough particles
     IF (nsmax>3.5*nstar) EXIT
!     WRITE(6,*) 'done', delta, nsmax
   END DO
   WRITE(6,*) 'Made initial cubic distribution'
!
! cut out sphere
   nsnow=nsmax
   DO is=1,nsmax
     rad=r(1,is)**2 + r(2,is)**2 + r(3,is)**2
     IF (rad>1.) THEN 
       r(1,is)=69.
       nsnow=nsnow - 1
     END IF
   END DO
   nswap=nsmax
   ndo=1
   DO
8    IF (r(1,ndo)==69.) THEN
       r(1:3,ndo)=r(1:3,nswap)
       nswap=nswap - 1
       GOTO 8
     END IF
     ndo=ndo + 1
     IF (ndo==nsnow) EXIT
   END DO
   nsmax=ndo - 1
   IF (nsmax<nstar) STOP 'nsmax<nstar - change random number seed'
!
! we now have a fractal of nsmax stars
! need to randomly remove stars until we get nstar
   DO
! star to remove
     is=INT(ran2(idum)*REAL(nsmax)) + 1
     r(1:3,is)=r(1:3,nsmax)
     nsmax=nsmax - 1
     IF (nsmax==nstar) EXIT
   END DO
   WRITE(6,*) 'Cut to size', nsmax
!
! now put r in rstar array
   DO is=1,nstar
     rstar(1:3,is)=r(1:3,is)
     vstar(1:3,is)=v(1:3,is)
   END DO
!
   IF (.NOT.vtype) THEN
     DO is=1,nstar
       DO k=1,3
         vstar(k,is)=(2.*ran2(idum) - 1.)
       END DO
     END DO
   END IF
!
! set coms to zero
   com=0.
   cov=0.
   DO i=1,nstar
     com(1:3)=com(1:3) + rstar(1:3,i)
     cov(1:3)=cov(1:3) + vstar(1:3,i)
   END DO
   com=com/REAL(nstar)
   cov=cov/REAL(nstar)
   WRITE(6,*) 'COM', com(1), com(2), com(3)
   DO i=1,nstar
     rstar(1:3,i)=rstar(1:3,i) - com(1:3)
     vstar(1:3,i)=vstar(1:3,i) - cov(1:3)
   END DO
!
   END SUBROUTINE makefractal
!
