 subroutine gravity(option)
  use pmgravMod
 include 'globals.h'
 character*4 option
 integer :: i,p

 call external_gravity(option)

 if(selfgrav) then
  if(adaptive_eps) call seteps('coll ')
 
  if(.NOT.directsum.AND.usepm) call pmaccgrav(option)
 
  if((.NOT.directsum).OR.adaptive_eps) call maketree
 
  call accgrav(option)
 endif
 
 if(fixthalo) call acchalo(option)

 if(option.NE.'pot ') THEN
  do i=1,npactive
    p=pactive(i)
    acc(p,4)=SQRT(acc(p,1)**2+acc(p,2)**2+acc(p,3)**2)
  enddo
 endif
 end subroutine

 subroutine accgrav(option)
  include 'globals.h'
  character*4 option
  integer :: p,i,j,nterms
  real :: lesoft,time1,time2,mintime,maxtime,tottime,utime1,utime2
  integer, parameter :: nbuf=32
  integer :: omp_get_max_threads,totalsearches
  integer :: maxthread, nchunk,k,imin,imax,ib, buf(nbuf),todo(nbuf),ntodo

  if(npactive.EQ.0) return
  nttot=0; ntmin=nbodies; ntmax=0
  mintime=1.e10; maxtime=0.; tottime=0

  maxthread=1
  nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,npactive/nbuf),maxthread)	
   totalsearches=0
!$omp parallel &
!$omp shared(root,nchunk) & 
!$omp private(p,nterms,lesoft,time1,time2, &
!$omp imin,imax,ib,i,k,buf,todo,ntodo) &
!$omp reduction( + : nttot, esofttot,tottime,totalsearches) & 
!$omp reduction( MIN : ntmin,mintime) &
!$omp reduction( MAX : ntmax,maxtime)
   call cpu_time(time1)
   ncalls=0;nsearches=0		
!$omp do schedule(guided,1)
   do k=1,nchunk
    buf=0
    imin=((k-1)*npactive)/nchunk+1
    imax=(k*npactive)/nchunk
    reuseflag=1; searchreuse=0
    do i=imin,imax
      call pretreewalk(i,imax,nbuf,buf,ntodo,todo)	   
      do ib=1,ntodo	 
       p=todo(ib)
       if(.NOT.directsum) then
        call pcond_treewalk(root,p,nterms)
       else
        nterms=nbodies
        do j=1,nbodies
         bodlist(j)=j
        enddo
       endif

       if(.not.periodic) then
        lesoft=0.
        call pgravsum(p,nterms,option,lesoft)
        if(npactive.EQ.nbodies) esofttot=esofttot+lesoft
       else
! in this case we will use pm gravity...
        call terror('periodic broken at the moment..')
       endif                 
       nttot=nttot+nterms
       ntmin=MIN(ntmin,nterms)
       ntmax=MAX(ntmax,nterms)
      enddo
    enddo  
   enddo	  
!$omp enddo nowait 
   call cpu_time(time2)
   mintime=MIN(mintime,time2-time1)
   maxtime=MAX(maxtime,time2-time1)
   tottime=tottime+time2-time1
   totalsearches=totalsearches+nsearches
!$omp end parallel
  ntavg=nttot/npactive
  if(verbosity.GT.0) then
   print*,'<accgrav> searches', npactive,totalsearches
   write(*,'(" <accgrav> time:", 2f8.2)') maxtime,mintime
   print*,'<accgrav> < a > t:',ntmin,ntavg,ntmax,nttot
  endif
 end subroutine

 subroutine pmaccgrav(option)
   use pmgravMod
  include 'globals.h'
  character*4 option
  integer :: p,i
  real :: lesoft,dtpm,vmax 
  real :: time1,time2,mintime,maxtime,tottime,utime1,utime2

  mintime=1.e10; maxtime=0.; tottime=0

! at the moment the box for pm grav is fixed, particles outside the 
! box are set to zero mass, and those are ignored in the pm CIC
! and in the acc. (hence they do not experience grav acc)
! this could be changed later to be more flexible)
 
   call cpu_time(utime1)
   if(tnow.GE.tpm) then 
    call pmgrav(nbodies,mass,pos)
    vmax=maxval(vel(1:nbodies,1:3))
    if(vmax.GT.0) dtpm=dmeshQ()/vmax/2.
    if(verbosity.GT.0) print*,'<pmaccgrav> timestep:',dtpm
    tpm=tpm+dtpm
   endif
   call cpu_time(utime2)
!$omp parallel  &
!$omp private(p,i,time1,time2) &
!$omp reduction( + : tottime) &
!$omp reduction( MIN : mintime) & 
!$omp reduction( MAX : maxtime)
   call cpu_time(time1)
!$omp do schedule(guided,200)
   do i=1,npactive
    p=pactive(i) 
    if(mass(p).GT.0) call pmgravsum(p,option)
   enddo
!$omp enddo nowait 
   call cpu_time(time2)
   mintime=MIN(mintime,time2-time1)
   maxtime=MAX(maxtime,time2-time1)
   tottime=tottime+time2-time1
!$omp end parallel

   if(verbosity.GT.0) then
    write(*,'(" <pmaccgrav> time:", 3f8.2)') utime2-utime1,maxtime,mintime
   endif
        
 end subroutine

subroutine acchalo(option)
 include 'globals.h'
 character*4 option
 real :: ppos(3),pacc(3),pphi
 integer :: i,p

 do i=1,npactive
  p=pactive(i)
  ppos(1:3)=pos(p,1:3)
  pacc(1:3)=0.
  pphi=0.
  call accpothalo(ppos,pacc,pphi,option)
  acc(p,1:3)=acc(p,1:3)+pacc
  phiext(p)=phiext(p)+pphi 
 enddo
end subroutine

subroutine system_gravity(peps,ppos,pphi,pacc,option,lesofttot)
 include 'globals.h'
 character*4 :: option
 real,intent(in) :: ppos(3),peps
 real,intent(inout) :: pphi,pacc(3),lesofttot
 integer i,nterms
! periodic - 
  if(.NOT.directsum) then
   nterms=0
   call treewalk(root,ppos,peps,0.,0.,nterms)
  else
   nterms=nbodies
   do i=1,nbodies
    bodlist(i)=i
   enddo
  endif 
  call gravsum(peps,ppos,pphi,pacc,nterms,option,lesofttot)
end subroutine

subroutine evaluate_gravity(peps,ppos,pphi,pacc,option)
 include 'globals.h'
 character*4 :: option
 real,intent(in) :: ppos(3),peps
 real,intent(inout) :: pphi,pacc(3)
 real :: dummy
 integer i
 
 pphi=0.
 pacc(1:3)=0
 if(selfgrav) call system_gravity(peps,ppos,pphi,pacc,option,dummy)
 if(fixthalo) call accpothalo(ppos,pacc,pphi,option)

end subroutine
