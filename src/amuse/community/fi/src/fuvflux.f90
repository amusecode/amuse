subroutine fuvflux
  include 'globals.h'
  integer p,i,j,nterms,nedge,nnonedge,totalsearches
  real :: time1,time2,mintime,maxtime,tottime,utime1,utime2
  integer, parameter :: nbuf=32
  integer :: omp_get_max_threads
  integer :: maxthread, nchunk,k,imin,imax,ib, buf(nbuf),todo(nbuf),ntodo
  integer :: niter
  
  if(nstar.LE.0.OR.nsphact.LE.0) return
  if(.NOT.directsum.AND.nstar.GT.ntreemin) call startree

  nttotfuv=0; ntminfuv=nbodies; ntmaxfuv=0
  mintime=1.e10; maxtime=0.; tottime=0
  totalsearches=0

  maxthread=1
  nchunk=1
  niter=0
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	

!$omp parallel &
!$omp shared(root,nchunk) & 
!$omp private(p,nterms,time1,time2, &
!$omp imin,imax,ib,i,k,buf,todo,ntodo) &
!$omp reduction( + : nttotfuv,tottime,totalsearches,niter) & 
!$omp reduction( MIN : ntminfuv, mintime) &
!$omp reduction( MAX : ntmaxfuv, maxtime)
  call wall_time(time1)
  ncalls=0;nsearches=0		
!$omp do schedule(guided,1)
  do k=1,nchunk
    buf=0
    imin=int(nsphact*float(k-1)/nchunk)+1
    imax=int(nsphact*float(k)/nchunk)
    reuseflag=1; searchreuse=0
    do i=imin,imax
      call prefuvwalk(i,imax,nbuf,buf,ntodo,todo)	   
      do ib=1,ntodo
        niter=niter+1
        p=todo(ib)
        if(.NOT.directsum.AND.nstar.GT.ntreemin) then
          call pcond_fuvwalk(root,p,nterms)
        else
          nterms=nstar
          do j=1,nstar
            bodlist(j)=j+nbodies-nstar
          enddo
        endif

        if(.not.periodic) then
          call pfuvsum(p,nterms)
        else
! check that periodic treewalk gives only nearest image
          call terror('periodic broken at the moment..')
        endif                 
        nttotfuv=nttotfuv+nterms
        ntminfuv=MIN(ntminfuv,nterms)
        ntmaxfuv=MAX(ntmaxfuv,nterms)
      enddo
    enddo  
  enddo	  
!$omp enddo nowait 
  call wall_time(time2)
  mintime=MIN(mintime,time2-time1)
  maxtime=MAX(maxtime,time2-time1)
  tottime=tottime+time2-time1
  totalsearches=totalsearches+nsearches
!$omp end parallel
  ntavgfuv=nttotfuv/nsphact
  if(verbosity.GT.0) then
    print*,'<fuvflux> parts,searches:', nsphact,totalsearches
    write(*,'(" <fuvflux> time:", 3f8.2)') maxtime,mintime,tottime
    print*,'<fuvflux> mn,av,mx:',ntminfuv,ntavgfuv,ntmaxfuv
  endif
  if(niter.NE.nsphact) call terror("fuvflux inconsistent iter count")
end subroutine
				
subroutine zerofuv
  include 'globals.h'  
  fuvheat(pactive(1:nsphact))=0.
end subroutine
