subroutine fuvflux
  include 'globals.h'
  integer p,i,j,nterms,nedge,nnonedge,totalsearches
  real :: time1,time2,mintime,maxtime,tottime,utime1,utime2
  integer, parameter :: nbuf=32
  integer :: omp_get_max_threads
  integer :: maxthread, nchunk,k,imin,imax,ib, buf(nbuf),todo(nbuf),ntodo
  
  if(nstar.LE.0.OR.nsphact.LE.0) return
  if(.NOT.directsum.AND.nstar.GT.ntreemin) call startree

  nttotfuv=0; ntminfuv=nbodies; ntmaxfuv=0
  mintime=1.e10; maxtime=0.; tottime=0
  totalsearches=0

  maxthread=1
  nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	

!$omp parallel &
!$omp shared(root,nchunk) & 
!$omp private(p,nterms,time1,time2, &
!$omp imin,imax,ib,i,k,buf,todo,ntodo) &
!$omp reduction( + : nttotfuv,tottime,totalsearches) & 
!$omp reduction( MIN : ntminfuv, mintime) &
!$omp reduction( MAX : ntmaxfuv, maxtime)
  call cpu_time(time1)
  ncalls=0;nsearches=0		
!$omp do schedule(guided,1)
  do k=1,nchunk
    buf=0
    imin=((k-1)*nsphact)/nchunk+1
    imax=(k*nsphact)/nchunk
    reuseflag=1; searchreuse=0
    do i=imin,imax
      call prefuvwalk(i,imax,nbuf,buf,ntodo,todo)	   
      do ib=1,ntodo	 
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
  call cpu_time(time2)
  mintime=MIN(mintime,time2-time1)
  maxtime=MAX(maxtime,time2-time1)
  tottime=tottime+time2-time1
  totalsearches=totalsearches+nsearches
!$omp end parallel
  ntavgfuv=nttotfuv/nsphact
  if(verbosity.GT.0) then
    print*,'<fuvflux> searches', nsphact,totalsearches
    write(*,'(" <fuvflux> time:", 3f8.2)') maxtime,mintime
    print*,'<fuvflux> < a > t:',ntminfuv,ntavgfuv,ntmaxfuv,nttotfuv
  endif
end subroutine
				
subroutine zerofuv
  include 'globals.h'  
  fuvheat(pactive(1:nsphact))=0.
end subroutine
