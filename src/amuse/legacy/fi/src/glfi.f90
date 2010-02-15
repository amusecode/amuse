program glfi
  include 'globals.h'
  integer omp_get_max_threads
  real time,rtime
  integer n1,n,ioerror
  
  print*,' --- GLFI ---'
  call initmem(nbodsmax,nsphmax,ncells)
!$  print*,' omp threads:',omp_get_max_threads()    
  time=rtime()
  print*,' '
  print*,' compiled with:'
  if(periodic) then
    print*,' > periodic boundary'
  else
    print*,' > vacuum boundary'
  endif
  if(sortpart) then
    print*,' > particle sorting'
  else
    print*,' > no particle sorting'
  endif
  if(H2cooling) then
    print*,' > h2 cooling'
  else
    print*,' > no h2 cooling'
  endif
  if(ionCool) then
    print*,' > variable ionization'
  else
    print*,' > fixed ionization'
  endif
  if(simpleeff) then
    print*,' > simple heating efficiency'
  else
    print*,' > full heating efficiency'
  endif
  if(cosmrayheat) then
    print*,' > cosmic ray heating'
  else
    print*,' > no cosmic ray heating'
  endif
  if(usepm) then
    print*,' > pm gravity'
  else
    print*,' > no pm gravity'
  endif
  print*,' > compiled constraints:'
  print*,'    no. of particles: ',nbodsmax
  print*,'    no. of sph part:  ',nsphmax
  print*,'    no. of cells:     ',ncells
  print*,' '
  open(unit=uboddump,file=dumpfile,status='OLD',iostat=ioerror)
  close(uboddump)
  if(ioerror.NE.0) then 
    call initsys
    n1=0
    print *,' ...starting simulation...'
  else
    call initdump(n1)
    print *,' ...restarting simulation...'
  endif 
  if(verbosity.GT.0) then
    print*,'  > verbosity level:', verbosity
  endif
  print*,' ...NO maps...'
  print *,' ...start loop...'
  call viewer
  do n=n1+1,nsteps
    call stepsystem(n)
  enddo
  print*,'*** GAME OVER ***'
  call viewer
end program
