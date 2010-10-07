subroutine report(level,message)
  include 'globals.h'
  integer level
  character(len=*) message
  
  if(verbosity.GT.level) print*, message 

end subroutine


subroutine diagnostics
  include 'globals.h'
  integer :: i,p,k,nblist,nbwarm
  logical :: firstc
  real :: temp

  mgas=0. 
  mstar=0.
  mtot=0.
  ektot=0.
  eptot=0.
  ethtot=0.

  mcold=0.
  mluke=0.
  mwarm=0.
  mhot=0.

  cmpos(1:3)=0.
  cmvel(1:3)=0.
  amvec(1:3)=0.
  mgas=sum(mass(1:nsph))
  mstar=sum(mass(nbodies-nstar+1:nbodies))

  do p=1,nbodies
    mtot=mtot+mass(p)
    eptot=eptot+mass(p)*phiext(p)
  enddo

  if(selfgrav) then
    if(periodic) call terror('check pot. for periodic')
    do p=1,nbodies
      eptot=eptot+.5*mass(p)*phi(p)
    enddo
  endif  

  if(.NOT.isotherm) then
    do p=1,nsph
      ethtot=ethtot+mass(p)*csound(p)**2/(gamma*gamma1)
    enddo
  else
    do p=1,nsph
      ethtot=ethtot+mass(p)*csound(p)**2
    enddo
  endif
  
  do k=1,3
    do p=1,nbodies
      ektot=ektot+.5*mass(p)*vel(p,k)*vel(p,k)
      cmpos(k)=cmpos(k)+mass(p)*pos(p,k)
      cmvel(k)=cmvel(k)+mass(p)*vel(p,k)
    enddo
    cmvel(k)=cmvel(k)/mtot
    cmpos(k)=cmpos(k)/mtot
  enddo
  etot=ektot+eptot+ethtot
  if(ndim.EQ.2) then
    do p=1,nbodies
      amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-pos(p,2)*vel(p,1))
    enddo
  endif
  if(ndim.EQ.3) then
    do p=1,nbodies
      amvec(1)=amvec(1)+mass(p)*(pos(p,2)*vel(p,3)-pos(p,3)*vel(p,2))
      amvec(2)=amvec(2)+mass(p)*(pos(p,3)*vel(p,1)-pos(p,1)*vel(p,3))
      amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-pos(p,2)*vel(p,1))
    enddo
  endif

  do p=1,nsph
    temp=meanmwt*mhboltz*csound(p)**2/gamma
    if(temp.LT.1000) mcold=mcold+mass(p)
    if(temp.GE.3.e4) mhot=mhot+mass(p)
    if(temp.LT.3.e4.AND.temp.GE.8000) mwarm=mwarm+mass(p)
  enddo
  mluke=mgas-mcold-mwarm-mhot
  if(verbosity.GT.0) print*,'<energy>',tnow,etot
end subroutine

subroutine outstate(n)
  include 'globals.h'
  logical :: testcrit
  integer :: n,ioerror,sstatus
  real :: dumpt,tend,tstep,tcurrent
  real :: rtime
  character(len=4) sstring
  real, save :: tlastdump=0.,tprevious=0.
  
  tcurrent=rtime()
  if(tcurrent.lt.tprevious) then
    print*,' timer error (overflow)'
    stop
  endif          
  tstep=tcurrent-tprevious
  tprevious=tcurrent

  sstatus=-1
  open(unit=upars,file='stop',status='OLD',iostat=ioerror)
  if(ioerror.EQ.0) then
    read(upars,'(a)',iostat=ioerror) sstring
    if(ioerror.EQ.0) then
      if(sstring.EQ.'exit') sstatus=3
      if(sstring.EQ.'stop') sstatus=2
      if(sstring.EQ.'snap') sstatus=1
      if(sstring.EQ.'dump') then
        sstatus=0
        read(upars,*,iostat=ioerror) dumpt
        if(ioerror.NE.0) then
          dumpt=12.
          print*,'...assuming tdump of 12 hrs...'
        endif 
        read(upars,*,iostat=ioerror) tend
        if(ioerror.NE.0) tend=8760.
      endif 
    endif
    if(sstatus.GT.0) then
      rewind(upars)
      write(upars,'(a)') 'ok..',sstring
    endif
    close(upars)
  endif

  if(sstatus.EQ.0.AND.tcurrent.GT.tend-2*tstep) sstatus=3

  if(sstatus.EQ.0) then
    if(tcurrent-tlastdump.GT.dumpt-2*tstep) then
      call writedump(n)
      tlastdump=tcurrent
    endif
  endif

  if(sstatus.GE.2) call writedump(n)

  if(MOD(n,steplog).EQ.0.OR.MOD(n,stepout).EQ.0.OR. &
        sstatus.EQ.1.OR.sstatus.EQ.2) then

    call activateparts

    testcrit=MOD(n,stepout).EQ.0.OR.n.EQ.0.OR.sstatus.EQ.1.OR.sstatus.EQ.2
    call corrpos(itimestp,'sync')
    if(n.NE.0) then
      call zeropot
      call gravity('pot ')
    endif

    if(MOD(n,steplog).EQ.0) then 
      call outlog(n)
      call diagnostics
      call outenrgy
    endif
    
    if(testcrit) call outbods('XXXXXX')
    call corrpos(itimestp,'desync')
  endif
   
  if(sstatus.GE.2) then
    print*,' stopped as requested'
    stop
  endif
end subroutine
