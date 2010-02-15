subroutine outhead
  include 'globals.h'
  character(len=4) type
  character(len=8) halocore,halodens
  character(len=10) halotype
  integer id,ierr
        
  open(UNIT=ulog,FILE=logfile,STATUS='OLD',position='APPEND',IOSTAT=ierr)
  if(ierr.NE.0) then
    if(verbosity.GT.0) print*,' logfile not found'
    return
  endif 

  write(ulog,'("initial nbodies, nsteps: ",4i9)') nbodies,nsteps
  write(ulog,'("output-interval bodies, energy: ",4i9)') stepout,steplog
  write(ulog,*) 'see parameter file'
  write(ulog,*)
  close(ulog)

end subroutine

subroutine startout
  include 'globals.h'
  integer ioerror
  character(len=4) sstring
  
  open(unit=upars,file='stop',status='UNKNOWN')
  read(upars,'(a)',iostat=ioerror) sstring 
  if(ioerror.EQ.0.AND.sstring.EQ.'stop') then
    print*,'stop -- stop file already present'
    print*,'     -- please remove or alter ' 
    stop
  endif 
  close(upars)

  print*,' ...making log...'
  open(UNIT=ulog,FILE=logfile,STATUS='NEW',IOSTAT=ioerror)
  if(ioerror.NE.0) then
    print*,'stop -- cannot make new log file: ',logfile
    stop
  endif
  write(ulog,*) ' Start of logfile output '
  close(ulog)
 
  open(unit=sfrfilenr,file=TRIM(outputfile)//'.sfh',status='unknown')
  write(sfrfilenr,*) ' star formation history '
  close(sfrfilenr)

  open(unit=erasenr,file=TRIM(outputfile)//'.del',status='unknown')
  write(erasenr,*) ' erased mass '
  write(erasenr,*) ' time, mass, com pos, com vel'
  close(erasenr)

  open(unit=bhfilenr,file=TRIM(outputfile)//'.bh',status='unknown')
  write(bhfilenr,*) ' blackhole accretion history'
  write(bhfilenr,'(a)',advance='NO') & 
    '1 origindex,2 time,3 reff,4 rho,5 csound,6 vbulk,7 drhodt,' 
  write(bhfilenr,'(a)',advance='NO') &
    '8 temp,9 mdot_bondi,10 mdot_divv,11 mdot_edd,12 mdot,'
  write(bhfilenr,'(a)') &
    '13 bhmass,14 h,15 x,16 y,17 z,18 drhodh,19 vx,20 vy,21 vz'
  close(bhfilenr)

end subroutine

subroutine outenrgy
  include 'globals.h'
  integer :: k,ioerr
  real :: cpunew,cpustep,etotrad,rtime
  real, save :: cpuold=0.0

  open(UNIT=ulog,FILE=logfile,STATUS='OLD',position='APPEND',IOSTAT=ioerr)
  if(ioerr.NE.0) then
    if(verbosity.GT.0) print*,' no log file'
    return
  endif  
  
  etotrad=etot-eradiate-snheat-enerror

  if(.NOT.starform) then
    write(ulog,'("mtot,mgas = ",2(1pe17.9))') mtot,mgas
  else
    write(ulog,'("mtot,mgas,mstar = ",3(1pe17.9))') mtot,mgas,mstar
    write(ulog,'("cold,luke,warm,hot = ",4(1pe17.9))') mcold,mluke,mwarm,mhot
  endif
  write(ulog,'("e, e_kin, e_pot = ",3(1pe17.9))') etotrad,ektot,eptot
  write(ulog,'("e_th, e_soft = ",2(1pe17.9))') ethtot,esofttot

  if(radiate) then
    write(ulog,'("e_radiated = ",1(1pe17.9))') eradiate
    write(ulog,'("e_fuv, e_cool, e_star = ",3(1pe17.9))') &
      efuvheat,eradcool,estar
  endif

  if(starform) then
    write(ulog,'("e_snheat = ",1(1pe17.9))')  snheat
  endif

  write(ulog,'("am_x, am_y, am_z = ",3(1pe17.9))') amvec(1),amvec(2),amvec(3)
  write(ulog,'("x_cm, y_cm, z_cm = ",3(1pe17.9))') cmpos(1),cmpos(2),cmpos(3)
  write(ulog,'("vx_cm, vy_cm, vz_cm = ",3(1pe17.9))') cmvel(1),cmvel(2),cmvel(3)

  cpunew=rtime()
  cpustep=cpunew-cpuold
  cpuold=cpunew
  write(ulog,'("wall clock time = ",1(1pe17.9))') cpustep
  write(ulog,*)
  close(ulog)

end subroutine

subroutine outlog(istep)
  include 'globals.h'
  integer :: istep,i,ibin,tstepbin(20),p,ioerr
  real :: tistep

  if(istep.EQ.0) call outhead
  open(UNIT=ulog,FILE=logfile,STATUS='OLD',position='APPEND',IOSTAT=ioerr)
  if(ioerr.NE.0) then
    if(verbosity.GT.0) print*,' no log file'
    return
  endif
	 
  if(cosmo) then
    write(ulog,'("time, a, Z, ncells = ",3(1pe17.9), i9)') tnow,1.0,0.,incellsg
  else
    write(ulog,'("time, ncells = ",(1pe17.9),i9)') tnow,incellsg
  endif
  write(ulog,'("nttot,ntmin,ntmax,ntavg = ",4i9)') nttot,ntmin,ntmax,ntavg
  write(ulog,'("nntot,ntmin,ntmax,ntavg = ",4i9)') nntot,nnmin,nnmax,nnavg
  if(adaptive_eps) &
    write(ulog,'("nstot,ntmin,ntmax,ntavg = ",4i9)') nstot,nsmin,nsmax,nsavg
  if(usesph) then
    tstepbin(1:20)=0
    do p=1,nsph
      tistep=itimestp(p)
      tistep=0.5+LOG(tistep)/LOG(2.)
      ibin=INT(tistep)+1
      tstepbin(ibin)=tstepbin(ibin)+1
    enddo
    write(ulog,'("timestep distribution SPH = ",20i9)') tstepbin(1:20)
  endif
  tstepbin(1:20)=0
  do p=nsph+1,nbodies
    tistep=itimestp(p)
    tistep=0.5+LOG(tistep)/LOG(2.)
    ibin=INT(tistep)+1
    tstepbin(ibin)=tstepbin(ibin)+1
  enddo
  write(ulog,'("timestep distribution non-SPH = ",20i9)') tstepbin(1:20)
  write(ulog,*)
  close(ulog)

end subroutine
