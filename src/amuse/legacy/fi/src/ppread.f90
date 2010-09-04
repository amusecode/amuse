! tempr has to be fixed?
subroutine postprocessread
  include 'globals.h'
  real,allocatable :: tempr(:)

! setup initial conditions for simulation according to input
 
  if(verbosity.GT.0) print*,' ...transforming input...'

  if(nsph.EQ.0) usesph=.FALSE.
 
! move tempvect to tempr to be certain
  if(nsph.GT.0.AND.usesph) then
    allocate(tempr(nsph))
    tempr(1:nsph)=tempvect(1:nsph)
  endif
 
  if(input(1).EQ.0) call terror(' input error: need masses')

  if(input(42).EQ.0) then
    call nbexistinit
  endif

  if(nsph.GT.0) then
    if(massres.le.0.) then
!  massres=sum(mass(1:nsph))/8./nsph*nsmooth
!  massres=MINVAL(mass(1:nsph))/8.*nsmooth
      massres=MAXVAL(mass(1:nsph))/8.*nsmooth
    endif
  endif

  if(input(2).EQ.0) call terror(' input error: need positions')

  if(input(3).EQ.0) call terror(' input error: need velocities')
 
  if(input(4).EQ.0.OR.any(epsgrav(1:nbodies).EQ.0)) then
    call initeps ! stars + collisionless
    if(usesph.and.nsph.gt.0) then
      if(eps_is_h) then  
        if(input(13).EQ.0.OR.any(hsmooth(1:nsph).EQ.0)) then 
          call seteps('sph ')  ! not so good in fact .....
        else
          epsgrav(1:nsph)=hsmooth(1:nsph)
        endif
      else
        call initepsgas
      endif
    endif
  endif
  if(input(4).EQ.1.AND.eps.EQ.0) eps=sum(epsgrav(1:nbodies))/nbodies
  
  if(input(5).EQ.0) then
    tform(1:nsph)=tnow
    tform(nsph+1:nbodies)=tnow-2*tbubble
  endif

  tvel(1:nbodies)=tnow

  call activateparts

  call zeroacc
  call zeropot
  acc(1:nbodies,4)=0.
  call gravity('both')
  if(.not.directsum.and.gdgop) then
    call zeroacc
    call zeropot
    call gravity('both')
  endif
 
  if(input(34).EQ.0) call starevolv

  if(input(35).EQ.0) snentropy(nbodies-nstar+1:nbodies)=0

  if(usesph.AND.nsph.GT.0) then
    veltpos(1:nsph,1:ndim)=vel(1:nsph,1:ndim) ! cosmofac veltpos is pec. vel.
    if(input(14).EQ.0.AND.radiate) then
      call zerofuv
      call fuvflux
    endif
    call makesphtree
    if(input(13).EQ.0) then
      hsmooth(1:nsph)=0. 
    endif
    if(hupdatemethod.NE.'mass')then
      call terror('TBD mass update')
!   call inithsm ! reimplement
!   if(input(13).EQ.0.OR.sphinit) then
!    CALL omp_neighbor('correct') 
!   endif
!   if(input(10).EQ.0.OR.sphinit) then
!    CALL initdens ! inclding symm and no divv
!   endif    
!  if(input(24).EQ.0) call veldisp
    else
      if(input(10).EQ.0.OR.input(13).EQ.0.OR.input(21).EQ.0.OR. &
       input(23).EQ.0.OR.sphinit) then
        call densnhsmooth
      endif
    endif
! quick fix:
    call tree_reduction(root,incells,'sph ') 
 
    if(input(20).NE.0) csound(1:nsph)=sqrt(gamma*csound(1:nsph)/rho(1:nsph))

    if( input(19).NE.0 .OR. input(20).NE.0 .OR. &
       input(11).NE.0 .OR. input(12).NE.0 ) then
      if(uentropy) then
        if(input(12).EQ.0) then
          if(input(11).NE.0) then
            entropy(1:nsph)=gamma1*tempr(1:nsph)/rho(1:nsph)**gamma1
          else
            entropy(1:nsph)=csound(1:nsph)**2/gamma/rho(1:nsph)**gamma1
          endif   
        endif
        if(input(11).NE.0.OR.input(12).NE.0) &
         csound(1:nsph)=sqrt(gamma*rho(1:nsph)**gamma1*entropy(1:nsph)) !cosmof3
      else
        if(.not.isotherm) then
          if(input(11).EQ.0) then
            if(input(12).NE.0) then
              ethermal(1:nsph)=tempr(1:nsph)/gamma1*rho(1:nsph)**gamma1
            else
              ethermal(1:nsph)=csound(1:nsph)**2/gamma/gamma1
            endif
          endif
        if(input(11).NE.0.OR.input(12).NE.0) &
         csound(1:nsph)=sqrt(gamma*gamma1*ethermal(1:nsph))
        else
          if(input(19).EQ.0.AND.input(20).EQ.0) then
            if(input(11).NE.0) then
              csound(1:nsph)=sqrt(ethermal(1:nsph))
            else
              csound(1:nsph)=sqrt(tempr(1:nsph))
            endif   
          else
	    ethermal(1:nsph)=csound(1:nsph)**2
	  endif
        endif 
      endif
    else
      call terror('no u,A,P or Csound for gas')
    endif
   
    ethold(1:nsph)=ethermal(1:nsph) 
   
    if(smoothinput) then
      call terror('checksmoothall')
      call smoothall
    endif

    dethdt(1:nsph)=0
    derad(1:nsph)=0

!  IF(smoothinput.OR.(hupdatemethod.NE.'mass'.AND.(input(21).EQ.0.OR. &
!       input(23).EQ.0))) CALL initdivv ! check why initvet does not change vel
! if ok can be integrated with density

    if(input(22).EQ.0) mumaxdvh(1:nsph)=0.
  
    if((input(17).EQ.0.OR.input(18).EQ.0).AND.radiate) call temperature
    call ethdotinit
   
    if(input(15).EQ.0.OR.feedback.NE.'solh') esnthdt(1:nsph)=0 !check
    if(input(16).EQ.0) tcollaps(1:nsph)=tnow   ! should be improved(tnow-.5*tff)!
!  if(input(24).EQ.0) call veldisp checkcheck

    if(input(25).EQ.0.AND.radiate) then
      h2frac(1:nsph)=0.
      h2time=tnow-1.e8/timescale*year
      call molecules
    else
      h2time=tnow
    endif

! hsmooth for stars+bh is not read..
    hsmooth(nbodies-nstar+1:nbodies)=0

    if(nbh.GT.0) then
      tform(nbodies-nbh+1:nbodies)=tnow-2*tbubble
      call blackholes
    endif

    if(nbh+nstar.GT.0) CALL mech_feedback
    
  endif

  if(input(43).EQ.0) then
    itimestp(1:nbodies)=1
    call refresh_itimestp
  endif
   
  if(sortpart) call mortonsort

  if(allocated(tempr)) deallocate(tempr)
  
end subroutine
