subroutine corrpos(ctimestp,rc)
  include 'globals.h'
  character*7 rc
  integer :: i,p,k,ctimestp(*)
  real dt2,rcsign,acceff,dt

  if(rc(1:4).NE.'sync'.AND.rc(1:6).NE.'desync') &
    call terror('unknown sync option in corrpos')
  if(rc(1:4).EQ.'sync') then
    if(syncflag.EQ.0) return
    syncflag=0
    rcsign=-1.
    if(verbosity.GT.0) print*,'<corrpos> sync'
  endif
  if(rc(1:6).EQ.'desync') then	
    if(syncflag.EQ.1) return
    syncflag=1
    rcsign=1.
    if(verbosity.GT.0) print*,'<corrpos> desync'
  endif

  if(.not.periodic) then
    do k=1,ndim
      do i=1,npactive
        p=pactive(i)
        dt2=(dtime/ctimestp(p))**2
        pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
      enddo
    enddo
  else            
    do k=1,ndim
      do i=1,npactive
        p=pactive(i)
        dt2=(dtime/ctimestp(p))**2
        pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
        if(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
        if(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
      enddo
    enddo
  endif
end subroutine

subroutine steppos
  include 'globals.h'
  integer i,ib,p,k,nkeep
  real acceff,distance,csdtime

  if(.not.periodic) then

    do k=1,ndim
      do p=1,nbodies
        pos(p,k)=pos(p,k)+vel(p,k)*tsteppos
      enddo
    enddo
  else
    do k=1,ndim
      do p=1,nbodies
        pos(p,k)=pos(p,k)+vel(p,k)*tsteppos
        if(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
        if(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
      enddo
    enddo  
  endif

  tpos=tpos+tsteppos
  tnow=tpos
end subroutine

subroutine stepvel
  include 'globals.h'
  integer p,k,i
  real acceff,vfac1now,vfac2now,dt,go,maxacc
! adhoc acc limiter
  k=0
  if(usesph) then
    do i=1,nsphact
      p=pactive(i)
      acceff=sqrt(sum(acc(p,1:3)**2))
      maxacc=tstepcrit**2*hsmooth(p)*(itimestp(p)/dtime)**2
      if(acceff.GT.maxacc) then
        k=k+1
        acc(p,1:3)=acc(p,1:3)/acceff*maxacc
      endif	
    enddo
  endif
  if(k.GT.0.AND.(verbosity.GT.0).OR.k.GT.0.001*nsph) &
    print*,' > acc limiter',tnow,k
!

  do k=1,ndim
    do i=1,npactive
      p=pactive(i)
      vel(p,k)=vel(p,k)+acc(p,k)*dtime/itimestp(p)
    enddo
  enddo  
  do i=1,npactive
    p=pactive(i)
    tvel(p)=tvel(p)+dtime/itimestp(p)
    if(ABS(tvel(p)-tnow).gt.dtime/itimestp(p)) then
      print*,p,tvel(p),tnow,itimestp(p),npactive
      call terror(' stepvel error: tvel mismatch') 
    endif
  enddo
end subroutine

subroutine zeroacc
  include 'globals.h'
  acc(pactive(1:npactive),1:3)=0.
end subroutine

subroutine zeropot
  include 'globals.h'
  phi(pactive(1:npactive))=0.
  phiext(pactive(1:npactive))=0.
  if(npactive.EQ.nbodies) esofttot=0.0
end subroutine

subroutine vextrap
  include 'globals.h'
  integer p,k
  do k=1,ndim
    do p=1,nsph
      veltpos(p,k)=vel(p,k)+acc(p,k)*(tnow-tvel(p))
    enddo
  enddo  
end subroutine

subroutine allethdot
  include 'globals.h'

  if(.NOT.consph) then
    call terror(' non conservative TBD')
!    if(sph_visc.EQ.'sphv') call omp_ethdotcv(pc)
!    if(sph_visc.EQ.'bulk') call omp_ethdotbv(pc)
!    if(sph_visc.EQ.'sph ') call omp_ethdot(pc)
  else
!    if(sph_visc.EQ.'sph ') call omp_ethdotco(pc)
    if(sph_visc.EQ.'bulk') call terror(' allethdot error')
    if(sph_visc.EQ.'sphv') call terror(' allethdot error')
  endif

end subroutine

subroutine allentdot
  include 'globals.h'

  if(sph_visc.EQ.'sph ') call omp_entdot
  if(sph_visc.EQ.'bulk') call terror(' allentdot error')
  if(sph_visc.EQ.'sphv') call terror(' allentdot error')	

end subroutine
 	 	 
subroutine allaccsph
  include 'globals.h'

  if(.NOT.consph) then
    call terror(' non conservative TBD')
!    if(sph_visc.EQ.'sphv') call omp_accsphcv
!    if(sph_visc.EQ.'bulk') call omp_accsphbv
!    if(sph_visc.EQ.'sph ') call omp_accsph
  else
    if(sph_visc.EQ.'sph ') call omp_accsphco
    if(sph_visc.EQ.'bulk') call terror(' allaccsph error')
    if(sph_visc.EQ.'sphv') call terror(' allaccsph error')
  endif

end subroutine

subroutine alldentacc
  include 'globals.h'

  if(.NOT.consph) then
    call terror(' non conservative TBD')
!    if(sph_visc.EQ.'sphv') call omp_accsphcv
!    if(sph_visc.EQ.'bulk') call omp_accsphbv
!    if(sph_visc.EQ.'sph ') call omp_accsph
  else
    if(sph_visc.EQ.'sph ') call omp_entdotaccsphco
    if(sph_visc.EQ.'bulk') call terror(' allaccsph error')
    if(sph_visc.EQ.'sphv') call terror(' allaccsph error')
  endif
 	 	 
end subroutine

subroutine alldethacc
  include 'globals.h'
  if(.NOT.consph) then
    call terror(' dethacc TBD ')
!    if(sph_visc.EQ.'sphv') call omp_accsphcv
!    if(sph_visc.EQ.'bulk') call omp_accsphbv
!    if(sph_visc.EQ.'sph ') call omp_accsph
  else
    if(sph_visc.EQ.'sph ') call omp_ethdotaccsphco
    if(sph_visc.EQ.'bulk') call terror(' allaccsph error')
    if(sph_visc.EQ.'sphv') call terror(' allaccsph error')
  endif	 

end subroutine

subroutine stepsph
  include 'globals.h'
  integer i,j,npnear,p

  call makesphtree

  if(hupdatemethod.NE.'mass') then
    call terror('TBD: hupdate=mass')
!    call omp_stepnear !change to hsmcal
!    call omp_density  !easychange 
!    call veldisp ! lesseasy
  else
    call densnhsmooth
  endif

! call update_reduction !!!!!! not so urgent
! quick fix: (not parallelized)
  call tree_reduction(root,incells,'sph ') 

  if(nbh.GT.0) call blackholes
  if(nstar+nbh.GT.0) call mech_feedback
        
  if(endstep.OR.nsphact.EQ.0) return

  if(.NOT.isotherm) then 
    if(.NOT.uentropy) then
      call allethdot
      call exstep2(1)
      call alldethacc
      call exstep2(2)
    else
      call allentdot   ! entdot is 2x faster than accsph 
      call exstep2(1)
      call alldentacc !  call allentdot
      call exstep2(2)
    endif
      call exstep2(3)
  else
    call alldethacc
  endif

end subroutine

subroutine step
  include 'globals.h'

  do
    call setrnd()

    if(starform.and.nsph.gt.0) then 
      if(verbosity.GT.0) print*,'<stepsys> starform...'
      call newstar
    endif

    if(nbh.gt.0) then
      if(verbosity.GT.0) print*,'<stepsys> bh mergers...'
      call bhmergers
    endif

    if(verbosity.GT.0) print*,'<stepsys> partremoval1...'
    call partremoval

    if(verbosity.GT.0) print*,'<stepsys> timestep..'
    call timestep

    if(verbosity.GT.0) print*,'<stepsys> steppos..'
    call steppos
    
! note that particles are not removed in clean, just set to zero mass
! - we have moved this from before partremoval. all routines should be
! 'zero-mass-safe'
    if(verbosity.GT.0) print*,'<stepsys> clean...'
    call clean

!    if(verbosity.GT.0) print*,'<stepsys> cosmo...'

!TBD    if(usesph.and.radiate) call cosmicray

    if(usesph) then
      if(verbosity.GT.0) print*,'<stepsys> extrapolate..'
      call vextrap
      call extrapeth
    endif
    
    if(npactive.gt.0.and.(.not.endstep)) then
      if(verbosity.GT.0) print*,'<stepsys> gravity..'
      call zeroacc
      call gravity('acc ')
    endif 

    if((starform.or.radiate)) then
      if(verbosity.GT.0) print*,'<stepsys> starevolv..'
      call starevolv
    endif
    
    if(usesph.and.nsphact.gt.0.and.(.not.endstep).and.radiate) then
      if(verbosity.GT.0) print*,'<stepsys> fuvflux..'
      call zerofuv
      call fuvflux
    endif
 
    if(usesph.and.radiate) then
      if(verbosity.GT.0) print*,'<stepsys> molecules...'
      call molecules
    endif
  
    if(usesph) then
      if(verbosity.GT.0) print*,'<stepsys> stepsph..'
      call stepsph
    endif 

    if(.not.endstep) then
      if(verbosity.GT.0) print*,'<stepsys> stepvel..'
!    if(sphfiltr.and.usesph.and.nsphact.gt.0) call vfilter tbd
      call stepvel
    endif
    
    if(endstep) exit
  enddo

  if(verbosity.GT.0) print*,'<stepsys> partremoval2...'
  call partremoval

  if(verbosity.GT.0) print*,'<stepsys> nbodies,nsph,nstar:',nbodies,nsph,nstar     
  
  if(sortpart) then
    if(verbosity.GT.0) print*,'<stepsys> sorting...'
    if(sortpart) call mortonsort
  endif
  
end subroutine

subroutine stepsystem(n)
  include 'globals.h'
  integer n
  
  call step
  print*,'<stepsystem> step completed:',n
  call outstate(n)

end subroutine
