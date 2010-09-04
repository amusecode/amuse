subroutine mech_feedback
 include 'globals.h'
 integer i
 real :: dt
 
 dt=tnow-tsnfeedback
  
 if(dt.lt.0) call terror(' feedback time error')
 
 select case (feedback)
  case('fuv ')
  case('pres')
   call pp_feedback(dt,tnow)
  case('kine')
   call randomvelfeedback(dt,tnow)
  case('solo')
   call terror('currently dangerous')
   call terror('solo not isotherm ')
   call energyfeedback(dt,tnow)
  case('solh')
   call heatingfeedback(dt,tnow)
  case default
   call terror(' feedback error 2')
 end select 
 
 tsnfeedback=tnow
 
 end subroutine mech_feedback
 
subroutine pp_feedback(dt,tnu)
 include 'globals.h'
 
   integer i,p,npp
   integer mythread,totalthread,low,up
!$ integer omp_get_thread_num,omp_get_num_threads 
   real tnu,dt
   real starpos(1:3),searchrad,sn_energy,starent
   real sn_activity,lsnheat,mc2fac,ts
   real, allocatable :: p_acc(:,:)
   real time1,time2,mintime,maxtime,tottime
   npp=0
   mc2fac=(lightspeed/velscale)**2
! test  
   ts=salpetertime*bh_rad_eff*year*1.e6/timescale
! stars
   do p=nbodies-nstar+1,nbodies   
    if(tnu-tform(p)-dt.le.tbubble.OR.snentropy(p).GT.0) then
     npp=npp+1
     templist(npp)=p
    endif
   enddo
! bh (tform=time of last mdot(=tfeedb) !=0 )
     
   if(verbosity.GT.0)print*,'<pressureparticle> particles:', npp  
   mintime=1.e10; maxtime=0.; tottime=0
   lsnheat=0.
!$omp  parallel &
!$omp private(p,i,sn_energy,starpos,searchrad, &
!$omp  mythread,totalthread,starent,low,up,time1,time2,p_acc) &
!$omp shared( npp,tnu,dt,lsnheat,ts) if(npp.gt.50) &
!$omp reduction( max : maxtime) &
!$omp reduction(+ : tottime) &
!$omp reduction( min : mintime)
        mythread=0
        totalthread=1
!$      mythread=omp_get_thread_num()
!$      totalthread=omp_get_num_threads()
  allocate(p_acc(nsph,3))
  call cpu_time(time1)
  p_acc(1:nsph,1:3)=0.
!$omp do reduction(+ : lsnheat)
 do i=1,npp   
  p=templist(i)
  if(p.LE.nbodies-nbh) then
   sn_energy=sne_eff*mass(p)*(sn_activity(tnu-tform(p))-sn_activity(tnu-tform(p)-dt))  
  else
! test  
   sn_energy=bh_rad_eff*bh_fdbk_eff*(tfeedb(p)-0.25*mass(p)/ts)*mc2fac*MIN(dt,tbubble)
  endif 
  lsnheat=lsnheat+sn_energy
  if(sn_energy.GT.0.OR.snentropy(p).GT.0) then
   starpos=pos(p,1:3)
   searchrad=hsmooth(p)
   starent=snentropy(p)
   if(searchrad.LT.0.OR.searchrad.GT.hboxsize) then
    print*,searchrad
    call terror('pres. searchrad')
   endif 
   call pressureparticle(starpos,starent,sn_energy,searchrad,p_acc)
   hsmooth(p)=searchrad
   snentropy(p)=starent
   if(snentropy(p).lt.0) snentropy(p)=0.
   if(tnu-tform(p).GT.tbubble) then
    snentropy(p)=0. 
    hsmooth(p)=0.
   endif 
  endif
 enddo
!$omp enddo nowait
 call cpu_time(time2)

if(nsphact.GT.0) then
  p_acc(1:pactive(1)-1,1:3)=0.
  do i=2,nsphact
   if(pactive(i-1).GE.pactive(i)) then
!$omp critical   
    call terror(' serious feedback error')
!$omp end critical
   endif
   p_acc(pactive(i-1)+1:pactive(i)-1,1:3)=0.
  enddo
 p_acc(pactive(nsphact)+1:nsph,1:3)=0.

 do i=0,totalthread-1
  low=mod(i+mythread,totalthread)*nsph/totalthread+1
  up=min(nsph,(mod(i+mythread,totalthread)+1)*nsph/totalthread)
  do p=low,up
   acc(p,1)=acc(p,1)+p_acc(p,1)
   acc(p,2)=acc(p,2)+p_acc(p,2)
   acc(p,3)=acc(p,3)+p_acc(p,3)
  enddo
!$omp barrier
 enddo
endif 
 mintime=MIN(mintime,time2-time1)
 maxtime=MAX(maxtime,time2-time1)
 tottime=tottime+time2-time1
 deallocate(p_acc)
!$omp end parallel
 snheat=snheat+lsnheat
 if(verbosity.GT.0) then
   write(*,'(" <pressureparticle> time:", 3f8.2)') maxtime,mintime,tottime
 endif
 
end subroutine pp_feedback


subroutine pressureparticle(spos,starent,dstarent,hsearch,p_acc)
 include 'globals.h'
 real,intent(in) :: spos(3),dstarent
 real, intent(inout) :: hsearch,starent,p_acc(nsph,3)
 real :: dens,ddens
 integer :: i,p,nneigh
  searchreuse=0
  nneigh=0
  call hsmdenspos2(spos,hsearch,dens,ddens,nneigh)
  if(nneigh.EQ.0.OR.dens.eq.0) return
  starent=starent+dstarent*gamma1/dens**gamma1
  if(starent.GT.0) call scatterfeedbackco(spos,hsearch,nneigh,dens,ddens,starent,p_acc)     
end

subroutine scatterfeedbackco(spos,hsearch,n,dens,ddens,starent,p_acc)
 include 'globals.h' 
 integer,intent(in) :: n
 real,intent(in) :: spos(3),dens,ddens,hsearch,starent
 real,intent(inout) :: p_acc(nsph,3)
 real :: rhog2,fj,dptot(3)
 real hsminv,dwnorm,distnorm,dx,dy,dz,dr2,dr2p,drw,dwsm,dwmass &
      , aij,drdp,mindrdp,wnorm,wsm
 integer i,p,iwsm,nearsph
 real mindist2
 integer, parameter :: mcorrect=1 ! 0=no, 1=smooth,2=closest,3=inproduct

 if(n.EQ.0.or.dens.le.0) then
!$omp critical
  call terror(' impossible feedback dens error')
!$omp end critical 
 endif

 rhog2=dens**(gamma-2)
 fj=1/(1+ddens*hsearch/3/(dens+rhomin))
 hsminv=1./hsearch
 wnorm=piinv*hsminv*hsminv*hsminv
 dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
 distnorm=hsminv**2*deldr2i
 dptot=0.
 
 mindist2=4*hsearch**2
 nearsph=0 
  
 do i=1,n
  p=srlist(i)
  tempvect(i)=0.
  dx=spos(1)-pos(p,1)
  dy=spos(2)-pos(p,2)
  dz=spos(3)-pos(p,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  dr2p=dr2*distnorm
  if(ninterp.GT.dr2p) then
  if(dr2.LE.mindist2.and.mass(p).GT.0) then
   nearsph=p
   mindist2=dr2
  endif 
  iwsm=INT(dr2p)
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm

  dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
  dwmass=dwnorm*dwsm
  
  aij=fj*starent*rhog2*dwmass
  p_acc(p,1)=p_acc(p,1)+aij*dx
  p_acc(p,2)=p_acc(p,2)+aij*dy
  p_acc(p,3)=p_acc(p,3)+aij*dz

  dptot(1)=dptot(1)+aij*mass(p)*dx
  dptot(2)=dptot(2)+aij*mass(p)*dy
  dptot(3)=dptot(3)+aij*mass(p)*dz
  
  endif
 enddo   
 
 if(mcorrect.EQ.3) then
  nearsph=0
  mindrdp=2*hsearch*sqrt(dptot(1)**2+dptot(2)**2+dptot(3)**2)
  do i=1,n
   p=srlist(i)
   dx=pos(p,1)-spos(1)
   dy=pos(p,2)-spos(2)
   dz=pos(p,3)-spos(3)
   if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
   if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
   if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
   if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
   if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
   if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
   drdp=dx*dptot(1)+dy*dptot(2)+dz*dptot(3)  
   if(drdp.LE.mindrdp.AND.tempvect(i).GT.0.AND.mass(p).GT.0) then
    nearsph=p
    mindrdp=drdp
   endif 
  enddo
 endif
 
 if(mcorrect.GT.1) then
  if(nearsph.eq.0) then 
!$omp critical
  call terror(' feedback oeps')
!$omp end critical  
  endif
  p_acc(nearsph,1)=p_acc(nearsph,1)-dptot(1)/mass(nearsph) 
  p_acc(nearsph,2)=p_acc(nearsph,2)-dptot(2)/mass(nearsph) 
  p_acc(nearsph,3)=p_acc(nearsph,3)-dptot(3)/mass(nearsph) 
 endif
 
 if(mcorrect.EQ.1) then
  do i=1,n
   p=srlist(i)
   p_acc(p,1)=p_acc(p,1)-dptot(1)*tempvect(i)/dens
   p_acc(p,2)=p_acc(p,2)-dptot(2)*tempvect(i)/dens
   p_acc(p,3)=p_acc(p,3)-dptot(3)*tempvect(i)/dens
  enddo
 endif
 
end subroutine
 
! ******** sn_activity ************************************************
! ***** returns cumulative sn energy/ unit mass as function of age ****
! (timestep independence)

function sn_activity(age)
 include 'globals.h'
 real, intent(in) :: age
 real :: sn_activity
 
 if(age.lt.tsnbeg) then
  sn_activity=0.
  return  
 endif 
 if(age.gt.tbubble) then
  sn_activity=snenergy
  return
 endif
 sn_activity=snenergy*(age-tsnbeg)/(tbubble-tsnbeg)
 
end function sn_activity

subroutine randomvelfeedback(dt,tnu)
 include 'globals.h'
 
   integer i,p,npp
   integer mythread,totalthread,low,up
!$ integer omp_get_thread_num,omp_get_num_threads 
   real tnu,dt
   real starpos(ndim), searchrad,sn_energy
   real sn_activity,lsnheat,mc2fac
   real, allocatable :: p_vel(:,:)
   
   if(nsphact.EQ.0.OR.dt.eq.0)return
   npp=0   
   mc2fac=(lightspeed/velscale)**2

   do p=nbodies-nstar+1,nbodies
    if(tnu-tform(p)-dt.LE.tbubble) then
      npp=npp+1
      templist(npp)=p
    endif
   enddo
   if(verbosity.GT.0) print*,'<randomvelfeedback> sn particles:', npp   
   
   lsnheat=0.
!$omp parallel &
!$omp private(p,i,sn_energy,starpos,searchrad,mythread,totalthread,low,up,p_vel) &
!$omp shared( npp,tnu,dt,lsnheat) if(npp.gt.50)
        mythread=0
        totalthread=1
!$      mythread=omp_get_thread_num()
!$      totalthread=omp_get_num_threads()
  allocate(p_vel(1:nsph,3))
  p_vel(1:nsph,1:3)=0.

!$omp do reduction(+ : lsnheat)
 do i=1,npp   
  p=templist(i)
! note: sn_energy is really a rate here !!
! (to be able to take account of individual timesteps)  
  if(p.LE.nbodies-nbh) then
   sn_energy=sne_eff*mass(p)*(sn_activity(tnu-tform(p))-sn_activity(tnu-tform(p)-dt))/dt  
  else
   sn_energy=bh_rad_eff*bh_fdbk_eff*tfeedb(p)*mc2fac
  endif 

  lsnheat=lsnheat+sn_energy*dt
  if(sn_energy.GT.0) then
   starpos=pos(p,1:3)
   searchrad=hsmooth(p)
   if(searchrad.LT.0.OR.searchrad.GT.hboxsize) then
    print*,searchrad
    call terror('kine searchrad')
   endif    
   call randomvelo(starpos,sn_energy,searchrad,p_vel,p)
   hsmooth(p)=searchrad
  endif
 enddo
!$omp enddo nowait

! p_vel is calculated with sn_energy rate, scales with sqrt(t)
! (maintain independence of timestep)
 if(nsphact.GT.0) then
  p_vel(1:pactive(1)-1,1:3)=0.
  p_vel(pactive(1),1:3)=p_vel(pactive(1),1:3)*SQRT(dtime/2**(itimestp(pactive(1))-1))
  do i=2,nsphact
   do p=pactive(i-1)+1,pactive(i)-1
    p_vel(p,1:3)=0.
   enddo
  p_vel(pactive(i),1:3)=p_vel(pactive(i),1:3)*SQRT(dtime/2**(itimestp(pactive(i))-1))
  enddo
  p_vel(pactive(nsphact)+1:nsph,1:3)=0.

 do i=0,totalthread-1
  low=mod(i+mythread,totalthread)*nsph/totalthread+1
  up=min(nsph,(mod(i+mythread,totalthread)+1)*nsph/totalthread)
  do p=low,up
   vel(p,1)=vel(p,1)+p_vel(p,1)
   vel(p,2)=vel(p,2)+p_vel(p,2)
   vel(p,3)=vel(p,3)+p_vel(p,3)
  enddo
!$omp barrier
 enddo
 endif
 deallocate(p_vel)
!$omp end parallel
 snheat=snheat+lsnheat

end subroutine

subroutine randomvelo(spos,energy,hsearch,p_vel,pid)
 include 'globals.h'
 real,intent(in) :: spos(3),energy
 real,intent(inout) :: hsearch,p_vel(nsph,3)
 real :: dens,ddens
 integer :: nneigh,i,p,pid,j
 real :: dx,dy,dz,dvx,dvy,px,py,pz,pRND
 real :: hsminv,wnorm,distnorm, dr2p,wsm,dr2,drw
 integer :: iwsm 
 
  searchreuse=0
  nneigh=0
  call hsmdenspos2(spos,hsearch,dens,ddens,nneigh)
  if(dens.EQ.0.OR.nneigh.EQ.0) return
  px=0;py=0;pz=0
  j=0
 hsminv=1./hsearch
 wnorm=piinv*hsminv*hsminv*hsminv
 distnorm=hsminv**2*deldr2i  
 do i=1,nneigh
  p=srlist(i)
  tempvect(i)=0.
  dx=spos(1)-pos(p,1)
  dy=spos(2)-pos(p,2)
  dz=spos(3)-pos(p,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  dr2p=dr2*distnorm
  if(ninterp.GT.dr2p) then
  iwsm=INT(dr2p)
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm
  endif
 enddo   
  do i=1,nneigh
   p=srlist(i)   
100     continue
         dx=2*pRND(p+3*(pid+j))-1.
         dy=2*pRND(p+3*(pid+j)+1)-1.
         dz=2*pRND(p+3*(pid+j)+2)-1.
         dvx=dx**2+dy**2+dz**2
         j=j+1
	 if(dvx.GE.1.) goto 100 
    dvy=4./3.*SQRT(2*energy*tempvect(i)/dens*sne_eff)
    p_vel(p,1)=p_vel(p,1)+dx*dvy
    p_vel(p,2)=p_vel(p,2)+dy*dvy
    p_vel(p,3)=p_vel(p,3)+dz*dvy
    px=px+dx*dvy*mass(p)
    py=py+dy*dvy*mass(p)
    pz=pz+dz*dvy*mass(p)
  enddo

  do i=1,nneigh
   p=srlist(i)
    p_vel(p,1)=p_vel(p,1)-px*tempvect(i)/dens ! tempvect already divided by mass
    p_vel(p,2)=p_vel(p,2)-py*tempvect(i)/dens
    p_vel(p,3)=p_vel(p,3)-pz*tempvect(i)/dens
  enddo
end

subroutine energyfeedback(dt,tnu)
 include 'globals.h'
 
   integer i,p,npp
   integer mythread,totalthread,low,up
!$ integer omp_get_thread_num,omp_get_num_threads 
   real tnu,dt
   real starpos(ndim), searchrad,sn_energy
   real sn_activity,ethtoent,lsnheat,mc2fac
   real, allocatable :: p_eth(:)
   
   npp=0
   mc2fac=(lightspeed/velscale)**2

   do p=nbodies-nstar+1,nbodies
    if(tnu-tform(p)-dt.LE.tbubble) then
     npp=npp+1
     templist(npp)=p
    endif
   enddo
   lsnheat=0.
!$omp parallel &
!$omp private(p,i,sn_energy,starpos,searchrad,mythread,totalthread,low,up,p_eth) &
!$omp shared( npp,tnu,dt,lsnheat) if(npp.gt.50)
        mythread=0
        totalthread=1
!$      mythread=omp_get_thread_num()
!$      totalthread=omp_get_num_threads()
  allocate(p_eth(nsph))
  p_eth(1:nsph)=0.

!$omp do reduction(+ : lsnheat)
 do i=1,npp   
  p=templist(i)
  if(p.LE.nbodies-nbh) then
   sn_energy=sne_eff*mass(p)*(sn_activity(tnu-tform(p))-sn_activity(tnu-tform(p)-dt))  
  else
   sn_energy=bh_rad_eff*bh_fdbk_eff*tfeedb(p)*mc2fac*MIN(dt,tbubble)
  endif 
  
  lsnheat=lsnheat+sn_energy
  if(sn_energy.GT.0) then
   starpos=pos(p,1:3)
   searchrad=epsgrav(p)
   call energyreturn(starpos,sn_energy,searchrad,p_eth)
  endif
 enddo
!$omp enddo nowait

 if(nsphact.GT.0) then
  p_eth(1:pactive(1)-1)=0.
  do i=2,nsphact
   do p=pactive(i-1)+1,pactive(i)-1
    p_eth(p)=0.
   enddo
  enddo
  p_eth(pactive(nsphact)+1:nsph)=0.

  do i=0,totalthread-1
  low=mod(i+mythread,totalthread)*nsph/totalthread+1
  up=min(nsph,(mod(i+mythread,totalthread)+1)*nsph/totalthread)
  do p=low,up
   ethermal(p)=ethermal(p)+p_eth(p)
   ethold(p)=ethold(p)+p_eth(p)
  enddo
!$omp barrier
  enddo
 endif
 deallocate(p_eth)
!$omp end parallel
 snheat=snheat+lsnheat

 if(uentropy) then
  do i=1,nsphact
   p=pactive(i)
   ethtoent=gamma1/rho(p)**gamma1
   csound(p)=SQRT(gamma*gamma1*ethermal(p)/ethtoent)
  enddo
 else
  do i=1,nsphact
   p=pactive(i)
   ethtoent=1
   csound(p)=SQRT(gamma*gamma1*ethermal(p)/ethtoent)
  enddo
 endif
 
end subroutine

subroutine energyreturn(spos,energy,hsearch,p_eth)
 include 'globals.h'
 real,intent(in) :: spos(3),energy
 real,intent(inout) :: hsearch,p_eth(nsph)
 integer :: nneigh,i,p,nearsph
 real :: dist2,dist2a
 
 nneigh=0
 do while(nneigh.eq.0)
 hsearch=2*hsearch
 call search_d2(root,hsearch,spos,nneigh,srlist,tempvect)
 enddo
 if(nneigh.gt.0) then
  dist2=hsearch**2
  nearsph=0
  do i=1,nneigh
      dist2=tempvect(i)
      p=srlist(i)
      if(dist2.gt.dist2a.and.mass(p).GT.0) then
        dist2=dist2a
        nearsph=p
      endif 
  enddo
  if(nearsph.eq.0) call terror(' weird solo error')
  if(uentropy) then
   p_eth(nearsph)=p_eth(nearsph)+energy/mass(nearsph)*gamma1/rho(nearsph)**gamma1
  else
   p_eth(nearsph)=p_eth(nearsph)+energy/mass(nearsph)  
  endif
 else
  print*,'*** feedback does not find particle ??'
 endif 
end

subroutine heatingfeedback(dt,tnu)
 include 'globals.h'
 
   integer i,p,npp
   integer mythread,totalthread,low,up
!$ integer omp_get_thread_num,omp_get_num_threads 
   real tnu,dt
   real starpos(ndim), searchrad,sn_energy
   real sn_activity,ethtoent,lsnheat,mc2fac
   real, allocatable :: pesn(:)
   
   npp=0
   mc2fac=(lightspeed/velscale)**2
   if(nsphact.GT.0.OR.dt.eq.0)return

! be careful: 
   do i=1,nsphact
    esnthdt(pactive(i))=0.
   enddo
   
   do p=nbodies-nstar+1,nbodies
    if(tnu-tform(p)-dt.LE.tbubble) then
     npp=npp+1
     templist(npp)=p
    endif
   enddo
   lsnheat=0.
!$omp parallel &
!$omp private(p,i,sn_energy,starpos,searchrad,mythread,totalthread,low,up,pesn) &
!$omp shared( npp,tnu,dt,lsnheat) if(npp.gt.50)
        mythread=0
        totalthread=1
!$      mythread=omp_get_thread_num()
!$      totalthread=omp_get_num_threads()
  allocate(pesn(nsph))
  pesn(1:nsph)=0.

!$omp do reduction(+ : lsnheat)
 do i=1,npp   
  p=templist(i)
  if(p.LE.nbodies-nbh) then
   sn_energy=sne_eff*mass(p)*(sn_activity(tnu-tform(p))-sn_activity(tnu-tform(p)-dt))  
  else
   sn_energy=bh_rad_eff*bh_fdbk_eff*tfeedb(p)*mc2fac*MIN(dt,tbubble)
  endif 
  lsnheat=lsnheat+sn_energy*dt
  sn_energy=sn_energy/dt
  if(sn_energy.GT.0) then
   starpos=pos(p,1:3)
   searchrad=epsgrav(p)
   call singlepheating(starpos,sn_energy,searchrad,pesn)
  endif
 enddo
!$omp enddo nowait

 if(nsphact.GT.0) then
  pesn(1:pactive(1)-1)=0.
  do i=2,nsphact
   do p=pactive(i-1)+1,pactive(i)-1
    pesn(p)=0.
   enddo
  enddo
  pesn(pactive(nsphact)+1:nsph)=0.

  do i=0,totalthread-1
  low=mod(i+mythread,totalthread)*nsph/totalthread+1
  up=min(nsph,(mod(i+mythread,totalthread)+1)*nsph/totalthread)
  do p=low,up
   esnthdt(p)=esnthdt(p)+pesn(p)
  enddo
!$omp barrier
  enddo
 endif
 deallocate(pesn)
!$omp end parallel
! snheat=snheat+lsnheat
 
end subroutine

subroutine singlepheating(spos,energy,hsearch,pesn)
 include 'globals.h'
 real,intent(in) :: spos(3),energy
 real,intent(inout) :: hsearch,pesn(nsph)
 integer :: nneigh,i,p,nearsph
 real :: dist2,dist2a
 
 nneigh=0
 do while(nneigh.eq.0)
 hsearch=2*hsearch
 call search_d2(root,hsearch,spos,nneigh,srlist,tempvect)
 enddo
 if(nneigh.gt.0) then
  dist2=hsearch**2
  nearsph=0
  do i=1,nneigh
      dist2=tempvect(i)
      p=srlist(i)
      if(dist2.gt.dist2a.and.mass(p).GT.0) then
        dist2=dist2a
        nearsph=p
      endif 
  enddo
  if(nearsph.eq.0) call terror(' weird solh error')
  pesn(nearsph)=pesn(nearsph)+energy/mass(nearsph)  
 else
  print*,'*** feedback does not find particle ??'
 endif 
end

