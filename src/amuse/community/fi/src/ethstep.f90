subroutine exstepp(p,totaldt,eth,deth,drad,lerad,lhe,imax,jmax)
 include 'globals.h'
 integer p,i,j,jmax,imax
 real eth,totaldt,deth,drad
 real eth0,eth1,eth2,deth0,deth1,deth2
 real tpart,dt,he,co,ethmax,ethmin
 logical tnowpassed,laststep
 real ttolerance,tmin,tmax,safefac
 real lerad,lhe
 parameter( ttolerance=0.01,safefac=.9,tmin=10.,tmax=1.e9)
 
 if(totaldt.LE.0) return
 i=0  
 jmax=0
 lhe=0
 lerad=0
 ethmax=tmax/mumhkg1
 ethmin=tmin/mumhkg1
 if(.NOT.radiate) then
  i=i+1
  eth0=eth+totaldt*deth	
  if(eth0.LT.ethmin) eth0=ethmin
  if(eth0.GT.ethmax) eth0=ethmax
 else
  tpart=0	  
  eth0=eth
  call heco(he,co,eth0,p)	   
  deth1=deth+(he-co)
  if(eth0.GT.ethmax) eth0=ethmax
  if(eth0.LT.ethmin) eth0=ethmin
  dt=totaldt
  laststep=.TRUE.
  if(deth1.NE.0) dt=min(dt,SQRT(ABS(eth0*dt/deth1)))
  if(dt.LT.totaldt) laststep=.FALSE.
  tnowpassed=.FALSE.
  do while((.NOT.tnowpassed).AND.i.LE.250)
   i=i+1	   	   	   
   if(ethmax-eth0.LE.dt*deth1) then
    laststep=.FALSE.
    if(deth1.NE.0) dt=(ethmax-eth0)/deth1
   endif
   if(ethmin-eth0.GE.dt*deth1) then
    laststep=.FALSE.
    if(deth1.NE.0) dt=(ethmin-eth0)/deth1
   endif 
   j=0
10 continue
  j=j+1
  if(j.GT.1000000) call terror('j overflow')
  eth1=eth0+dt*deth1	  
  eth2=eth0+0.5*dt*deth1	  
  call heco(he,co,eth2,p)	   
  deth2=deth+(he-co)
  eth2=eth2+0.5*dt*deth2
  if(ABS(eth1-eth2).LE.ttolerance*eth0) then	    	    
   tpart=tpart+dt
   tnowpassed=laststep

   lerad=lerad+dt*(he-co)
   lhe=lhe+dt*he

   laststep=.TRUE.
   if(eth1.NE.eth2) then
    dt=safefac*dt*SQRT(ABS((ttolerance*eth0)/(eth1-eth2)))
   else
    dt=totaldt-tpart
   endif
   if(dt.GE.totaldt-tpart) then
    dt=totaldt-tpart
   else
    laststep=.FALSE.
   endif
   eth0=eth2
   if(eth0.LE.ethmin) then
      eth0=ethmin
      tnowpassed=.TRUE.
   endif  
   if(eth0.GE.ethmax) then
      eth0=ethmax
      tnowpassed=.TRUE.
   endif  
  else	    
   laststep=.FALSE.
   dt=safefac*dt*SQRT(ABS((ttolerance*eth0)/(eth1-eth2)))
   goto 10	   
  endif
  jmax=MAX(jmax,j)
  call heco(he,co,eth0,p)	   
  deth1=deth+(he-co)
 enddo
 endif
 drad=he-co
 deth=(eth0-eth)/totaldt
 eth=eth0
 imax=max(i,imax)
end subroutine

subroutine extrapethrho   ! extrapolate rho, etc??
 include 'globals.h'
 real dt
 dt=tnow-teth
 if(.NOT.isotherm) then
  call extrapeth(dt)
  if(radiate) call temperature
 else
  call extraprho(dt)
 endif
 estar=estar+sum(starfuv(nbodies-nstar+1:nbodies))*dt
 teth=tnow
end subroutine

subroutine extraprho(dt)
  include 'globals.h'
  real dt
  integer p
!$omp parallel do shared(dt) private(p)
  do p=1,nsph
    rho(p)=rho(p)*exp(-hsmdivv(p)/hsmooth(p)*dt)
    hsmooth(p)=hsmooth(p)*exp(hsmdivv(p)/hsmooth(p)*dt/3)
  enddo
end subroutine

subroutine extrapeth(dt)
!  extrapolate eth to current tnow
!  using old dethdt
!  take eth/ ent step (isochoric approx.) by explicit integration
!  stabilized by use of adaptive int. 
!  note: interpolation in density ?
!  note2: control variables (eradiate etc) ok?
!  note3: cosmo ok?		
 include 'globals.h'
 integer p,i,imax,j,jmax
 real dt,eth,deth,drad,lerad,lhe
 real l_eradiate,l_snheat,l_efuvheat,l_eradcool
 real ethtoent
 real time1,time2,mintime,maxtime,tottime
				
  imax=0
  jmax=0
  mintime=1.e10
  maxtime=0.
  tottime=0.

!$omp parallel private(deth,eth,p,lhe,lerad,drad, &
!$omp ethtoent,time1,time2) &
!$omp shared(dt) &
!$omp reduction(+ : tottime) & 
!$omp reduction(max : jmax,imax,maxtime)	&
!$omp reduction(min : mintime)
  call cpu_time(time1)
!$omp do
  do p=1,nsph
   if(rho(p).EQ.0) cycle
   dethdt(p)=dethdt(p) ! note that here extrapolation should be done (only cosmo??)	         
   if(uentropy) then
    ethtoent=gamma1/rho(p)**gamma1
   else
    ethtoent=1
   endif  
   eth=ethermal(p)/ethtoent
   deth=dethdt(p)/ethtoent
   drad=derad(p)
   call exstepp(p,dt,eth,deth,drad,lerad,lhe,imax,jmax)

   rho(p)=rho(p)*exp(-hsmdivv(p)/hsmooth(p)*dt)
   hsmooth(p)=hsmooth(p)*exp(hsmdivv(p)/hsmooth(p)*dt/3)
   if(uentropy) then
    ethtoent=gamma1/rho(p)**gamma1
   else
    ethtoent=1
   endif  

   ethermal(p)=eth*ethtoent
   csound(p)=SQRT(gamma*gamma1*eth)
   derad(p)=drad  
  enddo
!$omp enddo nowait 
  call cpu_time(time2)
  mintime=MIN(mintime,time2-time1)
  maxtime=MAX(maxtime,time2-time1)
  tottime=tottime+time2-time1
!$omp end parallel
  if(verbosity.GT.0) then
   write(*,'(" <exstep> time:", 3f8.2)') mintime,maxtime,tottime
   print*,'<exstep> max iterations:',imax,jmax
  endif	
 end


subroutine exstep2(pc)
!  take half an ethermal step
!  using ethold
!  take eth/ ent step (isochoric approx.) by explicit integration
!  stabilized by use of adaptive int. 
!  note: interpolation in density ?
!  note2: control variables (eradiate etc) ok?
!  note3: cosmo ok?		
 include 'globals.h'
 integer pc
 integer k,p,i,imax,j,jmax
 real dt,eth,deth,drad,lerad,lhe
 real l_eradiate,l_snheat,l_efuvheat,l_eradcool
 real ethtoent
 real time1,time2,mintime,maxtime,tottime
				
  imax=0
  jmax=0
  mintime=1.e10
  maxtime=0.
  tottime=0.
!$omp parallel private(deth,eth,p,k,drad,dt, &
!$omp l_eradiate,l_snheat,l_efuvheat,l_eradcool, &
!$omp ethtoent,time1,time2,lhe,lerad) &
!$omp reduction(+ : eradiate,snheat,efuvheat,eradcool,tottime) & 
!$omp reduction(max : jmax,imax,maxtime)	 &
!$omp reduction(min : mintime) 
        call cpu_time(time1)
!$omp do
  do k=1,nsphact
   p=pactive(k)
   if(rho(p).EQ.0) cycle
   if(uentropy) then
    ethtoent=gamma1/rho(p)**gamma1
   else
    ethtoent=1
   endif  
   eth=ethold(p)/ethtoent
   deth=dethdt(p)/ethtoent ! note: no cosmo correction
   drad=derad(p)
   dt=dtime/2**itimestp(p)
!   if(dt.NE.tnow-tvel(p)) print*,'TE',p,dt,tnow-tvel(p) ! should be the same 
   call exstepp(p,dt,eth,deth,drad,lerad,lhe,imax,jmax)
   if(pc.EQ.1.OR.pc.EQ.2) then
    csound(p)=SQRT(gamma*gamma1*eth)
    ethermal(p)=eth*ethtoent
    derad(p)=drad  
   endif
   if(pc.EQ.2.OR.pc.EQ.3) then
    ethold(p)=eth*ethtoent
    l_eradiate=mass(p)*lerad
    l_snheat=dt*mass(p)*esnthdt(p)
    l_efuvheat=mass(p)*(lhe/fhydrogn-2*dt*esnthdt(p))
    l_eradcool=l_eradiate/fhydrogn
    eradiate=eradiate+l_eradiate
    snheat=snheat+l_snheat
    efuvheat=efuvheat+l_efuvheat
    eradcool=eradcool+l_eradcool
   endif 
  enddo
!$omp enddo nowait 
  call cpu_time(time2)
  mintime=MIN(mintime,time2-time1)
  maxtime=MAX(maxtime,time2-time1)
  tottime=tottime+time2-time1
!$omp end parallel
  if(verbosity.GT.0) then
   write(*,'(" <exstep2> time:", 3f8.2)') mintime,maxtime,tottime
   print*,'<exstep2> max iterations:',imax,jmax
  endif	
 end
