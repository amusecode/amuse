subroutine omp_entdotaccsphco
 include 'globals.h'
 integer p,k,nneigh,kmin,kmax
 integer, parameter :: nbuf=32
 integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
 integer omp_get_max_threads,nchunk,chunk,maxthread
 real time1,time2,mintime,maxtime,tottime,utime1,utime2

 if(nsphact.EQ.0) return
 nnmin=nbodies; nnmax=0; nntot=0
 mintime=1.e10; maxtime=0; tottime=0

 maxthread=1
 nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	
 totalsearches=0
!$omp parallel private(p,k,nneigh,time1,time2,kmin,kmax, &
!$omp  buf,todo,ntodo,ib,chunk) shared(nchunk)&
!$omp reduction( + : tottime,nntot,totalsearches) & 
!$omp reduction( MIN : mintime,nnmin) & 
!$omp reduction( MAX : maxtime,nnmax)
 call cpu_time(time1)
 ncalls=0;nsearches=0
!$omp do schedule(guided,1)
 do chunk=1,nchunk
  kmin=((chunk-1)*nsphact)/nchunk+1
  kmax=(chunk*nsphact)/nchunk
  buf=0
  reuseflag=1
  searchreuse=0
  do k=kmin,kmax
   call precomsearch(k,kmax,nbuf,buf,ntodo,todo)
   do ib=1,ntodo
    p=todo(ib)
    call pcond_comsrch(root,p,nneigh,srlist)
    call pentdotaccsphco(p,nneigh)
    nnmin=MIN(nnmin,nneigh)
    nnmax=MAX(nnmax,nneigh)
    nntot=nntot+nneigh
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
 nnavg=nntot/nsphact
 if(verbosity.GT.0) print*,'<dentacc> parts,searches', nsphact,totalsearches 
 if(verbosity.GT.0) print*,'<dentacc> < a > t',nnmin,nnavg,nnmax,nntot
 if(verbosity.GT.0) write(*,'(" <dentacc> time:", 3f8.2)') mintime,maxtime,tottime
end

subroutine omp_accsphco
 include 'globals.h'
 integer p,k,nneigh,kmin,kmax
 integer, parameter :: nbuf=32
 integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
 integer omp_get_max_threads,nchunk,chunk,maxthread
 real time1,time2,mintime,maxtime,tottime,utime1,utime2

 if(nsphact.EQ.0) return
 nnmin=nbodies; nnmax=0; nntot=0
 mintime=1.e10; maxtime=0; tottime=0

 maxthread=1
 nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	
totalsearches=0
!$omp parallel private(p,k,nneigh,time1,time2,kmin,kmax, &
!$omp  buf,todo,ntodo,ib,chunk) shared(nchunk)&
!$omp reduction( + : tottime,nntot) & 
!$omp reduction( MIN : mintime,nnmin) & 
!$omp reduction( MAX : maxtime,nnmax)
 call cpu_time(time1)
 ncalls=0;nsearches=0
!$omp do schedule(guided,1)
 do chunk=1,nchunk
  kmin=((chunk-1)*nsphact)/nchunk+1
  kmax=(chunk*nsphact)/nchunk
  buf=0
  reuseflag=1
  searchreuse=0
  do k=kmin,kmax
   call precomsearch(k,kmax,nbuf,buf,ntodo,todo)
   do ib=1,ntodo
    p=todo(ib)
    call pcond_comsrch(root,p,nneigh,srlist)
    call paccsphco(p,nneigh)
    nnmin=MIN(nnmin,nneigh)
    nnmax=MAX(nnmax,nneigh)
    nntot=nntot+nneigh
   enddo
  enddo
 enddo
!$omp enddo nowait  
 mintime=MIN(mintime,time2-time1)
 maxtime=MAX(maxtime,time2-time1)
 tottime=tottime+time2-time1
 totalsearches=totalsearches+nsearches
!$omp end parallel
 nnavg=nntot/nsphact
 if(verbosity.GT.0) print*,'<accsph> parts,searches', nsphact,totalsearches 
 if(verbosity.GT.0) print*,'<accsph> < a > t',nnmin,nnavg,nnmax,nntot
 if(verbosity.GT.0) write(*,'(" <accsph> time:", 3f8.2)') mintime,maxtime,tottime
end


subroutine omp_entdot
 include 'globals.h'
 integer p,k,nneigh,kmin,kmax
 integer, parameter :: nbuf=32
 integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
 integer omp_get_max_threads,nchunk,chunk,maxthread
 real time1,time2,mintime,maxtime,tottime,utime1,utime2

 if(nsphact.EQ.0) return
 nnmin=nbodies; nnmax=0; nntot=0
 mintime=1.e10; maxtime=0; tottime=0

 maxthread=1
 nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	
 totalsearches=0 
!$omp parallel private(p,k,nneigh,time1,time2,kmin,kmax, &
!$omp  buf,todo,ntodo,ib,chunk) shared(nchunk)&
!$omp reduction( + : tottime,nntot) & 
!$omp reduction( MIN : mintime,nnmin) & 
!$omp reduction( MAX : maxtime,nnmax)
 call cpu_time(time1)
 ncalls=0;nsearches=0
!$omp do schedule(guided,1)
 do chunk=1,nchunk
  kmin=((chunk-1)*nsphact)/nchunk+1
  kmax=(chunk*nsphact)/nchunk
  buf=0
  reuseflag=1
  searchreuse=0
  do k=kmin,kmax
   call precomsearch(k,kmax,nbuf,buf,ntodo,todo)
   do ib=1,ntodo
    p=todo(ib)
    call pcond_comsrch(root,p,nneigh,srlist)
    call pentdot(p,nneigh)
    nnmin=MIN(nnmin,nneigh)
    nnmax=MAX(nnmax,nneigh)
    nntot=nntot+nneigh
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
 call cpu_time(utime2)
 nnavg=nntot/nsphact
 if(verbosity.GT.0) print*,'<entdot> parts,searches', nsphact,totalsearches
 if(verbosity.GT.0) print*,'<entdot> < a > t',nnmin,nnavg,nnmax,nntot
 if(verbosity.GT.0) write(*,'(" <entdot> time:", 3f8.2)') mintime,maxtime,tottime
end

subroutine pentdot(p,n)
 include 'globals.h'
 integer n,p
 real h,ppos(3),pvel(3),lrho,lcsound,ldentdt,imumax
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  lrho=rho(p)
  lcsound=csound(p)
  call entdot(n,h,ppos,pvel,lrho,lcsound,ldentdt,imumax)
  dentdt(p)=ldentdt
  mumaxdvh(p)=imumax
end

subroutine paccsphco(p,n)
 include 'globals.h'
 integer n,p
 real h,ppos(3),pvel(3),pacc(3),lrho,lcsound,ldrhodh
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  lrho=rho(p)
  ldrhodh=drhodh(p)
  lcsound=csound(p)
  call accsphco(n,h,ppos,pvel,pacc,lrho,lcsound,ldrhodh)
  acc(p,1:3)=acc(p,1:3)+pacc(1:3)
end


subroutine entdot(n,hsm,ppos,pvel,lrho,lcsound,dent,tmuij)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: ppos,pvel
 real,intent(in) :: hsm,lrho,lcsound
 real,intent(inout) :: dent,tmuij
 integer :: i,nb,iwsm,iwsm1
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,drw1,tmass,dr2,dr2i,wsm
 real :: dwmass,dwmass1,dwsm,dwsm1,dwmnbi,dwnorm,wmass,muij,eij,piij
 real :: vdotdr,hsmavg
 dent=0
 tmuij=0
 if(lrho.EQ.0) return
 if(symmetry.EQ.'hk') then
  hsminv=1./hsm
  dwnorm=piinv*hsminv**5
  distnorm=hsminv*hsminv*deldr2i
  do i=1,n
  nb=srlist(i)
  dx=ppos(1)-pos(nb,1)
  dy=ppos(2)-pos(nb,2)
  dz=ppos(3)-pos(nb,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  if(dr2.GT.4*hsm**2.AND.dr2.GT.4*hsmooth(nb)**2) cycle
  dr2p=dr2*distnorm
  if(dr2p.LT.ninterp) then
   iwsm=INT(dr2p)
   drw=dr2p-iwsm
   dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
   dwmass=dwnorm*dwsm
  else
   dwmass=0.
  endif
  dr2i=dr2*deldr2i/(hsmooth(nb)*hsmooth(nb))
  if(dr2i.LT.ninterp) then
   iwsm1=INT(dr2i)
   drw1=dr2i-iwsm1
   dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
   dwmass1=piinv/(hsmooth(nb)**5)*dwsm1
  else
   dwmass1=0.
  endif
  vdotdr=(pvel(1)-veltpos(nb,1))*dx+ &
  	 (pvel(2)-veltpos(nb,2))*dy+ &
  	 (pvel(3)-veltpos(nb,3))*dz
  hsmavg=.5*(hsm+hsmooth(nb))
  muij=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
  if(vdotdr.GT.0.) muij=0.
  piij=(-alpha*muij*(max(lcsound,csound(nb)))+beta*muij**2)/ &
                                                 (min(lrho,rho(nb)))
! gadget visc
!  if(dr2.NE.0) muij=vdotdr/sqrt(dr2)
!  if(vdotdr.GT.0..OR.dr2.EQ.0) muij=0.
!  piij=-alpha*muij*(lcsound+csound(nb)-3*muij)/(lrho+rho(nb))
  eij=vdotdr*0.5*piij
  dent=dent+eij*mass(nb)*0.5*(dwmass+dwmass1)
  tmuij=MAX(tmuij,ABS(muij))
! gadget visc  
!  tmuij=MAX(tmuij, lcsound+csound(nb)-3*muij)  
  enddo
 else
  do i=1,n
  nb=srlist(i)
  dx=ppos(1)-pos(nb,1)
  dy=ppos(2)-pos(nb,2)
  dz=ppos(3)-pos(nb,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  if(dr2.GT.(hsm+hsmooth(nb))**2) cycle
  hsminv=2./(hsm+hsmooth(nb))
  dwnorm=piinv*hsminv**5
  distnorm=hsminv**2*deldr2i
  dr2p=dr2*distnorm
  if(dr2p.LT.ninterp) then
   iwsm=INT(dr2p)
   drw=dr2p-iwsm
   dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
   dwmass=dwnorm*dwsm
  else
   dwmass=0.
  endif
   vdotdr=(pvel(1)-veltpos(nb,1))*dx+ &
	  (pvel(2)-veltpos(nb,2))*dy+ &
	  (pvel(3)-veltpos(nb,3))*dz
   hsmavg=0.5*(hsm+hsmooth(nb))
   muij=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
   if(vdotdr.GT.0.) muij=0.
   piij=(-alpha*muij*(max(lcsound,csound(nb)))+beta*muij**2)/ &
                                                 (min(lrho,rho(nb)))
   eij=vdotdr*0.5*piij
   dent=dent+eij*mass(nb)*dwmass
   tmuij=MAX(tmuij,ABS(muij))
  enddo
 endif
 dent=dent*gamma1/(lrho**gamma1)
end

subroutine accsphco(n,hsm,ppos,pvel,pacc,lrho,lcsound,ldrhodh)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: ppos,pvel
 real,dimension(3),intent(inout) :: pacc
 real,intent(in) :: hsm,lrho,lcsound,ldrhodh
 integer :: i,nb,iwsm,iwsm1
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,drw1,tmass,dr2,dr2i,wsm
 real :: dwmass,dwmass1,dwsm,dwsm1,dwmnbi,dwnorm,wmass,muij,eij,piij
 real :: vdotdr,hsmavg,aij,aijmass,fi,fnbi

 pacc=0
 if(lrho.EQ.0) return
 hsminv=1./hsm
 dwnorm=piinv*hsminv**5
 distnorm=hsminv*hsminv*deldr2i
 do i=1,n
  nb=srlist(i)
  dx=ppos(1)-pos(nb,1)
  dy=ppos(2)-pos(nb,2)
  dz=ppos(3)-pos(nb,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  if(dr2.GT.4*hsm**2.AND.dr2.GT.4*hsmooth(nb)**2) cycle
  dr2p=dr2*distnorm
  if(dr2p.LT.ninterp) then
   iwsm=INT(dr2p)
   drw=dr2p-iwsm
   dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
   dwmass=dwnorm*dwsm
  else
   dwmass=0.
  endif
  dr2i=dr2*deldr2i/(hsmooth(nb)*hsmooth(nb))		    
  if(dr2i.LT.ninterp) then
   iwsm1=INT(dr2i)
   drw1=dr2i-iwsm1
   dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
   dwmass1=piinv/(hsmooth(nb)**5)*dwsm1
  else
   dwmass1=0.
  endif
  vdotdr=(pvel(1)-veltpos(nb,1))*dx+ &
	 (pvel(2)-veltpos(nb,2))*dy+ &
	 (pvel(3)-veltpos(nb,3))*dz
  hsmavg=.5*(hsm+hsmooth(nb))
  muij=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
  if(vdotdr.GT.0.) muij=0.
  piij=(-alpha*muij*(max(lcsound,csound(nb)))+beta*muij**2)/ &
                                                  (min(lrho,rho(nb)))
!  gadget visc
!  if(dr2.NE.0) muij=vdotdr/sqrt(dr2)
!  if(vdotdr.GT.0..OR.dr2.EQ.0) muij=0.
!  piij=-alpha*muij*(lcsound+csound(nb)-3*muij)/(lrho+rho(nb))
 		  
  fi=1/(1+ldrhodh*hsm/3/(rhomin+lrho))
  fnbi=1/(1+drhodh(nb)*hsmooth(nb)/3/(rhomin+rho(nb)))
  
  aij=fi*lcsound**2/(gamma*lrho)*dwmass+ &
      fnbi*csound(nb)**2/(gamma*rho(nb))*dwmass1
  
  aij=aij+piij*0.5*(dwmass+dwmass1)
  aijmass=aij*mass(nb)
  pacc(1)=pacc(1)-aijmass*dx
  pacc(2)=pacc(2)-aijmass*dy
  pacc(3)=pacc(3)-aijmass*dz
 enddo
end

subroutine pentdotaccsphco(p,n)
 include 'globals.h'
 integer n,p
 real h,ppos(3),pvel(3),pacc(3),lrho,lcsound,ldrhodh,ldentdt,imumax
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  lrho=rho(p)
  ldrhodh=drhodh(p)
  lcsound=csound(p)
  call entdotaccsphco(n,h,ppos,pvel,pacc,lrho,lcsound,ldrhodh,ldentdt,imumax)
  acc(p,1:3)=acc(p,1:3)+pacc(1:3)
  dentdt(p)=ldentdt
  mumaxdvh(p)=imumax
end


subroutine entdotaccsphco(n,hsm,ppos,pvel,pacc,lrho,lcsound,ldrhodh,dent,tmuij)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: ppos,pvel
 real,intent(inout) :: pacc(3),dent,tmuij
 real,intent(in) :: hsm,lrho,lcsound,ldrhodh
 integer :: i,nb,iwsm,iwsm1
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,drw1,tmass,dr2,dr2i,wsm
 real :: dwmass,dwmass1,dwsm,dwsm1,dwmnbi,dwnorm,wmass,muij,eij,piij
 real :: vdotdr,hsmavg,aij,aijmass,fi,fnbi

 pacc=0
 dent=0
 tmuij=0
 if(lrho.EQ.0) return
 hsminv=1./hsm
 dwnorm=piinv*hsminv**5
 distnorm=hsminv*hsminv*deldr2i
 do i=1,n
  nb=srlist(i)
  dx=ppos(1)-pos(nb,1)
  dy=ppos(2)-pos(nb,2)
  dz=ppos(3)-pos(nb,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  if(dr2.GT.4*hsm**2.AND.dr2.GT.4*hsmooth(nb)**2) cycle
  dr2p=dr2*distnorm
  if(dr2p.LT.ninterp) then
   iwsm=INT(dr2p)
   drw=dr2p-iwsm
   dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
   dwmass=dwnorm*dwsm
  else
   dwmass=0.
  endif
  dr2i=dr2*deldr2i/(hsmooth(nb)*hsmooth(nb))		    
  if(dr2i.LT.ninterp) then
   iwsm1=INT(dr2i)
   drw1=dr2i-iwsm1
   dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
   dwmass1=piinv/(hsmooth(nb)**5)*dwsm1
  else
   dwmass1=0.
  endif
  
  vdotdr=(pvel(1)-veltpos(nb,1))*dx+ &
	 (pvel(2)-veltpos(nb,2))*dy+ &
	 (pvel(3)-veltpos(nb,3))*dz
  hsmavg=.5*(hsm+hsmooth(nb))
  muij=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
  if(vdotdr.GT.0.) muij=0.
  piij=(-alpha*muij*(max(lcsound,csound(nb)))+beta*muij**2)/ &
                                                  (min(lrho,rho(nb)))
!  gadget visc
!  if(dr2.NE.0) muij=vdotdr/sqrt(dr2)
!  if(vdotdr.GT.0..OR.dr2.EQ.0) muij=0.
!  piij=-alpha*muij*(lcsound+csound(nb)-3*muij)/(lrho+rho(nb))


 		  
!  fi=1/(1+ldrhodh*hsm/3/(rhomin+lrho))
!  fnbi=1/(1+drhodh(nb)*hsmooth(nb)/3/(rhomin+rho(nb)))
  fi=(rhomin+lrho)/((rhomin+lrho)+ldrhodh*hsm/3.)
  fnbi=(rhomin+rho(nb))/((rhomin+rho(nb))+drhodh(nb)*hsmooth(nb)/3.)
  
  aij=fi*lcsound**2/(gamma*lrho)*dwmass+ &
      fnbi*csound(nb)**2/(gamma*rho(nb))*dwmass1
  
  aij=aij+piij*0.5*(dwmass+dwmass1)
  aijmass=aij*mass(nb)
  pacc(1)=pacc(1)-aijmass*dx
  pacc(2)=pacc(2)-aijmass*dy
  pacc(3)=pacc(3)-aijmass*dz
  
  eij=vdotdr*0.5*piij
  dent=dent+eij*mass(nb)*0.5*(dwmass+dwmass1)
  tmuij=MAX(tmuij,ABS(muij))
!  gadget visc
!  tmuij=MAX(tmuij, lcsound+csound(nb)-3*muij)
 enddo
 dent=dent*gamma1/(lrho**gamma1)
end


subroutine omp_ethdotaccsphco
 include 'globals.h'
 integer p,k,nneigh,kmin,kmax
 integer, parameter :: nbuf=32
 integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
 integer omp_get_max_threads,nchunk,chunk,maxthread
 real time1,time2,mintime,maxtime,tottime,utime1,utime2

 if(nsphact.EQ.0) return
 nnmin=nbodies; nnmax=0; nntot=0
 mintime=1.e10; maxtime=0; tottime=0

 maxthread=1
 nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	
 totalsearches=0
!$omp parallel private(p,k,nneigh,time1,time2,kmin,kmax, &
!$omp  buf,todo,ntodo,ib,chunk) shared(nchunk)&
!$omp reduction( + : tottime,nntot,totalsearches) & 
!$omp reduction( MIN : mintime,nnmin) & 
!$omp reduction( MAX : maxtime,nnmax)
 call cpu_time(time1)
 ncalls=0;nsearches=0
!$omp do schedule(guided,1)
 do chunk=1,nchunk
  kmin=((chunk-1)*nsphact)/nchunk+1
  kmax=(chunk*nsphact)/nchunk
  buf=0
  reuseflag=1
  searchreuse=0
  do k=kmin,kmax
   call precomsearch(k,kmax,nbuf,buf,ntodo,todo)
   do ib=1,ntodo
    p=todo(ib)
    call pcond_comsrch(root,p,nneigh,srlist)
    call pethdotaccsphco(p,nneigh)
    nnmin=MIN(nnmin,nneigh)
    nnmax=MAX(nnmax,nneigh)
    nntot=nntot+nneigh
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
 nnavg=nntot/nsphact
 if(verbosity.GT.0) print*,'<dethacc> parts,searches', nsphact,totalsearches 
 if(verbosity.GT.0) print*,'<dethacc> < a > t',nnmin,nnavg,nnmax,nntot
 if(verbosity.GT.0) write(*,'(" <dethacc> time:", 3f8.2)') mintime,maxtime,tottime
end

subroutine omp_ethdotco
 include 'globals.h'
 integer p,k,nneigh,kmin,kmax
 integer, parameter :: nbuf=32
 integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
 integer omp_get_max_threads,nchunk,chunk,maxthread
 real time1,time2,mintime,maxtime,tottime,utime1,utime2

 if(nsphact.EQ.0) return
 nnmin=nbodies; nnmax=0; nntot=0
 mintime=1.e10; maxtime=0; tottime=0

 maxthread=1
 nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)	
 totalsearches=0
!$omp parallel private(p,k,nneigh,time1,time2,kmin,kmax, &
!$omp  buf,todo,ntodo,ib,chunk) shared(nchunk)&
!$omp reduction( + : tottime,nntot,totalsearches) & 
!$omp reduction( MIN : mintime,nnmin) & 
!$omp reduction( MAX : maxtime,nnmax)
 call cpu_time(time1)
 ncalls=0;nsearches=0
!$omp do schedule(guided,1)
 do chunk=1,nchunk
  kmin=((chunk-1)*nsphact)/nchunk+1
  kmax=(chunk*nsphact)/nchunk
  buf=0
  reuseflag=1
  searchreuse=0
  do k=kmin,kmax
   call precomsearch(k,kmax,nbuf,buf,ntodo,todo)
   do ib=1,ntodo
    p=todo(ib)
    call pcond_comsrch(root,p,nneigh,srlist)
    call pethdotco(p,nneigh)
    nnmin=MIN(nnmin,nneigh)
    nnmax=MAX(nnmax,nneigh)
    nntot=nntot+nneigh
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
 nnavg=nntot/nsphact
 if(verbosity.GT.0) print*,'<deth> parts,searches', nsphact,totalsearches 
 if(verbosity.GT.0) print*,'<deth> < a > t',nnmin,nnavg,nnmax,nntot
 if(verbosity.GT.0) write(*,'(" <deth> time:", 3f8.2)') mintime,maxtime,tottime
end


subroutine pethdotaccsphco(p,n)
 include 'globals.h'
 integer n,p
 real h,ppos(3),pvel(3),pacc(3),lrho,lcsound,ldrhodh,ldethdt,imumax
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  lrho=rho(p)
  ldrhodh=drhodh(p)
  lcsound=csound(p)
  call ethdotaccsphco(n,h,ppos,pvel,pacc,lrho,lcsound,ldrhodh,ldethdt,imumax)
  acc(p,1:3)=acc(p,1:3)+pacc(1:3)
!  if(.NOT.isotherm) dethdt(p)=ldethdt
  dethdt(p)=ldethdt
  mumaxdvh(p)=imumax
end

subroutine pethdotco(p,n)
 include 'globals.h'
 integer n,p
 real h,ppos(3),pvel(3),pacc(3),lrho,lcsound,ldrhodh,ldethdt,imumax
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  lrho=rho(p)
  ldrhodh=drhodh(p)
  lcsound=csound(p)
  call ethdotaccsphco(n,h,ppos,pvel,pacc,lrho,lcsound,ldrhodh,ldethdt,imumax)
!  acc(p,1:3)=acc(p,1:3)+pacc(1:3)
  dethdt(p)=ldethdt
  mumaxdvh(p)=imumax
end


subroutine ethdotaccsphco(n,hsm,ppos,pvel,pacc,lrho,lcsound,ldrhodh,deth,tmuij)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: ppos,pvel
 real,intent(inout) :: pacc(3),deth,tmuij
 real,intent(in) :: hsm,lrho,lcsound,ldrhodh
 integer :: i,nb,iwsm,iwsm1
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,drw1,tmass,dr2,dr2i,wsm
 real :: dwmass,dwmass1,dwsm,dwsm1,dwmnbi,dwnorm,wmass,muij,eij,piij
 real :: vdotdr,hsmavg,aij,aijmass,fi,fnbi

 pacc=0
 deth=0
 tmuij=0
 if(lrho.EQ.0) return
 hsminv=1./hsm
 dwnorm=piinv*hsminv**5
 distnorm=hsminv*hsminv*deldr2i
 do i=1,n
  nb=srlist(i)
  dx=ppos(1)-pos(nb,1)
  dy=ppos(2)-pos(nb,2)
  dz=ppos(3)-pos(nb,3)
  if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
  if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
  if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
  if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
  if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
  if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
  dr2=dx*dx+dy*dy+dz*dz
  if(dr2.GT.4*hsm**2.AND.dr2.GT.4*hsmooth(nb)**2) cycle
  dr2p=dr2*distnorm
  if(dr2p.LT.ninterp) then
   iwsm=INT(dr2p)
   drw=dr2p-iwsm
   dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
   dwmass=dwnorm*dwsm
  else
   dwmass=0.
  endif
  dr2i=dr2*deldr2i/(hsmooth(nb)*hsmooth(nb))		    
  if(dr2i.LT.ninterp) then
   iwsm1=INT(dr2i)
   drw1=dr2i-iwsm1
   dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
   dwmass1=piinv/(hsmooth(nb)**5)*dwsm1
  else
   dwmass1=0.
  endif
  
  vdotdr=(pvel(1)-veltpos(nb,1))*dx+ &
	 (pvel(2)-veltpos(nb,2))*dy+ &
	 (pvel(3)-veltpos(nb,3))*dz
  hsmavg=.5*(hsm+hsmooth(nb))
  muij=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
  if(vdotdr.GT.0.) muij=0.
  piij=(-alpha*muij*(max(lcsound,csound(nb)))+beta*muij**2)/ &
                                                  (min(lrho,rho(nb)))
 		  
!  fi=1/(1+ldrhodh*hsm/3/(rhomin+lrho))
!  fnbi=1/(1+drhodh(nb)*hsmooth(nb)/3/(rhomin+rho(nb)))
  fi=(rhomin+lrho)/((rhomin+lrho)+ldrhodh*hsm/3.)
  fnbi=(rhomin+rho(nb))/((rhomin+rho(nb))+drhodh(nb)*hsmooth(nb)/3.)
  
  aij=fi*lcsound**2/(gamma*lrho)*dwmass+ &
      fnbi*csound(nb)**2/(gamma*rho(nb))*dwmass1
  
  aij=aij+piij*0.5*(dwmass+dwmass1)
  aijmass=aij*mass(nb)
  pacc(1)=pacc(1)-aijmass*dx
  pacc(2)=pacc(2)-aijmass*dy
  pacc(3)=pacc(3)-aijmass*dz
  
  eij=0.5*vdotdr*aij
  deth=deth+eij*mass(nb)
  tmuij=MAX(tmuij,ABS(muij))
 enddo
end
