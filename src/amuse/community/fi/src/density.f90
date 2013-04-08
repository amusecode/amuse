! preempt rho=0 case ******

subroutine minigatter(n,ppos,hsearch,dens,ddensdh)
  include 'globals.h' 
  integer,intent(in) :: n
  real,intent(in) :: hsearch,ppos(3)
  real,intent(inout) :: dens, ddensdh
  integer :: i,p,iwsm
  real :: wnorm,distnorm,hsminv,dr2p,dr2,drw,wsm
  real :: dwsm,dwmnbi,dwnorm,dx,dy,dz
  dens=0.
  ddensdh=0.    
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i      
  do i=1,n
    p=srlist(i)
    dx=ppos(1)-pos(p,1)
    dy=ppos(2)-pos(p,2)
    dz=ppos(3)-pos(p,3)
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
      dens=dens+wnorm*wsm
      ddensdh=ddensdh-3*hsminv*wnorm*wsm
      dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
      dwmnbi=dwnorm*dwsm
      ddensdh=ddensdh-dr2*hsminv*dwmnbi  
    endif
  enddo   
end subroutine

subroutine hsmfunc(ppos,h,fi,dfi)
  include 'globals.h'
  integer nneigh
  real fi,dfi,h,ppos(3),dens,ddensdh
  call cond_srch(ppos,h,root,nneigh,srlist)
  call minigatter(nneigh,ppos,h,dens,ddensdh)
  fi=4*pi/3*h**3*(dens+rhomin/(massres*8./nsmooth))-nsmooth/8.
  dfi=4*pi*h**2*(dens+rhomin/(massres*8./nsmooth))+4*pi/3*h**3*ddensdh
end

subroutine newtonhsm(ppos,h,success,i)
  include 'globals.h'
  integer success,nmax,nneigh,i
  real hsmnw,hmin,hmax,prec,fi,dfi,h,ppos(3)
  parameter(nmax=10,prec=0.001) 
  success=1
  hmin=.5*h;hmax=2*h 
  hsmnw=h
  i=0
10 if(hsmnw.LT.hmax.AND.hsmnw.GT.hmin.AND.i.LT.nmax) then
      i=i+1
      call hsmfunc(ppos,hsmnw,fi,dfi)
      if(dfi.NE.0) then 
        hsmnw=hsmnw-fi/dfi
      else
        hsmnw=0
      endif 
      if(ABS(fi).LT.ABS(hsmnw*prec*dfi)) then
        h=hsmnw
        return
      endif
      goto 10
   else
      success=0
   endif
end subroutine

subroutine brackethsm(ppos,hmin,hmax,success)
  include 'globals.h'
  integer p,success,ntry,j
  real hmin,hmax,ppos(3),fac,f1,f2,df,getlocalscale
  real dum1,dum2
  parameter(fac=1.5,ntry=40)
  external getlocalscale
  success=1
  hmax=fac*getlocalscale(ppos)
  hmin=hmax/fac
  call hsmfunc(ppos,hmax,f2,df)
  call hsmfunc(ppos,hmin,f1,df)    
  do j=1,ntry
    if(f1*f2.LT.0) return
    if(abs(f1).lt.abs(f2)) then
      hmin=hmin/fac
      call hsmfunc(ppos,hmin,f1,df)
    else
      hmax=hmax*fac
      call hsmfunc(ppos,hmax,f2,df)
    endif
  enddo 
  success=0
end subroutine

subroutine safehsm(ppos,hguess,hmin,hmax)
  include 'globals.h'
  integer p,maxit,j
  real hmin,hmax,prec,hguess,dum1,dum2
  real df,dx,dxold,f,h,fh,fl,temp,xl,xh,ppos(3)
  parameter(maxit=25,prec=0.001)
  
  call hsmfunc(ppos,hmax,fh,df)
  call hsmfunc(ppos,hmin,fl,df)
  
  if(fl.eq.0) then
    hguess=hmin
    return
  endif 
  if(fh.eq.0) then
    hguess=hmax
    return
  endif
  if(hmin.GE.hmax.OR.fl.GE.fh.OR.fl*fh.GT.0) call terror('hsafepos error')
  xl=hmin
  xh=hmax  
  
  hguess=.5*(xl+xh)
  dxold=(xh-xl)
  dx=dxold
  call hsmfunc(ppos,hguess,f,df)
  do j=1,maxit
    if(((hguess-xh)*df-f)*((hguess-xl)*df-f).GT.0.OR. &
          abs(2*f).GT.ABS(dxold*df)) then
      dxold=dx
      dx=0.5*(xh-xl)
      hguess=xl+dx
      if(xl.EQ.hguess)return
    else
      dxold=dx
      dx=f/df
      temp=hguess
      hguess=hguess-dx
      if(temp.EQ.hguess)return
    endif
    if(abs(dx).LT.prec*hmax) return
    call hsmfunc(ppos,hguess,f,df)      
    if(f.lt.0) then
      xl=hguess 
    else
      xh=hguess
    endif  
  enddo
end subroutine

subroutine gethsm(ppos,h,i,j)
  include 'globals.h'
  integer i,j,success
  real hsmn,fi,dfi,hmin,hmax
  real h,ppos(3),getlocalscale
  i=0;j=0
  if(h.LE.0) h=1.75*getlocalscale(ppos)
  hsmn=h
  call newtonhsm(ppos,hsmn,success,i)
  if(success.EQ.0) then    
    j=1
    call brackethsm(ppos,hmin,hmax,success) 
    if(success.EQ.0.OR.hmin.LE.0.OR.hmax.LE.0) then
!$omp critical
      print*,getlocalscale(ppos)
      print*, ppos
      print*,rsize,root
      print*,bottom(root,1:3)
      print*,success,hmin,hmax
      call terror(' hsmpos error 1')
!$omp end critical
    endif
    call safehsm(ppos,hsmn,hmin,hmax)
    if(.NOT.(hsmn.GE.hmin.AND.hsmn.LE.hmax)) then
!$omp critical
      print*,hmin,hsmn,hmax
      call terror(' hsmpos error 2')
!$omp end critical
    endif  
  endif
  h=hsmn
end subroutine

subroutine pdensity(p,n)
  include 'globals.h'
  integer n,p
  real h,ppos(3),pvel(3),lrho,ldrhodh,lhdivv,lhcurlv,lvdisp
  h=hsmooth(p)
  ppos=pos(p,1:3)
  pvel=veltpos(p,1:3)
  call skinnydensity(n,h,ppos,pvel,lrho,ldrhodh,lhdivv,lhcurlv,lvdisp)
  rho(p)=lrho
  drhodh(p)=ldrhodh
  hsmdivv(p)=lhdivv
  hsmcurlv(p)=lhcurlv
  vdisp(p)=lvdisp
end

subroutine skinnydensity(n,h,ppos,pvel,lrho,ldrhodh,lhdivv,lhcurlv,lvdisp)
  include 'globals.h'
  integer i,iwsm,n,p
  real hsminv,wnorm,distnorm,dr2p,drw,wsm,wmass,dr2, &
       wmass1,dwnorm,vdotdr,dwsm,dwmass,dx,dy,dz, &
       dvx,dvy,dvz,dwmnbi,curlvx,curlvy,curlvz,mv(3),mv2
  real h,ppos(3),pvel(3),lrho,ldrhodh,lhdivv,lhcurlv,lvdisp
  lrho=0.
  ldrhodh=0.
  lhdivv=0.
  lhcurlv=0.
  lvdisp=0.
  curlvx=0.;curlvy=0.;curlvz=0.
  mv=0.;mv2=0.
  hsminv=1./h
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i
  do i=1,n
    p=srlist(i)
    dx=ppos(1)-pos(p,1)
    dy=ppos(2)-pos(p,2)
    dz=ppos(3)-pos(p,3)   
    if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
    if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
    if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
    if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
    if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
    if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
    dr2=dx*dx+dy*dy+dz*dz
    dr2p=dr2*distnorm
    if(ninterp.GE.dr2p) then
      iwsm=INT(dr2p)   
      drw=dr2p-iwsm
      wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
      wmass=wnorm*wsm
      lrho=lrho+mass(p)*wmass
      ldrhodh=ldrhodh-3*hsminv*mass(p)*wmass

      mv(1)=mv(1)+veltpos(p,1)*mass(p)*wmass
      mv(2)=mv(2)+veltpos(p,2)*mass(p)*wmass
      mv(3)=mv(3)+veltpos(p,3)*mass(p)*wmass
      mv2=mv2+(veltpos(p,1)**2+veltpos(p,2)**2+ &
        veltpos(p,3)**2)*mass(p)*wmass

      dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
      dwmass=dwnorm*dwsm
      dwmnbi=mass(p)*dwmass
      dvx=pvel(1)-veltpos(p,1)
      dvy=pvel(2)-veltpos(p,2)
      dvz=pvel(3)-veltpos(p,3)
      vdotdr=dvx*dx+dvy*dy+dvz*dz     
      lhdivv=lhdivv-dwmnbi*vdotdr
      ldrhodh=ldrhodh-dr2*hsminv*dwmnbi

      curlvx=curlvx+dwmnbi*(dz*dvy-dy*dvz)
      curlvy=curlvy+dwmnbi*(dx*dvz-dz*dvx)
      curlvz=curlvz+dwmnbi*(dy*dvx-dx*dvy)
    endif
  enddo   
  if(lrho.GT.0) then
    lhdivv=h*lhdivv/lrho
    lhdivv=lhdivv/(1+h*ldrhodh/3/(rhomin+lrho)) ! f_i correction 
    lhcurlv=h*SQRT(curlvx**2+curlvy**2+curlvz**2)/lrho
   
    mv=mv/lrho
    mv2=mv2/lrho
  endif
   
  lvdisp=MAX(0.,(mv2-mv(1)**2-mv(2)**2-mv(3)**2))
  lvdisp=SQRT(lvdisp)
end subroutine

subroutine phsm(p,i,j)
  include 'globals.h'
  real hneigh,h,ppos(3)
  integer p,i,j,nneigh
  ppos=pos(p,1:3)
  h=hsmooth(p)
  call gethsm(ppos,h,i,j)
  hsmooth(p)=h
end subroutine

subroutine densnhsmooth
  include 'globals.h'
  integer p,i,j,k,imax,jtot,nneigh,kmin,kmax
  integer omp_get_max_threads,nchunk,chunk,maxthread
  integer, parameter :: nbuf=32
  integer :: ib,buf(nbuf),todo(nbuf),ntodo,totalsearches
  real :: time1,time2,mintime,maxtime,tottime,utime1,utime2
  real oldrho, drhosum
 
  if(nsphact.EQ.0) return
  imax=0
  jtot=0
  nnmin=nbodies;nnmax=0;nntot=0
  mintime=1.e10; maxtime=0.; tottime=0
  drhosum=0

  maxthread=1
  nchunk=1
!$  maxthread=omp_get_max_threads()
!$  nchunk=MAX(MIN(10*maxthread,nsphact/nbuf),maxthread)
  totalsearches=0
!$omp parallel shared(nchunk) &
!$omp private(i,j,p,k,nneigh,kmin,kmax, buf,todo,ntodo,ib,chunk,oldrho) &
!$omp reduction( max : imax,nnmax,maxtime) &
!$omp reduction(+ : jtot,nntot,tottime,totalsearches,drhosum) &
!$omp reduction( min : nnmin,mintime) 
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
      call presearch(k,kmax,nbuf,buf,ntodo,todo)
      do ib=1,ntodo
        p=todo(ib)
        call phsm(p,i,j)  
        call pcond_srch(root,p,nneigh,srlist)
        oldrho=rho(p)
        call pdensity(p,nneigh)
        if(oldrho.NE.0) drhosum=drhosum+abs(rho(p)-oldrho)/oldrho
        nnmin=MIN(nnmin,nneigh)
        nnmax=MAX(nnmax,nneigh)
        nntot=nntot+nneigh
        imax=MAX(i,imax)
        jtot=jtot+j 
        if(eps_is_h) epsgrav(p)=hsmooth(p)
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
  if(verbosity.GT.0) print*,'<densnhsmooth> parts,searches', nsphact,totalsearches 
  if(verbosity.GT.0) print*,'<densnhsmooth> < a > t',nnmin,nnavg,nnmax,nntot
  if(verbosity.GT.0) then
    write(*,'(" <densnhsmooth> time:", 3f8.2)') maxtime,mintime,tottime
    print*,'<densnhsmooth> max iter, fails:',imax,jtot
  endif
  
  if(verbosity.GT.0.AND.drhosum/npactive.GT.0.01) &
      print*,' *** drho warning *** ',drhosum/npactive, npactive
end subroutine

subroutine gatterdens(n,spos,hsearch,dens,ddensdh)
  include 'globals.h' 
  integer,intent(in) :: n
  real,dimension(3),intent(in) :: spos
  real,intent(in) :: hsearch
  real,intent(inout) :: dens, ddensdh
  integer :: i,p,iwsm
  real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm
  real :: dwmass,dwsm,dwmnbi,dwnorm,wmass
  dens=0.
  ddensdh=0.    
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i      
  do i=1,n
    p=srlist(i)
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
      wmass=wnorm*wsm
      dens=dens+mass(p)*wmass
      ddensdh=ddensdh-3*hsminv*mass(p)*wmass
      dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
      dwmass=dwnorm*dwsm
      dwmnbi=mass(p)*dwmass
      ddensdh=ddensdh-dr2*hsminv*dwmnbi  
    endif
  enddo   
end subroutine 

subroutine gatter_hydro_state(n,spos,svel,hsearch,dens,rhov,rhov2,rhoe)
  include 'globals.h' 
  integer,intent(in) :: n
  real,dimension(3),intent(in) :: spos, svel
  real,intent(in) :: hsearch
  real,intent(inout) :: dens,rhov(3),rhoe,rhov2
  integer :: i,p,iwsm
  real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm
  real :: dwmass,dwsm,dwmnbi,dwnorm,wmass
  dens=0.
  rhov=0.
  rhoe=0.
  rhov2=0.
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  distnorm=hsminv**2*deldr2i      
  do i=1,n
    p=srlist(i)
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
      wmass=wnorm*wsm
      dens=dens+mass(p)*wmass
      rhov=rhov+mass(p)*(vel(p,1:3)-svel(1:3))*wmass
      rhov2=rhov2+mass(p)*sum((vel(p,1:3)-svel(1:3))**2)*wmass
      rhoe=rhoe+mass(p)*ethermal(p)*wmass
    endif
  enddo   
end subroutine 

subroutine hsmdenspos2(ppos,h,dens,ddensdh,nneigh)
  include 'globals.h'
  integer nneigh,i,j
  real dens,ddensdh,h,ppos(3)
  call gethsm(ppos,h,i,j)
  call cond_srch(ppos,h,root,nneigh,srlist)
  call gatterdens(nneigh,ppos,h,dens,ddensdh)
end subroutine
