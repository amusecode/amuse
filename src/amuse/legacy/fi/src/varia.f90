! weights: returns normalized weights through tempvect for position spos 

subroutine weights(n,spos,secondsmooth)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: secondsmooth
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm
 
 tmass=0.
 if(n.EQ.0) return
  
 do i=1,n
  p=bodlist(i)
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
  hsminv=1./max(hsmooth(p),secondsmooth)
  wnorm=piinv*hsminv*hsminv*hsminv
  distnorm=hsminv**2*deldr2i
  dr2p=dr2*distnorm
  iwsm=INT(dr2p)
  if(ninterp.LT.iwsm) iwsm=ninterp
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm
  tmass=tmass+mass(p)*tempvect(i)
 enddo
 
 if(tmass.eq.0) call terror(' weights error ')
 
 do i=1,n
  tempvect(i)=tempvect(i)/tmass
 enddo
   
end subroutine weights


!  gatter: returns density through tempvect for position spos

function gatter(n,spos,secondsmooth)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: secondsmooth
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatter
 
 tmass=0.
 if(n.EQ.0) return
    
 do i=1,n
  p=bodlist(i)
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
  hsminv=1./max(hsmooth(p),secondsmooth)
  wnorm=piinv*hsminv*hsminv*hsminv
  distnorm=hsminv**2*deldr2i
  dr2p=dr2*distnorm
  iwsm=INT(dr2p)
  if(ninterp.LT.iwsm) iwsm=ninterp
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm
  tmass=tmass+mass(p)*tempvect(i)
 enddo
 if(tmass.eq.0) call terror(' weights error ')
 gatter=tmass
   
end function gatter


function gatterdensity(position)
 include 'globals.h'
 real,dimension(3),intent(in) :: position
 real :: gatterdensity
 real :: gatter
 integer :: inear
 inear=0;
 call comsearch(root,0.,0.,position,inear,bodlist)
 if(inear.LE.0) print*,' gatterdens inear<1'
 gatterdensity=gatter(inear,position,0.)
 if(gatterdensity.LE.0) print*,' gatterdensity:',inear, gatterdensity,position
end function gatterdensity



function gatscatdensity(position)
 include 'globals.h'
 real,dimension(3),intent(in) :: position
 real :: gatscatdensity
 real :: gatscat,hsearch,sethsmooth
 integer :: inear
 hsearch=sethsmooth(position)
 inear=0;
 call comsearch(root,2*hsearch,0.,position,inear,bodlist)
 if(inear.LT.0) print*,' gatscatdens inear<nsmooth!'
 gatscatdensity=gatscat(inear,position,hsearch)
 if(gatscatdensity.LE.0) print*,' gatscatdensity:',inear, gatscatdensity,position
end function gatscatdensity

function gatscat(n,spos,secondsmooth)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: secondsmooth
 integer :: i,p,iwsm,iwsm1
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatscat
 real :: wnorm1,dr2i,wsm1,wmass,wmass1,drw1
 
 
 tmass=0.
 gatscat=0.
 if(n.EQ.0) return
    
 hsminv=1./secondsmooth
 wnorm=piinv*hsminv*hsminv*hsminv
 distnorm=hsminv*hsminv*deldr2i
    
 do i=1,n
  p=bodlist(i)
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
  iwsm=INT(dr2p)
  if(ninterp.LT.iwsm) iwsm=ninterp
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  wmass=wnorm*wsm
  dr2i=dr2*deldr2i/(hsmooth(p)*hsmooth(p))
  iwsm1=INT(dr2i)
  if(ninterp.LT.iwsm1) iwsm1=ninterp
  drw1=dr2i-iwsm1
  wsm1=(1.-drw1)*wsmooth(iwsm1)+drw1*wsmooth(1+iwsm1)
  wmass1=piinv/(hsmooth(p)*hsmooth(p)*hsmooth(p))*wsm1
  tmass=tmass+.5*mass(p)*(wmass+wmass1)
 enddo
 if(tmass.le.0) call terror(' weights error ')
 gatscat=tmass
   
end function gatscat


subroutine denscomp
include 'globals.h'
integer p
real x(3),gatterdensity
do p=1,nsph
x(1)=pos(p,1);x(2)=pos(p,2);x(3)=pos(p,3)
print*,'denscomp',rho(p),gatterdensity(x)
enddo
end subroutine denscomp


 
subroutine gatterweight(n,spos,hsearch,dens,ddensdh,weights)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: hsearch,weights(*)
 real,intent(inout) :: dens, ddensdh
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatter
 real :: dwmass,dwsm,dwmnbi,dwnorm
 dens=0.
 ddensdh=0.    
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i 
 do i=1,n
  p=bodlist(i)
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
  iwsm=INT(dr2p)
  if(ninterp.LT.iwsm) iwsm=ninterp
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm
  dens=dens+weights(p)*tempvect(i)
  ddensdh=ddensdh-3*hsminv*weights(p)*tempvect(i)
  dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
  dwmass=dwnorm*dwsm
  dwmnbi=weights(p)*dwmass
  ddensdh=ddensdh-dr2*hsminv*dwmnbi  
 enddo   
    
end subroutine gatterweight


subroutine gatter2(n,spos,hsearch,dens,ddensdh)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: hsearch
 real,intent(inout) :: dens, ddensdh
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatter
 real :: dwmass,dwsm,dwmnbi,dwnorm
 dens=0.
 ddensdh=0.    
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i 
 do i=1,n
  tempvect(i)=0.
  p=bodlist(i)
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
  dens=dens+mass(p)*tempvect(i)
  ddensdh=ddensdh-3*hsminv*mass(p)*tempvect(i)
  dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
  dwmass=dwnorm*dwsm
  dwmnbi=mass(p)*dwmass
  ddensdh=ddensdh-dr2*hsminv*dwmnbi  
  endif
 enddo   

end subroutine gatter2

subroutine gatter3(n,spos,hsearch,dens,ddensdh,meana,ain)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 real,intent(in) :: hsearch, ain(*)
 real,intent(inout) :: dens, ddensdh,meana
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatter
 real :: dwmass,dwsm,dwmnbi,dwnorm
 dens=0.
 ddensdh=0.
 meana=0.    
  hsminv=1./hsearch
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i 
 do i=1,n
  tempvect(i)=0.
  p=bodlist(i)
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
  dens=dens+mass(p)*tempvect(i)
  meana=meana+ain(p)*tempvect(i)
  ddensdh=ddensdh-3*hsminv*mass(p)*tempvect(i)
  dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
  dwmass=dwnorm*dwsm
  dwmnbi=mass(p)*dwmass
  ddensdh=ddensdh-dr2*hsminv*dwmnbi  
  endif
 enddo   
end subroutine gatter3


subroutine hfuncpos(spos,h,dens,ddensdh,fi,dfi,nneigh,hneigh)
 include 'globals.h'
  integer p,nneigh,j
  real fi,dfi,hneigh,h,spos(3),dens,ddensdh
  if(nneigh.EQ.0.OR.h.GT.hneigh.OR.h.LT.hneigh/2) then
   nneigh=0
   hneigh=h
   call search(root,2*h,spos,nneigh,bodlist)
  endif
  call gatter2(nneigh,spos,h,dens,ddensdh)
   fi=4*pi/3*h**3*(dens+rhomin)-massres
   dfi=4*pi*h**2*(dens+rhomin)+4*pi/3*h**3*ddensdh
  return
end

subroutine hnewtonpos(spos,h,success,nneigh,hneigh)
include 'globals.h'
integer p,success,nmax,nneigh,i
real hsmnw,hmin,hmax,prec,hneigh,fi,dfi,h,dens,ddensdh
real spos(3)
parameter(nmax=10,prec=0.001) 
 success=1
 hmin=.5*h;hmax=2*h 
 hsmnw=h
 i=0
10 if(hsmnw.LT.hmax.AND.hsmnw.GT.hmin.AND.i.LT.nmax) then
  i=i+1
  call hfuncpos(spos,hsmnw,dens,ddensdh,fi,dfi,nneigh,hneigh)
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
end 


subroutine brackethpos(spos,hmin,hmax,success,nneigh,hneigh)
 include 'globals.h'
 integer p,success,nneigh,ntry,j
 real hmin,hmax,hneigh,spos(3),fac,f1,f2,df,sethsmooth
 real dum1,dum2
 parameter(fac=2.,ntry=20)
 external sethsmooth
 success=1
 hmax=fac*sethsmooth(spos)
 hmin=hmax/fac**2
 nneigh=0
 hneigh=0
 call hfuncpos(spos,hmax,dum1,dum2,f2,df,nneigh,hneigh)
 call hfuncpos(spos,hmin,dum1,dum2,f1,df,nneigh,hneigh) 
 do j=1,ntry
  if(f1*f2.LT.0) return
  if(abs(f1).lt.abs(f2)) then
   hmin=hmin/fac
   call hfuncpos(spos,hmin,dum1,dum2,f1,df,nneigh,hneigh)
  else
   hmax=hmax*fac
   call hfuncpos(spos,hmax,dum1,dum2,f2,df,nneigh,hneigh)
  endif
 enddo 
success=0
return
end 

subroutine hsafepos(spos,hguess,hmin,hmax,nneigh,hneigh)
include 'globals.h'
integer p,maxit,j,nneigh
real hmin,hmax,prec,hneigh,hguess,dum1,dum2
real df,dx,dxold,f,h,fh,fl,temp,xl,xh,spos(3)
parameter(maxit=25,prec=0.001)
 
nneigh=0
hneigh=0
call hfuncpos(spos,hmax,dum1,dum2,fh,df,nneigh,hneigh)
call hfuncpos(spos,hmin,dum1,dum2,fl,df,nneigh,hneigh)
 
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
 call hfuncpos(spos,hguess,dum1,dum2,f,df,nneigh,hneigh)
 do j=1,maxit
  if(((hguess-xh)*df-f)*((hguess-xl)*df-f).GT.0.OR. &
&  abs(2*f).GT.ABS(dxold*df)) then
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
  call hfuncpos(spos,hguess,dum1,dum2,f,df,nneigh,hneigh) 
  if(f.lt.0) then
   xl=hguess    
  else
   xh=hguess
  endif  
 enddo
 return 
end 

   subroutine hsmdenspos(spos,h,dens,ddensdh,nneigh,hneigh)
   include 'globals.h'
   integer p,i,j,success,nneigh
   real hsmn,hsmorg,fi,dfi,hmin,hmax,hneigh
   real dens,ddensdh,h,spos(3),getlocalscale
   i=0;j=0
   if(h.LE.0) h=getlocalscale(spos)
   hsmn=h
   nneigh=0
   hneigh=0
   call hnewtonpos(spos,hsmn,success,nneigh,hneigh)
   if(success.EQ.0) then 
   call brackethpos(spos,hmin,hmax,success,nneigh,hneigh) 
if(success.EQ.0) then
!$omp critical
call terror(' hsmdenspos error')
!$omp end critical
    endif
    call hsafepos(spos,hsmn,hmin,hmax,nneigh,hneigh)
   endif
   h=hsmn
   call hfuncpos(spos,h,dens,ddensdh,fi,dfi,nneigh,hneigh)
   end

  subroutine hsmpos(spos,h,nneigh,hneigh)
   include 'globals.h'
   integer p,i,success,nneigh
   real hsmn,fi,dfi,hmin,hmax,hneigh
   real dens,ddensdh,h,spos(3)
i=0
if(h.LE.0) h=consthsm
hsmn=h
nneigh=0
hneigh=0
   call hnewtonpos(spos,hsmn,success,nneigh,hneigh)
if(success.EQ.0) then
 call brackethpos(spos,hmin,hmax,success,nneigh,hneigh)
 call hsafepos(spos,hsmn,hmin,hmax,nneigh,hneigh)
endif
h=hsmn
end

subroutine denspos(spos,h,dens,ddensdh,newroot,newres)
  include 'globals.h'
  real spos(3),h,dens,ddensdh,newres,tmp2,hneigh,dum2,dum3
  integer newroot,tmp,nneigh
  
  tmp=root
  tmp2=massres
  root=newroot
  massres=newres
  call hsmpos(spos,h,nneigh,hneigh)
  call hfuncpos(spos,h,dens,ddensdh,dum2,dum3,nneigh,hneigh)
  root=tmp
  massres=tmp2
 end subroutine denspos

function setlocallength(x,newroot)
 include 'globals.h'
 real,dimension(3),intent(in) :: x
 integer,intent(in) :: newroot
 integer :: tmp,tmp2
 real :: sethsmooth,setlocallength
 
 tmp=root
 root=newroot
 setlocallength=sethsmooth(x)
 root=tmp
 return
 end function setlocallength

subroutine fuvfuncpos(spos,h,dens,ddensdh,nneigh,hneigh)
include 'globals.h'
integer p,nneigh,j
real fi,dfi,hneigh,h,spos(3),dens,ddensdh
if(nneigh.EQ.0.OR.h.GT.hneigh.OR.h.LT.hneigh/2) then
 nneigh=0
 hneigh=h
 call search(root,2*h,spos,nneigh,bodlist)
endif
call gatterweight(nneigh,spos,h,dens,ddensdh,starfuv)
return
end

 subroutine fuvpos(spos,h,dens,ddensdh,newroot,newres)
  include 'globals.h'
  real spos(3),h,dens,ddensdh,newres,tmp2,hneigh
  integer newroot,tmp,nneigh
  
  tmp=root
  tmp2=massres
  root=newroot
  massres=newres
  call hsmpos(spos,h,nneigh,hneigh)
  call fuvfuncpos(spos,h,dens,ddensdh,nneigh,hneigh)
  root=tmp
  massres=tmp2
 end subroutine fuvpos
  
 subroutine scattereps(n,spos)
 include 'globals.h' 
 integer,intent(in) :: n
 real,dimension(3),intent(in) :: spos
 integer :: i,p,iwsm
 real :: dx,dy,dz,wnorm,distnorm,hsminv,dr2p,drw,tmass,dr2,wsm,gatter
 tmass=0.
 if(n.EQ.0) return
    
 do i=1,n
  p=bodlist(i)
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
  hsminv=1./epsgrav(p)
  wnorm=piinv*hsminv*hsminv*hsminv
  distnorm=hsminv**2*deldr2i
  dr2p=dr2*distnorm
  iwsm=INT(dr2p)
  if(ninterp.LT.iwsm) iwsm=ninterp
  drw=dr2p-iwsm
  wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
  tempvect(i)=wnorm*wsm
  tmass=tmass+mass(p)*tempvect(i)
 enddo
 gatter=tmass
   
end subroutine

 
 subroutine emissdenspos(spos,h,dens,ddensdh,emiss,newroot,newres)
  use CoolingMod
  use ElementsMod
  include 'globals.h'
  real spos(3),h,dens,ddensdh,newres,tmp2,hneigh,dum2,dum3,emiss,t
  integer newroot,tmp,nneigh,i,p
  
  tmp=root
  tmp2=massres
  root=newroot
  massres=newres
  call hsmpos(spos,h,nneigh,hneigh)
  call hfuncpos(spos,h,dens,ddensdh,dum2,dum3,nneigh,hneigh)
  emiss=0
  t=0
  do i=1,nneigh
  p=bodlist(i)
   t=t+tempvect(i)*mass(p)*ethermal(p)
  enddo
   t=t/dens
   if(t.lt.20000) then
   if(t.gt.10000) t=10000
   emiss=abundanceQ(3)*elementCool(0.1,t,3)*dens**2*cool_par 
   endif
  root=tmp
  massres=tmp2
 end subroutine emissdenspos


function sethsmooth(position)
include 'globals.h'
real,dimension(3),intent(in) :: position
real :: hsearch,getlocalscale,hs
real :: hmin,hmax,sethsmooth
integer :: inear,i

 i=0
 hsearch=getlocalscale(position)
 hmin=0.;hmax=hsearch/2.;inear=0
 do while(inear.LE.nsmooth)
  i=i+1
  inear=0
  call search(root,2*hmax,position,inear,bodlist)
  hmax=hmax*1.3
 enddo
 hsearch=hmax/1.3
1 continue
 inear=0
 i=i+1
 call search(root,2*hsearch,position,inear,bodlist)
 if(inear.LT.(1-nsmtol)*nsmooth.OR.inear.GT.(1+nsmtol)*nsmooth) then
  if(inear.LT.nsmooth) then
   hmin=hsearch
   hsearch=0.5*(hsearch+hmax)
  endif
  if(inear.GT.nsmooth) then
   hmax=hsearch
   hsearch=0.5*(hsearch+hmin)
  endif
  goto 1
 endif
 sethsmooth=hsearch
 return
end function sethsmooth

