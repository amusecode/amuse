!> \file kerneldensity.F90

!> \brief Module for calculating the kernel density.
!<

module kernel_density_mod
use particle_system_mod, only: particle_system_type, particle_type
use sphpar_mod, only: sphpar_type 
use nblist_mod
implicit none

 private 
! public :: localpart_type
 public :: initdensity,gradrho,density,hsmdens,fromsphpar

 public :: massres,rhomin,spherevol

 integer,parameter :: ndim=3

 real :: massres
 real :: rhomin
 integer, parameter :: TARGETNN=32
 integer, parameter :: Ninterpolation=1000
 integer, parameter :: HSM_NEWTON_MAX_ITERATION=10
 real,    parameter :: HSM_NEWTON_PREC=0.0001 
 integer, parameter :: HSM_BISECTION_MAX_ITERATION=50
 real,    parameter :: HSM_BISECTION_PREC=0.001 
 integer, parameter :: HSM_BRACKET_TRY=50
 real,    parameter :: HSM_BRACKET_FAC=2. 
 real,    parameter :: Pi=3.141592653589793238
 real,    parameter :: kfac(3)=(/ 2.6666666666666666667,   &
                                  1.8189136353359466945,   &
                                  2.5464790894703253723 /)
 real :: spherevol,dspherevol

 real :: wsmooth(0:Ninterpolation+1),dwsmooth(0:Ninterpolation+1) 
 real, parameter :: deldr2=1./Ninterpolation

 integer :: nnewton=0,nbisect=0

 contains

 subroutine initdensity(mass)
   real :: mass
     massres=TARGETNN*mass
     rhomin=0.0
     call initkernel
     spherevol=UnitSphereVolume(ndim)
     dspherevol=spherevol*ndim
 end subroutine initdensity


 subroutine initkernel
   integer :: i
   real :: xw2,xw
     do i=0,1+Ninterpolation
        xw2=i*deldr2*4.
        xw=SQRT(xw2)
        if(xw2.LE.1.) then
           wsmooth(i)=1.-1.5*xw2+0.75*xw*xw2
           dwsmooth(i)=-3.+2.25*xw
        else
           wsmooth(i)=0.25*(2.-xw)**3
           dwsmooth(i)=-0.75*(2.-xw)**2/xw
        endif
        if(xw2.GE.4.) then
           wsmooth(i)=0.
           dwsmooth(i)=0.
        endif
        dwsmooth(i)=dwsmooth(i)*4 ! geloofik
     enddo
 end subroutine initkernel


 subroutine hsmdens(sphpar,psys,nblist)
 use oct_tree_mod, only: getlocalscale
   type(particle_system_type) :: psys
   type(sphpar_type) :: sphpar
   type(nblist_type) :: nblist
!    try if part%hsmooth=0.0 then get local scale
     if (sphpar%hsmooth==0.0) then
        sphpar%hsmooth = getlocalscale(sphpar%pos,nblist%stree)
     end if
     call resetnblist(nblist,sphpar%hsmooth,sphpar%pos)
     call hnewton(sphpar,psys,nblist)
     call hfunc(sphpar,psys,nblist)
 end subroutine hsmdens


 subroutine hnewton(sphpar,psys,nblist)
   type(particle_system_type) :: psys
   type(sphpar_type) :: sphpar
   type(nblist_type) :: nblist
   real :: hsmold,hmin,hmax
   integer :: i
     hmin=0.5*sphpar%hsmooth
     hmax=2*sphpar%hsmooth	 
     hsmold=sphpar%hsmooth
     nnewton=nnewton+1
     i=0
10   if(sphpar%hsmooth.LT.hmax .AND. & 
        sphpar%hsmooth.GT.hmin .AND. & 
        i.LT.HSM_NEWTON_MAX_ITERATION) then
        i=i+1
        call hfunc(sphpar,psys,nblist)
        if(sphpar%dfi.NE.0) then
           sphpar%hsmooth=sphpar%hsmooth-sphpar%fi/sphpar%dfi
!   print*,'nw',sphpar%hsmooth
        else
           sphpar%hsmooth=0
        endif
        if(ABS(sphpar%fi).LT.ABS(sphpar%hsmooth*HSM_NEWTON_PREC*sphpar%dfi)) then
           return
        endif
        goto 10
     else
        sphpar%hsmooth=hsmold
        call bisect(sphpar,psys,nblist)
     endif
 end subroutine hnewton



 subroutine hfunc(sphpar,psys,nblist)
   type(particle_system_type) :: psys
   type(sphpar_type) :: sphpar
   type(nblist_type) :: nblist
   
     if(nblist%nnb.LE.0.OR. &
        sphpar%hsmooth.GT.nblist%searchrange.OR. &
        sphpar%hsmooth.LT.0.5*nblist%searchrange) then
        call resetnblist(nblist,range=sphpar%hsmooth)
     endif
     call fidfi(sphpar,psys,nblist)
     return
 end subroutine hfunc

 subroutine do_nbsum(sphpar,psys,nblist,sumroutine)
   type(sphpar_type) :: sphpar
   type(particle_system_type) :: psys
   type(nblist_type) :: nblist

   interface
      subroutine sumroutine(pa,nb)
        use sphpar_mod, only: sphpar_type
        use nblist_mod, only: nblist_type
        type(sphpar_type) :: pa
        type(nblist_type) :: nb
      end subroutine sumroutine
   end interface
   
   if(nblist%reuseable) then
      call sumroutine(sphpar,nblist)  
      return
   endif
   call resetnblist(nblist)
   do while(nblist%searchcell.NE.0)
      nblist%nnb=0
      call fullsearch(psys,nblist)
      call sumroutine(sphpar,nblist)
   enddo
 end subroutine do_nbsum
 
subroutine gradrhosum(sphpar,nblist)
  type(sphpar_type) :: sphpar
  type(nblist_type)  :: nblist
  integer iwsm,i
  real :: dx(ndim),hinv,dknorm,dist2norm,dr2p,dist2,drw
  real :: dkernel
  real :: particle_mass
  hinv=1./sphpar%hsmooth
  dknorm=kfac(ndim)*hinv**(ndim+2)
  dist2norm=hinv**2/deldr2
  do i=1,nblist%nnb
   dx=sphpar%pos-nblist%sphpar(i)%pos
   dist2=sum(dx**2)
   dr2p=dist2*dist2norm
   if(Ninterpolation.GE.dr2p) then
    iwsm=INT(dr2p)                 
    drw=dr2p-iwsm
    dkernel=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
    dkernel=dknorm*dkernel

    particle_mass = nblist%sphpar(i)%mass

    sphpar%gradrho=sphpar%gradrho+dx*particle_mass*dkernel                 
   endif
  enddo
end subroutine gradrhosum


 subroutine gradrho(sphpar,nblist)
  type(sphpar_type) :: sphpar
  type(nblist_type) :: nblist

  type(particle_system_type) psys

  sphpar%gradrho=0
  call do_nbsum(sphpar,psys,nblist,gradrhosum)
 end subroutine gradrho


 subroutine densitysum(sphpar,nblist)

  type(sphpar_type) sphpar
  type(nblist_type) nblist
  integer iwsm, i
  real :: dx(ndim), hinv, knorm, dknorm, dist2norm, dr2p, dist2, drw
  real :: kernel
  real :: particle_mass

!  real :: dkernel

  hinv=1./sphpar%hsmooth
  knorm=kfac(ndim)*hinv**ndim
  dknorm=kfac(ndim)*hinv**(ndim+2)
  dist2norm=hinv**2/deldr2
  do i=1,nblist%nnb
     dx=sphpar%pos-nblist%sphpar(i)%pos
     dist2=sum(dx**2)
     dr2p=dist2*dist2norm
     if(Ninterpolation.GE.dr2p) then
        sphpar%nnb=sphpar%nnb+1
        iwsm=INT(dr2p)                 
        drw=dr2p-iwsm
        kernel=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
        kernel=knorm*kernel
        
        particle_mass = nblist%sphpar(i)%mass
    
        sphpar%rho=sphpar%rho+particle_mass*kernel
        sphpar%xHII = sphpar%xHII + particle_mass * nblist%sphpar(i)%xHII * kernel
        sphpar%T    = sphpar%T    + particle_mass * nblist%sphpar(i)%T    * kernel

     endif
  enddo

 end subroutine densitysum

 subroutine density(sphpar,psys,nblist)
  type(particle_system_type) psys
  type(sphpar_type) sphpar
  type(nblist_type) nblist
  sphpar%nnb=0
  sphpar%rho=0
  call do_nbsum(sphpar,psys,nblist,densitysum)  
  if(sphpar%rho.NE.0) sphpar%xHII = sphpar%xHII / sphpar%rho
  if(sphpar%rho.NE.0) sphpar%T    = sphpar%T    / sphpar%rho

  ! I added this to see about the vertical stripes
!  if(sphpar%rho.LT.rhomin) then
!     sphpar%rho = rhomin
!     sphpar%xHII = 1.00
!     sphpar%T = 1.0e4
!  end if


 end subroutine density


 subroutine fidfisum(sphpar,nblist)
  type(sphpar_type) sphpar
  type(nblist_type) nblist
  integer iwsm,i
  real :: dx(ndim),hinv, knorm,dknorm,dist2norm,dr2p,dist2,drw
  real :: kernel,dkernel
  real :: particle_mass
  hinv=1./sphpar%hsmooth
  knorm=kfac(ndim)*hinv**ndim
  dknorm=kfac(ndim)*hinv**(ndim+2)
  dist2norm=hinv**2/deldr2
  do i=1,nblist%nnb
     dx=sphpar%pos-nblist%sphpar(i)%pos
     dist2=sum(dx**2)
     dr2p=dist2*dist2norm
     if(dr2p.LT.0) then
        print*,'!!!!',dr2p
        stop
     endif
     if(Ninterpolation.GE.dr2p) then
        sphpar%nnb=sphpar%nnb+1
        iwsm=INT(dr2p)                 
        drw=dr2p-iwsm
        kernel=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
        dkernel=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
        kernel=knorm*kernel
        dkernel=dknorm*dkernel
        
        particle_mass = nblist%sphpar(i)%mass

        sphpar%rho=sphpar%rho+particle_mass*kernel
        sphpar%drhodh=sphpar%drhodh-ndim*hinv*particle_mass*kernel
        sphpar%drhodh=sphpar%drhodh-dist2*hinv*particle_mass*dkernel 
     endif
  enddo
 end subroutine fidfisum

 subroutine fidfi(sphpar,psys,nblist)
   type(particle_system_type) psys
   type(sphpar_type) sphpar
   type(nblist_type) nblist
     sphpar%nnb=0
     sphpar%rho=0
     sphpar%drhodh=0
     call do_nbsum(sphpar,psys,nblist,fidfisum)
     sphpar%fi=spherevol*sphpar%hsmooth**ndim*(sphpar%rho+rhomin)-massres
     sphpar%dfi=spherevol*ndim*sphpar%hsmooth**(ndim-1)*(sphpar%rho+rhomin)+ &
              spherevol*sphpar%hsmooth**ndim*sphpar%drhodh
 end subroutine fidfi



 subroutine bisect(sphpar,psys,nblist)
   type(particle_system_type) :: psys
   type(sphpar_type) :: sphpar
   type(nblist_type) :: nblist
   real :: hmin,hmax
  
     nbisect=nbisect+1
     call bracketh(sphpar,psys,nblist,hmin,hmax)
     call hsafe(sphpar,psys,nblist,hmin,hmax)  
     call hfunc(sphpar,psys,nblist)
 end subroutine bisect

 subroutine hsafe(sphpar,psys,nblist,hmin,hmax)
  type(particle_system_type) :: psys
  type(sphpar_type) :: sphpar
  type(nblist_type) :: nblist
  integer :: j
  real :: hmin, hmax
  real :: dx, dxold, fh, fl, temp, xl, xh  

!  integer :: p, maxit, nneigh
!  real :: prec, hneigh

  sphpar%hsmooth=hmax
  call hfunc(sphpar,psys,nblist)
  fh=sphpar%fi
  sphpar%hsmooth=hmin
  call hfunc(sphpar,psys,nblist)
  fl=sphpar%fi
  if(hmin.GE.hmax.OR.fl.GE.fh) call kernelError('hsafe error')
  if(fl.eq.0) then
   sphpar%hsmooth=hmin
   return
   if(fh.eq.0) then
    sphpar%hsmooth=hmax
   return
   endif
  endif 
  xl=hmin
  xh=hmax	 
  sphpar%hsmooth=.5*(xl+xh)
  dxold=(xh-xl)
  dx=dxold
  call hfunc(sphpar,psys,nblist)
  do j=1,HSM_BISECTION_MAX_ITERATION
   if( ((sphpar%hsmooth-xh)*sphpar%dfi - sphpar%fi) * &
       ((sphpar%hsmooth-xl)*sphpar%dfi - sphpar%fi) .GT. 0 .OR. &
       abs(2*sphpar%fi) .GT. ABS(dxold*sphpar%dfi) ) then
	    dxold=dx
	    dx=0.5*(xh-xl)
	    sphpar%hsmooth=xl+dx
	    if(xl.EQ.sphpar%hsmooth)return
   else
	    dxold=dx
	    dx=sphpar%fi/sphpar%dfi
	    temp=sphpar%hsmooth
	    sphpar%hsmooth=temp-dx
	    if(temp.EQ.sphpar%hsmooth)return
   endif
   if(abs(dx).LT.HSM_BISECTION_PREC*hmax) return
   call hfunc(sphpar,psys,nblist)
   if(sphpar%fi.lt.0) then 
    xl=sphpar%hsmooth
   else    
    xh=sphpar%hsmooth
   endif
  enddo
  return	 
 end subroutine hsafe
	 
 subroutine bracketh(sphpar,psys,nblist,hmin,hmax)
 type(particle_system_type) psys
 type(sphpar_type) sphpar
 type(nblist_type) nblist
 integer j
 real hmin,hmax,f1,f2
  hmax=sphpar%hsmooth
  hmin=hmax/HSM_BRACKET_FAC
  call hfunc(sphpar,psys,nblist)
  f2=sphpar%fi
  sphpar%hsmooth=hmin
  call hfunc(sphpar,psys,nblist)
  f1=sphpar%fi
  do j=1,HSM_BRACKET_TRY
   if(f1*f2.LT.0) return
   if(abs(f1).lt.abs(f2)) then
    hmin=hmin/HSM_BRACKET_FAC
    sphpar%hsmooth=hmin
    call hfunc(sphpar,psys,nblist)
    f1=sphpar%fi    
   else
    hmax=hmax*HSM_BRACKET_FAC
    sphpar%hsmooth=hmax
    call hfunc(sphpar,psys,nblist)
    f2=sphpar%fi
   endif
  enddo 
  print*,sphpar,hmin,hmax
  call kernelError('hsm bracketing')
  return
 end subroutine bracketh
  
 subroutine kernelError(string,i)
  character(len=*) :: string
  integer, optional :: i
  
  print*,'error detected'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif
  stop
 end subroutine kernelError

 function UnitSphereVolume(m) result(vol)
 real :: vol
 integer,intent(in) :: m
  if(mod(m,2).eq.0) then
   vol=Pi**(m/2)/fac(m/2)
  else
   vol=2**((m+1)/2)*Pi**((m-1)/2)/fac2(m) 
  endif
 end function UnitSphereVolume

 function fac(n)
 integer :: fac
 integer,intent(in) :: n
 integer :: i
 fac=1
 do i=1,n
    fac=fac*i
 enddo
 end function fac
 
 recursive function fac2(n) result(f2)
   integer :: f2
   integer,intent(in) :: n
     f2=1
     if(n.le.0) return
     if(mod(n,2).eq.0) then 
        f2=fac(n/2)*2**(n/2) 
     else 
        f2=fac(n)/(fac2(n-1))
     endif
 end function fac2

 function fromsphpar(part)
   type(particle_type) :: part
   type(sphpar_type) :: fromsphpar
     fromsphpar%pos=part%pos
     fromsphpar%rho=massres/spherevol/part%hsml**ndim-rhomin
     fromsphpar%nnb=0
     fromsphpar%gradrho=0.
     fromsphpar%drhodh=0.
     fromsphpar%hsmooth=part%hsml
     fromsphpar%fi=0.
     fromsphpar%dfi=0.
 end function fromsphpar



end module kernel_density_mod


