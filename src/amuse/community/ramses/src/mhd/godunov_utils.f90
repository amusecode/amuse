!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector),save::a2,B2,rho,ctot
  
  real(dp)::dtcell,smallp,cf,cc,bc,bn
  integer::k,idim
  
  smallp = smallr*smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
     rho(k)=uu(k,1)
  end do
  do idim = 1,3
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/rho(k)
     end do
  end do

  do k = 1,ncell
     B2(k)=zero
  end do
  do idim = 1,3
     do k = 1, ncell
        Bc = half*(uu(k,5+idim)+uu(k,nvar+idim))
        B2(k)=B2(k)+Bc**2
        uu(k,5) = uu(k,5)-half*uu(k,1)*uu(k,idim+1)**2-half*Bc**2
     end do
  end do

  ! Compute thermal sound speed
  do k = 1, ncell
     uu(k,5) = max((gamma-one)*uu(k,5),smallp)
     a2(k)=gamma*uu(k,5)/uu(k,1)
  end do

  ! Compute maximum wave speed (fast magnetosonic)
  do k = 1, ncell
     ctot(k)=zero
  end do
  if(scheme.eq.'induction')then
     do idim = 1,ndim   ! WARNING: ndim instead of 3  
        do k = 1, ncell
           ctot(k)=ctot(k)+abs(uu(k,idim+1))
        end do
     end do
  else
     do idim = 1,ndim   ! WARNING: ndim instead of 3  
        do k = 1, ncell
           cc=half*(B2(k)/rho(k)+a2(k))
           BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
           cf=sqrt(cc+sqrt(cc**2-a2(k)*BN**2/rho(k)))
           ctot(k)=ctot(k)+abs(uu(k,idim+1))+cf
        end do
     end do
  endif

  ! Compute gravity strength ratio
  do k = 1, ncell
     rho(k)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        rho(k)=rho(k)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     rho(k)=rho(k)*dx/ctot(k)**2
     rho(k)=MAX(rho(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell=dx/ctot(k)*(sqrt(one+two*courant_factor*rho(k))-one)/rho(k)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,nn,ilevel)
  use amr_parameters
  use hydro_parameters
  use amr_commons, ONLY: emag_tot
  use const
  implicit none
  ! dummy arguments
  integer nn,ilevel
  real(dp)::ug(1:nvector,1:nvar+3)
  real(dp)::um(1:nvector,1:nvar+3)
  real(dp)::ud(1:nvector,1:nvar+3)
  logical ::ok(1:nvector)
  
  integer::k,idim
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp),dimension(1:nvector),save::emagg,emagm,emagd
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error,emag_loc,ethres
  
  ! Convert to primitive variables
  do k = 1,nn
     ug(k,1) = max(ug(k,1),smallr)
     um(k,1) = max(um(k,1),smallr)
     ud(k,1) = max(ud(k,1),smallr)
  end do
  ! Velocity
  do idim = 1,3
     do k = 1,nn
        ug(k,idim+1) = ug(k,idim+1)/ug(k,1)
        um(k,idim+1) = um(k,idim+1)/um(k,1)
        ud(k,idim+1) = ud(k,idim+1)/ud(k,1)
     end do
  end do
  ! Pressure
  do k = 1,nn
     eking(k) = zero
     ekinm(k) = zero
     ekind(k) = zero
  end do
  do idim = 1,3
     do k = 1,nn
        eking(k) = eking(k) + half*ug(k,1)*ug(k,idim+1)**2
        ekinm(k) = ekinm(k) + half*um(k,1)*um(k,idim+1)**2
        ekind(k) = ekind(k) + half*ud(k,1)*ud(k,idim+1)**2
     end do
  end do
  do k = 1,nn
     emagg(k) = zero
     emagm(k) = zero
     emagd(k) = zero
  end do
  do idim = 1,3
     do k = 1,nn
        emagg(k) = emagg(k) + half*(half*(ug(k,idim+5)+ug(k,idim+nvar)))**2
        emagm(k) = emagm(k) + half*(half*(um(k,idim+5)+um(k,idim+nvar)))**2
        emagd(k) = emagd(k) + half*(half*(ud(k,idim+5)+ud(k,idim+nvar)))**2
     end do
  end do
  do k = 1,nn
     ug(k,5) = (gamma-one)*(ug(k,5)-eking(k)-emagg(k))
     um(k,5) = (gamma-one)*(um(k,5)-ekinm(k)-emagm(k))
     ud(k,5) = (gamma-one)*(ud(k,5)-ekind(k)-emagd(k))
  end do  
  ! Passive scalars
#if NVAR > 8
  do idim = 9,nvar
     do k = 1,nn
        ug(k,idim) = ug(k,idim)/ug(k,1)
        um(k,idim) = um(k,idim)/um(k,1)
        ud(k,idim) = ud(k,idim)/ud(k,1)
     end do
  end do
#endif

  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,nn
        dg=ug(k,1); dm=um(k,1); dd=ud(k,1)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm+floor_d)) , &
             & ABS((dm-dg)/(dm+dg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
     do k=1,nn
     end do
  end if

  if(err_grad_p >= 0.)then
     do k=1,nn
        pg=ug(k,5); pm=um(k,5); pd=ud(k,5)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_p)), &
             & ABS((pm-pg)/(pm+pg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if

  if(err_grad_b2 >= 0.)then
     do k=1,nn
        pg=emagg(k); pm=emagm(k); pd=emagd(k)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_b2)), &
             & ABS((pm-pg)/(pm+pg+floor_b2)) )
        ok(k) = ok(k) .or. error > err_grad_b2
     end do
  end if

  if(err_grad_A >= 0.)then
     idim = 1
     do k=1,nn
        vg=0.5*(ug(k,5+idim)+ug(k,nvar+idim))
        vm=0.5*(um(k,5+idim)+um(k,nvar+idim))
        vd=0.5*(ud(k,5+idim)+ud(k,nvar+idim))
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=2.0d0*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_A)) , &
             & ABS((vm-vg)/(cm+cg+floor_A)) )
        ok(k) = ok(k) .or. error > err_grad_A
     end do
  end if

  if(err_grad_B >= 0.)then
     idim = 2
     do k=1,nn
        vg=0.5*(ug(k,5+idim)+ug(k,nvar+idim))
        vm=0.5*(um(k,5+idim)+um(k,nvar+idim))
        vd=0.5*(ud(k,5+idim)+ud(k,nvar+idim))
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=2.0d0*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_B)) , &
             & ABS((vm-vg)/(cm+cg+floor_B)) )
        ok(k) = ok(k) .or. error > err_grad_B
     end do
  end if

  if(err_grad_C >= 0.)then
     idim = 3
     do k=1,nn
        vg=0.5*(ug(k,5+idim)+ug(k,nvar+idim))
        vm=0.5*(um(k,5+idim)+um(k,nvar+idim))
        vd=0.5*(ud(k,5+idim)+ud(k,nvar+idim))
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=2.0d0*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_C)) , &
             & ABS((vm-vg)/(cm+cg+floor_C)) )
        ok(k) = ok(k) .or. error > err_grad_C
     end do
  end if

  if(err_grad_u >= 0.)then
     do idim = 1,3
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
           cg=sqrt(max(gamma*ug(k,5)/ug(k,1),floor_u**2))
           cm=sqrt(max(gamma*um(k,5)/um(k,1),floor_u**2))
           cd=sqrt(max(gamma*ud(k,5)/ud(k,1),floor_u**2))
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

  if(scheme.eq.'induction')then
     if(m_refine(ilevel) >= 0.)then
        if(emag_tot==0.0)then
           emag_loc=1.5
        else
           emag_loc=emag_tot
        endif
        ethres=m_refine(ilevel)*emag_loc/boxlen**3
        do k=1,nn
           ok(k)=ok(k).or. emagm(k) > ethres
        end do
     endif
  endif

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE upwind(qleft,qright,fgdnv,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D Upwind Riemann solver
  IMPLICIT NONE  
  REAL(dp)::zero_flux
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv

  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean
  INTEGER ::l
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright)
  
  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux
  
  ! find the mean normal velocity
  vleft = half * ( qleft(3) + qright(3) )
  
  ! difference between the 2 states
  udiff = half * ( uright - uleft )
  
  ! the Upwind flux
  fgdnv = fmean - ABS(vleft) * udiff

END SUBROUTINE upwind
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE lax_friedrich(qleft,qright,fgdnv,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D local Lax-Friedrich Riemann solver
  IMPLICIT NONE  
  REAL(dp)::zero_flux
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv

  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean
  INTEGER ::l
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright) 
  
  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux
  
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_info(qleft ,vleft )
  CALL find_speed_info(qright,vright)
  
  ! difference between the 2 states
  udiff  = half * ( uright - uleft )
  
  ! the local Lax-Friedrich flux
  fgdnv = fmean - MAX(vleft,vright) * udiff

END SUBROUTINE lax_friedrich
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hll(qleft,qright,fgdnv)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D HLL Riemann solver
  IMPLICIT NONE  
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv
  REAL(dp),DIMENSION(1:nvar)::fleft,fright
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean,cfleft,cfright,SL,SR
  INTEGER ::l
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright) 
  
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_fast(qleft ,cfleft )
  CALL find_speed_fast(qright,cfright)
  vleft =qleft (3)
  vright=qright(3)
  SL=min(min(vleft,vright)-max(cfleft,cfright),zero)
  SR=max(max(vleft,vright)+max(cfleft,cfright),zero)

  ! the HLL flux
  fgdnv = (SR*fleft-SL*fright+SR*SL*(uright-uleft))/(SR-SL)

END SUBROUTINE hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hlld(qleft,qright,fgdnv)
  USE hydro_parameters, ONLY: nvar, gamma
  USE const
  ! HLLD Riemann solver (Miyoshi & Kusano, 2005, JCP, 208, 315)
  IMPLICIT NONE
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv
  REAL(dp)::SL,SR,SAL,SAR
  REAL(dp)::entho,A,sgnm
  REAL(dp)::rl,pl,ul,vl,wl,cl,ecinl,emagl,etotl,ptotl,vdotbl,bl,el
  REAL(dp)::rr,pr,ur,vr,wr,cr,ecinr,emagr,etotr,ptotr,vdotbr,br,er
  REAL(dp)::cfastl,calfvenl,rcl,rstarl,vstarl,wstarl,bstarl,cstarl,vdotbstarl
  REAL(dp)::cfastr,calfvenr,rcr,rstarr,vstarr,wstarr,bstarr,cstarr,vdotbstarr
  REAL(dp)::sqrrstarl,etotstarl,etotstarstarl
  REAL(dp)::sqrrstarr,etotstarr,etotstarstarr
  REAL(dp)::ustar,ptotstar,estar,vstarstar,wstarstar,bstarstar,cstarstar,vdotbstarstar
  REAL(dp)::ro,uo,vo,wo,bo,co,ptoto,etoto,vdotbo

  INTEGER ::ivar
  
  entho = one/(gamma-one)

  ! Enforce continuity of normal component
  A=half*(qleft(4)+qright(4))
  sgnm=sign(one,A)
  qleft(4)=A; qright(4)=A

  ! Left variables
  rl=qleft(1); Pl=qleft(2); ul=qleft(3)
  vl=qleft(5); Bl=qleft(6); wl=qleft(7); Cl=qleft(8)
  ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
  etotl = Pl*entho+ecinl+emagl
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*cl

  ! Right variables
  rr=qright(1); Pr=qright(2); ur=qright(3)
  vr=qright(5); Br=qright(6); wr=qright(7); Cr=qright(8)
  ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
  etotr = Pr*entho+ecinr+emagr
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr

  ! Find the largest eigenvalues in the normal direction to the interface
  CALL find_speed_fast(qleft ,cfastl)
  CALL find_speed_fast(qright,cfastr)

  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)
!  SL=ul-cfastl
!  SR=ur+cfastr
  
  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  ! Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  estar =rl*(SL-ul)*(SL-ustar)-A**2
  el    =rl*(SL-ul)*(SL-ul   )-A**2
  if(estar==0)then
     vstarl=vl
     Bstarl=Bl
     wstarl=wl
     Cstarl=Cl
  else
     vstarl=vl-A*Bl*(ustar-ul)/estar
     Bstarl=Bl*el/estar
     wstarl=wl-A*Cl*(ustar-ul)/estar
     Cstarl=Cl*el/estar
  endif
  vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
  sqrrstarl=sqrt(rstarl)
  calfvenl=abs(A)/sqrrstarl
  SAL=ustar-calfvenl

  ! Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  estar =rr*(SR-ur)*(SR-ustar)-A**2
  er    =rr*(SR-ur)*(SR-ur   )-A**2
  if(estar==0)then
     vstarr=vr
     Bstarr=Br
     wstarr=wr
     Cstarr=Cr
  else
     vstarr=vr-A*Br*(ustar-ur)/estar
     Bstarr=Br*er/estar
     wstarr=wr-A*Cr*(ustar-ur)/estar
     Cstarr=Cr*er/estar
  endif
  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  ! Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
  
  ! Sample the solution at x/t=0
  if(SL>0d0)then
     ro=rl
     uo=ul
     vo=vl
     wo=wl
     Bo=Bl
     Co=Cl
     Ptoto=Ptotl
     etoto=etotl
     vdotBo=vdotBl
  else if(SAL>0d0)then
     ro=rstarl
     uo=ustar
     vo=vstarl
     wo=wstarl
     Bo=Bstarl
     Co=Cstarl
     Ptoto=Ptotstar
     etoto=etotstarl
     vdotBo=vdotBstarl
  else if(ustar>0d0)then
     ro=rstarl
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarl
     vdotBo=vdotBstarstar
  else if(SAR>0d0)then
     ro=rstarr
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarr
     vdotBo=vdotBstarstar
  else if (SR>0d0)then
     ro=rstarr
     uo=ustar
     vo=vstarr
     wo=wstarr
     Bo=Bstarr
     Co=Cstarr
     Ptoto=Ptotstar
     etoto=etotstarr
     vdotBo=vdotBstarr
  else
     ro=rr
     uo=ur
     vo=vr
     wo=wr
     Bo=Br
     Co=Cr
     Ptoto=Ptotr
     etoto=etotr
     vdotBo=vdotBr
  end if

  ! Compute the Godunov flux
  fgdnv(1) = ro*uo
  fgdnv(2) = (etoto+Ptoto)*uo-A*vdotBo
  fgdnv(3) = ro*uo*uo+Ptoto-A*A
  fgdnv(4) = zero
  fgdnv(5) = ro*uo*vo-A*Bo
  fgdnv(6) = Bo*uo-A*vo
  fgdnv(7) = ro*uo*wo-A*Co
  fgdnv(8) = Co*uo-A*wo
#if NVAR > 8
  do ivar = 9,nvar
     if(fgdnv(1)>0)then
        fgdnv(ivar) = fgdnv(1)*qleft (ivar)
     else
        fgdnv(ivar) = fgdnv(1)*qright(ivar)
     endif
  end do
#endif

END SUBROUTINE hlld
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_mhd_flux(qvar,cvar,ff) 
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! compute the 1D MHD fluxes from the conservative variables
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1, Btransverse1, Vtransverse2, Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i , ivar, ib
  REAL(dp),DIMENSION(1:nvar):: qvar, cvar, ff
  REAL(dp) :: ecin,emag,etot,d,u,v,w,A,B,C,P,Ptot,entho

  ! Local variables
  entho = one/(gamma-one)
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  ecin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
  etot = P*entho+ecin+emag
  Ptot = P + emag
  
  ! Compute conservative variables
  cvar(1) = d
  cvar(2) = etot
  cvar(3) = d*u
  cvar(4) = A
  cvar(5) = d*v
  cvar(6) = B
  cvar(7) = d*w
  cvar(8) = C
#if NVAR > 8
  do ivar = 9,nvar
     cvar(ivar) = d*qvar(ivar)
  end do
#endif

  ! Compute fluxes
  ff(1) = d*u
  ff(2) = (etot+Ptot)*u-A*(A*u+B*v+C*w)
  ff(3) = d*u*u+Ptot-A*A
  ff(4) = zero
  ff(5) = d*u*v-A*B
  ff(6) = B*u-A*v
  ff(7) = d*u*w-A*C
  ff(8) = C*u-A*w
#if NVAR > 8
  do ivar = 9,nvar
     ff(ivar) = d*u*qvar(ivar)
  end do
#endif

END SUBROUTINE find_mhd_flux 
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_info(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fastest velocity at which information is exchanged 
  !! at the interface
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf+abs(u)

END SUBROUTINE find_speed_info
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_fast(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf

END SUBROUTINE find_speed_fast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_alfven(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the alfven velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  d=qvar(1); A=qvar(4)
  vel_info = sqrt(A*A/d)

END SUBROUTINE find_speed_alfven
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE athena_roe(qleft,qright,fmean,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  IMPLICIT NONE

  REAL(dp),DIMENSION(1:nvar)::qleft,qright
  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean,qmean,uleft,uright,vleft,vright,udiff

  REAL(dp), DIMENSION(7,7) :: lem, rem
  REAL(dp), DIMENSION(7)        :: lambda, lambdal, lambdar, a

  REAL(dp) :: droe, vxroe, vyroe, vzroe, sqrtdl, sqrtdr
  REAL(dp) :: hroe, byroe, bzroe, Xfactor, Yfactor
  REAL(dp) :: pbr, pbl

  REAL(dp) :: fluxd, fluxmx, fluxmy, fluxmz, fluxe
  REAL(dp) :: fluxby, fluxbz, byl, byr, bzl, bzr, bx
  REAL(dp) :: dl, vxl, vyl, vzl, pl, hl
  REAL(dp) :: dr, vxr, vyr, vzr, pr, hr
  REAL(dp) :: mxl, myl, mzl, el
  REAL(dp) :: mxr, myr, mzr, er
  REAL(dp) :: coef, bx_mean

  REAL(dp) :: zero_flux,dim,pim,eim,mxm,mym,mzm,bym,bzm,etm,l1,l2

  INTEGER :: i, n, m
  LOGICAL :: llf

  ! Enforce continuity of normal component
  bx_mean=0.5d0*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  ! compute the fluxes and the conserved variables
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright)

  ! define the primitive quantities explicitly 
  dl  = qleft (1)
  dr  = qright(1)

  pl  = qleft (2)
  pr  = qright(2)

  vxl = qleft (3)
  vxr = qright(3)

  vyl = qleft (5)
  vyr = qright(5)

  byl = qleft (6)
  byr = qright(6)

  vzl = qleft (7)
  vzr = qright(7)

  bzl = qleft (8)
  bzr = qright(8)

  ! finally attribute bx value (left and right are identical)
  bx  = 0.5d0*(qleft(4)+qright(4))

  ! define explicitly the conserved quantities
  el  = uleft (2)
  er  = uright(2)

  mxl = uleft (3)
  mxr = uright(3)

  myl = uleft (5)
  myr = uright(5)

  mzl = uleft (7)
  mzr = uright(7)

  ! the magnetic pressure
  pbl = half * (bx*bx + byl*byl + bzl*bzl)
  pbr = half * (bx*bx + byr*byr + bzr*bzr)
  
  ! the total (specific) enthalpy
  hl = (el+pl+pbl)/dl
  hr = (er+pr+pbr)/dr
  !
  ! Step 2 : Compute Roe-averaged data from left and right states
  !
  sqrtdl  = sqrt(dl)
  sqrtdr  = sqrt(dr)
  droe    = sqrtdl*sqrtdr
  vxroe   = (sqrtdl*vxl + sqrtdr*vxr)/(sqrtdl+sqrtdr)
  vyroe   = (sqrtdl*vyl + sqrtdr*vyr)/(sqrtdl+sqrtdr)
  vzroe   = (sqrtdl*vzl + sqrtdr*vzr)/(sqrtdl+sqrtdr)
  byroe   = (sqrtdr*byl + sqrtdl*byr)/(sqrtdl+sqrtdr)
  bzroe   = (sqrtdr*bzl + sqrtdl*bzr)/(sqrtdl+sqrtdr)
  hroe    = (sqrtdl*hl  + sqrtdr*hr )/(sqrtdl+sqrtdr)
  Xfactor = ((byroe*byroe-byl*byr)+(bzroe*bzroe-bzl*bzr))/(2*droe)
  Yfactor = (dl+dr)/(2*droe)
  !
  ! Step 3 : Compute eigenvalues and eigenmatrices from Roe-averaged values
  !
  call eigen_cons(droe,vxroe,vyroe,vzroe,hroe,bx,byroe,bzroe,Xfactor,Yfactor,lambda,rem,lem)
  !
  ! Step 4: Compute eigenvalues from left and right states
  !  
  call eigenvalues(dl,vxl,vyl,vzl,pl,bx,byl,bzl,lambdal)
  call eigenvalues(dr,vxr,vyr,vzr,pr,bx,byr,bzr,lambdar)
  !
  ! Step 5 : Create intermediate states from eigenmatrices
  !
  do n = 1, 7
     a(n)  = 0.0
     a(n)  = a(n)  + (dr -dl ) * lem(1,n)
     a(n)  = a(n)  + (mxr-mxl) * lem(2,n)
     a(n)  = a(n)  + (myr-myl) * lem(3,n)
     a(n)  = a(n)  + (mzr-mzl) * lem(4,n)
     a(n)  = a(n)  + (er - el) * lem(5,n)
     a(n)  = a(n)  + (byr-byl) * lem(6,n)
     a(n)  = a(n)  + (bzr-bzl) * lem(7,n)
  end do

  llf = .false.
  dim = dl 
  mxm = mxl
  mym = myl
  mzm = mzl
  eim = el 
  bym = byl
  bzm = bzl
  do n = 1, 7
     dim = dim + a(n) * rem(n,1)
     mxm = mxm + a(n) * rem(n,2)
     mym = mym + a(n) * rem(n,3)
     mzm = mzm + a(n) * rem(n,4)
     eim = eim + a(n) * rem(n,5)
     bym = bym + a(n) * rem(n,6)
     bzm = bzm + a(n) * rem(n,7)
     etm = eim-0.5*(mxm*mxm+mym*mym+mzm*mzm)/dim-0.5*(bx*bx+bym*bym+bzm*bzm)
     if(dim.le.zero.or.etm.le.zero)then
        llf=.true.
     endif
  end do

  IF( llf ) THEN
     fmean = half * ( fright + fleft ) * zero_flux
     CALL find_speed_info(qleft ,vleft )
     CALL find_speed_info(qright,vright)
     udiff = half * ( uright - uleft )
     fmean = fmean - MAX(vleft,vright) * udiff
     RETURN
  END IF
  !
  ! Step 6 : Entropy fix for genuinely non linear waves
  !
  do n = 1, 7, 2
     l1 = min(lambdal(n),lambda(n))
     l2 = max(lambdar(n),lambda(n))
     if(l1.lt.zero.and.l2.gt.zero)then
        lambda(n)=(lambda(n)*(l2+l1)-two*l2*l1)/(l2-l1)
     endif
  end do
  !
  ! Step 6 : Compute fluxes at interface using  Roe  solver
  !
  ! add the left and right fluxes
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  fleft  = fleft  * zero_flux
  fright = fright * zero_flux
  fluxd  = fleft(1) + fright(1)
  fluxe  = fleft(2) + fright(2)
  fluxmx = fleft(3) + fright(3)
  fluxmy = fleft(5) + fright(5)
  fluxby = fleft(6) + fright(6)
  fluxmz = fleft(7) + fright(7)
  fluxbz = fleft(8) + fright(8)

  ! now compute the Roe fluxes
  ! remember that the convention of athena's eigenvalues is : 
  ! rho,rhovx,rhovy,rhovz,E,by,bz
  DO n=1,7
      coef = ABS(lambda(n))*a(n)
      fluxd   = fluxd  - coef*rem(n,1)
      fluxe   = fluxe  - coef*rem(n,5)
      fluxmx  = fluxmx - coef*rem(n,2)
      fluxmy  = fluxmy - coef*rem(n,3)
      fluxby  = fluxby - coef*rem(n,6)
      fluxmz  = fluxmz - coef*rem(n,4)
      fluxbz  = fluxbz - coef*rem(n,7)
  ENDDO
  ! take half and put into the fmean variables 
  fmean(1) = half * fluxd 
  fmean(2) = half * fluxe 
  fmean(3) = half * fluxmx
  fmean(4) = zero
  fmean(5) = half * fluxmy
  fmean(6) = half * fluxby
  fmean(7) = half * fluxmz
  fmean(8) = half * fluxbz
#if NVAR > 8
  DO n=9,nvar
     if(fmean(1)>0)then
        fmean(n)=qleft (n)*fmean(1)
     else
        fmean(n)=qright(n)*fmean(1)
     endif
  END DO
#endif

END SUBROUTINE athena_roe
!###########################################################
!###########################################################
!###########################################################
!########################################################### 
SUBROUTINE hydro_acoustic(qleft,qright,fgdnv)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  real(dp),dimension(1:nvar)::qleft,qright,fgdnv,qgdnv,ugdnv

  ! local variables
  integer::i,n
  real(dp)::smallp, bx_mean
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   
  real(dp)::ro   ,uo   ,po   ,co   
  real(dp)::rstar,ustar,pstar,cstar
  real(dp)::wl   ,wr   ,wo   
  real(dp)::sgnm ,spin ,spout,ushock
  real(dp)::frac

  ! constants
  smallp = smallr*smallc**2

  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  ! Initial states pressure, density and velocity
  rl=max(qleft (1),smallr)
  rr=max(qright(1),smallr)

  pl=max(qleft (2),smallp)
  pr=max(qright(2),smallp)

  ul = qleft (3)
  ur = qright(3)

  ! Acoustic star state
  cl = sqrt(gamma*pl/rl)
  cr = sqrt(gamma*pr/rr)
  wl = cl*rl
  wr = cr*rr
  pstar = ((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr)
  ustar = ((wr*ur+wl*ul)+(pl-pr))/(wl+wr)

  ! Left going or right going contact wave
  sgnm = sign(one,ustar)

  ! Left or right unperturbed state
  if(sgnm==one)then
     ro = rl
     uo = ul
     po = pl
     wo = wl
     co = cl
  else
     ro = rr
     uo = ur
     po = pr
     wo = wr
     co = cr
  end if

  ! Star region density and sound speed
  rstar = ro+(pstar-po)/co**2
  rstar = max(rstar,smallr)
  cstar = sqrt(abs(gamma*pstar/rstar))
  cstar = max(cstar,smallc)

  ! Head and tail speed of rarefaction
  spout = co   -sgnm*uo   
  spin  = cstar-sgnm*ustar

  ! Shock speed
  ushock = half*(spin+spout)
  ushock = max(ushock,-sgnm*ustar)
  if(pstar>=po)then
     spout=ushock
     spin =spout 
  end if

  ! Sample the solution at x/t=0
  if(spout<zero)then      ! Initial state
     qgdnv(1) = ro
     qgdnv(2) = po
     qgdnv(3) = uo
  else if(spin>=zero)then  ! Star region
     qgdnv(1) = rstar
     qgdnv(2) = pstar
     qgdnv(3) = ustar
  else                        ! Rarefaction
     frac = spout/(spout-spin)
     qgdnv(1) = frac*rstar + (one - frac)*ro
     qgdnv(2) = frac*pstar + (one - frac)*po
     qgdnv(3) = frac*ustar + (one - frac)*uo
  end if

  ! Passive scalars
  do n = 4,nvar
     if(sgnm==one)then
        qgdnv(n) = qleft (n)
     else
        qgdnv(n) = qright(n)
     end if
  end do

  CALL find_mhd_flux(qgdnv,ugdnv,fgdnv)

end subroutine hydro_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eigenvalues(d,vx,vy,vz,p,bx,by,bz,lambda)
!
! MHD adiabatic eigenvalues
!
! Input Arguments:
!   Bx    = magnetic field in sweep direction
!   d     = density
!   vx    = X velocity 
!   vy    = Y velocity 
!   vz    = Z velocity 
!   by    = Y magnetic field
!   bz    = Z magnetic field
!   p     = thermal pressure
!
! Output Arguments:
!
!   lambda  = eigenvalues
! 
  USE amr_parameters
  USE const
  USE hydro_parameters
  IMPLICIT NONE

  real (dp), intent(IN) :: d, vx, vy, vz, p
  real (dp), intent(IN) :: bx, by, bz
  real (dp), dimension(1:7), intent(OUT) :: lambda
  ! local variables
  real (dp) :: btsq, bt, vsq, vax,  vaxsq, hp
  real (dp) :: asq, astarsq, cfsq, cfast, cssq, cslow

  vsq = vx**2+vy**2+vz**2
  btsq = by**2+bz**2
  bt = sqrt(btsq)
  vaxsq = Bx**2/d
  vax = sqrt(vaxsq)
  asq = gamma*p/d
  asq = MAX(asq,smallc**2)
  astarsq = asq+vaxsq+btsq/d
  
  cfsq = .5*(astarsq + sqrt(astarsq**2-4.0*asq*vaxsq))
  cfast = sqrt(cfsq)
  
  cssq = .5*(astarsq - sqrt(astarsq**2-4.0*asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)
  
  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast

end subroutine eigenvalues
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE eigen_cons(d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac,lambda,rem,lem)
!
! Input Arguments:
!   d     = Roe density
!   vx    = Roe X velocity 
!   vy    = Roe Y velocity 
!   vz    = Roe Z velocity 
!   Bx    = magnetic field in sweep direction
!   by    = Roe Y magnetic field
!   bz    = Roe Z magnetic field
!   h     = Roe enthalpy
!   Xfac  = ((by^2-byl*byr)+bz^2-bzl*bzr))/(2*rho)
!   Yfac  = (rho_l+rho_r)/(2*rho)
!
! Output Arguments:
!
!   lambda = eigenvalues
!   rem     = right eigenmatrix
!   lem     = left  eigenmatrix
!
  USE amr_parameters
  USE const
  USE hydro_parameters
  IMPLICIT NONE
  
  REAL (dp) :: d, vx, vy, vz, h
  REAL (dp) :: Bx, by, bz, Xfac, Yfac
  
  REAL(dp), DIMENSION(7,7) :: lem, rem
  REAL(dp), DIMENSION(7)        :: lambda
  
  REAL(dp) :: btsq, bt_starsq, bt, bt_star
  REAL(dp) :: vsq, vax,  vaxsq, hp, twid_asq, q_starsq 
  REAL(dp) :: cfsq, cfast, cssq, cslow
  REAL(dp) :: beta_y, beta_z, beta_ystar, beta_zstar, beta_starsq, vbeta
  REAL(dp) :: alpha_f, alpha_s, droot, s, twid_a
  REAL(dp) :: Qfast, Qslow, af_prime, as_prime, Afpbb, Aspbb, na
  REAL(dp) :: cff, css, af, as, Afpb, Aspb, vqstr, norm
  REAL(dp) :: Q_ystar, Q_zstar
  
  REAL(dp) :: smalle
  
  real(dp),dimension(7)::toto
  integer :: i
  
  vsq = vx*vx+vy*vy+vz*vz
  btsq = by*by+bz*bz
  bt_starsq = (gamma-1. - (gamma-2.)*Yfac)*btsq
  bt=sqrt(btsq)
  bt_star = sqrt(bt_starsq)
  
  vaxsq = Bx*Bx/d
  vax = sqrt(vaxsq)
  
  hp = h - (vaxsq + btsq/d)
  twid_asq = ((gamma-1.)*(hp-.5*vsq)-(gamma-2.)*Xfac)
  twid_asq = MAX(twid_asq,smallc*smallc)
  q_starsq = twid_asq+(vaxsq+bt_starsq/d)
  
  cfsq = .5*(q_starsq + sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  cfast = sqrt(cfsq)
  
  cssq = .5*(q_starsq - sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)
  
  if (bt .eq. 0) then
     beta_y = .5*sqrt(2.)
     beta_z = .5*sqrt(2.)
     beta_ystar = .5*sqrt(2.)
     beta_zstar = .5*sqrt(2.)
  else
     beta_y = by/bt
     beta_z = bz/bt
     beta_ystar = by/bt_star
     beta_zstar = bz/bt_star
  endif
  beta_starsq = beta_ystar*beta_ystar + beta_zstar*beta_zstar
  vbeta = vy*beta_ystar + vz*beta_zstar
  
  if ( (cfsq - cssq) .eq. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else if ( (twid_asq - cssq) .le. 0.) then
     alpha_f = 0.0
     alpha_s = 1.0
  else if ( (cfsq - twid_asq) .le. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else
     alpha_f = sqrt((twid_asq-cssq)/(cfsq-cssq))
     alpha_s = sqrt((cfsq-twid_asq)/(cfsq-cssq))
  endif
  !
  ! compute Qs and As for eigenmatrices
  !
  droot = sqrt(d)
  s  = SIGN(one,bx)
  twid_a = sqrt(twid_asq)
  Qfast = s*cfast*alpha_f
  Qslow = s*cslow*alpha_s
  af_prime = twid_a*alpha_f/droot
  as_prime = twid_a*alpha_s/droot
  Afpbb = af_prime*bt_star*beta_starsq
  Aspbb = as_prime*bt_star*beta_starsq
  !
  ! eigenvalues
  !
  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast
  !
  ! eigenmatrix
  !
  rem(1,1) = alpha_f
  rem(1,2) = alpha_f*(vx-cfast)
  rem(1,3) = alpha_f*vy + Qslow*beta_ystar
  rem(1,4) = alpha_f*vz + Qslow*beta_zstar
  rem(1,5) = alpha_f*(hp-vx*cfast) + Qslow*vbeta + Aspbb
  rem(1,6) = as_prime*beta_ystar
  rem(1,7) = as_prime*beta_zstar

  rem(2,1) = 0.
  rem(2,2) = 0.
  rem(2,3) = -beta_z
  rem(2,4) =  beta_y
  rem(2,5) = -(vy*beta_z - vz*beta_y)
  rem(2,6) = -s*beta_z/droot
  rem(2,7) =  s*beta_y/droot

  rem(3,1) = alpha_s
  rem(3,2) = alpha_s*(vx-cslow)
  rem(3,3) = alpha_s*vy - Qfast*beta_ystar
  rem(3,4) = alpha_s*vz - Qfast*beta_zstar
  rem(3,5) = alpha_s*(hp - vx*cslow) - Qfast*vbeta - Afpbb
  rem(3,6) = -af_prime*beta_ystar
  rem(3,7) = -af_prime*beta_zstar

  rem(4,1) = 1.0
  rem(4,2) = vx
  rem(4,3) = vy
  rem(4,4) = vz
  rem(4,5) = 0.5*vsq + (gamma-2.)*Xfac/(gamma-1.)
  rem(4,6) = 0.
  rem(4,7) = 0.

  rem(5,1) =  alpha_s
  rem(5,2) =  alpha_s*(vx+cslow)
  rem(5,3) =  alpha_s*vy + Qfast*beta_ystar
  rem(5,4) =  alpha_s*vz + Qfast*beta_zstar
  rem(5,5) =  alpha_s*(hp+vx*cslow) + Qfast*vbeta - Afpbb
  rem(5,6) =  rem(3,6)
  rem(5,7) =  rem(3,7)

  rem(6,1) = 0.
  rem(6,2) = 0.
  rem(6,3) =  beta_z
  rem(6,4) = -beta_y
  rem(6,5) = -rem(2,5)
  rem(6,6) =  rem(2,6)
  rem(6,7) =  rem(2,7)

  rem(7,1) =  alpha_f
  rem(7,2) = alpha_f*(vx+cfast)
  rem(7,3) = alpha_f*vy - Qslow*beta_ystar
  rem(7,4) = alpha_f*vz - Qslow*beta_zstar
  rem(7,5) = alpha_f*(hp + vx*cfast) - Qslow*vbeta + Aspbb
  rem(7,6) =  rem(1,6)
  rem(7,7) =  rem(1,7)
  ! 
  ! Left eignematrix
  !
  ! normalize some of the quantities by 1/(2a^2)
  ! some by (gamma-1)/2a^2
  !
  na  = 0.5/twid_asq
  cff = na*alpha_f*cfast
  css = na*alpha_s*cslow
  Qfast = Qfast*na
  Qslow = Qslow*na
  af = na*af_prime*d
  as = na*as_prime*d
  Afpb = na*af_prime*bt_star
  Aspb = na*as_prime*bt_star
  !
  alpha_f = (gamma-1.)*na*alpha_f
  alpha_s = (gamma-1.)*na*alpha_s
  Q_ystar = beta_ystar/beta_starsq
  Q_zstar = beta_zstar/beta_starsq
  vqstr   = (vy*Q_ystar+vz*Q_zstar)
  norm = (gamma-1.)*2.*na

  lem(1,1) = alpha_f*(vsq-hp) + cff*(cfast+vx) - Qslow*vqstr - Aspb
  lem(2,1) = -alpha_f*vx - cff
  lem(3,1) = -alpha_f*vy + Qslow*Q_ystar
  lem(4,1) = -alpha_f*vz + Qslow*Q_zstar
  lem(5,1) =  alpha_f
  lem(6,1) = as*Q_ystar - alpha_f*by
  lem(7,1) = as*Q_zstar - alpha_f*bz

  lem(1,2) =  0.5*(vy*beta_z - vz*beta_y)
  lem(2,2) =  0.
  lem(3,2) = -0.5*beta_z
  lem(4,2) =  0.5*beta_y
  lem(5,2) =  0.
  lem(6,2) = -0.5*droot*beta_z*s
  lem(7,2) =  0.5*droot*beta_y*s

  lem(1,3) =  alpha_s*(vsq-hp) + css*(cslow+vx) + Qfast*vqstr + Afpb
  lem(2,3) = -alpha_s*vx - css
  lem(3,3) = -alpha_s*vy - Qfast*Q_ystar
  lem(4,3) = -alpha_s*vz - Qfast*Q_zstar
  lem(5,3) =  alpha_s
  lem(6,3) = -af*Q_ystar - alpha_s*by
  lem(7,3) = -af*Q_zstar - alpha_s*bz

  ! Attention, il y a une difference de signe avec ramses_mhd...
  lem(1,4) =  1. - norm*(.5*vsq - (gamma-2.)*Xfac/(gamma-1.))
  ! This is the old version...
  !  lem(1,4) =  1. - norm*(.5*vsq + (gamma-2.)*Xfac/(gamma-1.))
  lem(2,4) =  norm*vx
  lem(3,4) =  norm*vy
  lem(4,4) =  norm*vz
  lem(5,4) = -norm
  lem(6,4) =  norm*by
  lem(7,4) =  norm*bz

  lem(1,5) =  alpha_s*(vsq-hp) + css*(cslow-vx) - Qfast*vqstr + Afpb
  lem(2,5) = -alpha_s*vx + css 
  lem(3,5) = -alpha_s*vy + Qfast*Q_ystar 
  lem(4,5) = -alpha_s*vz + Qfast*Q_zstar 
  lem(5,5) =  alpha_s
  lem(6,5) =  lem(6,3)
  lem(7,5) =  lem(7,3)

  lem(1,6) = -lem(1,2)
  lem(2,6) =  0.
  lem(3,6) = -lem(3,2)
  lem(4,6) = -lem(4,2)
  lem(5,6) =  0.
  lem(6,6) =  lem(6,2)
  lem(7,6) =  lem(7,2)

  lem(1,7) =  alpha_f*(vsq-hp) + cff*(cfast-vx) + Qslow*vqstr -Aspb
  lem(2,7) = -alpha_f*vx + cff
  lem(3,7) = -alpha_f*vy - Qslow*Q_ystar
  lem(4,7) = -alpha_f*vz - Qslow*Q_zstar
  lem(5,7) =  alpha_f
  lem(6,7) =  lem(6,1)
  lem(7,7) =  lem(7,1)

END SUBROUTINE eigen_cons
!###########################################################
!###########################################################
!###########################################################
!###########################################################
