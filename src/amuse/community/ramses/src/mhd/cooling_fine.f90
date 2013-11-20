subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef RT
  use rt_parameters, only: rt_UV_hom,rt_isDiffuseUVsrc
  use rt_cooling_module, only: update_UVrates
  use UV_module
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  endif
#ifdef RT
  if(neq_chem.and.ilevel==levelmin) then
     if(cosmo)call update_rt_c
     if(cosmo .and. rt_UV_hom)call update_UVrates
     if(cosmo .and. rt_isDiffuseUVsrc)call update_UVsrc
     if(ilevel==levelmin) call output_rt_stats
  endif
#endif

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroup
  use rt_hydro_commons
  use rt_cooling_module, only: n_U,iNpU,iFpU,rt_solve_cooling
#endif
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant,Fpnew,Npnew
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk
  real(kind=8),dimension(1:nvector),save::emag
#ifdef RT
  real(dp)::scale_Np,scale_Fp
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp),dimension(1:nvector,n_U),save::U,U_old
  real(dp),dimension(1:nvector,nGroups),save::Fp, Fp_precool
  real(dp),dimension(1:nvector,nGroups),save::dNpdt=0., dFpdt=0.
#endif
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf.eq.0)cycle

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),5)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        emag(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+5)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=exp(-nH(i)/0.01)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     if(cooling)then
        ! Compute thermal temperature by substracting polytrope
        do i=1,nleaf
           T2(i) = max(T2(i)-T2min(i),T2_min_fix)
        end do
     endif

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

#ifdef RT
     if(neq_chem) then
        ! Get gas thermal temperature
        do i=1,nleaf
           U(i,1) = T2(i)
        end do

        ! Get the ionization fractions
        do ivar=0,nIons-1
           do i=1,nleaf
              U(i,2+ivar) = uold(ind_leaf(i),iIons+ivar)/uold(ind_leaf(i),1)
           end do
        end do

        ! Get photon densities and flux magnitudes
        do ivar=1,nGroups
           do i=1,nleaf
              U(i,iNpU(ivar)) = scale_Np * rtuold(ind_leaf(i),iGroups(ivar))
              U(i,iFpU(ivar)) = scale_Fp &
                   * sqrt(sum((rtuold(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim))**2))
           enddo
           if(rt_smooth) then                           ! Smooth RT update
              do i=1,nleaf !Calc addition per sec to Np, Fp for current dt
                 Npnew = scale_Np * rtunew(ind_leaf(i),iGroups(ivar))
                 Fpnew = scale_Fp &
                      * sqrt(sum((rtunew(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim))**2))
                 dNpdt(i,ivar) = (Npnew - U(i,iNpU(ivar))) / dtcool
                 dFpdt(i,ivar) = (Fpnew - U(i,iFpU(ivar))) / dtcool ! Change in magnitude
                 ! Update flux vector to get the right direction
                 rtuold(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim) = &
                      rtunew(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim)
                 Fp_precool(i,ivar)=Fpnew           ! For update after solve_cooling
              end do
           else
              do i=1,nleaf
                 Fp_precool(i,ivar)=U(i,iFpU(ivar)) ! For update after solve_cooling
              end do
           end if
        end do

        if(cooling .and. delayed_cooling) then
           cooling_on(1:nleaf)=.true.
           do i=1,nleaf
              if(uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1) .gt. 1d-3) &
                   cooling_on(i)=.false.
           end do
        end if
        if(isothermal)cooling_on(1:nleaf)=.false.
     endif
#endif

     ! Compute net cooling at constant nH
     if(cooling.and..not.neq_chem)then
        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
     endif
#ifdef RT
     if(neq_chem) then
        U_old=U
        call rt_solve_cooling(U, dNpdt, dFpdt, nH, cooling_on, Zsolar, dtcool, aexp, nleaf)
        do i=1,nleaf
           delta_T2(i) = U(i,1) - T2(i)
        end do
     endif
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0) + ekk(i) + emag(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),5)
     end do
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           T2(i) = T2(i)+delta_T2(i)
        end do
     endif
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),5) = T2min(i)
        end do
     else
        do i=1,nleaf
           uold(ind_leaf(i),5) = max(T2(i),T2min(i))
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=20d0*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=uold(ind_leaf(i),idelay)*damp_factor
        end do
     endif

#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
        do ivar=0,nIons-1
           do i=1,nleaf
              uold(ind_leaf(i),iIons+ivar) = U(i,2+ivar)*nH(i)
           end do
        end do
     endif
     if(rt) then
        ! Update photon densities and flux magnitudes
        do ivar=1,nGroups
           do i=1,nleaf
              rtuold(ind_leaf(i),iGroups(ivar)) = U(i,iNpU(ivar)) /scale_Np
              if(Fp_precool(i,ivar) .gt. 0.d0)then
                 rtuold(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim) = U(i,iFpU(ivar))/Fp_precool(i,ivar) &
                      & *rtuold(ind_leaf(i),iGroups(ivar)+1:iGroups(ivar)+ndim)
              endif
           enddo
        end do
     endif
#endif

  end do
  ! End loop over cells

end subroutine coolfine1



