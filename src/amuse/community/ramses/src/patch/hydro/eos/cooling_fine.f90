subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
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

  if(cooling.and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################

subroutine read_eos_params(eos_type, nH_H_cc_threshold)
  use amr_commons
  implicit none
  character(LEN=80)::infile
  character(len=32)::eos_type
  real(dp)::nH_H_cc_threshold
  
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/eos_params/eos_type, nH_H_cc_threshold
  
  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=eos_params,END=105)
105 continue
  close(1)

  !--------------------------------------------------
  ! Check equation-of-state type
  !--------------------------------------------------
  select case (eos_type)
      case ('isothermal')
          if(myid==1) write(*,*) "Chosen EOS type :'isothermal'"
      case ('pseudo_cooling')
          if(myid==1) write(*,*) "Chosen EOS type :'pseudo_cooling'"
      case ('gamma_support')
          if(myid==1) write(*,*) "Chosen EOS type :'gamma_support'"
      case default
          if(myid==1) write(*,*) "Chosen EOS type :'isothermal'"
          eos_type='isothermal'
  end select

end subroutine read_eos_params

!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  logical, save::init_nml=.false.
  character(len=32), save::eos_type='isothermal'
  real(dp), save::nH_H_cc_threshold=10.0D0

  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,T2min,Zsolar

  real(kind=8)::dx,dx_loc,scale,alpha_dx2
  real(kind=8),dimension(1:3)::skip_loc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined EOS parameters in the namelist
  if (.not. init_nml) then
      call read_eos_params(eos_type, nH_H_cc_threshold)
      init_nml = .true.
  end if
  
  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Mesh maximum resolution
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(icoarse_max-icoarse_min+1)*scale_l
  dx=0.5D0**(nlevelmax)
  dx_loc=dx*scale


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

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),ndim+3)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     !do i=1,nleaf
     !   T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
     !   if(cooling)T2min(i)=T2min(i)+T2_min_fix
     !end do
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================
     do i=1,nleaf
        if(nH(i) .LT. 1.0D-3) then ! Low-density gamma=5/3 polytropic EOS
            T2min(i) = 4.0D6*(nH(i) / 1.0D-3)**(gamma - 1.0D0)
        else
            select case (eos_type)
                case ('isothermal') ! Isothermal EOS
                    T2min(i) = T2_star
                case ('pseudo_cooling') ! Pseudo-coolong EOS
                    if(nH(i) .LT. 10.0D0**(-0.5D0)) then ! Isothermal T2_star EOS
                        T2min(i) = T2_star
                    else ! Cooling EOS
                        T2min(i) = T2_star * (nH(i)/10.0D0**(-0.5D0))**(-1.0D0/2.0D0)
                    end if
                case ('gamma_support') ! Adiabatic gas
                    if(nH(i) .LT. nH_H_cc_threshold) then ! Isothermal T2_star EOS
                        T2min(i) = T2_star
                    else ! (Gamma-1) polytropic EOS
                        T2min(i) = T2_star * (nH(i)/nH_H_cc_threshold)**(g_star-1.0D0)
                    end if
                case default ! Isothermal EOS
                    T2min(i) = T2_star
            end select
        end if
        ! high-density gamma=2 polytropic EOS (Jeans criterion)
        ! Spherical collapse
        !   alpha_dx2 = 16.0D0 * (32.0D0/3.0D0) * (mH/kB) / (scale_t**2 * scale_d * dacos(-1.0D0) * g_star) * dx_loc**2
        ! Infinite 1D sinusoidal perturbartion collapse
        alpha_dx2 = 16 * (mH/kB) / (scale_t**2 * scale_d * dacos(-1.0D0) * gamma) * dx_loc**2
        T2min(i) = max(T2min(i), alpha_dx2 * (nH(i)*mH))
        if(cooling)T2min(i)=T2min(i)+T2_min_fix
     end do


     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     ! Compute net cooling at constant nH
     if(cooling)then
        do i=1,nleaf
           T2(i)=MAX(T2(i),T2min(i))
        end do
        call solve_cooling(nH,T2,Zsolar,dtcool,delta_T2,nleaf)
     endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink
     if(cooling)then
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
     endif

     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0) + ekk(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),ndim+2)
     end do
     if(cooling)then
        do i=1,nleaf
           T2(i) = T2(i)+delta_T2(i)
        end do
     endif
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i)
        end do
     else
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = max(T2(i),T2min(i))
        end do
     endif

  end do
  ! End loop over cells

end subroutine coolfine1



