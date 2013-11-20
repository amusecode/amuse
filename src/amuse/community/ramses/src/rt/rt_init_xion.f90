!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_init_xion(ilevel)

! Initialize hydrogen ionization state in all cells at given level from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer:: ilevel
  integer:: ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save:: ind_grid
!-------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Do the initialization by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call rt_init_xion_vsweep(ind_grid, ngrid)
  end do

111 format('   Entering rt_init_xion for level',i2)

END SUBROUTINE rt_init_xion

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_init_xion_vsweep(ind_grid, ngrid)

! Vector sweep initialization of hydrogen ionization state
! ind_grid => Indexes of grids/octs to initialize
! ngrid    => Number of vaid indexes in ind_grid (i.e. number of grids)
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters,only:nIons,iIons
  use cooling_module,only:Y
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::ind_grid
  integer::i, ind, iskip, idim, nleaf
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  integer,dimension(1:nvector),save::ind_cell, ind_leaf
  real(dp)::nH, T2, ekk, x, mu
  real(dp),dimension(3)::phI_rates       ! Photoionization rates [s-1]
  real(dp),dimension(6)::nSpec           !          Species abundances
!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pHI_rates(:)=0.0                   ! No UV background for the time being
  ! Loop over cells in each oct
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
     if(nleaf .eq. 0) cycle

     do i=1,nleaf
        ! Compute rho
        nH = MAX(uold(ind_leaf(i),1),smallr)   !       Mass density of gas
        ! Compute pressure from energy density
        T2 = uold(ind_leaf(i),ndim+2)          ! Energy density (kin+heat)
        ekk = 0.0d0                            !            Kinetic energy
        do idim=1,ndim
           ekk = ekk+0.5*uold(ind_leaf(i),idim+1)**2/nH
        end do
        T2 = (gamma-1.0)*(T2-ekk)              !     Gamma is ad. exponent
        ! now T2 is pressure (in user units)   !    (relates p and energy)
        ! Compute T2=T/mu in Kelvin from pressure:
        T2 = T2/nH*scale_T2                    !        Ideal gas equation
        ! Compute nH in H/cc (number of H nuclei per cubic centimeter)
        nH = nH*scale_nH

        call cmp_Equilibrium_Abundances(T2, nH, pHI_rates, mu, nSpec)

        ! update ionization states
        x = nSpec(3)/(nSpec(2)+nSpec(3))                  !   HII fraction
        uold(ind_leaf(i),iIons) = x*uold(ind_leaf(i),1)
        if(Y .gt. 0.d0 .and. nIons .ge. 3) then
           x = nSpec(5)/(nSpec(4)+nSpec(5)+nSpec(6))      !  HeII fraction
           uold(ind_leaf(i),iIons+1) = x*uold(ind_leaf(i),1)         
           x = nSpec(6)/(nSpec(4)+nSpec(5)+nSpec(6))      ! HeIII fraction
           uold(ind_leaf(i),iIons+2) = x*uold(ind_leaf(i),1)         
        endif
      end do

  end do
  ! End loop over cells
END SUBROUTINE rt_init_xion_vsweep

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE calc_equilibrium_xion(vars, rtvars, xion)

! Calculate and return photoionization equilibrium abundance states for 
! a cell
! vars     => Cell variables (rho, v, u, w, etc)
! rtvars   => Cell RT variables (Np1, Fpx1, Fpy1, etc)
! xion     => Equilibrium ionization states of cell
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  use cooling_module,only:Y
  use rt_cooling_module,only:UVrates,signc
  implicit none
  real(dp),dimension(nvar)::vars
  real(dp),dimension(nrtvar)::rtvars
  real(dp),dimension(nIons)::xion
  integer::ip, iI, idim
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(dp)::scale_Np, scale_Fp, nH, T2, ekk, mu
  real(dp),dimension(nIons)::phI_rates       ! Photoionization rates [s-1]
  real(dp),dimension(6)::nSpec               !          Species abundances
!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np,scale_Fp)

  ! Calculate photoionization rates:
  phI_rates(:)=0.0
  do ip=1, nGroups
     do iI=1,nIons
        phI_rates(iI) = phI_rates(iI) &
                      + rtvars(iGroups(ip))*scale_Np*signc(ip,iI)
     end do
  end do


  nH = MAX(vars(1),smallr)                  !   Number density of gas [UU]

  ! Compute pressure from energy density
  T2 = vars(ndim+2)                         ! Energy dens. (kin+heat) [UU]
  ekk = 0.0d0                               !          Kinetic energy [UU]
  do idim=1,ndim
     ekk = ekk+0.5*vars(idim+1)**2/nH
  end do
  T2 = (gamma-1.0)*(T2-ekk)                 !        Gamma is ad. exponent
                                            !      now T2 is pressure [UU]
  T2 = T2/nH*scale_T2                       !                T/mu [Kelvin]
  nH = nH*scale_nH                          !        Number density [H/cc]

  if(rt_UV_hom .and. nH .lt. rt_UV_nHSS) &  !   UV backgr. photoionization
       phI_rates = phI_rates + UVrates(:,1)

  call cmp_Equilibrium_Abundances(T2, nH, pHI_rates, mu, nSpec)
  xion(1)=nSpec(3)/(nSpec(2)+nSpec(3))                    !   HII fraction
  if(Y .gt. 0.d0 .and. nIons .ge. 3) then
     xion(2) = nSpec(5)/(nSpec(4)+nSpec(5)+nSpec(6))      !  HeII fraction
     xion(3) = nSpec(6)/(nSpec(4)+nSpec(5)+nSpec(6))      ! HeIII fraction
  endif

END SUBROUTINE calc_equilibrium_xion

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_Equilibrium_Abundances(T2, nH, phI_rates, mu, nSpec)
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  use rt_cooling_module
  use rt_parameters,only:nIons
  implicit none
  real(dp) ::T2,nH
  real(dp),dimension(nIons)::phI_rates
  real(dp) ::mu
  real(dp),dimension(1:6)::nSpec!-----------------------------------------
  real(dp) ::mu_old, err_mu, mu_left, mu_right, T, nTot
  integer :: niter
!-------------------------------------------------------------------------
  ! Iteration to find mu                     ! n_E     = n_spec(1) ! e
  err_mu=1.                                  ! n_HI    = n_spec(2) ! H
  mu_left=0.5                                ! n_HII   = n_spec(3) ! H+
  mu_right=1.3                               ! n_HEI   = n_spec(4) ! He
  niter=0                                    ! n_HEII  = n_spec(5) ! He+
  do while (err_mu > 1.d-4 .and. niter <= 50)! n_HEIII = n_spec(6) ! He++
     mu_old=0.5*(mu_left+mu_right)
     T = T2*mu_old
     call cmp_chem_eq(T, nH, phI_rates, nSpec, nTot, mu)
     err_mu = (mu-mu_old)/mu_old
     if(err_mu>0.)then 
        mu_left =0.5*(mu_left+mu_right)
        mu_right=mu_right
     else
        mu_left =mu_left
        mu_right=0.5*(mu_left+mu_right)
     end if
     err_mu=ABS(err_mu)
     niter=niter+1
  end do
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_Equilibrium_Abundances : too many iterations.'
     STOP
  endif
    
END SUBROUTINE cmp_Equilibrium_Abundances
