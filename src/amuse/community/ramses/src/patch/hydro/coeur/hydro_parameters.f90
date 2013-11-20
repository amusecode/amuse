module hydro_parameters
  use amr_parameters

  ! Number of independant variables
  integer,parameter::nvar=ndim+2

  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:100,1:nvar)::boundary_var
  real(dp),dimension(1:100)::d_bound=0.0d0
  real(dp),dimension(1:100)::p_bound=0.0d0
  real(dp),dimension(1:100)::u_bound=0.0d0
  real(dp),dimension(1:100)::v_bound=0.0d0
  real(dp),dimension(1:100)::w_bound=0.0d0

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0                  ! Density gradient
  real(dp)::err_grad_u=-1.0                  ! Velocity gradient
  real(dp)::err_grad_p=-1.0                  ! Pressure gradient
  real(dp)::floor_d=1.d-10                   ! Density floor
  real(dp)::floor_u=1.d-10                   ! Velocity floor
  real(dp)::floor_p=1.d-10                   ! Pressure floor
  real(dp)::mass_sph=0.0D0                   ! mass_sph


  ! Initial conditions hydro variables
  real(dp),dimension(1:100)::d_region=0.
  real(dp),dimension(1:100)::u_region=0.
  real(dp),dimension(1:100)::v_region=0.
  real(dp),dimension(1:100)::w_region=0.
  real(dp),dimension(1:100)::p_region=0.

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::level_fix=0.08d0
  real(dp)::gamma=1.4d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  real(dp)::scaleL=1.0
  real(dp)::scaleT=1.0
  real(dp)::scaleM=1.0
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='acoustic'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  !parameter for the equation of state
  real(dp) :: f_rho_ad = 1.d0
  real(dp) :: rho_ad 
  real(dp) :: temp 

  !parameter for the initial conditions (Plummer sphere)
  real(dp) :: Mc       = 1.
  real(dp) :: fcut     = 0.8  
  real(dp) :: f0       = 0.2
  real(dp) :: alpha    = 0.5
  real(dp) :: beta_rot = 0.
  real(dp) :: omega
  real(dp) :: Pext

end module hydro_parameters
