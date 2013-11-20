subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  integer ::i,idim
  integer ,dimension(1:100)::bound_type
  real(dp)::scale,ek_bound
  real(dp)::mass_int,Etherm_int,Erot_int,ksi,pi,r0,rcut
  real(dp)::mass_plum,phi_plum,rot_plum,P_c,rho_c,E_grav,E_therm,E_rot
  integer ::np

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,nregion,region_type &
       & ,x_center,y_center,z_center &
       & ,length_x,length_y,length_z,exp_region &
       & ,d_region,u_region,v_region,w_region,p_region,Mc,alpha,beta_rot,f0,fcut
  namelist/hydro_params/gamma,courant_factor,smallr,smallc,c_star &
       & ,niter_riemann,slope_type,pressure_fix,cooling,scheme,riemann,f_rho_ad
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine,a_refine,b_refine,exp_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u,floor_d,floor_u,floor_p &
       & ,interpol_var,interpol_type
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max,kbound_min,kbound_max &
       & ,d_bound,u_bound,v_bound,w_bound,p_bound

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params)
  rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  if(simple_boundary)read(1,NML=boundary_params)



  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(c_star>0)star=.true.

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if

  ! Plummer sphere parameters
  ! assume r0 (pc), rcut(pc), Mc(in solar mass), alpha=Etherm/Egrav, beta_rot=Erot/Egrav are known
  pi = acos(-1.0d0)
  r0 = 0.5*boxlen*f0
  rcut = 0.5*boxlen*fcut
  mass_sph = Mc / (4.*pi/3.) * boxlen**3 / r0**3 / (twotondim)**levelmin 
  np  = 1000
  ksi = rcut / r0 / np
  ! calculate various integrals assuming a plummer sphere profile
  ! ie rho = rho_c / ( 1 + (r/r0)**2 )
  mass_plum = 0.
  phi_plum  = 0.
  rot_plum  = 0.
  DO i=1,np
     mass_plum = mass_plum + (i*ksi)**2*ksi / (1.d0 + (i*ksi)**2)
     phi_plum  = phi_plum  + (i*ksi)   *ksi / (1.d0 + (i*ksi)**2) * mass_plum
     ! assume solid body rotation
     rot_plum  = rot_plum  + (i*ksi)**4*ksi / (1.d0 + (i*ksi)**2)  
  ENDDO
  ! density in the center (in solar mass per pc^3)
  rho_c = Mc / (4.*pi) / r0**3 / mass_plum 
  ! from rho_c we can find the gravitational energy
  E_grav = (rho_c*4.*pi)**2*r0**5*phi_plum
  ! the ratio between Etherm and Egrav is an input parameter
  E_therm = alpha*E_grav
  ! the central pressure is given by Etherm
  P_c = E_therm / (4.d0*pi*r0**3) / mass_plum
  ! the external pressure 
  Pext = P_c / ( 1. + (rcut/r0)**2 ) 
  ! calculate the temperature. Will be used in barotrop.f90
  temp = P_c / rho_c
  ! calculate the density at which the gas becomes adiabatic (used in barotrop.f90)
  rho_ad = 1.d10 * 2.2 * 1.64d-24         !in g cm^-3
  rho_ad = rho_ad / 2.d33 * (3.08d18)**3  !normalise in code units (solar mass, pc)
  rho_ad = f_rho_ad * rho_ad              !f_rho_ad is readed in the .nml
  ! calculate the rotational energy, beta_rot is an input parameter
  E_rot = beta_rot*E_grav 
  ! calculate omega
  omega = sqrt(2.*E_rot/rho_c/r0**4/rot_plum)

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
  end do

end subroutine read_hydro_params

