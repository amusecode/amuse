function initialize_code() result(ret)
  use amuse_sphrayMod
  integer :: ret
  call sphray_init
  ret=0
end function

function cleanup_code() result(ret)
  use amuse_sphrayMod
  integer :: ret
  call sphray_end
  ret=0
end function

function commit_particles() result(ret)
  use amuse_sphrayMod
  integer :: ret
  call sphray_commit_particles
  ret=0
end function

function recommit_particles() result(ret)
  use amuse_sphrayMod
  integer :: ret
  call sphray_commit_particles
  ret=0
end function


function commit_parameters() result(ret)
  use amuse_sphrayMod
  integer :: ret
  call sphray_commit_parameters
  ret=0
end function

function recommit_parameters() result(ret)
  use amuse_sphrayMod
  integer :: ret
  ret=-1
end function

function get_number_of_gas_particles(n) result(ret)
  use amuse_sphrayMod
  integer n,ret
  ret=-1
end function

function new_gas_particle(id,mass,hsml,x,y,z,rho,xe,u,vx,vy,vz) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: mass,hsml,x,y,z,rho,xe,u,vx,vy,vz
  real(r4b) :: t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u,t_vx,t_vy,t_vz
  t_mass=mass;t_hsml=hsml;t_x=x;t_y=y;t_z=z;t_rho=rho;t_xe=xe;t_u=u
  t_vx=vx;t_vy=vy;t_vz=vz
  call sphray_add_gas_particle(id,t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u,t_vx,t_vy,t_vz)
  ret=0
end function

function new_src_particle(id,L,x,y,z,SpcType) result(ret)
  use amuse_sphrayMod
  integer :: ret,id
  real(r8b) :: L,x,y,z,SpcType
  real(r4b) :: t_L,t_x,t_y,t_z,t_SpcType
  t_L=L;t_x=x;t_y=y;t_z=z;t_SpcType=SpcType
  call sphray_add_src_particle(id,t_L,t_x,t_y,t_z,t_SpcType)
  ret=0
end function

function get_state_gas(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  use amuse_sphrayMod
  integer :: ret,id
  real(r8b) :: mass,hsml,x,y,z,rho,xe,u
  real(r4b) :: t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u
  ret= sphray_get_gas_particle_state(id,t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u)
  mass=t_mass;hsml=t_hsml;x=t_x;y=t_y;z=t_z;rho=t_rho;xe=t_xe;u=t_u;
end function

function set_state_gas(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: mass,hsml,x,y,z,rho,xe,u
  real(r4b) :: t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u
  t_mass=mass;t_hsml=hsml;t_x=x;t_y=y;t_z=z;t_rho=rho;t_xe=xe;t_u=u
  ret=sphray_set_gas_particle_state(id,t_mass,t_hsml,t_x,t_y,t_z,t_rho,t_xe,t_u)
end function

function set_pos_gas(id,x,y,z) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: x,y,z
  real(r4b) ::t_x,t_y,t_z
  t_x=x;t_y=y;t_z=z
  ret=sphray_set_gas_particle_pos(id,t_x,t_y,t_z)
end function
function set_hsml_gas(id,hsml) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: hsml
  real(r4b) :: t_hsml
  t_hsml=hsml
  ret=sphray_set_gas_particle_hsml(id,t_hsml)
end function
function set_rho_gas(id,rho) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: rho
  real(r4b) :: t_rho
  t_rho=rho
  ret=sphray_set_gas_particle_rho(id,t_rho)
end function
function set_u_gas(id,u) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: u
  real(r4b) :: t_u
  t_u=u
  ret=sphray_set_gas_particle_u(id,t_u)
end function

function set_dudt_gas(id,u) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: u
  real(r4b) :: t_u
  t_u=u
  ret=sphray_set_gas_particle_dudt(id,t_u)
end function

function set_vel_gas(id,vx,vy,vz) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: vx,vy,vz
  real(r4b) :: t_vx,t_vy,t_vz
  t_vx=vx;t_vy=vy;t_vz=vz
  ret=sphray_set_gas_particle_vel(id,t_vx,t_vy,t_vz)
end function

function get_vel_gas(id,vx,vy,vz) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: vx,vy,vz
  real(r4b) :: t_vx,t_vy,t_vz
  ret=sphray_get_gas_particle_vel(id,t_vx,t_vy,t_vz)
  vx=t_vx;vy=t_vy;vz=t_vz
end function

function get_dudt_gas(id,u) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r8b) :: u
  real(r4b) :: t_u
  ret=sphray_get_gas_particle_dudt(id,t_u)
  u=t_u
end function

function get_state_src(id,L,x,y,z,spctype) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: id,ret
  real(r8b) :: L,x,y,z,spctype
  real(r4b) :: t_L,t_x,t_y,t_z,t_spctype
  ret=sphray_get_src_particle_state(id,t_L,t_x,t_y,t_z,t_spctype)
  L=t_L;x=t_x;y=t_y;z=t_z;spctype=t_spctype
end function

function set_state_src(id,L,x,y,z,spctype) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: id,ret
  real(r8b) :: L,x,y,z,spctype
  real(r4b) :: t_L,t_x,t_y,t_z,t_SpcType
  t_L=L;t_x=x;t_y=y;t_z=z;t_SpcType=SpcType
  ret=sphray_set_src_particle_state(id,t_L,t_x,t_y,t_z,t_spctype)
end function

function evolve_model(tend) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: tend
  call sphray_evolve(tend)
  ret=0
end function

function set_spectra_file(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_set_spectra_file(x)
  ret=0
end function

function get_spectra_file(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_get_spectra_file(x)
  ret=0
end function

function set_sphray_data_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_set_data_directory(x)
  ret=0
end function

function get_sphray_data_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_get_data_directory(x)
  ret=0
end function

function set_sphray_output_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_set_output_directory(x)
  ret=0
end function

function get_sphray_output_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_get_output_directory(x)
  ret=0
end function

function remove_gas_particle(id) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  ret=sphray_remove_gas_particle(id)
end function

function remove_src_particle(id) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  ret=sphray_remove_src_particle(id)
end function

function set_isothermal(flag) result(ret)
  use amuse_sphrayMod
  integer flag,ret
  if(flag.NE.0) call sphray_set_isothermal(.TRUE.)
  if(flag.EQ.0) call sphray_set_isothermal(.FALSE.)
  ret=0
end function

function get_isothermal(flag) result(ret)
  use amuse_sphrayMod
  integer :: ret,flag
  logical :: x
  flag=0
  call sphray_get_isothermal(x)
  if(x) flag=1
  ret=0
end function

function set_momentum_kicks(flag) result(ret)
  use amuse_sphrayMod
  integer flag,ret
  if(flag.NE.0) call sphray_set_momentum_kicks(.TRUE.)
  if(flag.EQ.0) call sphray_set_momentum_kicks(.FALSE.)
  ret=0
end function

function get_momentum_kicks(flag) result(ret)
  use amuse_sphrayMod
  integer :: ret,flag
  logical :: x
  flag=0
  call sphray_get_momentum_kicks(x)
  if(x) flag=1
  ret=0
end function


function set_H_caseA(flag) result(ret)
  use amuse_sphrayMod
  integer flag,ret
  if(flag.NE.0) call sphray_set_H_caseA(.TRUE.)
  if(flag.EQ.0) call sphray_set_H_caseA(.FALSE.)
  ret=0
end function

function get_H_caseA(flag) result(ret)
  use amuse_sphrayMod
  integer :: ret,flag
  logical :: x
  flag=0
  call sphray_get_H_caseA(x)
  if(x) flag=1
  ret=0
end function

function set_He_caseA(flag) result(ret)
  use amuse_sphrayMod
  integer flag,ret
  if(flag.NE.0) call sphray_set_He_caseA(.TRUE.)
  if(flag.EQ.0) call sphray_set_He_caseA(.FALSE.)
  ret=0
end function

function get_He_caseA(flag) result(ret)
  use amuse_sphrayMod
  integer :: ret,flag
  logical :: x
  flag=0
  call sphray_get_He_caseA(x)
  if(x) flag=1
  ret=0
end function

function get_time(time) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: time
  time=sphray_model_time()
  ret=0
end function

function set_time(time) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: time
  call sphray_set_model_time(time)
  ret=-1
end function


function get_raynumber(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x  
  call sphray_get_raynumber(x)
  ret=0
end function

function set_raynumber(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x  
  call sphray_set_raynumber(x)
  ret=0
end function

function get_defaultspectype(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  call sphray_get_defaultspectype(x)
  ret=0
end function

function set_defaultspectype(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  call sphray_set_defaultspectype(x)
  ret=0
end function


function get_iontempsolver(N) result(ret)
  use amuse_sphrayMod
  integer :: ret
  integer(i4b) :: N
  call sphray_get_iontempsolver(N)
  ret=0
end function

function set_iontempsolver(N) result(ret)
  use amuse_sphrayMod
  integer :: ret
  integer(i4b) :: N
  call sphray_set_iontempsolver(N)
  ret=0
end function

function get_boundary(N) result(ret)
  use amuse_sphrayMod
  integer :: ret
  integer(i4b) :: N
  call sphray_get_boundary(N)
  ret=0
end function

function set_boundary(N) result(ret)
  use amuse_sphrayMod
  integer :: ret
  integer(i4b) :: N
  call sphray_set_boundary(N)
  ret=0
end function

function get_boxsize(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  call sphray_get_boxsize(x)
  ret=0
end function

function set_boxsize(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  call sphray_set_boxsize(x)
  ret=0
end function

function set_globalHefraction(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  ret=sphray_set_he_mass_frac(x)
end function

function get_globalHefraction(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: x
  ret=sphray_get_he_mass_frac(x)
end function


