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

function new_gas_particle(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: ret,id
  real(r4b) :: mass,hsml,x,y,z,rho,xe,u
  call sphray_add_gas_particle(id,mass,hsml,x,y,z,rho,xe,u)
  ret=0
end function

function new_src_particle(id,L,x,y,z,SpcType) result(ret)
  use amuse_sphrayMod
  integer :: ret,id
  real(r4b) :: L,x,y,z,SpcType
  call sphray_add_src_particle(id,L,x,y,z,SpcType)
  ret=0
end function

function get_state_gas(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  use amuse_sphrayMod
  integer :: ret,id
  real(r4b) :: mass,hsml,x,y,z,rho,xe,u
  ret= sphray_get_gas_particle_state(id,mass,hsml,x,y,z,rho,xe,u)
end function

function set_state_gas(id,mass,hsml,x,y,z,rho,xe,u) result(ret)
  use amuse_sphrayMod
  integer :: ret,id
  real(r4b) :: mass,hsml,x,y,z,rho,xe,u
  ret=sphray_set_gas_particle_state(id,mass,hsml,x,y,z,rho,xe,u)
end function

function get_state_src(id,L,x,y,z,spctype) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: id,ret
  real(r4b) :: L,x,y,z,spctype
  ret=sphray_get_src_particle_state(id,L,x,y,z,spctype)
end function

function set_state_src(id,L,x,y,z,spctype) result(ret)
  use amuse_sphrayMod
  integer(i4b) :: id,ret
  real(r4b) :: L,x,y,z,spctype
  ret=sphray_set_src_particle_state(id,L,x,y,z,spctype)
end function

function evolve_model(tend) result(ret)
  use amuse_sphrayMod
  integer :: ret
  real(r8b) :: tend
  call sphray_evolve(tend)
  ret=0
end function

function set_data_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_set_data_directory(x)
  ret=0
end function

function get_data_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_get_data_directory(x)
  ret=0
end function

function set_output_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_set_output_directory(x)
  ret=0
end function

function get_output_directory(x) result(ret)
  use amuse_sphrayMod
  integer :: ret
  character(len=512) :: x
  call sphray_get_output_directory(x)
  ret=0
end function
