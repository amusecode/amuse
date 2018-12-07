module AmuseInterface
contains

function initialize_code() result(ret)
  use amuse_mercuryMod
  use StoppingConditions
  integer :: ret
  integer :: error
  error = set_support_for_condition(TIMEOUT_DETECTION)
  error = set_support_for_condition(COLLISION_DETECTION)
  ret=mercury_init()
end function  

function cleanup_code() result(ret)
  use amuse_mercuryMod
  integer :: ret
  ret=mercury_end()
end function  

function commit_particles() result(ret)
  use amuse_mercuryMod
  integer :: ret
  
  ret=finish_init()
end function  

function commit_parameters() result(ret)
  use amuse_mercuryMod
  integer :: ret
  ret=mercury_commit_parameters()
end function  

function recommit_parameters() result(ret)
  use amuse_mercuryMod
  integer :: ret
  ret=0
end function  

function recommit_particles() result(ret)
  use amuse_mercuryMod
  integer :: ret

  ret=finish_init()
end function  

function get_time(time_out) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: time_out
  ret=mercury_time(time_out)
end function  

function set_initial_timestep(time_) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: time_
  ret=set_initial_timestep_src(time_)
end function  

function get_initial_timestep(time_) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: time_
  ret=get_initial_timestep_src(time_)
end function  

function evolve_model(tend) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: tend
  ret=evolve_mercury(tend)
end function  

function synchronize_model() result(ret)
  integer :: ret
  ret = 0
end function

function new_orbiter(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit
  ret=add_particle(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)
end function  

function new_central_particle(id,mass,radius,j2,j4,j6,Lx,Ly,Lz) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass, radius, oblateness(3), spin(3)
  real*8 :: j2, j4, j6, Lx, Lz, Ly
  oblateness(1)=j2;oblateness(2)=j4;oblateness(3)=j6
  spin(1)=Lx;spin(2)=Ly;spin(3)=Lz
  ret=set_central_body(mass=mass, radius=radius, oblateness=oblateness,spin=spin)
end function  

function get_orbiter_state(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit
  ret=get_particle_state(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)
  if(id.EQ.1) ret=1
end function

function set_orbiter_state(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit
  ret=set_particle_state(id,mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)
  write(6,*) "in set orbiter state", id
  if(id.EQ.1) ret=1
end function

function set_central_particle_state(id,mass,radius,j2,j4,j6,Lx,Ly,Lz) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass, radius, oblateness(3), spin(3),j2,j4,j6,Lx,Ly,Lz
  oblateness(1)=j2;oblateness(2)=j4;oblateness(3)=j6
  spin(1)=Lx;spin(2)=Ly;spin(3)=Lz
  ret=set_central_body(mass=mass, radius=radius, oblateness=oblateness,spin=spin)
end function

function get_central_particle_state(id,mass,radius,j2,j4,j6,Lx,Ly,Lz) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: mass, radius, oblateness(3), spin(3),j2,j4,j6,Lx,Ly,Lz
  ret=get_central_body(mass=mass, radius=radius, oblateness=oblateness,spin=spin)
  j2=oblateness(1);j4=oblateness(2);j6=oblateness(3)
  Lx=spin(1);Ly=spin(2);Lz=spin(3)
end function

function get_position(id, x, y, z) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: x,y,z
  ret=get_particle_state(id, x=x, y=y, z=z)
end function

function set_position(id, x, y, z) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: x,y,z
  ret=set_particle_state(id,x=x,y=y,z=z)
end function

function get_velocity(id, vx, vy, vz) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: vx,vy,vz
  ret=get_particle_state(id, vx=vx, vy=vy, vz=vz)
end function

function set_velocity(id, vx, vy, vz) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: vx,vy,vz
  ret=set_particle_state(id, vx=vx, vy=vy, vz=vz)
end function

function delete_particle(id) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  ret=remove_particle(id)
end function  

function set_density(id, density) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: density
  ret=set_particle_state(id, dens=density)
end function

function get_density(id, density) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: density
  ret=get_particle_state(id, dens=density)
end function

function set_mass(id, mass) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: mass
  ret=set_particle_state(id, mass=mass)
end function

function get_mass(id, mass) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: mass
  ret=get_particle_state(id, mass=mass)
end function

function set_central_mass(id, mass) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: mass
  ret=set_central_body(mass=mass)
end function

function get_central_mass(id, mass) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: mass
  ret=get_central_body(mass=mass)
end function

function set_central_radius(id, radius) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: radius
  ret=set_central_body(radius=radius)
end function

function get_central_radius(id, radius) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: radius
  ret=get_central_body(radius=radius)
end function

function set_celimit(id, celimit) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: celimit
  ret=set_particle_state(id, celimit=celimit)
end function

function get_celimit(id, celimit) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: celimit
  ret=get_particle_state(id, celimit=celimit)
end function

function set_central_oblateness(id, j2,j4,j6) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: j2,j4,j6,oblateness(3)
  oblateness(1)=j2;oblateness(2)=j4;oblateness(3)=j6
  ret=set_central_body(oblateness=oblateness)
end function

function get_central_oblateness(id, j2,j4,j6) result(ret)
  use amuse_mercuryMod
  integer :: ret,id
  real*8 :: j2,j4,j6,oblateness(3)
  ret=get_central_body(oblateness=oblateness)
  j2=oblateness(1);j4=oblateness(2);j6=oblateness(3)
end function

function set_central_spin(id, lx,ly,lz) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: lx,ly,lz,spin(3)
  spin(1)=lx;spin(2)=ly;spin(3)=lz
  ret=set_central_body(spin=spin)
end function

function get_central_spin(id, lx,ly,lz) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: lx,ly,lz,spin(3)
  ret=get_central_body(spin=spin)
  lx=spin(1);ly=spin(2);lz=spin(3)
end function

function set_angularmomentum(id, Lx,Ly,Lz) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: Lx,Ly,Lz
  ret=set_particle_state(id, sx=Lx, sy=Ly, sz=Lz)
end function

function get_angularmomentum(id, Lx,Ly,Lz) result(ret)
  use amuse_mercuryMod
  integer :: ret, id
  real*8 :: Lx,Ly,Lz
  ret=get_particle_state(id, sx=Lx, sy=Ly, sz=Lz)
end function

function get_kinetic_energy(ek) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: ek
  call energy_angular_momentum(0,ek=ek)
  ret=0
end function

function get_potential_energy(ep) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: ep
  call energy_angular_momentum(0,ep=ep)
  ret=0
end function

function get_total_energy(e_tot) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: e_tot
  call total_energy_angular_momentum(e_tot=e_tot)
  ret=0
end function

function get_total_angular_momentum(am_tot) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: am_tot
  call total_energy_angular_momentum(am_tot=am_tot)
  ret=0
end function

function get_energy_deviation(delta_e) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: delta_e
  call energy_angular_momentum_deviation(delta_e=delta_e)
  ret=0
end function
  
function get_number_of_orbiters(norbiters) result(ret)
  use amuse_mercuryMod
  integer :: ret,norbiters
  ret=get_number_of_particles(norbiters)
end function  

function get_begin_time(system_time) result(ret)
      use amuse_mercuryMod
      implicit none
      integer :: ret
      real*8, intent(out) :: system_time

      ret = mercury_get_begin_time(system_time)
end function  

function set_begin_time(system_time) result(ret)
      use amuse_mercuryMod
      implicit none
      integer :: ret
      real*8, intent(in) :: system_time

      ret = mercury_set_begin_time(system_time)
end function  

function set_integrator(t_) result(ret)
  use amuse_mercuryMod
  integer :: ret,t_
  ret=set_algor(t_)
end function  

function get_integrator(t_) result(ret)
  use amuse_mercuryMod
  integer :: ret,t_
  ret=get_algor(t_)
end function

function set_rmax(r_max) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: r_max
  ret = mercury_set_rmax(r_max)
end function
function get_rmax(r_max) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: r_max
  ret = mercury_get_rmax(r_max)
end function

function set_cefac(cefac_n1) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: cefac_n1
  ret = mercury_set_cefac(cefac_n1)
end function
function get_cefac(cefac_n1) result(ret)
  use amuse_mercuryMod
  integer :: ret
  real*8 :: cefac_n1
  ret = mercury_get_cefac(cefac_n1)
end function

function set_elements_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f1=s)
end function
function get_elements_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f1=s)
end function

function set_close_encounters_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f2=s)
end function
function get_close_encounters_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f2=s)
end function

function set_info_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f3=s)
end function
function get_info_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f3=s)
end function

function set_bigbody_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f4=s)
end function
function get_bigbody_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f4=s)
end function

function set_smallbody_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f5=s)
end function
function get_smallbody_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f5=s)
end function

function set_integration_parameters_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f6=s)
end function
function get_integration_parameters_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f6=s)
end function

function set_restart_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=set_outputfiles(f7=s)
end function
function get_restart_file(s) result(ret)
  use amuse_mercuryMod
  integer :: ret
  character*4096 :: s
  ret=get_outputfiles(f7=s)
end function

! get acceleration at (heliocentric) positions x1,y1,z1 smoothed with eps1
function get_gravity_at_point(eps1, x1, y1, z1, ax,ay,az, number_of_points) result(ret)
      use amuse_mercuryMod
      implicit none
      integer :: ret
      integer, intent(in) :: number_of_points
      real*8, intent(in) :: eps1(number_of_points)
      real*8, intent(in) :: x1(number_of_points), y1(number_of_points)
      real*8, intent(in) :: z1(number_of_points)
      real*8, intent(out) :: ax(number_of_points),ay(number_of_points),az(number_of_points)
      ret = amuse_get_gravity_at_point(eps1,x1,y1,z1,ax,ay,az,number_of_points)
end function


! get potential at (heliocentric) positions x1,y1,z1 smoothed with eps1
function get_potential_at_point(eps1, x1, y1, z1, phi, number_of_points) result(ret)
      use amuse_mercuryMod
      implicit none
      integer :: ret
      integer, intent(in) :: number_of_points
      real*8, intent(in) :: eps1(number_of_points)
      real*8, intent(in) :: x1(number_of_points), y1(number_of_points)
      real*8, intent(in) :: z1(number_of_points)
      real*8, intent(out) :: phi(number_of_points)
      ret = amuse_get_potential_at_point(eps1,x1,y1,z1,phi,number_of_points)
end function

end module
