MODULE AmuseInterface
    use MuseHelpers
    
CONTAINS

function initialize_code() result(ret)
  use StoppingConditions
  integer :: ret
  integer :: error
  call muse_start
  call muse_reset(0.0)
  call muse_set_time(0.0)
  error = set_support_for_condition(TIMEOUT_DETECTION)
  error = set_support_for_condition(NUMBER_OF_STEPS_DETECTION)
  error = set_support_for_condition(OUT_OF_BOX_DETECTION)
  error = set_support_for_condition(DENSITY_LIMIT_DETECTION)
  error = set_support_for_condition(INTERNAL_ENERGY_LIMIT_DETECTION)
  ret=0
end function

function cleanup_code() result(ret)
  integer :: ret
  call muse_end
  call muse_start
  call muse_reset(0.0)
  ret=0
end function


function commit_particles() result(ret)
  integer :: ret
  call muse_finalize_init 
  ret=0
end function

function recommit_particles() result(ret)
  integer :: ret
  ret=muse_reinitialize()
end function


function commit_parameters() result(ret)
  integer :: ret
  call muse_init
  ret=0
end function

function recommit_parameters() result(ret)
  integer :: ret
  ret=0
end function

function get_number_of_particles(n) result(ret)
  integer n,ret
  n=muse_get_nbodies() 
  ret=0
end function

function get_index_of_first_particle(id) result(ret)
  integer :: id,ret
  ret=muse_index_of_first_particle(id)
end function

function get_index_of_next_particle(id,id1) result(ret)
  integer :: id,id1,ret
  ret=muse_index_of_next_particle(id,id1)
end function

function get_time(t) result(ret)
  integer :: ret
  real*8 :: t
  t=muse_get_time()
  ret=0
end function

function set_begin_time(t) result(ret)
  integer :: ret
  real*8 :: t
  call muse_set_begin_time(t)
  ret = 0
end function

function get_begin_time(t) result(ret)
  integer :: ret
  real*8 :: t
  t = muse_get_begin_time()
  ret = 0
end function


function evolve_model(tend) result(ret)
  integer :: ret
  real*8 :: tend
  call muse_stepsys(tend,1)
  ret=0
end function

function synchronize_model() result(ret)
  integer :: ret
  real*8 :: dum1,dum2,dum3 
  ret=amuse_synchronize_model()
end function

function get_kinetic_energy(e) result(ret)
  integer :: ret
  real*8 :: e,ek,ep,eth
  call muse_energies(0,ek,ep,eth)
  e=ek
  ret=0
end function

function get_potential_energy(e) result(ret)
  integer :: ret
  real*8 :: e,ek,ep,eth
  call muse_energies(0,ek,ep,eth)
  e=ep
  ret=0
end function

function get_thermal_energy(e) result(ret)
  integer :: ret
  real*8 :: e,ek,ep,eth
  call muse_energies(0,ek,ep,eth)
  e=eth
  ret=0
end function

function get_total_energy(e) result(ret)
  integer :: ret
  real*8 :: e,ek,ep,eth
  call muse_energies(0,ek,ep,eth)
  e=ek+ep+eth
  ret=0
end function


function new_dm_particle(ids,mass,x,y,z,vx,vy,vz,eps) result(ret)
  integer :: ids,ret,oldnp
  real*8 :: mass,eps,x,y,z,vx,vy,vz
  integer :: idsa(1)
  real*8 :: massa(1),epsa(1),xa(1),ya(1),za(1),vxa(1),vya(1),vza(1)
  ids=new_id()
  oldnp=muse_get_nbodies()
  idsa = ids
  massa = mass
  xa = x
  ya = y
  za = z
  vxa = vx
  vya = vy
  vza = vz
  epsa = eps
  ret=add_dm_particle(idsa,massa,xa,ya,za,vxa,vya,vza,epsa,1)
  if(ret.EQ.oldnp+1) then
    ret=0
  else
    ret=-1
  endif 
end function

function new_sph_particle(ids,mass,x,y,z,vx,vy,vz,u,eps) result(ret)
  integer :: ids,ret,oldnp
  real*8 :: mass,eps,x,y,z,vx,vy,vz,u
  integer :: idsa(1)
  real*8 :: massa(1),epsa(1),xa(1),ya(1),za(1),vxa(1),vya(1),vza(1),ua(1)
  ids=new_id()
  oldnp=muse_get_nsph()
  
  idsa = ids
  massa = mass
  xa = x
  ya = y
  za = z
  vxa = vx
  vya = vy
  vza = vz
  epsa = eps
  ua = u
  
  ret=add_sph_particle(idsa,massa,xa,ya,za,vxa,vya,vza,epsa,ua,1)
  if(ret.EQ.oldnp+1) then
    ret=0
  else
    ret=-1
  endif 
end function

function new_star_particle(ids,mass,x,y,z,vx,vy,vz,tf,eps) result(ret)
  integer :: ids,ret,oldnp
  real*8 :: mass,eps,x,y,z,vx,vy,vz,tf
  integer :: idsa(1)
  real*8 :: massa(1),epsa(1),xa(1),ya(1),za(1),vxa(1),vya(1),vza(1),tfa(1)
  ids=new_id()
  oldnp=muse_get_nstar()
  idsa = ids
  massa = mass
  xa = x
  ya = y
  za = z
  vxa = vx
  vya = vy
  vza = vz
  epsa = eps
  tfa = tf
  
  ret=add_star_particle(idsa,massa,xa,ya,za,vxa,vya,vza,epsa,tfa,1)
  if(ret.EQ.oldnp+1) then
    ret=0
  else
    ret=-1
  endif 
end function

function add_dm_particle(ids,mass,x,y,z,vx,vy,vz,eps,npart) result(n)
  integer :: npart
  integer :: ids(npart),n
  real*8 :: mass(npart),x(npart),y(npart),z(npart), &
    vx(npart),vy(npart),vz(npart),eps(npart)
  call muse_add_particle_dm(ids,mass,x,y,z,vx,vy,vz,eps,npart)
  n=muse_get_nbodies()
end function

function add_star_particle(ids,mass,x,y,z,vx,vy,vz,eps,tf,npart) result(n)
  integer :: npart
  integer :: ids(npart),n
  real*8 :: mass(npart),x(npart),y(npart),z(npart), &
    vx(npart),vy(npart),vz(npart),eps(npart),tf(npart)
  call muse_add_particle_star(ids,mass,x,y,z,vx,vy,vz,eps,tf,npart)
  n=muse_get_nstar()
end function

function add_sph_particle(ids,mass,x,y,z,vx,vy,vz,eps,u,npart) result(n)
  integer :: npart
  integer :: ids(npart),n
  real*8 :: mass(npart),x(npart),y(npart),z(npart), &
    vx(npart),vy(npart),vz(npart),eps(npart),u(npart)
  call muse_add_particle_sph(ids,mass,x,y,z,vx,vy,vz,eps,u,npart)
  n=muse_get_nsph()
end function

function set_state(id,mass,x,y,z,vx,vy,vz,eps) result(ret)
  integer id,ret
  real*8 mass,eps,x,y,z,vx,vy,vz 
  ret=amuse_set_state(id,mass,x,y,z,vx,vy,vz,eps)
end function

function set_state_sph(id,mass,x,y,z,vx,vy,vz,u,eps) result(ret)
  integer id,ret
  real*8 mass,eps,x,y,z,vx,vy,vz,u 
  ret=amuse_set_state_sph(id,mass,x,y,z,vx,vy,vz,eps,u)
end function

function set_state_star(id,mass,x,y,z,vx,vy,vz,tf,eps) result(ret)
  integer id,ret
  real*8 mass,eps,x,y,z,vx,vy,vz,tf 
  ret=amuse_set_state_star(id,mass,x,y,z,vx,vy,vz,eps,tf)
end function

function get_state(id,mass,x,y,z,vx,vy,vz,eps) result(ret)
  integer :: id,ret
  real*8 :: mass,x,y,z,vx,vy,vz,eps
  ret=amuse_get_state(id,mass,x,y,z,vx,vy,vz,eps)
end function

function get_state_sph(id,mass,x,y,z,vx,vy,vz,u,eps) result(ret)
  integer :: id,ret
  real*8 :: mass,x,y,z,vx,vy,vz,eps,u
  ret=amuse_get_state_sph(id,mass,x,y,z,vx,vy,vz,eps,u)
end function

function get_state_star(id,mass,x,y,z,vx,vy,vz,tf,eps) result(ret)
  integer :: id,ret
  real*8 :: mass,x,y,z,vx,vy,vz,eps,tf
  ret=amuse_get_state_star(id,mass,x,y,z,vx,vy,vz,eps,tf)
end function

function set_mass(id,mass) result(ret)
  integer id,ret
  real*8 mass 
  ret=amuse_set_mass(id,mass)
end function
function set_radius(id,r) result(ret)
  integer id,ret
  real*8 :: r
  ret=amuse_set_epsgrav(id,r)
end function
function set_smoothing_length(id,eps) result(ret)
  integer id,ret
  real*8 eps
  ret=amuse_set_hsmooth(id,eps)
end function
function set_position(id,x,y,z) result(ret)
  integer id,ret
  real*8 x,y,z 
  ret=amuse_set_position(id,x,y,z)
end function
function set_velocity(id,vx,vy,vz) result(ret)
  integer id,ret
  real*8 vx,vy,vz 
  ret=amuse_set_velocity(id,vx,vy,vz)
end function
function set_internal_energy(id,u) result(ret)
  integer id,ret
  real*8 u 
  ret=amuse_set_internal_energy(id,u)
end function
function set_star_tform(id,tf) result(ret)
  integer id,ret
  real*8 tf 
  ret=amuse_set_star_tform(id,tf)
end function

function get_mass(id,mass) result(ret)
  integer :: id,ret
  real*8 :: mass
  ret=amuse_get_mass(id,mass)
end function
function get_radius(id,eps) result(ret)
  integer id,ret
  real*8 :: eps
  ret=amuse_get_epsgrav(id,eps)
end function
function get_smoothing_length(id,eps) result(ret)
  integer :: id,ret
  real*8 :: eps
  ret=amuse_get_hsmooth(id,eps)
end function
function get_density(id,density) result(ret)
  integer :: id,ret
  real*8 :: density
  ret=amuse_get_density(id,density)
end function

function get_pressure(id,pressure) result(ret)
  integer :: id,ret
  real*8 :: pressure
  ret=amuse_get_pressure(id,pressure)
end function

function get_position(id,x,y,z) result(ret)
  integer :: id,ret
  real*8 :: x,y,z
  ret=amuse_get_position(id,x,y,z)
end function
function get_velocity(id,vx,vy,vz) result(ret)
  integer :: id,ret
  real*8 :: vx,vy,vz
  ret=amuse_get_velocity(id,vx,vy,vz)
end function
function get_internal_energy(id,u) result(ret)
  integer :: id,ret
  real*8 :: u
  ret=amuse_get_internal_energy(id,u)
end function
function get_dinternal_energy_dt(id,dudt) result(ret)
  integer :: id,ret
  real*8 :: dudt
  ret=amuse_get_dinternal_energy_dt(id,dudt)
end function
function get_star_tform(id,tf) result(ret)
  integer :: id,ret
  real*8 :: tf
  ret=amuse_get_star_tform(id,tf)
end function

function delete_particle(id) result(ret)
  integer id,ret
  ret=muse_remove_particle(id)
end function

function get_gravity_at_point(eps, x, y, z, ax, ay, az, n) result(ret)
  integer :: ret,n  
  real*8 :: eps(n), x(n), y(n), z(n), ax(n), ay(n), az(n)
  ax=0;ay=0;az=0
  call muse_get_gravity(eps,x,y,z,ax,ay,az,n)
  ret=0
end function

function get_potential_at_point(eps, x, y, z, phi, n) result(ret)
  integer :: ret,n
  real*8 :: eps(n),x(n), y(n), z(n), phi(n)
  call muse_get_pot(eps,x,y,z,phi,n)
  ret=0  
end function

function get_potential(id, phi) result(ret)
  integer:: ret, id
  real*8 :: phi
  ret = amuse_get_potential(id, phi)
end function

function get_hydro_state_at_point(x, y, z, vx, vy, vz, rho, rhovx, rhovy, rhovz, rhoe,n) result(ret)
  integer :: ret,n
  real*8 :: x(n), y(n), z(n), vx(n), vy(n), vz(n), rho(n), rhovx(n), &
    rhovy(n), rhovz(n), rhoe(n)
  rho=0.
  rhovx=0.
  rhovy=0.
  rhovz=0.
  rhoe=0.  
  call muse_get_hydro_state(x,y,z,vx,vy,vz,rho,rhovx,rhovy,rhovz,rhoe,n)
  ret=0  
end function


! setting/ getting parameters

function set_use_hydro(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_usesph(.TRUE.)
  if(flag.EQ.0) call amuse_set_usesph(.FALSE.)
  ret=0
end function
function get_use_hydro(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_usesph(x) 
  if(x) flag=1
  ret=0
end function

function set_radiate(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_radiate(.TRUE.)
  if(flag.EQ.0) call amuse_set_radiate(.FALSE.)
  ret=0
end function
function get_radiate(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_radiate(x) 
  if(x) flag=1
  ret=0
end function

function set_starform(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_starform(.TRUE.)
  if(flag.EQ.0) call amuse_set_starform(.FALSE.)
  ret=0
end function
function get_starform(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_starform(x) 
  if(x) flag=1
  ret=0
end function

function set_cosmo(zeroiftrue) result(ret)
  integer :: zeroiftrue,ret
  if(zeroiftrue.EQ.0) call amuse_set_cosmo(.TRUE.)
  if(zeroiftrue.NE.0) call amuse_set_cosmo(.FALSE.)
  ret=0
end function
function get_cosmo(zeroiftrue) result(ret)
  integer :: ret,zeroiftrue
  logical :: x
  zeroiftrue=1
  call amuse_get_cosmo(x) 
  if(x) zeroiftrue=0
  ret=0
end function

function set_sqrttstp(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_sqrttstp(.TRUE.)
  if(flag.EQ.0) call amuse_set_sqrttstp(.FALSE.)
  ret=0
end function
function get_sqrttstp(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_sqrttstp(x) 
  if(x) flag=1
  ret=0
end function

function set_acc_tstp(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_acc_tstp(.TRUE.)
  if(flag.EQ.0) call amuse_set_acc_tstp(.FALSE.)
  ret=0
end function
function get_acc_tstp(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_acc_tstp(x) 
  if(x) flag=1
  ret=0
end function

function set_freetstp(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_freetstp(.TRUE.)
  if(flag.EQ.0) call amuse_set_freetstp(.FALSE.)
  ret=0
end function
function get_freetstp(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_freetstp(x) 
  if(x) flag=1
  ret=0
end function

function set_usequad(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_usequad(.TRUE.)
  if(flag.EQ.0) call amuse_set_usequad(.FALSE.)
  ret=0
end function
function get_usequad(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_usequad(x) 
  if(x) flag=1
  ret=0
end function

function set_directsum(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_directsum(.TRUE.)
  if(flag.EQ.0) call amuse_set_directsum(.FALSE.)
  ret=0
end function
function get_directsum(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_directsum(x) 
  if(x) flag=1
  ret=0
end function

function set_selfgrav(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_selfgrav(.TRUE.)
  if(flag.EQ.0) call amuse_set_selfgrav(.FALSE.)
  ret=0
end function
function get_selfgrav(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_selfgrav(x) 
  if(x) flag=1
  ret=0
end function

function set_fixthalo(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_fixthalo(.TRUE.)
  if(flag.EQ.0) call amuse_set_fixthalo(.FALSE.)
  ret=0
end function
function get_fixthalo(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_fixthalo(x) 
  if(x) flag=1
  ret=0
end function

function set_adaptive_eps(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_adaptive_eps(.TRUE.)
  if(flag.EQ.0) call amuse_set_adaptive_eps(.FALSE.)
  ret=0
end function
function get_adaptive_eps(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_adaptive_eps(x) 
  if(x) flag=1
  ret=0
end function

function set_gdgop(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_gdgop(.TRUE.)
  if(flag.EQ.0) call amuse_set_gdgop(.FALSE.)
  ret=0
end function
function get_gdgop(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_gdgop(x) 
  if(x) flag=1
  ret=0
end function

function set_smoothinput(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_smoothinput(.TRUE.)
  if(flag.EQ.0) call amuse_set_smoothinput(.FALSE.)
  ret=0
end function
function get_smoothinput(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_smoothinput(x) 
  if(x) flag=1
  ret=0
end function

function set_consph(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_consph(.TRUE.)
  if(flag.EQ.0) call amuse_set_consph(.FALSE.)
  ret=0
end function
function get_consph(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_consph(x) 
  if(x) flag=1
  ret=0
end function

function set_sphinit(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_sphinit(.TRUE.)
  if(flag.EQ.0) call amuse_set_sphinit(.FALSE.)
  ret=0
end function
function get_sphinit(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_sphinit(x) 
  if(x) flag=1
  ret=0
end function

function set_uentropy(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_uentropy(.TRUE.)
  if(flag.EQ.0) call amuse_set_uentropy(.FALSE.)
  ret=0
end function
function get_uentropy(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_uentropy(x) 
  if(x) flag=1
  ret=0
end function

function set_isotherm(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_isotherm(.TRUE.)
  if(flag.EQ.0) call amuse_set_isotherm(.FALSE.)
  ret=0
end function
function get_isotherm(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_isotherm(x) 
  if(x) flag=1
  ret=0
end function

function set_eps_is_h(flag) result(ret)
  integer :: flag,ret
  if(flag.NE.0) call amuse_set_eps_is_h(.TRUE.)
  if(flag.EQ.0) call amuse_set_eps_is_h(.FALSE.)
  ret=0
end function
function get_eps_is_h(flag) result(ret)
  integer :: ret,flag
  logical :: x
  flag=0
  call amuse_get_eps_is_h(x) 
  if(x) flag=1
  ret=0
end function


! integers

function set_firstsnap(i) result(ret)
  integer :: ret,i
  call amuse_set_firstsnap(i)
  ret=0
end function
function get_firstsnap(i) result(ret)
  integer :: ret,i
  call amuse_get_firstsnap(i) 
  ret=0
end function

function set_stepout(i) result(ret)
  integer :: ret,i
  call amuse_set_stepout(i)
  ret=0
end function
function get_stepout(i) result(ret)
  integer :: ret,i
  call amuse_get_stepout(i) 
  ret=0
end function

function set_steplog(i) result(ret)
  integer :: ret,i
  call amuse_set_steplog(i)
  ret=0
end function
function get_steplog(i) result(ret)
  integer :: ret,i
  call amuse_get_steplog(i) 
  ret=0
end function

function set_max_tbin(i) result(ret)
  integer :: ret,i
  call amuse_set_max_tbin(i)
  ret=0
end function
function get_max_tbin(i) result(ret)
  integer :: ret,i
  call amuse_get_max_tbin(i) 
  ret=0
end function

function set_minppbin(i) result(ret)
  integer :: ret,i
  call amuse_set_minppbin(i)
  ret=0
end function
function get_minppbin(i) result(ret)
  integer :: ret,i
  call amuse_get_minppbin(i) 
  ret=0
end function

function set_targetnn(i) result(ret)
  integer :: ret,i
  call amuse_set_targetnn(i)
  ret=0
end function
function get_targetnn(i) result(ret)
  integer :: ret,i
  call amuse_get_targetnn(i) 
  ret=0
end function

function set_verbosity(i) result(ret)
  integer :: ret,i
  call amuse_set_verbosity(i)
  ret=0
end function
function get_verbosity(i) result(ret)
  integer :: ret,i
  call amuse_get_verbosity(i) 
  ret=0
end function

function set_nsmooth(i) result(ret)
  integer :: ret,i
  call amuse_set_nsmooth(i)
  ret=0
end function
function get_nsmooth(i) result(ret)
  integer :: ret,i
  call amuse_get_nsmooth(i) 
  ret=0
end function

! reals

function set_pboxsize(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_pboxsize(x)
  ret=0
end function
function get_pboxsize(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_pboxsize(x) 
  ret=0
end function

function set_unitm_in_msun(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_unitm_in_msun(x)
  ret=0
end function
function get_unitm_in_msun(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_unitm_in_msun(x) 
  ret=0
end function

function set_unitl_in_kpc(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_unitl_in_kpc(x)
  ret=0
end function
function get_unitl_in_kpc(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_unitl_in_kpc(x) 
  ret=0
end function

function set_dtime(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_dtime(x)
  ret=0
end function

function set_time_step(time_step) result(ret)
  integer :: ret
  real*8 :: time_step
  ret=set_dtime(time_step)
end function

function get_dtime(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_dtime(x) 
  ret=0
end function

function get_time_step(time_step) result(ret)
  integer :: ret
  real*8 :: time_step
  ret=get_dtime(time_step) 
end function

function set_tstepcrit(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_tstepcrit(x)
  ret=0
end function

function get_tstepcrit(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_tstepcrit(x) 
  ret=0
end function

function set_tstpcr2(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_tstpcr2(x)
  ret=0
end function
function get_tstpcr2(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_tstpcr2(x) 
  ret=0
end function

function set_freev(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_freev(x)
  ret=0
end function

function get_freev(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_freev(x) 
  ret=0
end function

function set_freea(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_freea(x)
  ret=0
end function
function get_freea(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_freea(x) 
  ret=0
end function

function set_freevexp(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_freevexp(x)
  ret=0
end function
function get_freevexp(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_freevexp(x) 
  ret=0
end function

function set_freeaexp(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_freeaexp(x)
  ret=0
end function

function get_freeaexp(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_freeaexp(x) 
  ret=0
end function

function set_bh_tol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_bh_tol(x)
  ret=0
end function

function get_bh_tol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_bh_tol(x) 
  ret=0
end function

function set_eps(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_eps(x)
  ret=0
end function

function get_eps(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_eps(x) 
  ret=0
end function

function set_gdgtol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_gdgtol(x)
  ret=0
end function
function get_gdgtol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_gdgtol(x) 
  ret=0
end function

function set_nn_tol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_nn_tol(x)
  ret=0
end function

function get_nn_tol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_nn_tol(x) 
  ret=0
end function

function set_epsgas(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_epsgas(x)
  ret=0
end function
function get_epsgas(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_epsgas(x) 
  ret=0
end function

function set_gamma(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_gamma(x)
  ret=0
end function
function get_gamma(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_gamma(x) 
  ret=0
end function

function set_alpha(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_alpha(x)
  ret=0
end function

function get_alpha(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_alpha(x) 
  ret=0
end function

function set_beta(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_beta(x)
  ret=0
end function
function get_beta(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_beta(x) 
  ret=0
end function

function set_epssph(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_epssph(x)
  ret=0
end function
function get_epssph(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_epssph(x) 
  ret=0
end function

function set_courant(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_courant(x)
  ret=0
end function
function get_courant(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_courant(x) 
  ret=0
end function

function set_removgas(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_removgas(x)
  ret=0
end function
function get_removgas(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_removgas(x) 
  ret=0
end function

function set_consthsm(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_consthsm(x)
  ret=0
end function
function get_consthsm(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_consthsm(x) 
  ret=0
end function

function set_nsmtol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_nsmtol(x)
  ret=0
end function
function get_nsmtol(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_nsmtol(x) 
  ret=0
end function

function set_graineff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_graineff(x)
  ret=0
end function
function get_graineff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_graineff(x) 
  ret=0
end function

function set_crionrate(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_crionrate(x)
  ret=0
end function
function get_crionrate(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_crionrate(x) 
  ret=0
end function

function set_heat_par1(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_heat_par1(x)
  ret=0
end function
function get_heat_par1(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_heat_par1(x) 
  ret=0
end function

function set_heat_par2(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_heat_par2(x)
  ret=0
end function
function get_heat_par2(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_heat_par2(x) 
  ret=0
end function

function set_cool_par(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_cool_par(x)
  ret=0
end function
function get_cool_par(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_cool_par(x) 
  ret=0
end function

function set_optdepth(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_optdepth(x)
  ret=0
end function
function get_optdepth(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_optdepth(x) 
  ret=0
end function

function set_tcollfac(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_tcollfac(x)
  ret=0
end function
function get_tcollfac(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_tcollfac(x) 
  ret=0
end function

function set_masscrit(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_masscrit(x)
  ret=0
end function
function get_masscrit(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_masscrit(x) 
  ret=0
end function

function set_sfeff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_sfeff(x)
  ret=0
end function
function get_sfeff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_sfeff(x) 
  ret=0
end function

function set_tbubble(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_tbubble(x)
  ret=0
end function
function get_tbubble(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_tbubble(x) 
  ret=0
end function

function set_sne_eff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_sne_eff(x)
  ret=0
end function
function get_sne_eff(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_sne_eff(x) 
  ret=0
end function

function set_tsnbeg(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_tsnbeg(x)
  ret=0
end function
function get_tsnbeg(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_tsnbeg(x) 
  ret=0
end function

function set_rhomax(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_set_rhomax(x)
  ret=0
end function
function get_rhomax(x) result(ret)
  integer :: ret
  real*8 :: x
  call amuse_get_rhomax(x) 
  ret=0
end function


! character

function set_halofile(x) result(ret)
  integer :: ret
  character(len=30) :: x
  call amuse_set_halofile(x)
  ret=0
end function
function get_halofile(x) result(ret)
  integer :: ret
  character(len=30) :: x
  call amuse_get_halofile(x)
  ret=0
end function

function set_feedback(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_set_feedback(x)
  ret=0
end function
function get_feedback(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_get_feedback(x) 
  ret=0
end function

function set_sfmode(x) result(ret)
  integer :: ret
  character(len=10) :: x
  call amuse_set_sfmode(x)
  ret=0
end function
function get_sfmode(x) result(ret)
  integer :: ret
  character(len=10) :: x
  call amuse_get_sfmode(x)
  ret=0
end function

function set_hupdatemethod(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_set_hupdatemethod(x)
  ret=0
end function
function get_hupdatemethod(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_get_hupdatemethod(x) 
  ret=0
end function

function set_sph_visc(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_set_sph_visc(x)
  ret=0
end function
function get_sph_visc(x) result(ret)
  integer :: ret
  character(len=4) :: x
  call amuse_get_sph_visc(x)
  ret=0
end function

function set_fi_data_directory(x) result(ret)
  integer :: ret
  character(len=200) :: x
  call amuse_set_fi_data_directory(x)
  ret=0
end function
function get_fi_data_directory(x) result(ret)
  integer :: ret
  character(len=200) :: x
  call amuse_get_fi_data_directory(x)
  ret=0
end function


function get_number_of_sph_particles_removed(x)
  integer :: get_number_of_sph_particles_removed
  integer, intent(out) :: x
  get_number_of_sph_particles_removed = amuse_get_number_of_sph_particles_removed(x)
end function

function get_id_of_removed_sph_particle(x, id_of_removed_particle)
  integer :: get_id_of_removed_sph_particle
  integer, intent(in) :: x
  integer, intent(out) :: id_of_removed_particle
  get_id_of_removed_sph_particle = amuse_get_id_of_removed_sph_particle(x, id_of_removed_particle)
end function

END MODULE

