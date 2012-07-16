subroutine initparameters
        INTEGER input,output,firstsnap,stepout,steplog,nsteps,           &
     &    max_tbin, minppbin,targetnn,verbosity 
        REAL*8 pboxsize,unitl_in_kpc,unitm_in_msun,dtime, tstepcrit,     &
     &  tstpcr2,freev, freea, freevexp, freeaexp,bh_tol,eps,gdgtol,      &
     &  nn_tol
        LOGICAL usesph,radiate,starform,cosmo,sqrttstp,                  &
     &    acc_tstp,freetstp,usequad,directsum,selfgrav,fixthalo,         &
          adaptive_eps, gdgop
        CHARACTER*16 outputfile
        CHARACTER*30 halofile
        CHARACTER*16 inputfile
        CHARACTER*200 datadir
 
        COMMON/io_param/input(128),output(128),firstsnap,stepout,        &
     &	  steplog,verbosity,datadir,inputfile,outputfile,halofile 
        COMMON/sim_param/pboxsize,usesph,radiate,starform,cosmo
        COMMON/units_param/unitl_in_kpc, unitm_in_msun
        COMMON/step_param/dtime, tstepcrit, tstpcr2, freev, freea,       &
     &    freevexp, freeaexp, nsteps, max_tbin, minppbin,                &
     &    sqrttstp, acc_tstp, freetstp
        COMMON/grav_param/bh_tol,eps,gdgtol,nn_tol,targetnn,             &
     &    usequad,directsum,selfgrav,fixthalo,adaptive_eps, gdgop



        REAL*8 epsgas, gamma, alpha, beta, epssph, courant,              & 
     &    removgas,consthsm,nsmtol,graineff,crionrate,heat_par1,         &
     &    heat_par2,cool_par,optdepth,tcollfac,masscrit,sfeff,           &
     &    tbubble,sne_eff, tsnbeg, rhomax
        INTEGER nsmooth
	LOGICAL smoothinput, consph, sphinit, uentropy, isotherm,        &
     &    eps_is_h,balsara
        CHARACTER*4 feedback
        CHARACTER*10 sfmode
        CHARACTER*4 hupdatemethod
        CHARACTER*4 sph_visc

        COMMON/sph_param/epsgas, gamma, alpha, beta, epssph, courant,    & 
     &    removgas,consthsm, nsmtol, nsmooth, smoothinput, consph,       &
     &    sphinit, uentropy, isotherm,eps_is_h,balsara,                  &
     &    hupdatemethod,sph_visc
        COMMON/rad_param/graineff,crionrate,heat_par1,heat_par2,         &
     &    cool_par, optdepth
        COMMON/sf_param/tcollfac,masscrit,sfeff, tbubble,                &
     &    sne_eff, tsnbeg, rhomax, sfmode, feedback

 call muse_start

end subroutine

subroutine setup_module()
  call muse_reset
  call muse_init
end subroutine

subroutine cleanup_module
  call muse_end
  call muse_start
  call muse_reset
end subroutine

subroutine initialize_particles(time)
!f2py intent(callback,hide) external_acc
!f2py intent(callback,hide) external_pot
!f2py use call_external_acc__user__routines
!f2py use call_external_pot__user__routines
!f2py external external_acc
!f2py external external_pot
 real*8 :: time
  call muse_set_time(time)
  call muse_finalize_init 
end subroutine

subroutine reinitialize_particles
! ????????????
end subroutine

function add_dm_particle(ids,mass,x,y,z,vx,vy,vz,eps,npart) result(n)
integer :: npart
integer :: ids(npart),n,muse_get_nbodies
real*8 :: mass(npart),x(npart),y(npart),z(npart), &
 vx(npart),vy(npart),vz(npart),eps(npart)
  call muse_add_particle_dm(ids,mass,x,y,z,vx,vy,vz,eps,npart)
  n=muse_get_nbodies()
end function

function add_star_particle(ids,mass,x,y,z,vx,vy,vz,eps,tf,npart) result(n)
integer :: npart
integer :: ids(npart),n,muse_get_nbodies
real*8 :: mass(npart),x(npart),y(npart),z(npart), &
 vx(npart),vy(npart),vz(npart),eps(npart),tf(npart)
  call muse_add_particle_star(ids,mass,x,y,z,vx,vy,vz,eps,tf,npart)
  n=muse_get_nbodies()
end function

function add_sph_particle(ids,mass,x,y,z,vx,vy,vz,eps,u,npart) result(n)
integer :: npart
integer :: ids(npart),n,muse_get_nbodies
real*8 :: mass(npart),x(npart),y(npart),z(npart), &
 vx(npart),vy(npart),vz(npart),eps(npart),u(npart)
  call muse_add_particle_sph(ids,mass,x,y,z,vx,vy,vz,eps,u,npart)
  n=muse_get_nbodies()
end function

function set_particle(id1,id,mass,x,y,z,vx,vy,vz,eps,npart) result(n)
integer :: npart
integer :: id1(npart),id(npart),n
real*8 :: mass(npart),x(npart),y(npart),z(npart), &
 vx(npart),vy(npart),vz(npart),eps(npart)
! find id1, set it 
 n=1
end function

function set_mass(id,mass) result(error)
integer id,error
real*8 :: mass
 error=1
end function

function set_radius(id,r) result(error)
integer id,error
real*8 :: r
 error=1
end function

function set_pos(id,p) result(error)
integer id,error
real*8 :: p(3)
 error=1
end function

function set_vel(id,v) result(error)
integer id,error
real*8 :: v(3)
 error=1
end function

function remove_particle(id) result(n)
integer id,n
 n=-1
end function

function get_number() result(n)
 integer n,muse_get_nbodies
  n=muse_get_nbodies() 
end function

function get_dynamical_time_scale() result(td)
 real*8 :: td,muse_get_dynamical_time_scale
 td=muse_get_dynamical_time_scale()
end function

function get_time() result(t)
 real*8 :: t,muse_get_time
 t=muse_get_time()
end function

function get_time_step() result(dt)
 real*8 :: dt,muse_get_time_step
 dt=muse_get_time_step()
end function

function initialize_time_step() result(error)
integer :: error
 error=0
end function

function finalize_time_step() result(error)
integer :: error
 error=0
end function

function evolve(tend, sync) result(icoll)
!f2py intent(callback,hide) external_acc
!f2py intent(callback,hide) external_pot
!f2py use call_external_acc__user__routines
!f2py use call_external_pot__user__routines
!f2py external external_acc
!f2py external external_pot
 integer :: sync,icoll
 real*8 :: tend
 call muse_stepsys(tend,sync)
 icoll=-1
end function

function find_colliding_primary() result(ip)
 integer :: ip
 ip=0
end function 

function find_colliding_secondary(ip) result(is)
 integer :: ip,is
 is=0
end function 
 
subroutine get_state(id,mass,x,y,z,vx,vy,vz,eps)
integer :: id
real*8 :: mass,x,y,z,vx,vy,vz,eps
!f2py intent(in) id
!f2py intent(out) mass,x,y,z,vx,vy,vz,eps
  call muse_get_state(id,mass,x,y,z,vx,vy,vz,eps)
end subroutine

function get_mass(id) result(m)
 integer id
 real*8 :: m
 m=0
end function

function get_radius(id) result(r)
 integer id
 real*8 :: r
 r=0
end function

function get_kinetic_energy(mode) result(e)
!f2py intent(callback,hide) external_acc
!f2py intent(callback,hide) external_pot
!f2py use call_external_acc__user__routines
!f2py use call_external_pot__user__routines
!f2py external external_acc
!f2py external external_pot
 real*8 :: e,ek,ep,eth
 integer :: mode
 call muse_energies(mode,ek,ep,eth)
 e=ek
end function 

function get_potential_energy() result(e)
 real*8 :: e,ek,ep,eth
 call muse_energies(0,ek,ep,eth)
 e=ep
end function 

function get_thermal_energy() result(e)
 real*8 :: e,ek,ep,eth
 call muse_energies(0,ek,ep,eth)
 e=eth
end function 

function get_escaper() result(iesc)
 integer :: iesc
 iesc=-1
end function 

subroutine get_gravity(eps,x,y,z,ax,ay,az,n)
 integer :: n
 real*8 :: eps(n),x(n),y(n),z(n),ax(n),ay(n),az(n)
 call muse_get_gravity(eps,x,y,z,ax,ay,az,n) 
end subroutine

subroutine get_tidalfield(eps,x,y,z,tide,n)
 integer :: n
 real*8 :: eps(n),x(n),y(n),z(n),tide(6,n)
 call muse_get_tidalfield(eps,x,y,z,tide,n) 
end subroutine

subroutine get_pot(eps,x,y,z,pot,n)
 integer :: n
 real*8 :: eps(n),x(n),y(n),z(n),pot(n)
 call muse_get_pot(eps,x,y,z,pot,n) 
end subroutine

subroutine call_external_acc(eps,x,y,z,ax,ay,az,n)
  integer n
  double precision :: eps(n), x(n),y(n),z(n)
  double precision :: ax(n),ay(n),az(n)
!f2py intent(callback,hide) external_acc
!f2py external external_acc

 call external_acc(eps,x,y,z,ax,ay,az,n) 

end subroutine

subroutine call_external_pot(eps,x,y,z,phi,n)
  integer n
  double precision :: eps(n), x(n),y(n),z(n)
  double precision :: phi(n)
!f2py intent(callback,hide) external_pot
!f2py external external_pot

 call external_pot(eps,x,y,z,phi,n) 

end subroutine



subroutine set_nstar(nstar)
 integer :: nstar
 call muse_set_nstar(nstar) 
end subroutine
