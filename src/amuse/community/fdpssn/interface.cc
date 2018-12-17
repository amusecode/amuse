// For AMUSE
#include "worker_code.h"

// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif
// Include user-defined headers
#include "src/macro_defs.hpp"
#include "src/mathematical_constants.h"
#include "src/physical_constants.h"
#include "src/user_defined.hpp"
//#include "src/ic.hpp"
#include "src/leapfrog.hpp"

#include <map>

#define DMPARTICLE 0
#define SPHPARTICLE 1


std::map <int, int> type_id;
std::map <int, int> dm_id;
std::map <int, int> sph_id;

PS::ParticleSystem<FP_nbody> psys_nbody;  // Gravity-only
PS::ParticleSystem<FP_sph> psys_sph;  // SPH particles
PS::F64 system_time = 0.0;
PS::S32 nstep_dm = 0;
PS::S32 nstep_sph = 0;
PS::F64 dt;
PS::BOUNDARY_CONDITION bc;

PS::DomainInfo dinfo;
PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole tree_grav;
PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;

// Counter for all particles. This number must never decrease!
PS::U64 unique_particles = 0;

void calcDensity(PS::ParticleSystem<FP_sph> & psys,
                 PS::DomainInfo & dinfo, 
                 PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree) {
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Determine the density and the smoothing length so that Eq.(6) in Springel (2005) 
    // holds within a specified accuracy.
    SCF_smth = 1.25;
    PS::S32 iter = 0;
    for (;;) {
        iter++;
        // Compute density, etc.
        tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
        // Check convergence
        PS::S32 n_compl_loc = 0;
        for (PS::S32 i = 0; i < n_loc; i++) {
            if (psys[i].flag == 1) n_compl_loc++;
        }
        const PS::S64 n_compl = PS::Comm::getSum(n_compl_loc);
        if (n_compl == n_glb) break;
    }
    // Reset SCF_smth
    SCF_smth = 1.0;
#else
    SCF_smth = 1.0;
    tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
#endif
}

void setEntropy(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setEntropy();
    }
}

void setPressure(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setPressure();
    }
}

void setCircularVelocity(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        const PS::F64vec acc = psys[i].acc_grav + psys[i].acc_hydro;
        const PS::F64 r = std::sqrt(psys[i].pos * psys[i].pos);
        const PS::F64 R = std::sqrt(psys[i].pos.x * psys[i].pos.x 
                                   +psys[i].pos.y * psys[i].pos.y);
        const PS::F64 phi = std::atan2(psys[i].pos.y, psys[i].pos.x);
        const PS::F64 theta = std::atan2(R, psys[i].pos.z);
        PS::F64vec base_vect_r;
        base_vect_r.x = std::sin(theta) * std::cos(phi);
        base_vect_r.y = std::sin(theta) * std::sin(phi);
        base_vect_r.z = std::cos(theta);
        const PS::F64 vel_circ = std::sqrt(r * std::abs(acc * base_vect_r));
        psys[i].vel.x = - vel_circ * std::sin(phi);
        psys[i].vel.y =   vel_circ * std::cos(phi);
        psys[i].vel.z = 0.0;
    }
}

PS::F64 getTimeStep(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                    const PS::ParticleSystem<FP_sph>& psys_sph) {
   PS::F64 dt = DBL_MAX; 
   if (dt_max > 0.0) dt = dt_max;

   // Timescale for N-body system
#if defined(ENABLE_GRAVITY_INTERACT)
   for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
      const PS::F64 acc = std::sqrt(psys_nbody[i].acc * psys_nbody[i].acc);
      if (acc > 0.0)
          dt = std::min(dt, CFL_dyn * std::sqrt(eps_grav / acc));
   }
#endif

   // Timescale for SPH system
   for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
#if defined(ENABLE_GRAVITY_INTERACT)
      const PS::F64 acc = std::sqrt((psys_sph[i].acc_grav + psys_sph[i].acc_hydro)
                                   *(psys_sph[i].acc_grav + psys_sph[i].acc_hydro));
      if (acc > 0.0)
          dt = std::min(dt, CFL_dyn * std::sqrt(eps_grav / acc));
#endif
#if defined(ENABLE_HYDRO_INTERACT)
      dt = std::min(dt, psys_sph[i].dt);
#endif
   }
   return PS::Comm::getMinValue(dt);
}

void checkConservativeVariables(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                                const PS::ParticleSystem<FP_sph>& psys_sph) {
    PS::F64    ekin_loc = 0.0;
    PS::F64    epot_loc = 0.0;
    PS::F64    eth_loc  = 0.0; 
    PS::F64vec mom_loc  = 0.0; 
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        ekin_loc += 0.5 * psys_nbody[i].mass * psys_nbody[i].vel * psys_nbody[i].vel;
        epot_loc += 0.5 * psys_nbody[i].mass * (psys_nbody[i].pot + psys_nbody[i].mass / eps_grav);
        mom_loc  += psys_nbody[i].mass * psys_nbody[i].vel;
    }
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
        ekin_loc += 0.5 * psys_sph[i].mass * psys_sph[i].vel * psys_sph[i].vel;
        epot_loc += 0.5 * psys_sph[i].mass * (psys_sph[i].pot_grav + psys_sph[i].mass / eps_grav);
        eth_loc  += psys_sph[i].mass * psys_sph[i].eng;
        mom_loc  += psys_sph[i].mass * psys_sph[i].vel;
    }
    PS::F64 ekin    = PS::Comm::getSum(ekin_loc);
    PS::F64 epot    = PS::Comm::getSum(epot_loc);
    PS::F64 eth     = PS::Comm::getSum(eth_loc);
    PS::F64vec mom  = PS::Comm::getSum(mom_loc);

    static bool is_initialized = false;
    static PS::F64 emech_ini, etot_ini;
    if (is_initialized == false) {
        emech_ini = ekin + epot;
        etot_ini  = ekin + epot + eth;
        is_initialized = true;
    }

    if (PS::Comm::getRank() == 0){
        const PS::F64 emech = ekin + epot;
        const PS::F64 etot  = ekin + epot + eth;
        const PS::F64 relerr_mech = std::fabs((emech - emech_ini)/emech_ini);
        const PS::F64 relerr_tot  = std::fabs((etot  - etot_ini)/etot_ini);
        std::cout << "-------------------------" << std::endl;
        std::cout << "E_kin  = " << ekin  << std::endl;
        std::cout << "E_pot  = " << epot  << std::endl;
        std::cout << "E_th   = " << eth   << std::endl;
        std::cout << "E_mech = " << emech << " (" << relerr_mech << ")" << std::endl;
        std::cout << "E_tot  = " << etot  << " (" << relerr_tot  << ")" << std::endl;
        std::cout << "Mom (x) = " << mom.x << std::endl;
        std::cout << "Mom (y) = " << mom.y << std::endl;
        std::cout << "Mom (z) = " << mom.z << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}


int get_dm_id(int i){
  std::map<int, int>::iterator ii = dm_id.find(i);
  if (ii == dm_id.end()) return -1;
  else return ii->second;
}

int get_sph_id(int i){
  std::map<int, int>::iterator ii = sph_id.find(i);
  if (ii == sph_id.end()) return -1;
  else return ii->second;
}

int get_type(int i){
  std::map<int, int>::iterator ii = type_id.find(i);
  if (ii == type_id.end()) return -1;
  else return ii->second;
}

int initialize_code(){
  // // Initialize FDPS
  // PS::Initialize(argc, argv);

  // // Make a directory
  // char dir_name[1024];
  // sprintf(dir_name,"./result");
  // makeOutputDirectory(dir_name);

  // Display # of MPI processes and threads
  PS::S32 nprocs = PS::Comm::getNumberOfProc();
  PS::S32 nthrds = PS::Comm::getNumberOfThread();
  if (PS::Comm::getRank() == 0) {
      std::cout << "===========================================" << std::endl
                << " This is a sample program of "               << std::endl
                << " Nbody + SPH particle system on FDPS!"       << std::endl
                << " # of processes is " << nprocs               << std::endl
                << " # of thread is    " << nthrds               << std::endl
                << "===========================================" << std::endl;
  }

  // Make instances of ParticleSystem and initialize them
  psys_nbody.initialize();
  psys_sph.initialize();

  // Make an instance of DomainInfo and initialize it
  dinfo.initialize();

  // // Define local variables
  // PS::F64 dt, time_dump, dt_dump, time_end;

  // Defaults, set these with AMUSE setters
  eps_grav = 0.1; // gravitational softening
  bc = PS::BOUNDARY_CONDITION_OPEN;
  dt_max = 0.01; // Maximum allowed timestep, zero means no maximum?
  return 0;
}

int new_dm_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz){
  FP_nbody nbody_particle;
  PS::S64 pid = unique_particles;
  unique_particles++;

  nbody_particle.id = pid;
  nbody_particle.mass = mass;
  nbody_particle.pos.x = x;
  nbody_particle.pos.y = y;
  nbody_particle.pos.z = z;
  nbody_particle.vel.x = vx;
  nbody_particle.vel.y = vy;
  nbody_particle.vel.z = vz;
  nbody_particle.pot = 0.0;

  dm_id[pid] = psys_nbody.getNumberOfParticleGlobal();
  type_id[pid] = DMPARTICLE;
  psys_nbody.addOneParticle(nbody_particle);
  *index_of_the_particle = pid;
  return 0;
}

int new_sph_particle(int * index_of_the_particle, double mass, double x,
    double y, double z, double vx, double vy, double vz, double u, double h_smooth){
  FP_sph sph_particle;
  PS::S64 pid = unique_particles;
  unique_particles++;

  sph_particle.id = pid;
  sph_particle.mass = mass;
  sph_particle.pos.x = x;
  sph_particle.pos.y = y;
  sph_particle.pos.z = z;
  sph_particle.vel.x = vx;
  sph_particle.vel.y = vy;
  sph_particle.vel.z = vz;
  sph_particle.eng = u;
  sph_particle.smth = h_smooth;
  sph_particle.acc_grav = 0.0;
  sph_particle.pot_grav = 0.0;
  sph_particle.acc_hydro = 0.0;
  sph_particle.flag = 0;
  sph_particle.dens = 0.0;
  sph_particle.ent = 0.0;
  sph_particle.pres = 0.0;
  sph_particle.gradh = 0.0;
  sph_particle.divv = 0.0;
  sph_particle.rotv.x = 0.0;
  sph_particle.rotv.y = 0.0;
  sph_particle.rotv.z = 0.0;
  sph_particle.BalSW = 0.0;
  sph_particle.snds = 0.0;
  sph_particle.eng_dot = 0.0;
  sph_particle.ent_dot = 0.0;
  sph_particle.dt = 0.0;
  sph_particle.eng_half = 0.0;
  sph_particle.ent_half = 0.0;


  sph_id[pid] = psys_sph.getNumberOfParticleGlobal();
  type_id[pid] = SPHPARTICLE;
  psys_sph.addOneParticle(sph_particle);
  *index_of_the_particle = pid;
  return 0;
}

int commit_particles(){
  std::cout << "Commit particles called here" << std::endl;
  bool have_dm = false, have_sph = false;
  if (psys_nbody.getNumberOfParticleGlobal() > 0){
    have_dm = true;
    if (nstep_dm == 0){
      nstep_dm++;
    }
  }
  else {
    nstep_dm = 0;
  }
  if (psys_sph.getNumberOfParticleGlobal() > 0){
    have_sph = true;
    if (nstep_sph == 0){
      nstep_sph++;
    }
  }
  else {
    nstep_sph = 0;
  }
  if (nstep_dm == 1){
    // Perform domain decomposition 
    dinfo.collectSampleParticle(psys_nbody);
  }
  if (nstep_sph == 1){
    dinfo.collectSampleParticle(psys_sph,false);
  }

  if (nstep_sph == 1 or nstep_dm == 1){
    dinfo.decomposeDomain();
  
    // Perform particle exchange
    if (have_dm){
      psys_nbody.exchangeParticle(dinfo);
    }
    if (have_sph){
      psys_sph.exchangeParticle(dinfo);
    }

    // Make tree structures
    const PS::S64 numPtclSPH = std::max(psys_sph.getNumberOfParticleLocal(),1);
    const PS::S64 numPtclAll = psys_nbody.getNumberOfParticleLocal() + numPtclSPH;
  
    const PS::F32 theta_grav = 0.5;
    tree_grav.initialize(3 * numPtclAll, theta_grav);

    if (have_sph){
      tree_dens.initialize(3 * numPtclSPH);
  
      tree_hydro.initialize(3 * numPtclSPH);
    }

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_open();
    g5_set_eps_to_all(eps_grav);
#endif

    // Peform force calculations 
    //- Gravity calculations
#if defined(ENABLE_GRAVITY_INTERACT)
    if (have_dm){
      tree_grav.setParticleLocalTree(psys_nbody);
    }
    if (have_sph){
      tree_grav.setParticleLocalTree(psys_sph,false);
    }
    tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                  CalcGravity<PS::SPJMonopole>,
                                  dinfo);
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        psys_nbody[i].copyFromForce(tree_grav.getForce(i));
    }
    const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
        psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
    }
#endif

    if (have_sph){
      //- SPH calculations
#if defined(ENABLE_HYDRO_INTERACT)
      calcDensity(psys_sph, dinfo, tree_dens);
#if defined(USE_ENTROPY)
      setEntropy(psys_sph);
#endif
      setPressure(psys_sph);
      tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_sph, dinfo);
#endif

      // Set the initial velocity of gas particles
      // Skipping this in AMUSE because AMUSE handles initial conditions
// #if defined(SET_CIRCULAR_VELOCITY)
//       setCircularVelocity(psys_sph);
// #endif
    }

    // Get timestep
    dt = getTimeStep(psys_nbody, psys_sph);

    // Calculate energies
    if (have_sph and have_dm){
      checkConservativeVariables(psys_nbody, psys_sph);
    }
  }

  return 0;
}

int evolve_model(double time_end){
  std::cout << "evolve_model started" << std::endl;
  bool have_dm = false, have_sph = false;
  if (psys_nbody.getNumberOfParticleGlobal() > 0){
    have_dm = true;
    if (nstep_dm == 0){
      nstep_dm++;
    }
  }
  else {
    nstep_dm = 0;
  }
  if (psys_sph.getNumberOfParticleGlobal() > 0){
    have_sph = true;
    if (nstep_sph == 0){
      nstep_sph++;
    }
  }
  else {
    nstep_sph = 0;
  }
  // Main loop for time integration
  while (system_time < time_end){
  //for (PS::F64 time = 0.0; time < time_end; time += dt, nstep++){
    if (PS::Comm::getRank() == 0) {
        std::cout << "nstep_dm = " << nstep_dm
                  << " nstep_sph = " << nstep_sph
                  << " dt = " << dt 
                  << " time = " << system_time 
                  << " time_end = " << time_end
                  << std::endl;
    }

    // Leap frog: Initial Kick & Full Drift
    if (have_dm){
      InitialKick(psys_nbody, dt);
    }
    if (have_sph){
      InitialKick(psys_sph, dt);
    }
    if (have_dm){
      FullDrift(psys_nbody, dt);
    }
    if (have_sph){
      FullDrift(psys_sph, dt);
    }
    if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
      if (have_dm){
        psys_nbody.adjustPositionIntoRootDomain(dinfo);
      }
      if (have_sph){
        psys_sph.adjustPositionIntoRootDomain(dinfo);
      }
    }

    // Leap frog: Predict
    if (have_sph){
      Predict(psys_sph, dt);
    }

    // Perform domain decomposition again
    if (have_dm){
      dinfo.collectSampleParticle(psys_nbody);
    }
    if (have_sph){
      dinfo.collectSampleParticle(psys_sph,false);
    }
    dinfo.decomposeDomain();

    // Exchange the particles between the (MPI) processes
    if (have_dm){
      psys_nbody.exchangeParticle(dinfo);
    }
    if (have_sph){
      psys_sph.exchangeParticle(dinfo);
    }

    // Peform force calculations
    PS::F64 t_start; 
    //- Gravity calculations
    PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_GRAVITY_INTERACT)
    if (have_dm){
      tree_grav.setParticleLocalTree(psys_nbody);
    }
    if (have_sph){
      tree_grav.setParticleLocalTree(psys_sph,false);
    }
    tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                  CalcGravity<PS::SPJMonopole>,
                                  dinfo);
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        psys_nbody[i].copyFromForce(tree_grav.getForce(i));
    }
    const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
        psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
    }
#endif
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "t_grav = " << (PS::GetWtime() - t_start) << std::endl;
    //- SPH calculations
    PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_HYDRO_INTERACT)
    if (have_sph){
      calcDensity(psys_sph, dinfo, tree_dens);
      setPressure(psys_sph);
      tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_sph, dinfo);
    }
#endif
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "t_hydro = " << (PS::GetWtime() - t_start) << std::endl;

    // Get a new timestep
    dt = getTimeStep(psys_nbody, psys_sph);

    // Leap frog: Final Kick
    if (have_dm){
      FinalKick(psys_nbody, dt);
    }
    if (have_sph){
      FinalKick(psys_sph, dt);
    }

    // Calculate energies
    if (have_sph and have_dm){
      checkConservativeVariables(psys_nbody, psys_sph);
    }

    // Check the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
    if (nstep_sph % 100 == 0) 
        checkDensityFluctuation(psys_sph);
#endif
    system_time += dt;
    nstep_dm++;
    nstep_sph++;
  }

  return 0;
}

int get_mass(int index_of_the_particle, double * mass){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    *mass = psys_nbody[i].mass;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *mass = psys_sph[i].mass;
  }
  else{
    return -1;
  }
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    psys_nbody[i].mass = mass;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    psys_sph[i].mass = mass;
  }
  else{
    return -1;
  }
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    *x = psys_nbody[i].pos.x;
    *y = psys_nbody[i].pos.y;
    *z = psys_nbody[i].pos.z;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *x = psys_sph[i].pos.x;
    *y = psys_sph[i].pos.y;
    *z = psys_sph[i].pos.z;
  }
  else{
    return -1;
  }
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    psys_nbody[i].pos.x = x;
    psys_nbody[i].pos.y = y;
    psys_nbody[i].pos.z = z;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    psys_sph[i].pos.x = x;
    psys_sph[i].pos.y = y;
    psys_sph[i].pos.z = z;
  }
  else{
    return -1;
  }
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    *vx = psys_nbody[i].vel.x;
    *vy = psys_nbody[i].vel.y;
    *vz = psys_nbody[i].vel.z;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *vx = psys_sph[i].vel.x;
    *vy = psys_sph[i].vel.y;
    *vz = psys_sph[i].vel.z;
  }
  else{
    return -1;
  }
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    psys_nbody[i].vel.x = vx;
    psys_nbody[i].vel.y = vy;
    psys_nbody[i].vel.z = vz;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    psys_sph[i].vel.x = vx;
    psys_sph[i].vel.y = vy;
    psys_sph[i].vel.z = vz;
  }
  else{
    return -1;
  }
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    *ax = psys_nbody[i].acc.x;
    *ay = psys_nbody[i].acc.y;
    *az = psys_nbody[i].acc.z;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *ax = psys_sph[i].acc_grav.x;
    *ay = psys_sph[i].acc_grav.y;
    *az = psys_sph[i].acc_grav.z;
  }
  else{
    return -1;
  }
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
  return -1;
}

int get_radius(int index_of_the_particle, double * radius){
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  return 0;
}

int get_internal_energy(int index_of_the_particle, double * u){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    // Property not available for DM
    return -1;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *u = psys_sph[i].eng;
  }
  else{
    return -1;
  }
  return 0;
}

int set_internal_energy(int index_of_the_particle, double u){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    // Property not available for DM
    return -1;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    psys_sph[i].eng = u;
  }
  else{
    return -1;
  }
  return 0;
}


int get_smoothing_length(int index_of_the_particle, double * h_smooth){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    // Property not available for DM
    return -1;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *h_smooth = psys_sph[i].smth;
  }
  else{
    return -1;
  }
  return 0;
}

int set_smoothing_length(int index_of_the_particle, double h_smooth){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    // Property not available for DM
    return -1;
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    psys_sph[i].smth = h_smooth;
  }
  else{
    return -1;
  }
  return 0;
}

int get_pressure(int index_of_the_particle, double * pressure){
  int t = get_type(index_of_the_particle);
  if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *pressure = psys_sph[i].pres;
  }
  else{
    return -1;
  }
  return 0;
}

int get_density(int index_of_the_particle, double * rho){
  int t = get_type(index_of_the_particle);
  if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    *rho = psys_sph[i].dens;
  }
  else{
    return -1;
  }
  return 0;
}

int get_state_sph(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, double * u,
  double * h_smooth){
  int i = get_sph_id(index_of_the_particle);
  *mass = psys_sph[i].mass;
  *x  = psys_sph[i].pos.x;
  *y  = psys_sph[i].pos.y;
  *z  = psys_sph[i].pos.z;
  *vx = psys_sph[i].vel.x;
  *vy = psys_sph[i].vel.y;
  *vz = psys_sph[i].vel.z;
  *u  = psys_sph[i].eng;
  *h_smooth = psys_sph[i].smth;
  return 0;
}

int set_state_sph(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double u, double h_smooth){
  int i = get_sph_id(index_of_the_particle);
  psys_sph[i].mass = mass;
  psys_sph[i].pos.x = x;
  psys_sph[i].pos.y = y;
  psys_sph[i].pos.z = z;
  psys_sph[i].vel.x = vx;
  psys_sph[i].vel.y = vy;
  psys_sph[i].vel.z = vz;
  psys_sph[i].eng = u;
  psys_sph[i].smth = h_smooth;
  return 0;
}

int get_state_dm(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz){
  int i = get_dm_id(index_of_the_particle);
  *mass = psys_nbody[i].mass;
  *x = psys_nbody[i].pos.x;
  *y = psys_nbody[i].pos.y;
  *z = psys_nbody[i].pos.z;
  *vx = psys_nbody[i].vel.x;
  *vy = psys_nbody[i].vel.y;
  *vz = psys_nbody[i].vel.z;
  return 0;
}

int set_state_dm(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz){
  int i = get_dm_id(index_of_the_particle);
  psys_nbody[i].mass = mass;
  psys_nbody[i].pos.x = x;
  psys_nbody[i].pos.y = y;
  psys_nbody[i].pos.z = z;
  psys_nbody[i].vel.x = vx;
  psys_nbody[i].vel.y = vy;
  psys_nbody[i].vel.z = vz;
  return 0;
}

int delete_particle(int index_of_the_particle){
  int t = get_type(index_of_the_particle);
  if (t == DMPARTICLE){
    int i = get_dm_id(index_of_the_particle);
    FP_nbody p = psys_nbody[i];
    PS::S32 n_remove = 1;
    PS::S32 *id[n_remove] = { new PS::S32[n_remove] };
    *id[0] = p.id;
    psys_nbody.removeParticle(*id, n_remove);
    dm_id.erase(index_of_the_particle);
  }
  else if (t == SPHPARTICLE){
    int i = get_sph_id(index_of_the_particle);
    FP_sph p = psys_sph[i];
    PS::S32 n_remove = 1;
    PS::S32 *id[n_remove] = { new PS::S32[n_remove] };
    *id[0] = p.id;
    psys_sph.removeParticle(*id, n_remove);
    sph_id.erase(index_of_the_particle);
  }
  else{
    return -1;
  }
  return 0;
}

int get_time(double * time){
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int get_total_mass(double * mass){
  return 0;
}

int set_eps2(double epsilon_squared){
  return 0;
}

int get_begin_time(double * time){
  return 0;
}

int get_eps2(double * epsilon_squared){
  return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  return 0;
}

int synchronize_model(){
  return 0;
}

int get_time_step(double * time_step){
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  return 0;
}

int get_number_of_particles(int * number_of_particles){
  PS::S64 nbody_max = psys_nbody.getNumberOfParticleGlobal();
  PS::S64 sph_max = psys_sph.getNumberOfParticleGlobal();
  * number_of_particles = nbody_max + sph_max;
  return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int set_begin_time(double time){
  return 0;
}

int cleanup_code(){
  return 0;
}

int recommit_parameters(){
  return 0;
}

int get_potential_energy(double * potential_energy){
  return 0;
}

int commit_parameters(){
  return 0;
}

int set_epsilon_squared(double epsilon_squared){
  if (epsilon_squared <= 0.){
    return -2; // illegal value
  }
  eps_grav = sqrt(epsilon_squared);
  return 0;
}

int get_epsilon_squared(double * epsilon_squared){
  * epsilon_squared = eps_grav * eps_grav;
  return 0;
}
