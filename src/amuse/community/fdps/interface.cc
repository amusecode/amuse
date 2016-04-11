#include "worker_code.h"

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_gpu_cuda.hpp"
#endif
#include "user-defined.hpp"

#include <map>

PS::ParticleSystem<FPGrav> system_grav;
PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
PS::S32 n_loc = 0;
PS::S64 n_tot = 0;
PS::F32 time_sys = 0.0;
PS::F64 FPGrav::eps = sqrt(1.0/8.0);

// Not sure if these are needed, but they're used in the example so
// let's set them.
PS::F32 theta = 0.75;
PS::S32 n_leaf_limit = 8;
PS::S32 n_group_limit = 64;

PS::F32 time_end = 10.0;
PS::F32 dt = 1.0 / 64.0;
PS::S32 c;
PS::F32 coef_ema = 0.3;
PS::DomainInfo dinfo;

PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;

PS::S64 n_loop = 0;

std::map <int,int> inverse_id;

int get_inverse_id(int i)
{
    std::map<int,int>::iterator ii = inverse_id.find(i);
    if (ii == inverse_id.end()) return -1;
    else return ii->second;
}


template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void calcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
#else
    etot = etot_loc;
    epot = epot_loc;
    ekin = ekin_loc;
#endif
}

template<class Tpsys>
void calcPotentialEnergy(const Tpsys & system,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        epot = 0.0;
    }
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
    }
    epot_loc *= 0.5;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    epot = PS::Comm::getSum(epot_loc);
#else
    epot = epot_loc;
#endif
}

template<class Tpsys>
void calcKineticEnergy(const Tpsys & system,
                PS::F64 & ekin,
                const bool clear=true){
    if(clear){
        ekin = 0.0;
    }
    PS::F64 ekin_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
    }
    ekin_loc *= 0.5;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    ekin = PS::Comm::getSum(ekin_loc);
#else
    ekin = ekin_loc;
#endif
}


int get_mass(int index_of_the_particle, double * mass){
  int j = get_inverse_id(index_of_the_particle);
  
  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    *mass = system_grav[j].mass;
    return 0;
  }
  return -1;
}

int commit_particles(){
  dinfo.initialize(coef_ema);
  dinfo.collectSampleParticle(system_grav);
  dinfo.decomposeDomain();
  system_grav.exchangeParticle(dinfo);
  n_loc = system_grav.getNumberOfParticleLocal();

  tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

  tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                     CalcGravity<PS::SPJMonopole>,
                                     system_grav,
                                     dinfo);

  calcEnergy(system_grav, Etot0, Ekin0, Epot0);
  //PS::S64 n_loop = 0;
  //PS::S32 id_snap = 0;

  return 0;
}

int get_time(double * time){
  *time = time_sys;
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].mass = mass;
    return 0;
  }
  return -1;
}

int get_index_of_first_particle(int * index_of_the_particle){
  *index_of_the_particle = system_grav[0].id;
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int new_sph_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
  return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
  n_loc = system_grav.getNumberOfParticleLocal();

  int pid = n_loc;//FIXME: this must be a unique id, not to be reused!
  int new_id = n_loc;

  n_loc++;
  n_tot++;
  system_grav.setNumberOfParticleLocal(n_loc);

  system_grav[new_id].mass = mass; 
  system_grav[new_id].id = pid;
  system_grav[new_id].pos.x = x;
  system_grav[new_id].pos.y = y;
  system_grav[new_id].pos.z = z;
  system_grav[new_id].vel.x = vx;
  system_grav[new_id].vel.y = vy;
  system_grav[new_id].vel.z = vz;
  system_grav[new_id].radius = radius;

  //FIXME: need to create a "reverse id" lookup table and use that here
  inverse_id[pid] = new_id;
  *index_of_the_particle = pid;
    
  return 0;
}

int get_total_mass(double * mass){

  *mass = 0;

  for(int j = 0; j < system_grav.getNumberOfParticleLocal(); j++)
  {
    *mass += system_grav[j].mass;
  }
  
  return 0;
}

int evolve_model(double time){
  while(time_sys < time){
    calcEnergy(system_grav, Etot1, Ekin1, Epot1);
    kick(system_grav, dt * 0.5);
    time_sys += dt;
    drift(system_grav, dt);
    //if(n_loop % 4 == 0){
    //  dinfo.decomposeDomainAll(system_grav);
    //}
    //system_grav.exchangeParticle(dinfo);
    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav,
                                       dinfo);
    kick(system_grav, dt * 0.5);
          
    n_loop++;
  }

  return 0;
}

int set_eps2(double epsilon_squared){
  return -2;
}

int set_epsilon_squared(double epsilon_squared){
  FPGrav::eps = sqrt(epsilon_squared);
  return 0;
}

int get_begin_time(double * time){
  return 0;
}

int get_eps2(double * epsilon_squared){
  return -2;
}

int get_epsilon_squared(double * epsilon_squared){
  *epsilon_squared = FPGrav::eps * FPGrav::eps;
  return 0;
}

int get_theta_for_tree(double *theta_for_tree){
    *theta_for_tree = theta;
    return 0;
}
int set_theta_for_tree(double theta_for_tree){
    theta = theta_for_tree;
    return 0;
}

int get_group_limit_for_tree(int *group_limit){
    *group_limit = n_group_limit;
    return 0;
}
int set_group_limit_for_tree(int group_limit){
    n_group_limit = group_limit;
    return 0;
}

int get_leaf_limit_for_tree(int *leaf_limit){
    *leaf_limit = n_leaf_limit;
    return 0;
}
int set_leaf_limit_for_tree(int leaf_limit){
    n_leaf_limit = leaf_limit;
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
  int j = get_inverse_id(index_of_the_particle);
  n_loc = system_grav.getNumberOfParticleLocal();

  if (j < 0) return -1;//means "not found"
  else if (j >= n_loc-1) return 1;//means "at the end of the array"
  else *index_of_the_next_particle = system_grav[j+1].id;
  return 0;
}

int delete_particle(int index_of_the_particle){
  int j = get_inverse_id(index_of_the_particle);

  n_loc = system_grav.getNumberOfParticleLocal();
  for(PS::S32 i = j; i < n_loc-1; i++){
    system_grav[i].mass   = system_grav[i+1].mass;
    system_grav[i].pos    = system_grav[i+1].pos;
    system_grav[i].vel    = system_grav[i+1].vel;
    system_grav[i].acc    = system_grav[i+1].acc;
    system_grav[i].radius = system_grav[i+1].radius;
    system_grav[i].id     = system_grav[i+1].id;
    inverse_id[system_grav[i+1].id] = i;
  }
  n_loc--;
  n_tot--;
  system_grav.setNumberOfParticleLocal(n_loc);
  inverse_id.erase(index_of_the_particle);

  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  int j = get_inverse_id(index_of_the_particle);
  double pot;
  double r2_inv;
  double r_inv;
  r2_inv = FPGrav::eps * FPGrav::eps;
  r_inv = 1.0/sqrt(r2_inv);
  pot = system_grav[j].pot + system_grav[j].mass * r_inv;
  *potential = pot;
  return 0;
}

int synchronize_model(){
  return 0;
}

int set_internal_energy(int index_of_the_particle, double u, int length){
  return 0;
}

int get_internal_energy(int index_of_the_particle, double * u, int * length){
  return 0;
}

int set_state_sph(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double u){
  return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double radius){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].mass = mass;
    system_grav[j].pos.x = x;
    system_grav[j].pos.y = y;
    system_grav[j].pos.z = z;
    system_grav[j].vel.x = vx;
    system_grav[j].vel.y = vy;
    system_grav[j].vel.z = vz;
    system_grav[j].radius = radius;
    
    return 0;
  }
  return -1;
}

int get_state_sph(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, 
  double * u){
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, 
  double * radius){
  int j = get_inverse_id(index_of_the_particle);

  *mass   = system_grav[j].mass;
  *radius = system_grav[j].radius;
  *x      = system_grav[j].pos.x;
  *y      = system_grav[j].pos.y;
  *z      = system_grav[j].pos.z;
  *vx     = system_grav[j].vel.x;
  *vy     = system_grav[j].vel.y;
  *vz     = system_grav[j].vel.z;

  return 0;
}

int get_time_step(double * time_step){
  *time_step = dt;
  return 0;
}

int set_time_step(double time_step){
  dt = time_step;
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  PS::F64 ekin;
  calcKineticEnergy(system_grav, ekin);
  *kinetic_energy = ekin;
  return 0;
}

int get_number_of_particles(int * value){
  *value = system_grav.getNumberOfParticleLocal();
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].acc.x = ax;
    system_grav[j].acc.y = ay;
    system_grav[j].acc.z = az;
    return 0;
  }
  return -1;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  *x = 0; *y = 0; *z = 0;
  double M;
  double mass;

  get_total_mass(&M);

  for(int j = 0; j<system_grav.getNumberOfParticleLocal(); j++)
  {
    mass = system_grav[j].mass;
    *x += mass*system_grav[j].pos.x;
    *y += mass*system_grav[j].pos.y;
    *z += mass*system_grav[j].pos.z;
  }

  *x /= M;
  *y /= M;
  *z /= M;

  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  *vx = 0; *vy = 0; *vz = 0;
  double M;
  double mass;

  get_total_mass(&M);

  for(int j = 0; j<system_grav.getNumberOfParticleLocal(); j++)
  {
    mass = system_grav[j].mass;
    *vx += mass*system_grav[j].vel.x;
    *vy += mass*system_grav[j].vel.y;
    *vz += mass*system_grav[j].vel.z;
  }

  *vx /= M;
  *vy /= M;
  *vz /= M;

  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  int j = get_inverse_id(index_of_the_particle);
  *radius = system_grav[j].radius;
  return 0;
}

int set_begin_time(double time){
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].radius = radius;
    return 0;
  }
  return -1;
}

int cleanup_code(){
  return 0;
}

int recommit_parameters(){
  return 0;
}

int initialize_code(){

//    std::cout<<std::setprecision(15);
//    std::cerr<<std::setprecision(15);

//    PS::Initialize(argc, argv);
//    PS::F32 theta = 0.5;
//    PS::S32 n_leaf_limit = 8;
//    PS::S32 n_group_limit = 64;
//    PS::F32 time_end = 10.0;
//    PS::F32 dt = 1.0 / 128.0;
//    PS::F32 dt_diag = 1.0 / 8.0;
//    PS::F32 dt_snap = 1.0;
//    char dir_name[1024];
//    PS::S64 n_tot = 1024;
//    PS::S32 c;
//    sprintf(dir_name,"./result");
//    opterr = 0;
//    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
//        switch(c){
//        case 'o':
//            sprintf(dir_name,optarg);
//            break;
//        case 't':
//            theta = atof(optarg);
//            std::cerr << "theta =" << theta << std::endl;
//            break;
//        case 'T':
//            time_end = atof(optarg);
//            std::cerr << "time_end = " << time_end << std::endl;
//            break;
//        case 's':
//            dt = atof(optarg);
//            std::cerr << "time_step = " << dt << std::endl;
//            break;
//        case 'd':
//            dt_diag = atof(optarg);
//            std::cerr << "dt_diag = " << dt_diag << std::endl;
//            break;
//        case 'D':
//            dt_snap = atof(optarg);
//            std::cerr << "dt_snap = " << dt_snap << std::endl;
//            break;
//        case 'l':
//            n_leaf_limit = atoi(optarg);
//            std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
//            break;
//        case 'n':
//            n_group_limit = atoi(optarg);
//            std::cerr << "n_group_limit = " << n_group_limit << std::endl;
//            break;
//        case 'N':
//            n_tot = atoi(optarg);
//            std::cerr << "n_tot = " << n_tot << std::endl;
//            break;
//        case 'h':
//            if(PS::Comm::getRank() == 0) {
//                printHelp();
//            }
//            PS::Finalize();
//            return 0;
//        default:
//            if(PS::Comm::getRank() == 0) {
//                std::cerr<<"No such option! Available options are here."<<std::endl;
//                printHelp();
//            }
//            PS::Abort();
//        }
//    }

//    makeOutputDirectory(dir_name);

//    std::ofstream fout_eng;
//    char sout_de[1024];
//    sprintf(sout_de, "%s/t-de.dat", dir_name);
//    std::cerr << sout_de << std::endl;
//    fout_eng.open(sout_de);

//    if(PS::Comm::getRank() == 0) {
//        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
//        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
//    }

//    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
//    PS::S32 n_loc    = 0;
//    PS::F32 time_sys = 0.0;
//    if(PS::Comm::getRank() == 0) {
//        setParticlesColdUniformSphere(system_grav, n_tot, n_loc);
//    } else {
//        system_grav.setNumberOfParticleLocal(n_loc);
//    }
//
//    const PS::F32 coef_ema = 0.3;
//    PS::DomainInfo dinfo;
//    dinfo.initialize(coef_ema);
//    dinfo.collectSampleParticle(system_grav);
//    dinfo.decomposeDomain();
//    system_grav.exchangeParticle(dinfo);
//    n_loc = system_grav.getNumberOfParticleLocal();
//    
//#ifdef ENABLE_PHANTOM_GRAPE_X86
//    g5_open();
//    g5_set_eps_to_all(FPGrav::eps);
//#endif
//    
//    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
//    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
//#ifdef MULTI_WALK
//    const PS::S32 n_walk_limit = 200;
//    const PS::S32 tag_max = 1;
//    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
//                                                RetrieveKernel,
//                                                tag_max,
//                                                system_grav,
//                                                dinfo,
//                                                n_walk_limit);
//#else
//    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
//                                       CalcGravity<PS::SPJMonopole>,
//                                       system_grav,
//                                       dinfo);
//#endif
//    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
//    calcEnergy(system_grav, Etot0, Ekin0, Epot0);
//    PS::F64 time_diag = 0.0;
//    PS::F64 time_snap = 0.0;
//    PS::S64 n_loop = 0;
//    PS::S32 id_snap = 0;
//    while(time_sys < time_end){
//        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
//            char filename[256];
//            sprintf(filename, "%s/%04d.dat", dir_name, id_snap++);
//            FileHeader header;
//            header.time   = time_sys;
//            header.n_body = system_grav.getNumberOfParticleGlobal();
//            system_grav.writeParticleAscii(filename, header);
//            time_snap += dt_snap;
//        }
//
//        calcEnergy(system_grav, Etot1, Ekin1, Epot1);
//        
//        if(PS::Comm::getRank() == 0){
//            if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
//                fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 << std::endl;
//                fprintf(stderr, "time: %10.7f energy error: %+e\n",
//                        time_sys, (Etot1 - Etot0) / Etot0);
//                time_diag += dt_diag;
//            }            
//        }
//        
//        
//        kick(system_grav, dt * 0.5);
//        
//        time_sys += dt;
//        drift(system_grav, dt);
//        
//        if(n_loop % 4 == 0){
//            dinfo.decomposeDomainAll(system_grav);
//        }
//        
//        system_grav.exchangeParticle(dinfo);
//#ifdef MULTI_WALK
//        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
//                                                    RetrieveKernel,
//                                                    tag_max,
//                                                    system_grav,
//                                                    dinfo,
//                                                    n_walk_limit,
//                                                    true);
//#else
//        tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
//                                           CalcGravity<PS::SPJMonopole>,
//                                           system_grav,
//                                           dinfo);
//#endif
//        
//        kick(system_grav, dt * 0.5);
//        
//        n_loop++;
//    }
//    
//#ifdef ENABLE_PHANTOM_GRAPE_X86
//    g5_close();
//#endif
//
//    PS::Finalize();
//    return 0;

  return 0;
}

int get_potential_energy(double * potential_energy){
  //*potential_energy = Epot0;
  PS::F64 epot;
  calcPotentialEnergy(system_grav, epot);
  *potential_energy = epot;
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    *vx = system_grav[j].vel.x;
    *vy = system_grav[j].vel.y;
    *vz = system_grav[j].vel.z;
    return 0;
  }
  return -1;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    *x = system_grav[j].pos.x;
    *y = system_grav[j].pos.y;
    *z = system_grav[j].pos.z;
    return 0;
  }
  return -1;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].pos.x = x;
    system_grav[j].pos.y = y;
    system_grav[j].pos.z = z;
    return 0;
  }
  return -1;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    *ax = system_grav[j].acc.x;
    *ay = system_grav[j].acc.y;
    *az = system_grav[j].acc.z;
    return 0;
  }
  return -1;
}

int commit_parameters(){
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz){
  int j = get_inverse_id(index_of_the_particle);

  if (j >= 0 && j < system_grav.getNumberOfParticleLocal())
  {
    system_grav[j].vel.x = vx;
    system_grav[j].vel.y = vy;
    system_grav[j].vel.z = vz;
    return 0;
  }
  return -1;
}

