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


PS::ParticleSystem<FPGrav> system_grav;
PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
PS::S32 n_loc = 0;
PS::F32 time_sys = 0.0;
PS::F64 FPGrav::eps = 0.0;


int get_mass(int index_of_the_particle, double * mass){
  *mass = system_grav[index_of_the_particle].mass;
  return 0;
}

int commit_particles(){
  return 0;
}

int get_time(double * time){
  *time = time_sys;
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  system_grav[index_of_the_particle].mass = mass;
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
  *index_of_the_particle = system_grav[0].id;
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
  n_loc = system_grav.getNumberOfParticleLocal();

  int new_id = n_loc;

  n_loc++;
  system_grav.setNumberOfParticleLocal(n_loc);

  system_grav[new_id].mass = mass; 
  system_grav[new_id].id = new_id;
  system_grav[new_id].pos.x = x;
  system_grav[new_id].pos.y = y;
  system_grav[new_id].pos.z = z;
  system_grav[new_id].vel.x = vx;
  system_grav[new_id].vel.y = vy;
  system_grav[new_id].vel.z = vz;
  system_grav[new_id].radius = radius;

  *index_of_the_particle = new_id;
    
  return 0;
}

int get_total_mass(double * mass){
  return 0;
}

int evolve_model(double time){
  return 0;
}

int set_eps2(double epsilon_squared){
  FPGrav::eps = sqrt(epsilon_squared);
  return 0;
}

int get_begin_time(double * time){
  return 0;
}

int get_eps2(double * epsilon_squared){
  *epsilon_squared = FPGrav::eps * FPGrav::eps;
  return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
  return 0;
}

int delete_particle(int index_of_the_particle){
  n_loc = system_grav.getNumberOfParticleLocal();
  for(PS::S32 i = index_of_the_particle; i < n_loc-1; i++){
    system_grav[i].mass   = system_grav[i+1].mass;
    system_grav[i].pos    = system_grav[i+1].pos;
    system_grav[i].vel    = system_grav[i+1].vel;
    system_grav[i].acc    = system_grav[i+1].acc;
    system_grav[i].radius = system_grav[i+1].radius;
    system_grav[i].id     = system_grav[i+1].id;
  }
  n_loc--;
  system_grav.setNumberOfParticleLocal(n_loc);    

  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  *potential = system_grav[index_of_the_particle].pot;
  return 0;
}

int synchronize_model(){
  return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, 
  double * radius){

  *mass   = system_grav[index_of_the_particle].mass;
  *radius = system_grav[index_of_the_particle].radius;
  *x      = system_grav[index_of_the_particle].pos.x;
  *y      = system_grav[index_of_the_particle].pos.y;
  *z      = system_grav[index_of_the_particle].pos.z;
  *vx     = system_grav[index_of_the_particle].vel.x;
  *vy     = system_grav[index_of_the_particle].vel.y;
  *vz     = system_grav[index_of_the_particle].vel.z;

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

int get_number_of_particles(int * value){
  *value = system_grav.getNumberOfParticleLocal();
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
  system_grav[index_of_the_particle].acc.x = ax;
  system_grav[index_of_the_particle].acc.y = ay;
  system_grav[index_of_the_particle].acc.z = az;
  return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  *radius = system_grav[index_of_the_particle].radius;
  return 0;
}

int set_begin_time(double time){
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  system_grav[index_of_the_particle].radius = radius;
  return 0;
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
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz){
  *vx = system_grav[index_of_the_particle].vel.x;
  *vy = system_grav[index_of_the_particle].vel.y;
  *vz = system_grav[index_of_the_particle].vel.z;
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z){
  *x = system_grav[index_of_the_particle].pos.x;
  *y = system_grav[index_of_the_particle].pos.y;
  *z = system_grav[index_of_the_particle].pos.z;
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  system_grav[index_of_the_particle].pos.x = x;
  system_grav[index_of_the_particle].pos.y = y;
  system_grav[index_of_the_particle].pos.z = z;
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az){
  *ax = system_grav[index_of_the_particle].acc.x;
  *ay = system_grav[index_of_the_particle].acc.y;
  *az = system_grav[index_of_the_particle].acc.z;
  return 0;
}

int commit_parameters(){
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz){
  system_grav[index_of_the_particle].vel.x = vx;
  system_grav[index_of_the_particle].vel.y = vy;
  system_grav[index_of_the_particle].vel.z = vz;
  return 0;
}

