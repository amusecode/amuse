// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif

#include "leapfrog.hpp"
#include "mathematical_constants.h"

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


class Sui {
private:
    PS::S64 nstep;
    PS::F64 time;
    PS::S64 id_offset;
    PS::S64 N_sph_max;
    PS::F32 theta_gravity;
    PS::BOUNDARY_CONDITION boundary_condition;
    PS::F64ort pos_root_domain;
    PS::S32 nprocs;
    PS::S32 nthrds;
    //PS::ParticleSystem<FP_test> psys_test;
    PS::DomainInfo dinfo;
    PS::S64 numPtclSPH;
    PS::S64 numPtclAll;
    PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole tree_grav;
    PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
    PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;
    PS::F64 dt;

public:
    PS::ParticleSystem<FP_sph> psys_sph;
    PS::F64 epsilon_gravity;
    PS::F64 dt_max;
    Sui() {
        this->boundary_condition = PS::BOUNDARY_CONDITION_OPEN;
        this->nstep = 0;
        this->time = 0.;
        this->dt_max = 1./64.;
        this->id_offset = 0;
        this->epsilon_gravity = 1.0e-4;
        this->theta_gravity = 0.5;
        this->N_sph_max = 100000;
    }

    void set_defaults(){
        //PS::F32 theta_grav = 0.5;
        //this->epsilon_gravity = 1.0e-3;
    }

    void initialize(){
        int argc = 0;
        char **argv = NULL;
        PS::Initialize(argc, argv);
        this->psys_sph.initialize();
        //psys_test.initialize();
        this->psys_sph.setNumberOfParticleLocal(0);
        //psys_test.setNumberOfParticleLocal(0);
        this->dinfo.initialize();
        this->nprocs = PS::Comm::getNumberOfProc();
        this->nthrds = PS::Comm::getNumberOfThread();
        if (PS::Comm::getRank() == 0) {
            std::cout << "=================================" << std::endl
                      << " This is a test program of "       << std::endl
                      << " FDPS"                             << std::endl
                      << " # of processes is " << this->nprocs     << std::endl
                      << " # of thread is    " << this->nthrds     << std::endl
                      << "=================================" << std::endl;
        }
    }

    void getTimeStep(){	
        this->dt = DBL_MAX; 
        if (this->dt_max > 0.0) this->dt = this->dt_max;
     
        // Timescale for N-body system
	/*
//#if defined(ENABLE_GRAVITY_INTERACT)
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
           const PS::F64 acc = std::sqrt(psys_nbody[i].acc * psys_nbody[i].acc);
           if (acc > 0.0)
               this->dt = std::min(this->dt, CFL_dyn * std::sqrt(epsilon_gravity / acc));
        }
//#endif
        */
     
        // Timescale for SPH system
        for (PS::S32 i = 0; i < this->psys_sph.getNumberOfParticleLocal(); i++) {
//#if defined(ENABLE_GRAVITY_INTERACT)
           const PS::F64 acc = std::sqrt((psys_sph[i].acc_grav + psys_sph[i].acc_hydro)
                                        *(psys_sph[i].acc_grav + psys_sph[i].acc_hydro));
           if (acc > 0.0)
               this->dt = std::min(this->dt, CFL_dyn * std::sqrt(epsilon_gravity / acc));
//#endif
//#if defined(ENABLE_HYDRO_INTERACT)
           this->dt = std::min(this->dt, psys_sph[i].dt);
//#endif
        }
        this->dt = PS::Comm::getMinValue(this->dt);
    }

    PS::S64 get_number_of_particles_sph(){
        return id_offset;
    }

    PS::S64 add_sph_particle(
        PS::F64 mass,
        PS::F64 x,
        PS::F64 y,
        PS::F64 z,
        PS::F64 vx,
        PS::F64 vy,
        PS::F64 vz,
        PS::F64 eng,
        PS::F64 smth
        )
    {
        const PS::S64 N_sph = psys_sph.getNumberOfParticleLocal();
        FP_sph p;
        //psys_sph.setNumberOfParticleLocal(N_sph+1);
        this->id_offset += 1;
        p.id = id_offset;
        p.mass = mass;
        p.pos.x = x;
        p.pos.y = y;
        p.pos.z = z;
        //p.vel = 0.0;
        p.vel.x = vx;
        p.vel.y = vy;
        p.vel.z = vz;
        p.eng = eng;
        p.acc_grav = 0.0;
        //p.acc_grav.x = 0.;
        //p.acc_grav.y = 0.;
        //p.acc_grav.z = 0.;
        p.pot_grav = 0.0;
        p.acc_hydro = 0.0;
        //p.acc_hydro.x = 0.;
        //p.acc_hydro.y = 0.;
        //p.acc_hydro.z = 0.;
        p.flag = 0;
        p.dens = 0.0;
        p.ent = 0.0;
        p.pres = 0.0;
        p.smth = smth;
        p.gradh = 0.0;
        p.divv = 0.0;
        p.rotv = 0.0;
        //p.rotv = [0,0,0];
        p.BalSW = 0;
        p.snds = 0.0;
        p.eng_dot = 0.0;
        p.ent_dot = 0.0;
        p.dt = this->dt;
        p.vel_half = 0.0;
        //p.vel_half = [0,0,0];
        p.eng_half = 0;
        p.ent_half = 0;
        this->psys_sph.addOneParticle(p);
        //std::cout << id << " " << x << " " << y << " " << z << std::endl;
        return p.id;
    }

    PS::F64 get_mass(
        PS::S64 id
        ){
        return psys_sph[id].mass;
    }

    PS::F64vec get_position(
        PS::S64 id
        ){
        return psys_sph[id].pos;
    }

    PS::F64vec get_velocity(
        PS::S64 id
        ){
        return psys_sph[id].vel;
    }

    void commit_particles(){
        // Broadcast 
        PS::Comm::broadcast(&this->boundary_condition, 1, 0);
        PS::Comm::broadcast(&this->pos_root_domain, 1, 0);
        PS::Comm::broadcast(&this->epsilon_gravity, 1, 0);
        //PS::Comm::broadcast(&dt_dump, 1, 0);
        //PS::Comm::broadcast(&time_dump, 1, 0);
        //PS::Comm::broadcast(&time_end, 1, 0);
        PS::Comm::broadcast(&dt_max, 1, 0);

        // Set the boundary condition and the size of the computational domain if needed.
        this->dinfo.setBoundaryCondition(boundary_condition);
        if (boundary_condition != PS::BOUNDARY_CONDITION_OPEN) {
            this->dinfo.setPosRootDomain(pos_root_domain.low_,
                                         pos_root_domain.high_);
        }

        // Compute the average mass of SPH particles
        PS::F64 m_sum_loc = 0.0; 
        for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++)
            m_sum_loc += psys_sph[i].getCharge();
        mass_avg = PS::Comm::getSum(m_sum_loc) / psys_sph.getNumberOfParticleGlobal();


        // Perform domain decomposition 
        // TODO: what does 'false' mean here?
        this->dinfo.collectSampleParticle(this->psys_sph,false);
        this->dinfo.decomposeDomain();
        
        // Perform particle exchange
        this->psys_sph.exchangeParticle(this->dinfo);
    
        
        // Make tree structures
        numPtclSPH = std::max(psys_sph.getNumberOfParticleLocal(),1);
        //numPtclAll = psys_test.getNumberOfParticleLocal() + numPtclSPH;
        numPtclAll = std::max(psys_sph.getNumberOfParticleLocal(),1);
        std::cout << "N_sph: " << psys_sph.getNumberOfParticleLocal() << std::endl;

        //PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole this->tree_grav;
        this->tree_grav.initialize(3 * numPtclAll, theta_gravity);
    
        //PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather this->tree_dens;
        this->tree_dens.initialize(3 * numPtclSPH);
    
        //PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry this->tree_hydro;
        this->tree_hydro.initialize(3 * numPtclSPH);
    
#if defined(ENABLE_PHANTOM_GRAPE_X86)
        g5_open();
        g5_set_eps_to_all(epsilon_gravity);
#endif
    
        // Peform force calculations 
        //- Gravity calculations
//#if defined(ENABLE_GRAVITY_INTERACT)
        //tree_grav.setParticleLocalTree(psys_nbody);
        tree_grav.setParticleLocalTree(psys_sph,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                      CalcGravity<PS::SPJMonopole>,
                                      dinfo);
        //for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        //    psys_nbody[i].copyFromForce(tree_grav.getForce(i));
        //}
        //const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
        std::cout << "N_SPH: " << psys_sph.getNumberOfParticleLocal() << std::endl;
        for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
            //psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
            psys_sph[i].copyFromForce(tree_grav.getForce(i));
        }
//#endif

    //- SPH calculations
//#if defined(ENABLE_HYDRO_INTERACT)
		calcDensity(psys_sph, dinfo, tree_dens);
#if defined(USE_ENTROPY)
		setEntropy(psys_sph);
#endif
		setPressure(this->psys_sph);
		this->tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), this->psys_sph, this->dinfo);
//#endif

        
	// Set the timestep
	getTimeStep();
    //checkConservativeVariables(psys_nbody, psys_sph);
        std::cout << "finished commit particles" << std::endl;

    }

    void evolve_model(PS::F64 time_end){
    for (
                    this->time;
                    this->time < time_end;
                    this->time += this->dt,
                    this->nstep++
        ){
        if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << this->nstep 
                      << " dt = " << this->dt 
                      << " time = " << this->time 
                      << " time_end = " << time_end
                      << std::endl;
        }
        //std::cout << "***** a: " << psys_sph[0].acc_grav << " " << psys_sph[0].acc_hydro << std::endl;
        // Leap frog: Initial Kick & Full Drift
        //InitialKick(psys_nbody, dt);
        InitialKick(this->psys_sph, this->dt);
        //FullDrift(psys_nbody, this->dt);
        FullDrift(this->psys_sph, this->dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            //psys_nbody.adjustPositionIntoRootDomain(dinfo);
            psys_sph.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(this->psys_sph, this->dt);

        // Perform domain decomposition again
        //dinfo.collectSampleParticle(psys_nbody);
        dinfo.collectSampleParticle(psys_sph,false);
        dinfo.decomposeDomain();

        // Exchange the particles between the (MPI) processes
        //psys_nbody.exchangeParticle(dinfo);
        psys_sph.exchangeParticle(dinfo);

        // Peform force calculations
        PS::F64 t_start; 
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
//#if defined(ENABLE_GRAVITY_INTERACT)
        //tree_grav.setParticleLocalTree(psys_nbody);
        //tree_grav.setParticleLocalTree(psys_sph);
        tree_grav.setParticleLocalTree(psys_sph,true);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                      CalcGravity<PS::SPJMonopole>,
                                      dinfo);
        //for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        //    psys_nbody[i].copyFromForce(tree_grav.getForce(i));
        //}
        //const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
            //psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
            psys_sph[i].copyFromForce(tree_grav.getForce(i));
        }
//#endif
        PS::Comm::barrier();
        //if (PS::Comm::getRank() == 0) std::cout << "t_grav = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
//#if defined(ENABLE_HYDRO_INTERACT)

        calcDensity(psys_sph, dinfo, tree_dens);
        setPressure(psys_sph);
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_sph, dinfo);

//#endif
        PS::Comm::barrier();
        //if (PS::Comm::getRank() == 0) std::cout << "t_hydro = " << (PS::GetWtime() - t_start) << std::endl;

        // Get a new timestep
        //this->dt = getTimeStep(psys_nbody, psys_sph);
        getTimeStep();

        // Leap frog: Final Kick
        //FinalKick(psys_nbody, dt);
        FinalKick(psys_sph, this->dt);
        }
    }

    void finalize(){
#if defined(ENABLE_PHANTOM_GRAPE_X86)
        g5_close();
#endif
        PS::Finalize(); 
    }
};
