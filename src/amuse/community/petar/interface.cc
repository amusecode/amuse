#include "petar.hpp"
#include "interface.h"

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Interface code
 */

    static PeTar* ptr=NULL;
    static double time_start = 0.0;
#ifdef USE_SIMD
    static CalcForcePPSimd<ParticleBase,FPSoft> fcalc; // Force calculator
#else
    static CalcForcePPNoSimd<ParticleBase,FPSoft> fcalc;
#endif
    static int n_particle_in_interrupt_connected_cluster_glb; // 

    // flags
    static bool particle_list_change_flag=true;
    static bool initial_particle_flag=false;

    // common

    int initialize_code() {
        PeTar::initial_fdps_flag = true;
        ptr = new PeTar;

        int argc = 0;
        char **argv=NULL;

        n_particle_in_interrupt_connected_cluster_glb = 0;

        //No second MPI init
        //ptr->initialFDPS(argc,argv);
        //ptr->initial_fdps_flag = true;
        //ptr->my_rank= PS::Comm::getRank();
        //ptr->n_proc = PS::Comm::getNumberOfProc();

#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: Initialize_code start\n";
#endif

        // set print flag to rank 0
        ptr->input_parameters.print_flag = (ptr->my_rank==0) ? true: false;
        // set writing flag to false
        ptr->input_parameters.write_style.value = 0;

        // default input
        int flag= ptr->readParameters(argc,argv);

        // set id_offset
        ptr->input_parameters.id_offset.value = 10000000;

        // set restart flat to false
        ptr->file_header.nfile = 0; 

#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: Initialize_code end\n";
#endif
        return flag;
    }

    int cleanup_code() {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: cleanup_code start\n";
#endif
        delete ptr;
        ptr=NULL;
//#ifdef INTERFACE_DEBUG_PRINT
//        if(ptr->my_rank==0) std::cout<<"PETAR: cleanup_code end\n";
//#endif
        return 0;
    }

    int commit_parameters() {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: commit_parameters start\n";
#endif
        if (!ptr->read_parameters_flag) return -1;

        // set stopping condtions support
        set_support_for_condition(COLLISION_DETECTION);
        set_support_for_condition(PAIR_DETECTION);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        mpi_setup_stopping_conditions();
#endif
        
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: commit_parameters end\n";
#endif
        return 0;
    }

    int recommit_parameters() {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: recommit_parameters start\n";
#endif
        ptr->input_parameters.n_glb.value = ptr->stat.n_real_glb;
        ptr->input_parameters.update_changeover_flag = true;
        ptr->input_parameters.update_rsearch_flag = true;
        ptr->initialParameters();
        ptr->initial_step_flag = false;

        // set stopping condtions support
        set_support_for_condition(COLLISION_DETECTION);
        set_support_for_condition(PAIR_DETECTION);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        mpi_setup_stopping_conditions();
#endif

#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: recommit_parameters end\n";
#endif
        return 0;
    }

    // GravitationalDynamicsInterface

    int new_particle(int* index_of_the_particle,  double mass, double x, double y, double z, double vx, double vy, double vz, double radius) {
        // if not yet initial the system
#ifdef INTERFACE_DEBUG_PRINT
        if (!ptr->read_data_flag && ptr->my_rank==0) {
            std::cout<<"New particle, rank "<<ptr->my_rank<<std::endl;
            std::cout<<std::setw(20)<<"ID";
            ParticleBase::printColumnTitle(std::cout);
            std::cout<<std::endl;
        }
#endif

        ptr->read_data_flag = true;
        
        PS::S64 id_offset = ptr->input_parameters.id_offset.value;
        PS::S64 n_glb = ptr->stat.n_real_glb;
        PS::S64 n_loc = ptr->stat.n_real_loc;
#ifdef INTERFACE_DEBUG
        assert(n_loc == ptr->system_soft.getNumberOfParticleLocal());
        assert(n_glb == ptr->system_soft.getNumberOfParticleGlobal());
#endif

        if(ptr->my_rank==0) {
            FPSoft p;
            p.mass = mass;
            p.pos.x = x;
            p.pos.y = y;
            p.pos.z = z;
            p.vel.x = vx;
            p.vel.y = vy;
            p.vel.z = vz;
            p.radius = radius;

            if (ptr->initial_parameters_flag) p.calcRSearch(ptr->input_parameters.dt_soft.value);
            p.id = n_glb+1;
            if (p.id>=id_offset) return -1;
            p.group_data.artificial.setParticleTypeToSingle();

            p.rank_org = ptr->my_rank;
            p.adr = n_loc;

            if (initial_particle_flag) {
                PS::F64 m_fac = p.mass*Ptcl::mean_mass_inv;
                PS::F64 r_out = ptr->input_parameters.r_out.value;
                PS::F64 r_in = r_out*ptr->input_parameters.ratio_r_cut.value;
                p.changeover.setR(m_fac, r_in, r_out);
                p.calcRSearch(ptr->input_parameters.dt_soft.value);
            }
            
            ptr->system_soft.addOneParticle(p);

            ptr->stat.n_real_loc++;
            *index_of_the_particle = p.id;
#ifdef INTERFACE_DEBUG_PRINT
            std::cout<<std::setprecision(14);
            std::cout<<std::setw(20)<<p.id;
            p.ParticleBase::printColumn(std::cout);
            std::cout<<std::endl;
#endif

        }
        ptr->stat.n_real_glb++;

        particle_list_change_flag = true;

        return 0;
    }

    int delete_particle(int index_of_the_particle) {
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            ptr->remove_list.push_back(index);
            //ptr->stat.n_real_loc--;
#ifdef INTERFACE_DEBUG_PRINT
            std::cout<<"Remove particle index "<<index<<" id "<<index_of_the_particle<<" rank "<<ptr->my_rank<<std::endl;
#endif
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        particle_list_change_flag = true;

        return 0;
    }

    int get_state(int index_of_the_particle,
                  double * mass, 
                  double * x, double * y, double * z,
                  double * vx, double * vy, double * vz, double * radius){
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int rank_mask = index==-1 ? 0 : ptr->my_rank;
        int particle_rank = PS::Comm::getSum(rank_mask);
        if (particle_rank==0) {
            if (index==-1) return -1;
            FPSoft* p = &(ptr->system_soft[index]);
            *mass = p->mass;
            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
            *radius = p->radius;
        }
        else {
            if (ptr->my_rank==particle_rank) { // sender
                FPSoft* p = &(ptr->system_soft[index]);
                MPI_Send(p, 7, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (ptr->my_rank==0) { // receiver
                ParticleBase p;
                MPI_Recv(&p, 7, MPI_DOUBLE, particle_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                *mass = p.mass;
                *x = p.pos.x;
                *y = p.pos.y;
                *z = p.pos.z;
                *vx = p.vel.x;
                *vy = p.vel.y;
                *vz = p.vel.z;
                *radius = p.radius;
            }
        }
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *mass = p->mass;
            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
            *radius = p->radius;
        }
        else return -1;
#endif
        return 0;
    }

    int set_state(int index_of_the_particle,
                  double mass, 
                  double x, double y, double z,
                  double vx, double vy, double vz, double radius) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->mass = mass;
            p->pos.x = x;
            p->pos.y = y;
            p->pos.z = z;
            p->vel.x = vx;
            p->vel.y = vy;
            p->vel.z = vz;
            p->radius= radius;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int get_mass(int index_of_the_particle, double * mass) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        double mass_local = 0.0;
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            mass_local = p->mass;
        }    
        *mass = PS::Comm::getSum(mass_local);
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *mass = p->mass;
        }    
        else return -1;
#endif
        return 0;
    }

    int set_mass(int index_of_the_particle, double mass) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->mass = mass;
        }    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int get_radius(int index_of_the_particle, double * radius) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        double radius_local = 0.0;
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            radius_local = p->radius;
        }    
        *radius = PS::Comm::getSum(radius_local);
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *radius = p->radius;
        }    
        else return -1;
#endif
        return 0;
    }

    int set_radius(int index_of_the_particle, double radius) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->radius = radius;
        }    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int set_position(int index_of_the_particle,
                     double x, double y, double z) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->pos.x = x;
            p->pos.y = y;
            p->pos.z = z;
        }    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int get_position(int index_of_the_particle,
                     double * x, double * y, double * z) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int rank_mask = index==-1 ? 0 : ptr->my_rank;
        int particle_rank = PS::Comm::getSum(rank_mask);
        if (particle_rank==0) {
            if (index==-1) return -1;
            FPSoft* p = &(ptr->system_soft[index]);

            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
        }
        else {
            if (ptr->my_rank==particle_rank) { // sender
                FPSoft* p = &(ptr->system_soft[index]);
                double* pos = &(p->pos.x);
                MPI_Send(pos, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (ptr->my_rank==0) { // receiver
                double pos[3];
                MPI_Recv(pos, 3, MPI_DOUBLE, particle_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                *x = pos[0];
                *y = pos[1];
                *z = pos[2];
            }
        }
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
        }    
        else return -1;
#endif
        return 0;
    }

    int set_velocity(int index_of_the_particle,
                     double vx, double vy, double vz) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->vel.x = vx;
            p->vel.y = vy;
            p->vel.z = vz;
        }    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int get_velocity(int index_of_the_particle,
                     double * vx, double * vy, double * vz) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int rank_mask = index==-1 ? 0 : ptr->my_rank;
        int particle_rank = PS::Comm::getSum(rank_mask);
        if (particle_rank==0) {
            if (index==-1) return -1;
            FPSoft* p = &(ptr->system_soft[index]);

            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
        }
        else {
            if (ptr->my_rank==particle_rank) { // sender
                FPSoft* p = &(ptr->system_soft[index]);
                double* vel = &(p->vel.x);
                MPI_Send(vel, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (ptr->my_rank==0) { // receiver
                double vel[3];
                MPI_Recv(vel, 3, MPI_DOUBLE, particle_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                *vx = vel[0];
                *vy = vel[1];
                *vz = vel[2];
            }
        }
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
        }    
        else return -1;
#endif
        return 0;
    }

    int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int rank_mask = index==-1 ? 0 : ptr->my_rank;
        int particle_rank = PS::Comm::getSum(rank_mask);
        if (particle_rank==0) {
            if (index==-1) return -1;
            FPSoft* p = &(ptr->system_soft[index]);

            *ax = p->acc.x;
            *ay = p->acc.y;
            *az = p->acc.z;
        }
        else {
            if (ptr->my_rank==particle_rank) { // sender
                FPSoft* p = &(ptr->system_soft[index]);
                double* acc = &(p->acc.x);
                MPI_Send(acc, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (ptr->my_rank==0) { // receiver
                double acc[3];
                MPI_Recv(acc, 3, MPI_DOUBLE, particle_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                *ax = acc[0];
                *ay = acc[1];
                *az = acc[2];
            }
        }
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *ax = p->acc.x;
            *ay = p->acc.y;
            *az = p->acc.z;
        }    
        else return -1;
#endif
        return 0;
    }

    int set_acceleration(int index_of_the_particle, double ax, double ay, double az) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->acc.x = ax; 
            p->acc.y = ay; 
            p->acc.z = az; 
        }    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        else return -1;
#endif
        return 0;
    }

    int get_potential(int index_of_the_particle, double * potential) {
        reconstruct_particle_list();
        int index = ptr->getParticleAdrFromID(index_of_the_particle);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        double pot_local = 0.0;
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            pot_local = p->pot_tot;
        }    
        *potential = PS::Comm::getSum(pot_local);
        int check = PS::Comm::getMaxValue(index);
        if (check==-1) return -1;
#else
        if (index>=0) {
            FPSoft* p = &(ptr->system_soft[index]);
            *potential = p->pot_tot;
        }    
        else return -1;
#endif
        return 0;
    }

    int evolve_model(double time_next) {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: evolve models to "<<time_next<< "start\n";
#endif

        if (ptr->stat.n_real_glb==0) {// escape if no particle
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"PETAR: evolve models end\n";
#endif
            return 0;
        }

        // check whether interrupted cases, exist, if so, copy back data to local particles
        int n_interrupt_isolated = ptr->system_hard_isolated.getNumberOfInterruptClusters();
        for (int i=0; i<n_interrupt_isolated; i++) {
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"interrupt isolated: "<<i<<"\n";
#endif
            auto interrupt_hard_int = ptr->system_hard_isolated.getInterruptHardIntegrator(i);
            for (int k=0; k<2; k++) {
                auto pk = interrupt_hard_int->interrupt_binary.adr->getMember(k);
                
                // copy data from global particle array
                assert(pk->id==ptr->system_soft[pk->adr_org].id);
                pk->DataCopy(ptr->system_soft[pk->adr_org]);
            }
        }            

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // check connected cluster case
        if (n_particle_in_interrupt_connected_cluster_glb>0) {
            ptr->search_cluster.SendPtcl(ptr->system_soft, ptr->system_hard_connected.getPtcl());
            int n_interrupt_connected = ptr->system_hard_connected.getNumberOfInterruptClusters();
        
            for (int i=0; i<n_interrupt_connected; i++) {
#ifdef INTERFACE_DEBUG_PRINT
                if(ptr->my_rank==0) std::cout<<"interrupt connected: "<<i<<"\n";
#endif
                auto interrupt_hard_int = ptr->system_hard_connected.getInterruptHardIntegrator(i);
                for (int k=0; k<2; k++) {
                    auto pk = interrupt_hard_int->interrupt_binary.adr->getMember(k);
                    int pk_index = interrupt_hard_int->interrupt_binary.adr->getMemberIndex(k);
                    // copy data from ptcl hard or globall array
                    if (pk->adr_org>=0) {
                        assert(pk->id==ptr->system_soft[pk->adr_org].id);
                        pk->DataCopy(ptr->system_soft[pk->adr_org]);
                    }
                    else {
                        assert(pk->id==interrupt_hard_int->ptcl_origin[pk_index].id);
                        pk->DataCopy(interrupt_hard_int->ptcl_origin[pk_index]);
                    }
                }
            }
        }

        int mpi_distribute_stopping_conditions();
#endif

        if (ptr->n_interrupt_glb==0) {
            reconstruct_particle_list();
            if (!ptr->initial_step_flag) ptr->initialStep();
            ptr->input_parameters.time_end.value = time_next*2;
        }

        int is_collision_detection_enabled;
        is_stopping_condition_enabled(COLLISION_DETECTION, &is_collision_detection_enabled);
        int is_pair_detection_enabled;
        is_stopping_condition_enabled(PAIR_DETECTION, &is_pair_detection_enabled);
        if (is_collision_detection_enabled||is_pair_detection_enabled) {
            ptr->input_parameters.interrupt_detection_option.value = 1;
            ptr->hard_manager.ar_manager.interrupt_detection_option = 1;
        }
        else {
            ptr->input_parameters.interrupt_detection_option.value = 0;
            ptr->hard_manager.ar_manager.interrupt_detection_option = 0;
        }

        // record interrupt binaries in stopping condition container.
        int n_interrupt = ptr->integrateToTime(time_next);

        reset_stopping_conditions();    

        if (n_interrupt>0) {
            // isolate clusters
            int n_interrupt_isolated = ptr->system_hard_isolated.getNumberOfInterruptClusters();
            for (int i=0; i<n_interrupt_isolated; i++) {
                int stopping_index  = next_index_for_stopping_condition();
                auto interrupt_hard_int = ptr->system_hard_isolated.getInterruptHardIntegrator(i);
                auto interrupt_state = interrupt_hard_int->interrupt_binary.adr->getLeftMember()->getBinaryInterruptState();
                switch (interrupt_state) {
                case BinaryInterruptState::form:
                    set_stopping_condition_info(stopping_index, PAIR_DETECTION);
                    break;
                case BinaryInterruptState::exchange:
                    set_stopping_condition_info(stopping_index, PAIR_DETECTION);
                    break;
                case BinaryInterruptState::collision:
                    set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                    break;
                case BinaryInterruptState::none:
                    continue;
                default:
                    return -1;
                }

                for (int k=0; k<2; k++) {
                    auto pk = interrupt_hard_int->interrupt_binary.adr->getMember(k);
                    set_stopping_condition_particle_index(stopping_index, k, pk->id);
                
                    // copy back data to global particle array
                    ptr->system_soft[pk->adr_org].DataCopy(*pk);
                }
            }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            int n_particle_in_interrupt_connected_cluster=0;
            int n_interrupt_connected = ptr->system_hard_connected.getNumberOfInterruptClusters();
            for (int i=0; i<n_interrupt_connected; i++) {
                int stopping_index  = next_index_for_stopping_condition();
                auto interrupt_hard_int = ptr->system_hard_connected.getInterruptHardIntegrator(i);
                auto interrupt_state = interrupt_hard_int->interrupt_binary.adr->getLeftMember()->getBinaryInterruptState();
                switch (interrupt_state) {
                case BinaryInterruptState::form:
                    set_stopping_condition_info(stopping_index, PAIR_DETECTION);
                    break;
                case BinaryInterruptState::exchange:
                    set_stopping_condition_info(stopping_index, PAIR_DETECTION);
                    break;
                case BinaryInterruptState::collision:
                    set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                    break;
                case BinaryInterruptState::none:
                    continue;
                default:
                    return -1;
                }
                for (int k=0; k<2; k++) {
                    auto pk = interrupt_hard_int->interrupt_binary.adr->getMember(k);
                    int pk_index = interrupt_hard_int->interrupt_binary.adr->getMemberIndex(k);
                    set_stopping_condition_particle_index(stopping_index, k, pk->id);
                
                    // copy back data to global particle array
                    if (pk->adr_org>=0) {
                        assert(ptr->system_soft[pk->adr_org].id==pk->id);
                        ptr->system_soft[pk->adr_org].DataCopy(*pk);
                    }
                    else {
                        // if particle is in remote node, copy back to ptcl_hard and wait for MPI_send/recv
                        n_particle_in_interrupt_connected_cluster++;
                        assert(pk->id == interrupt_hard_int->ptcl_origin[pk_index].id);
                        interrupt_hard_int->ptcl_origin[pk_index].DataCopy(*pk);
                    }
                }
            }
            // if particle in remote node need update, call MPI send/recv
            n_particle_in_interrupt_connected_cluster_glb = PS::Comm::getSum(n_particle_in_interrupt_connected_cluster);
            if (n_particle_in_interrupt_connected_cluster_glb>0)
                ptr->search_cluster.writeAndSendBackPtcl(ptr->system_soft, ptr->system_hard_connected.getPtcl(), ptr->remove_list);

            mpi_collect_stopping_conditions();
#endif      
        }
//#ifdef PROFILE
        //ptr->printProfile();
        //ptr->clearProfile();
//#endif        
        ptr->reconstructIdAdrMap();
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: evolve models end\n";
#endif
        return 0;
    }

    int commit_particles() {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: commit_particles start\n";
#endif
        
        if (!ptr->read_parameters_flag) return -1;
        if (!ptr->read_data_flag) return -1;
        ptr->input_parameters.n_glb.value = ptr->stat.n_real_glb;
        ptr->initialParameters();
        ptr->initialStep();
        ptr->reconstructIdAdrMap();
        particle_list_change_flag = false;
        initial_particle_flag = true;
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: commit_particles end\n";
#endif
        return 0;
    }

    int synchronize_model() {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: synchronize_model\n";
#endif
        return 0;
    }

    int reconstruct_particle_list() {
        if (particle_list_change_flag) {
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"PETAR: reconstruct particle list start\n";
#endif
            ptr->removeParticles();
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"PETAR: remove particles end\n";
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            ptr->domainDecompose(true);
            ptr->exchangeParticle();
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"PETAR: exchange particles end\n";
#endif
#endif
            ptr->reconstructIdAdrMap();
#ifdef INTERFACE_DEBUG_PRINT
            if(ptr->my_rank==0) std::cout<<"PETAR: reconstruct particle list end\n";
#endif
            particle_list_change_flag = false;
        }
        return 0;
    }

    int recommit_particles() {
        // this function is called too frequent (every time when set_xx is used).
        // thus only register the flag and do update at once in the begining of evolve_model
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: recommit_particles start\n";
#endif
        if (ptr->n_interrupt_glb==0) ptr->initial_step_flag = false;
        reconstruct_particle_list();
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: recommit_particles end\n";
#endif
        return 0;
    }

    int get_eps2(double * epsilon_squared) {
        *epsilon_squared = ptr->input_parameters.eps.value ;
        return 0;
    }

    int set_eps2(double epsilon_squared) {
        ptr->input_parameters.eps.value = epsilon_squared;
        return 0;
    }

    // set changeover radius outer boundary (if zero, auto-determine)
    int set_changeover_rout(double r_out) {
        ptr->input_parameters.r_out.value = r_out;
        return 0;
    }

    int get_changeover_rout(double* r_out) {
        *r_out = ptr->input_parameters.r_out.value;
        return 0;
    }

    // set changeover ratio of inner / outer boundary (if zero, auto-determine)
    int set_changeover_ratio(double ratio_r_cut) {
        ptr->input_parameters.ratio_r_cut.value = ratio_r_cut;
        return 0;
    }

    int get_changeover_ratio(double* ratio_r_cut) {
        *ratio_r_cut = ptr->input_parameters.ratio_r_cut.value;
        return 0;
    }

    // set group detection maximum radius to switch on AR (if zero, auto-determine)
    int set_group_radius(double r_bin) {
        ptr->input_parameters.r_bin.value = r_bin;
        return 0;
    }

    int get_group_radius(double* r_bin) {
        *r_bin = ptr->input_parameters.r_bin.value;
        return 0;
    }

    // set neighbor search radius minimum (if zero, auto-determine)
    int set_rsearch_min(double r_search_min) {
        ptr->input_parameters.r_search_min.value = r_search_min;
        return 0;
    }

    int get_rsearch_min(double* r_search_min) {
        *r_search_min = ptr->input_parameters.r_search_min.value;
        return 0;
    }

    // set tree opening angle
    int set_theta(double theta) {
        ptr->input_parameters.theta.value = theta;
        return 0;
    }

    int get_theta(double* theta) {
        *theta = ptr->input_parameters.theta.value;
        return 0;
    }

    // set tree time step
    int set_tree_step(double dt_soft) {
        ptr->input_parameters.dt_soft.value = dt_soft;
        return 0;
    }

    int get_tree_step(double* dt_soft) {
        *dt_soft = ptr->input_parameters.dt_soft.value;
        return 0;
    }

    int set_output_step(double dt_snap) {
        ptr->input_parameters.dt_snap.value = dt_snap;
        return 0;
    }

    int get_output_step(double* dt_snap) {
        *dt_snap = ptr->input_parameters.dt_snap.value;
        return 0;
    }


    //// set gravitational constant
    //int set_gravitational_constant(double G) {
    //    ptr->input_parameters.unit_set.value=-1;
    //    ptr->input_parameters.gravitational_constant.value = G;
    //    return 0;
    //}
    // 
    //// get gravitational constant
    //int get_gravitational_constant(double* G) {
    //    *G = ptr->input_parameters.gravitational_constant.value;
    //    return 0;
    //}

    int get_kinetic_energy(double * kinetic_energy) {
        // update particle array first if necessary
        reconstruct_particle_list();
        if (!ptr->initial_step_flag) ptr->initialStep();
        *kinetic_energy = ptr->stat.energy.ekin;
        return 0;
    }

    int get_potential_energy(double * potential_energy) {
        // update particle array first if necessary
        reconstruct_particle_list();
        if (!ptr->initial_step_flag) ptr->initialStep();
        *potential_energy = ptr->stat.energy.epot;
        return 0;
    }

    int get_time(double * time) {
        * time = ptr->stat.time;
        return 0;
    }

    int get_begin_time(double * time) {
        * time = time_start;
        return 0;
    }

    int set_begin_time(double time) {
#ifdef INTERFACE_DEBUG_PRINT
        if(ptr->my_rank==0) std::cout<<"PETAR: set begin time from "<<ptr->stat.time<<" to "<<time<<std::endl;
#endif 
        time_start = time;
        ptr->stat.time = time_start;
        return 0;
    }

    int get_time_step(double * time_step) {
        * time_step = ptr->input_parameters.dt_soft.value;
        return 0;
    }

    int get_total_mass(double * mass) {
        // update particle array first if necessary
        reconstruct_particle_list();
        * mass = ptr->stat.pcm.mass;
        return 0;
    }

    int get_center_of_mass_position(double * x, double * y, double * z) {
        // update particle array first if necessary
        reconstruct_particle_list();
        * x = ptr->stat.pcm.pos.x;
        * y = ptr->stat.pcm.pos.y;
        * z = ptr->stat.pcm.pos.z;
        return 0;
    }

    int get_center_of_mass_velocity(double * x, double * y, double * z) {
        // update particle array first if necessary
        reconstruct_particle_list();
        * x = ptr->stat.pcm.pos.x;
        * y = ptr->stat.pcm.pos.y;
        * z = ptr->stat.pcm.pos.z;
        return 0;
    }

    int get_total_radius(double * radius) {
        // update particle array first if necessary
        reconstruct_particle_list();
        * radius = ptr->stat.half_mass_radius;
        return 0;
    }

    int get_number_of_particles(int * number_of_particles) {
        // update particle array first if necessary
        reconstruct_particle_list();
        * number_of_particles = ptr->stat.n_real_glb;
        return 0;
    }

    int get_index_of_first_particle(int * index_of_the_particle) {
        * index_of_the_particle = 1;
        return 0;
    }

    int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle) {
        * index_of_the_next_particle = index_of_the_particle + 1;
        return 0;
    }

    // GravityFieldInterface

    int get_gravity_at_point(double * eps, double * x, double * y, double * z, 
                             double * forcex, double * forcey, double * forcez, int n)  {
        // update particle array first if necessary
        reconstruct_particle_list();

        // transform data
        ParticleBase ptmp[n];
        ForceSoft force[n];
        for (int i=0; i<n; i++) {
            ptmp[i].pos.x = x[i];
            ptmp[i].pos.y = y[i];
            ptmp[i].pos.z = z[i];
            force[i].acc.x = 0;
            force[i].acc.y = 0;
            force[i].acc.z = 0;
        };

        
        fcalc(ptmp, n, &(ptr->system_soft[0]), ptr->system_soft.getNumberOfParticleLocal(), force);

        for (int i=0; i<n; i++) {
            forcex[i] = force[i].acc.x;
            forcey[i] = force[i].acc.y;
            forcez[i] = force[i].acc.z;
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (ptr->my_rank==0) {
            MPI_Reduce(MPI_IN_PLACE, forcex, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, forcey, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(MPI_IN_PLACE, forcez, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        else {
            MPI_Reduce(forcex, NULL, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(forcey, NULL, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(forcez, NULL, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
#endif
        return 0;
    }


    int get_potential_at_point(double * eps,
                               double * x, double * y, double * z, 
                               double * phi, int n)  {
        // update particle array first if necessary
        reconstruct_particle_list();

        // transform data
        ParticleBase ptmp[n];
        ForceSoft force[n];
        for (int i=0; i<n; i++) {
            ptmp[i].pos.x = x[i];
            ptmp[i].pos.y = y[i];
            ptmp[i].pos.z = z[i];
            force[i].pot = 0;
        };

        fcalc(ptmp, n, &(ptr->system_soft[0]), ptr->system_soft.getNumberOfParticleLocal(), force);

        for (int i=0; i<n; i++) phi[i] = force[i].pot;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (ptr->my_rank==0) MPI_Reduce(MPI_IN_PLACE, phi,  n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        else                 MPI_Reduce(phi,          NULL, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        return 0;
    }    

#ifdef __cplusplus
}
#endif
