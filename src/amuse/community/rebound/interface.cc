#include "worker_code.h"

extern "C" {
#define restrict 

#include "rebound.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "boundary.h"
#include "gravity.h"
//#include "output.h"
#include "particle.h"
#include "collision.h"
//#include "simulationarchive.h"
}

#ifdef OPENMP_ENABLED
#include <openmp.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include <stopcond.h>

typedef struct _particle_location {
    struct reb_simulation * code; //= NULL;
    struct reb_particle * p;// = NULL;
} particle_location;

typedef struct _code_state {
    reb_simulation * code;
    bool has_removal;
    bool has_unsorted_massless_particles;
    double time_offset;
    int subset;
    _code_state(reb_simulation * c):code(c),has_removal(false), has_unsorted_massless_particles(false), time_offset(0), subset(0) {}
    _code_state(reb_simulation * c, double time_offset, int s):code(c),has_removal(false), has_unsorted_massless_particles(false), time_offset(time_offset),subset(s) {}
} code_state;


typedef struct _particle_sort {
    int ref_index;
    reb_particle* p;
    _particle_sort(int i, reb_particle* p):ref_index(i),p(p){}
} particle_sort;
typedef std::vector<_particle_sort>  ParticleSortVector;

static int max_id = 0;

typedef std::vector<code_state> ReboundSimulationVector;
static ReboundSimulationVector codes;
static double _time;
static double timestep = 0.0001;

static inline particle_location get_particle_from_identity(int index_of_the_particle)
{
    particle_location particle = {NULL,NULL};

    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        particle.code = cs.code;
        particle.p = reb_get_particle_by_hash(particle.code, index_of_the_particle);
        if (particle.p != NULL) break;
        //*i = cs;
    }
    return particle;
}

int get_mass(int index_of_the_particle, double * mass){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    *mass = p->m;
    return 0;
}


bool sort_particles (particle_sort i,particle_sort j) { 
    return (i.p->m > j.p->m); 
}

int commit_particles(){
    return 0;
}

int get_time(int code_index, double * time){
    if(code_index < 0) {
        *time = _time;
        return 0;
    }
    if(code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    *time = codes[code_index].code->t;
    return 0;
}

int set_mass(int index_of_the_particle, double mass){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}

    if(p->m==0){
        if(mass>0){
            int index_old=reb_get_particle_index(p);
            if(index_old!=code->N_active){
                struct reb_particle tmp = code->particles[index_old];
                for(int j=index_old; j>code->N_active; j--){
                    code->particles[j] = code->particles[j-1];
                }
                code->particles[code->N_active] = tmp;
            }
            code->N_active++;
        }
    }
    else {
        if(mass==0){
            int index_old=reb_get_particle_index(p);
            code->N_active--;

            if(index_old!=code->N_active){
                struct reb_particle tmp = code->particles[index_old];
                for( int j = index_old; j<code->N_active; j++){
                    code->particles[j] = code->particles[j+1];
                }
                code->particles[code->N_active] = tmp;
            }
        }
    }
    p->m = mass;
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
    return 0;
}

int get_total_radius(double * radius){
    return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius, int code_index){
      if(code_index < 0 || code_index >= (signed) codes.size()){
        *index_of_the_particle=0;
        return -10;
      }
      uint32_t new_hash = max_id++;
      struct reb_particle pt;
      pt.x = x;
      pt.y = y;
      pt.z = z;
      pt.vx = vx;
      pt.vy = vy;
      pt.vz = vz;
      pt.m = mass;
      pt.r = radius; 
      pt.hash = new_hash;
      reb_add(codes[code_index].code, pt);
      //std::cout<<"new particle :"<<pt.id<< " << "<<code_index<<" << "<<pt.x<<std::endl;
      *index_of_the_particle = new_hash;

      //make sure massless particles are last and that N_active is equal to massive particles
      if(pt.m != 0.0) {
          int N_active = codes[code_index].code->N_active;
          for(int j=codes[code_index].code->N-1;j>N_active;j--){
              codes[code_index].code->particles[j] = codes[code_index].code->particles[j-1];
          }
          codes[code_index].code->particles[N_active] = pt;
          codes[code_index].code->N_active++;
      }
      return 0;
}

int get_total_mass(double * mass){
    return 0;
}
int _evolve_code(double _tmax, code_state * cs);
int evolve_model(double _tmax, int code_index){
    int result = 0;
    
    reset_stopping_conditions();
    
    if(code_index == -1){
        for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
            code_state cs = *i;
            if(cs.code) {
                _evolve_code(_tmax, &cs);
            }
        }
        _time = _tmax;
        return result;
    }else if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    } else {
        result = _evolve_code(_tmax, &codes[code_index]);
        return result;
    }   
}


int _evolve_code(double _tmax, code_state * cs){
    reb_simulation * code = cs->code;
    if(!code) {return -1;}
    
    int error = 0;
    int is_collision_detection_enabled = 0;
    int is_timeout_detection_enabled = 0;
    int is_out_of_box_detection_enabled = 0;
    double center_of_mass[3];
    double box_size_squared = out_of_box_parameter*out_of_box_parameter;    
  
 
    error = is_stopping_condition_enabled(
        COLLISION_DETECTION, 
        &is_collision_detection_enabled
    );
    error = is_stopping_condition_enabled(
        TIMEOUT_DETECTION, 
        &is_timeout_detection_enabled
    );
    error = is_stopping_condition_enabled(
        OUT_OF_BOX_DETECTION,
        &is_out_of_box_detection_enabled
    );
    int is_condition_set = 0;
    
    
    // original : rebound_integrate
    int exact_finish_time = 1;
    double maxR = 0;
    //double minD = 0;
    //struct timeval tim;
    //gettimeofday(&tim, NULL);
    //double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
    double tmax = _tmax;
    code-> dt_last_done = code->dt;
    int last_step = 0;
    int ret_value = 0;
    double ke = 0.0, ke1 = 0.0;
    const double dtsign = copysign(1.,code->dt);                // Used to determine integration direction
    
    time_t starttime, currenttime;
    time(&starttime);
    get_kinetic_energy(cs->subset, &ke);
    //printf("Code time: %d ,  %f -> %f (%f)\n",cs->subset , code->t, tmax, ke);
    while(code->t*dtsign<tmax*dtsign && last_step<2 && ret_value==0){
        if (code->N<=0){
            fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
            return(1);
        }
       
        reb_step(code);                                 // 0 to not do timing within step 
        
        if ((code->t+code->dt)*dtsign>=tmax*dtsign && exact_finish_time==1){
            reb_integrator_synchronize(code);
            code->dt = tmax-code->t;
            last_step++;
        }else{
            code->dt_last_done = code->dt;
        }
        if (maxR){
            // Check for escaping particles
            const double maxR2 = maxR*maxR;
            for (int i=0;i<code->N-code->N_var;i++){
                struct reb_particle p = code->particles[i];
                double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
                if (r2>maxR2){
                    ret_value = 2;
                }
            }
        }
        if (is_collision_detection_enabled){
            // Check for close encounters
            for (int i=0;i<code->N-code->N_var;i++){
                struct reb_particle pi = code->particles[i];
                for (int j=0;j<i;j++){
                    struct reb_particle pj = code->particles[j];
                    const double x = pi.x-pj.x;
                    const double y = pi.y-pj.y;
                    const double z = pi.z-pj.z;
                    const double r2 = x*x + y*y + z*z;
                    
                    const double rsum = pi.r+pj.r;
                    if (r2<(rsum*rsum)){
                        int stopping_index  = next_index_for_stopping_condition();
                        if(stopping_index < 0)
                        {
                            
                        }
                        else
                        {
                            set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                            set_stopping_condition_particle_index(stopping_index, 0, pi.hash);
                            set_stopping_condition_particle_index(stopping_index, 1, pj.hash);
                        }
                        is_condition_set = 1;
                    }
                }
            }
            if(is_condition_set) {
                break;
            }
        }
        // AMUSE STOPPING CONDITIONS
        if (is_out_of_box_detection_enabled) {
            int i,k;
            for (k = 0; k < 3; k++) {
                center_of_mass[k] = 0.0;
            }
            if(use_center_of_mass_parameter) {
                
                double total = 0.0;
                for (i = 0; i < code->N-code->N_var; i++) {
                    double mass = code->particles[i].m;
                    center_of_mass[0] +=  mass * code->particles[i].x;
                    center_of_mass[1] +=  mass * code->particles[i].y;
                    center_of_mass[2] +=  mass * code->particles[i].z;
                    
                    total += mass;
                }
                for (k = 0; k < 3; k++) {
                    center_of_mass[k] /= total;
                }
            } 
            for (i = 0; i < code->N-code->N_var; i++) {
               
                double dx = (code->particles[i].x - center_of_mass[0]);
                double dy = (code->particles[i].y - center_of_mass[1]);
                double dz = (code->particles[i].z - center_of_mass[2]);
                double sqr_distance_wrt_origin = 0.0;
                sqr_distance_wrt_origin += dx*dx;
                sqr_distance_wrt_origin += dy*dy;
                sqr_distance_wrt_origin += dz*dz;
                
                if (sqr_distance_wrt_origin > box_size_squared) {
                    is_condition_set = 1;
                    int stopping_index = next_index_for_stopping_condition();
                    if(stopping_index >= 0){
                        set_stopping_condition_info(stopping_index, OUT_OF_BOX_DETECTION);
                        set_stopping_condition_particle_index(stopping_index, 0, code->particles[i].hash);
                    } else {
                        printf("Run out of storable out of box events\n");
                    }
                }
            }
            if(is_condition_set) {
                break;
            }
        }
        // AMUSE STOPPING CONDITIONS
        if(is_timeout_detection_enabled) {
            time(&currenttime);
            if((currenttime - starttime) > timeout_parameter) {
                int stopping_index  = next_index_for_stopping_condition();
                set_stopping_condition_info(stopping_index, TIMEOUT_DETECTION);
                break;
            }
        }
    }
    reb_integrator_synchronize(code);
    code->dt = code->dt_last_done;
    get_kinetic_energy(cs->subset, &ke1);
    //printf("Code time: %d ,  %f -> %f (%f,%f)\n",cs->subset , code->t, tmax, ke1, (ke1-ke)/ke);
    //gettimeofday(&tim, NULL);
    //double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
    //double timing = timing_final-timing_initial;
    return ret_value;
}

int set_eps2(double epsilon_squared, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;

    code->softening = sqrt(epsilon_squared);
    return 0;
}

int get_begin_time(double * time){
    return 0;
}

int get_eps2(int code_index, double * epsilon_squared){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;

    *epsilon_squared = code->softening * code->softening;
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
    return 0;
}

int delete_particle(int index_of_the_particle){
    return 0;    
}

int _delete_particle(int index_of_the_particle, int code_index){
    int keepSorted = 1;
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    reb_remove_by_hash(codes[code_index].code, index_of_the_particle, keepSorted);
    return 0;
}

int get_potential(int index_of_the_particle, double * potential){
    return 0;
}

int synchronize_model(){
    return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
    double z, double vx, double vy, double vz, double radius){
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    
    p->x = x;
    p->y = y;
    p->z = z;
    p->vx = vx;
    p->vy = vy;
    p->vz = vz;
    
    p->m = mass;
    p->r = radius;
    if (code->integrator == reb_simulation::REB_INTEGRATOR_JANUS){
        code->ri_janus.recalculate_integer_coordinates_this_timestep = 1;
    }
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
    double * y, double * z, double * vx, double * vy, double * vz, 
    double * radius, int * subset){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    
    *x    = p->x;
    *y    = p->y;
    *z    = p->z;
    *vx   = p->vx;
    *vy   = p->vy;
    *vz   = p->vz;
    *mass = p->m;

#ifndef COLLISIONS_NONE
    *radius = p->r;
#else
    *radius = 0;
#endif // COLLISIONS_NONE
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        p = reb_get_particle_by_hash(cs.code, index_of_the_particle);
        if (p != NULL) {
            *subset = cs.subset;
            break;
        }
    }
    return 0;
}

int get_time_step(int code_index, double * time_step){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(code_index == 0) {
        *time_step = timestep;
        return 0;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    *time_step = code->dt;
    return 0;
}

int set_time_step(double time_step, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(code_index == 0) {
        timestep = time_step;
    } else if(!codes[code_index].code) {
        return -11;
    }
    if(codes[code_index].code) {
        reb_simulation * code = codes[code_index].code;
        code->dt = time_step;
    }
    return 0;
}

int recommit_particles(){
    return commit_particles();
}

int get_kinetic_energy(int code_index, double * kinetic_energy){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    double e_kin = 0.;
    for (int i=0;i<code->N-code->N_var;i++){
        struct reb_particle pi = code->particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    *kinetic_energy = e_kin;
    return 0;
}

int get_number_of_particles(int * number_of_particles){
    int N = 0;
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        N += cs.code->N;
    }

    *number_of_particles = N;
    return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    p->ax = ax;
    p->ay = ay;
    p->az = az;
    if (code->integrator == reb_simulation::REB_INTEGRATOR_JANUS){
        code->ri_janus.recalculate_integer_coordinates_this_timestep = 1;
    }
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    return 0;
}

int get_radius(int index_of_the_particle, double * radius){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    //if(loc.index < 0) {*radius = 0; return -1;}
    *radius = p->r;
    return 0;
}

int get_subset(int index_of_the_particle, int * subset){
    //FIXME
    struct reb_particle* p=NULL;
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        p = reb_get_particle_by_hash(cs.code, index_of_the_particle);
        if (p != NULL) {
            *subset = cs.subset;
            break;
        }
    }
    if(p == NULL) {return -1;}
    //if(loc.index < 0) {*subset = -2; return -1;}
    return 0;
}

int set_begin_time(double time){
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    p->r = radius;
    return 0;
}

int set_subset(int index_of_the_particle, int subset){
    
    struct reb_particle* p=NULL;
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        p = reb_get_particle_by_hash(cs.code, index_of_the_particle);
        if (p != NULL) {
            if(cs.subset != subset) {return -2;}
            break;
        }
    }
    if(p == NULL) {return -1;}

    return 0;
}

int cleanup_code() {
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        if(cs.code){
            reb_remove_all(cs.code);
            reb_free_simulation(cs.code);
            cs.code = 0;
            *i = cs;
        }
    }
    codes.clear();
    max_id = 0;
    timestep = 0.0001;
    _time=0;
    return 0;
}

int recommit_parameters(){
    return 0;
}

int initialize_code(){
    initialize_stopping_conditions();
    max_id = 0;
    timestep = 0.0001;
    _time=0;
#ifdef OPENMP_ENABLED
    int nt = omp_get_max_threads();
    omp_set_num_threads(nt);
#endif
    reb_simulation * code = reb_create_simulation();
    codes.push_back(code_state(code));
    code->integrator = reb_simulation::REB_INTEGRATOR_IAS15;
    code->N_active = 0;
    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);
    set_support_for_condition(TIMEOUT_DETECTION);
    //set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    set_support_for_condition(OUT_OF_BOX_DETECTION);
    
    return 0;
}

int get_potential_energy(int code_index, double * potential_energy){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    double e_pot = 0.;
    for (int i=0;i<code->N-code->N_var;i++){
        struct reb_particle pi = code->particles[i];
        for (int j=i+1;j<code->N-code->N_var;j++){
            struct reb_particle pj = code->particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= code->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    *potential_energy = e_pot;
    return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
    double * vz){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    *vx = p->vx;
    *vy = p->vy;
    *vz = p->vz;
    
    return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
      double * z){
    
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    *x = p->x;
    *y = p->y;
    *z = p->z;
    return 0;
}


int set_position(int index_of_the_particle, double x, double y, double z){
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    p->x = x;
    p->y = y;
    p->z = z;
    if (code->integrator == reb_simulation::REB_INTEGRATOR_JANUS){
        code->ri_janus.recalculate_integer_coordinates_this_timestep = 1;
    }
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
      double * az){
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    //reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    *ax = p->ax;
    *ay = p->ay;
    *az = p->az;
    return 0;
}

int commit_parameters(){
    return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
    double vz){
    particle_location particle = get_particle_from_identity(index_of_the_particle);
    struct reb_particle* p = particle.p;
    reb_simulation * code = particle.code;
    if(p == NULL) {return -1;}
    p->vx = vx;
    p->vy = vy;
    p->vz = vz;
    if (code->integrator == reb_simulation::REB_INTEGRATOR_JANUS){
        code->ri_janus.recalculate_integer_coordinates_this_timestep = 1;
    }
    return 0;
}

int new_subset(int * index, double time_offset) {
    reb_simulation * code = reb_create_simulation();
    reb_integrator_reset(code);
    code->dt = timestep;
    if(time_offset < 0) {time_offset = _time;}
    code->integrator = reb_simulation::REB_INTEGRATOR_IAS15;
    code->N_active = 0;
    code->t = time_offset;
    codes.push_back(code_state(code, time_offset, codes.size()));
    *index = codes.size() - 1;
    //printf("Code time: %d ,  %f\n",*index , code->t);
    return 0;
}
int stop_subset(int code_index) {
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    code_state cs = codes[code_index];
    if(cs.code) {
        reb_simulation * code = cs.code;
        reb_remove_all(code);
        reb_free_simulation(code);
        cs.code = 0;
        codes[code_index] = cs;
    }
    return 0;
}
int _set_integrator(int value, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    switch(value){
        case 0:
             code->integrator = reb_simulation::REB_INTEGRATOR_IAS15;
            break;
        case 1:
            code->integrator = reb_simulation::REB_INTEGRATOR_WHFAST;
            break;
        case 2:
            code->integrator = reb_simulation::REB_INTEGRATOR_SEI;
            break;
        case 3:
            // This integrator was removed
            return -1;
            break;
        case 4:
            code->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
            break;
        case 5:
            code->integrator = reb_simulation::REB_INTEGRATOR_HERMES;
            break;
        case 6:
            code->integrator = reb_simulation::REB_INTEGRATOR_WHFASTHELIO;
            break;
        case 7:
            code->integrator = reb_simulation::REB_INTEGRATOR_NONE;
            break;
        case 8:
            code->integrator = reb_simulation::REB_INTEGRATOR_JANUS;
            break;
        default:
            code->integrator = reb_simulation::REB_INTEGRATOR_NONE;
            return -1;
            break;
            
    }
    return 0;
}

int _get_integrator(int code_index, int * value){
    
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    *value = code->integrator; 
    return 0;
}

int _set_solver(int value, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    switch(value){
        case 0:
            code->gravity = reb_simulation::REB_GRAVITY_NONE;
            break;
        case 1:
            code->gravity = reb_simulation::REB_GRAVITY_BASIC;
            break;
        case 2:
            code->gravity = reb_simulation::REB_GRAVITY_COMPENSATED;
            break;
        case 3:
            code->gravity = reb_simulation::REB_GRAVITY_TREE;
            break;
        default:
            code->gravity = reb_simulation::REB_GRAVITY_NONE;
            return -1;
            break;
            
    }
    return 0;
}

int _get_solver(int code_index, int * value){
    
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    *value = code->gravity; 
    return 0;
}

int get_opening_angle2(int code_index, double * opening_angle2){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    *opening_angle2 = code->opening_angle2;
    return 0;
}

int set_opening_angle2(double opening_angle2, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    code->opening_angle2 = opening_angle2;
    return 0;
}

int _set_boundary(int value, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    switch(value){
        case 0:
            code->boundary = reb_simulation::REB_BOUNDARY_NONE;
            break;
        case 1:
            code->boundary = reb_simulation::REB_BOUNDARY_OPEN;
            break;
        case 2:
            code->boundary = reb_simulation::REB_BOUNDARY_PERIODIC;
            break;
        case 3:
            code->boundary = reb_simulation::REB_BOUNDARY_SHEAR;
            break;
        default:
            code->boundary = reb_simulation::REB_BOUNDARY_NONE;
            return -1;
            break;
            
    }
    return 0;
}

int _get_boundary(int code_index, int * value){
    
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    *value = code->boundary; 
    return 0;
}

int get_boundary_size(int code_index, double * boundary_size){
        if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    
    *boundary_size = code->root_size; 
    return 0;
}

int set_boundary_size(double boundary_size, int code_index){
    if(code_index < 0 || code_index >= (signed) codes.size()){
        return -10;
    }
    if(!codes[code_index].code) {
        return -11;
    }
    reb_simulation * code = codes[code_index].code;
    reb_configure_box(code,boundary_size,1,1,1);
    return 0;
}

