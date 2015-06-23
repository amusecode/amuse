#include "worker_code.h"

extern "C" {

#include "integrator.h"
#include "integrator_whfast.h"
#include "boundaries.h"
#include "gravity.h"
#include "output.h"
#include "particle.h"
#include "librebound.h"

//missing from librebound.h
void reset(void); 
int rebound_integrate(double _tmax, int exact_finish_time, double maxR, double minD);
//missing from librebound.h

}

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
#include <iostream>

#include <stopcond.h>


typedef std::map<int, int> IndexMap;
static IndexMap indexMap;

static int max_id = 0;
static bool has_removal = false;

static inline int get_index_from_identity(int id)
{
    IndexMap::iterator i = indexMap.find(id);
    if(i == indexMap.end()) {
        return -1;
    }
    return (*i).second;
}

static inline int get_identity_from_index(int id)
{
    for( IndexMap::iterator i = indexMap.begin(); i != indexMap.end(); i++) {
        if(indexMap [(*i).first] == id) {
            return (*i).first;
        }
    }
    return -1;
}

int get_mass(int index_of_the_particle, double * mass){
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {*mass = 0; return -1;}
    *mass = particles[index].m;
    return 0;
}

int commit_particles(){
    if(has_removal) {
        /*
         * after deletion of one or more particles, clear out
         * the original data and rebuild without these particles.
         */
        struct particle * previous = (struct particle *) malloc(sizeof(struct particle) * N);
        memcpy(previous, particles, sizeof(struct particle) * N);
        particles_remove_all();
        IndexMap newMap;
        for( IndexMap::iterator i = indexMap.begin(); i != indexMap.end(); i++) {
            struct particle * p = previous + (*i).second;
            particles_add(*p);
            newMap[(*i).first] = N - 1;
        }
        indexMap = newMap;
        has_removal = false;
    }
    
    return 0;
}

int get_time(double * time){
    return 0;
}

int set_mass(int index_of_the_particle, double mass){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    particles[index].m = mass;
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
    return 0;
}

int get_total_radius(double * radius){
    return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
      int new_id = max_id++;
      struct particle pt;
      pt.x = x;
      pt.y = y;
      pt.z = z;
      pt.vx = vx;
      pt.vy = vy;
      pt.vz = vz;
      pt.m = mass;
      pt.ax = 0;
      pt.ay = 0;
#ifndef COLLISIONS_NONE
      pt.r = radius; 
      pt.lastcollision = 0;
#endif // COLLISIONS_NONE
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
      pt.cell = NULL;
#endif // TREE
      particles_add(pt);
      indexMap[new_id] = N - 1;
      *index_of_the_particle = new_id;
      return 0;
}

int get_total_mass(double * mass){
    return 0;
}

int evolve_model(double _tmax){
    int error = 0;
    int is_collision_detection_enabled = 0;
    error = is_stopping_condition_enabled(
                COLLISION_DETECTION, 
                &is_collision_detection_enabled
    );
    int is_condition_set = 0;
    
    // original : rebound_integrate
    int exact_finish_time = 1;
    double maxR = 0;
    double minD = 0;
    struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	tmax = _tmax;
	double dt_last_done = dt;
	int last_step = 0;
	int ret_value = 0;
	const double dtsign = copysign(1.,dt); 				// Used to determine integration direction
	while(t*dtsign<tmax*dtsign && last_step<2 && ret_value==0){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return(1);
		}
		rebound_step(0); 								// 0 to not do timing within step
		if ((t+dt)*dtsign>=tmax*dtsign && exact_finish_time==1){
			integrator_synchronize();
			dt = tmax-t;
			last_step++;
		}else{
			dt_last_done = dt;
		}
		if (maxR){
			// Check for escaping particles
			const double maxR2 = maxR*maxR;
			for (int i=0;i<N-N_megno;i++){
				struct particle p = particles[i];
				double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
				if (r2>maxR2){
					ret_value = 2;
				}
			}
		}
		if (is_collision_detection_enabled){
			// Check for close encounters
			for (int i=0;i<N-N_megno;i++){
				struct particle pi = particles[i];
				for (int j=0;j<i;j++){
					struct particle pj = particles[j];
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
                            set_stopping_condition_particle_index(stopping_index, 0, get_identity_from_index(i));
                            set_stopping_condition_particle_index(stopping_index, 1, get_identity_from_index(j));
                        }
                        is_condition_set = 1;
					}
				}
			}
            if(is_condition_set) {
                break;
            }
		}
	}
	integrator_synchronize();
	dt = dt_last_done;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	timing = timing_final-timing_initial;
	return ret_value;
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

int delete_particle(int index_of_the_particle){
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    
    indexMap.erase(index_of_the_particle);
    
    has_removal = true;
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
     int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    
    
    particles[index].x = x;
    particles[index].y = y;
    particles[index].z = z;
    particles[index].vx = vx;
    particles[index].vy = vy;
    particles[index].vz = vz;
    
    particles[index].m = mass;
    
#ifndef COLLISIONS_NONE
    particles[index].r = radius;
#endif // COLLISIONS_NONE
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
    double * y, double * z, double * vx, double * vy, double * vz, 
    double * radius){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    
    *x = particles[index].x;
    *y = particles[index].y;
    *z = particles[index].z;
    *vx = particles[index].vx;
    *vy = particles[index].vy;
    *vz = particles[index].vz;
    *mass = particles[index].m;
    
#ifndef COLLISIONS_NONE
    *radius = particles[index].r;
#else
    *radius = 0;
#endif // COLLISIONS_NONE
    return 0;
}

int get_time_step(double * time_step){
  return 0;
}

int recommit_particles(){
    return commit_particles();
}

int get_kinetic_energy(double * kinetic_energy){
    return 0;
}

int get_number_of_particles(int * number_of_particles){
    *number_of_particles = indexMap.size();
    return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    particles[index].ax = ax;
    particles[index].ay = ay;
    particles[index].az = az;
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    return 0;
}

int get_radius(int index_of_the_particle, double * radius){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {*radius = 0; return -1;}
    *radius = particles[index].r;
    return 0;
}

int set_begin_time(double time){
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    particles[index].r = radius;
    return 0;
}

int cleanup_code(){
    if(particles != NULL) {
        particles_remove_all();
    }
    indexMap.clear();
    max_id = 0;
    has_removal = false;
    reset();
    return 0;
}

int recommit_parameters(){
    return 0;
}

int initialize_code(){
    max_id = 0;
    
    reset();
    
    integrator = WHFAST;
    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);
    //set_support_for_condition(TIMEOUT_DETECTION);
    //set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    //set_support_for_condition(OUT_OF_BOX_DETECTION);
    
    return 0;
}

int get_potential_energy(double * potential_energy){
    return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
    double * vz){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) { return -1;}
    
    *vx = particles[index].vx;
    *vy = particles[index].vy;
    *vz = particles[index].vz;
    
    return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
      double * z){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    *x = particles[index].x;
    *y = particles[index].y;
    *z = particles[index].z;
    return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    particles[index].x = x;
    particles[index].y = y;
    particles[index].z = z;
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
      double * az){
    
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    
    *ax = particles[index].ax;
    *ay = particles[index].ay;
    *az = particles[index].az;
    
    return 0;
}

int commit_parameters(){
    return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
    double vz){
    
    int index = get_index_from_identity(index_of_the_particle);
    if(index < 0) {return -1;}
    particles[index].vx = vx;
    particles[index].vy = vy;
    particles[index].vz = vz;
    return 0;
}


int _set_integrator(int value){
    
    switch(value){
        case 0:
            integrator = IAS15;
            break;
        case 1:
            integrator = WHFAST;
            break;
        case 2:
            integrator = SEI;
            break;
        case 3:
            integrator = WH;
            break;
        case 4:
            integrator = LEAPFROG;
            break;
        case 5:
            integrator = HYBRID;
            break;
        case 6:
            integrator = NONE;
            break;
        default:
            integrator = NONE;
            return -1;
            break;
            
    }
    return 0;
}

int _get_integrator(int * value){
    *value = integrator; 
    return 0;
}

