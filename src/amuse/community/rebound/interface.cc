#include "worker_code.h"

extern "C" {
#define restrict 

#include "rebound.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "boundary.h"
#include "gravity.h"
#include "output.h"
#include "particle.h"
#include "collision.h"

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
#include <vector>
#include <iostream>
#include <algorithm>

#include <stopcond.h>

typedef struct _particle_location {
    reb_simulation * code;
    int index;
    int subset;
    _particle_location():code(0),index(-1),subset(-1){}
    _particle_location(reb_simulation * c, int i, int s):code(c),index(i),subset(s) {}
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

typedef std::map<int, particle_location> IndexMap;
static IndexMap indexMap;

typedef struct _particle_sort {
    int ref_index;
    reb_particle * p;
    _particle_sort(int i, reb_particle * p):ref_index(i),p(p){}
} particle_sort;
typedef std::vector<_particle_sort>  ParticleSortVector;

static int max_id = 0;
//m ,m mnstatic reb_simulation * code;

typedef std::vector<code_state> ReboundSimulationVector;
static ReboundSimulationVector codes;
static particle_location sentinel = particle_location();
static double _time;
static double timestep = 0.0001;

static inline particle_location get_index_from_identity(int id)
{
    IndexMap::iterator i = indexMap.find(id);
    if(i == indexMap.end()) {
        return sentinel;
    }
    return (*i).second;
}


static inline int get_identity_from_index(particle_location location)
{
    for( IndexMap::iterator i = indexMap.begin(); i != indexMap.end(); i++) {
        if((*i).second.index == location.index &&  indexMap [(*i).first].code == location.code) {
            return (*i).first;
        }
    }
    return -1;
}

int get_mass(int index_of_the_particle, double * mass){
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {*mass = 0;  return -1;}
    *mass = loc.code->particles[loc.index].m;
    return 0;
}


bool sort_particles (particle_sort i,particle_sort j) { 
    return (i.p->m > j.p->m); 
}

int commit_particles(){
    /*
     * after deletion of one or more particles, clear out
     * the original data and rebuild without these particles.
     */
     
    for( ReboundSimulationVector::iterator j = codes.begin(); j != codes.end(); j++) {
        code_state cs = *j;
        bool has_removal = cs.has_removal;
        if(has_removal) {
            if(cs.code) {
                struct reb_particle * previous = (struct reb_particle *) malloc(sizeof(struct reb_particle) * cs.code->N);
                memcpy(previous, cs.code->particles, sizeof(struct reb_particle) * cs.code->N);
                reb_remove_all(cs.code);
                for( IndexMap::iterator i = indexMap.begin(); i != indexMap.end(); i++) {
                    if( (*i).second.code == cs.code ) {
                        struct reb_particle * p = previous + (*i).second.index;
                        reb_add(cs.code, *p);
                        indexMap[(*i).first] = particle_location(cs.code, cs.code->N - 1, cs.subset);
                    }
                }
            }
            cs.has_removal = false;
            *j = cs;
        }
        if(cs.has_unsorted_massless_particles || has_removal) {
            if(cs.code) {
                //std::cout<<"HAS MASSLESS, WILL SORT"<<std::endl;
                ParticleSortVector sortvector;
                
                //std::cout<<"NA:"<<(cs.code)->N_active<<std::endl;
                struct reb_particle * previous = (struct reb_particle *) malloc(sizeof(struct reb_particle) * cs.code->N);
                memcpy(previous, cs.code->particles, sizeof(struct reb_particle) * cs.code->N);
                reb_remove_all(cs.code);
                for( IndexMap::iterator i = indexMap.begin(); i != indexMap.end(); i++) {
                    if( (*i).second.code == cs.code ) {
                        struct reb_particle * p = previous + (*i).second.index;
                        sortvector.push_back(particle_sort((*i).first, p));
                        
                        /*reb_add(cs.code, *p);
                        indexMap[(*i).first] = particle_location(cs.code, cs.code->N - 1, cs.subset);*/
                        
                    }
                }
                
                std::sort (sortvector.begin(), sortvector.end(), sort_particles);
                cs.code->N_active = 0;
                for( ParticleSortVector::iterator i = sortvector.begin(); i != sortvector.end(); i++) {
                    reb_add(cs.code,*(*i).p);
                    
                    //std::cout<<"p:"<<(*i).ref_index<<", m:"<<(*i).p->m<<std::endl;
                    if((*i).p->m > 0) {
                        cs.code->N_active++;
                    }
                    indexMap[(*i).ref_index] = particle_location(cs.code, cs.code->N - 1, cs.subset);
                }
                //std::cout<<"NA:"<<(cs.code)->N_active<<std::endl;
                sortvector.clear();
            }
            cs.has_unsorted_massless_particles = false;
            *j = cs;
        }
    }
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
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    loc.code->particles[loc.index].m = mass;
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
        return -10;
      }
      int new_id = max_id++;
      struct reb_particle pt;
      pt.x = x;
      pt.y = y;
      pt.z = z;
      pt.vx = vx;
      pt.vy = vy;
      pt.vz = vz;
      pt.ax = 0;
      pt.ay = 0;
      pt.az = 0;
      pt.m = mass;
      pt.ax = 0;
      pt.ay = 0;
      pt.r = radius; 
      pt.lastcollision = 0;
      pt.c = NULL;
      pt.id = new_id;
      if(pt.m == 0.0) {
          codes[code_index].has_unsorted_massless_particles = true;
      }
      reb_add(codes[code_index].code, pt);
      //std::cout<<"new particle :"<<pt.id<< " << "<<code_index<<" << "<<pt.x<<std::endl;
      indexMap[new_id] = particle_location(codes[code_index].code, codes[code_index].code->N - 1, code_index);
      *index_of_the_particle = new_id;
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
    double minD = 0;
    struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	double tmax = _tmax;
	code-> dt_last_done = code->dt;
	int last_step = 0;
	int ret_value = 0;
	const double dtsign = copysign(1.,code->dt); 				// Used to determine integration direction
    
    time_t starttime, currenttime;
    time(&starttime);
    
    while(code->t*dtsign<tmax*dtsign && last_step<2 && ret_value==0){
		if (code->N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return(1);
		}
       
		reb_step(code); 								// 0 to not do timing within step 
        
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
                            set_stopping_condition_particle_index(stopping_index, 0, get_identity_from_index(particle_location(code,i, -1)));
                            set_stopping_condition_particle_index(stopping_index, 1, get_identity_from_index(particle_location(code,j, -1)));
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
                        set_stopping_condition_particle_index(stopping_index, 0, get_identity_from_index(particle_location(code,i, -1)));
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
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double timing = timing_final-timing_initial;
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
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    
    indexMap.erase(index_of_the_particle);
    for( ReboundSimulationVector::iterator i = codes.begin(); i != codes.end(); i++) {
        code_state cs = *i;
        if(cs.code == loc.code){
            cs.has_removal = true;
        }
        *i = cs;
    }
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
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    
    
    loc.code->particles[loc.index].x = x;
    loc.code->particles[loc.index].y = y;
    loc.code->particles[loc.index].z = z;
    loc.code->particles[loc.index].vx = vx;
    loc.code->particles[loc.index].vy = vy;
    loc.code->particles[loc.index].vz = vz;
    
    loc.code->particles[loc.index].m = mass;
    loc.code->particles[loc.index].r = radius;
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
    double * y, double * z, double * vx, double * vy, double * vz, 
    double * radius, int * subset){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    
    *x = loc.code->particles[loc.index].x;
    *y = loc.code->particles[loc.index].y;
    *z = loc.code->particles[loc.index].z;
    *vx = loc.code->particles[loc.index].vx;
    *vy = loc.code->particles[loc.index].vy;
    *vz = loc.code->particles[loc.index].vz;
    *mass = loc.code->particles[loc.index].m;
    *subset = loc.subset;
#ifndef COLLISIONS_NONE
    *radius = loc.code->particles[loc.index].r;
#else
    *radius = 0;
#endif // COLLISIONS_NONE
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
    *number_of_particles = indexMap.size();
    return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    loc.code->particles[loc.index].ax = ax;
    loc.code->particles[loc.index].ay = ay;
    loc.code->particles[loc.index].az = az;
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    return 0;
}

int get_radius(int index_of_the_particle, double * radius){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {*radius = 0; return -1;}
    *radius = loc.code->particles[loc.index].r;
    return 0;
}
int get_subset(int index_of_the_particle, int * subset){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {*subset = 0; return -1;}
    *subset = loc.subset;
    return 0;
}

int set_begin_time(double time){
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    loc.code->particles[loc.index].r = radius;
    return 0;
}
int set_subset(int index_of_the_particle, int subset){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    if(loc.subset != subset) {return -2;}
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
    indexMap.clear();
    codes.clear();
    max_id = 0;
    return 0;
}

int recommit_parameters(){
    return 0;
}

int initialize_code(){
    max_id = 0;
    _time=0;
    reb_simulation * code = reb_create_simulation();
    codes.push_back(code_state(code));
    code->integrator = reb_simulation::REB_INTEGRATOR_WHFAST;
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
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    *vx = loc.code->particles[loc.index].vx;
    *vy = loc.code->particles[loc.index].vy;
    *vz = loc.code->particles[loc.index].vz;
    
    return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
      double * z){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    *x = loc.code->particles[loc.index].x;
    *y = loc.code->particles[loc.index].y;
    *z = loc.code->particles[loc.index].z;
    return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    loc.code->particles[loc.index].x = x;
    loc.code->particles[loc.index].y = y;
    loc.code->particles[loc.index].z = z;
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
      double * az){
    
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    
    *ax = loc.code->particles[loc.index].ax;
    *ay = loc.code->particles[loc.index].ay;
    *az = loc.code->particles[loc.index].az;
    
    return 0;
}

int commit_parameters(){
    return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
    double vz){
    
    particle_location loc = get_index_from_identity(index_of_the_particle);
    if(loc.index < 0) {return -1;}
    loc.code->particles[loc.index].vx = vx;
    loc.code->particles[loc.index].vy = vy;
    loc.code->particles[loc.index].vz = vz;
    return 0;
}

int new_subset(int * index, double time_offset) {
    reb_simulation * code = reb_create_simulation();
    reb_integrator_reset(code);
    code->dt = timestep;
    if(time_offset < 0) {time_offset = _time;}
    codes.push_back(code_state(code, time_offset, codes.size()));
    code->integrator = reb_simulation::REB_INTEGRATOR_WHFAST;
    code->t = time_offset;
    *index = codes.size() - 1;
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
            code->integrator = reb_simulation::REB_INTEGRATOR_WH;
            break;
        case 4:
            code->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
            break;
        case 5:
            code->integrator = reb_simulation::REB_INTEGRATOR_HYBRID;
            break;
        case 6:
            code->integrator = reb_simulation::REB_INTEGRATOR_NONE;
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

