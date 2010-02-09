#include "g6lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_NUMBER_OF_PARTICLES 100000

typedef double v4df[3];
 
struct g6_j_particle_tag {
    int    id;
    
    double tj;
    double dtj;
    double mass;
    
    v4df a2by18; 
    v4df a1by6;
    v4df aby2; 
    
    v4df v; 
    v4df x;
};
typedef struct g6_j_particle_tag g6_j_particle;

struct g6_i_particle_tag{
    int id;
    v4df x;
    v4df v;
    double eps2;
    double h2;  
    
    v4df acc;
    v4df jerk;
    double pot; 
    double nearest_r_squared;
    int    nearest_j;
};
typedef struct g6_i_particle_tag g6_i_particle;


struct g6_j_predict_tag {
    int    id;
    v4df v; 
    v4df x;
};
typedef struct g6_j_predict_tag g6_j_predicted_particle;

struct g6_unit_tag{
    g6_j_particle j_particles[MAX_NUMBER_OF_PARTICLES];
    g6_i_particle i_particles[MAX_NUMBER_OF_PARTICLES];
    g6_j_predicted_particle jp_particles[MAX_NUMBER_OF_PARTICLES];
    
    double ti;
    double nj;
    double ni;
    
};
typedef struct g6_unit_tag g6_unit;

g6_unit * unit = NULL;

void predict_positions_and_velocities_for_j_particles() {
    int n,k;
    for( n = 0 ; n < unit->nj; n++) {
        g6_j_particle  * current_j  = unit->j_particles + n;
        g6_j_predicted_particle * current_jp = unit->jp_particles + n;
        
        current_jp->id = current_j->id;
        
        double delta_t = unit->ti - current_j->tj;
        double delta_t2 = delta_t * delta_t;
        double delta_t3 = delta_t2 * delta_t;
        double delta_t4 = delta_t2 * delta_t2;
        
        for( k = 0; k < 3; k++) {
            double delta_x = (current_j->v[k] * delta_t);
            delta_x += (current_j->aby2[k] * delta_t2);
            delta_x += (current_j->a1by6[k] * delta_t3);
            delta_x += (current_j->a2by18[k] * 3.0 / 4.0 * delta_t4);
            
            current_jp->x[k] = delta_x + current_j->x[k];
            
            double delta_v = ((current_j->aby2[k]) * 2.0 * delta_t);
            delta_v += ((current_j->a1by6[k]) * 3.0 * delta_t3);
            delta_v += ((current_j->a2by18[k]) * 6.0  * delta_t4);
            
            current_jp->v[k] = delta_v + current_j->v[k];
        }
    }
}

static inline double vec_squared(v4df v){
    int k;
    double result = 0.0;
    for(k = 0; k < 3; k++) {
        result += v[k] * v[k];
    }
    return result;
}

static inline double vec_dot(v4df va, v4df vb){
    return (va[0] * vb[0]) + (va[1] * vb[1]) + (va[2] * vb[2]);
}

inline clear_i_particle(g6_i_particle * particle) {
    int k;
    for(k = 0; k < 3; k++) {
        particle->acc[k] = 0.0;
        particle->jerk[k] = 0.0;
    }
    particle->pot = 0.0;
    particle->nearest_r_squared = 0.0;
    particle->nearest_j = -1;     
}

void calculate_acceleration_jerk_and_potential_for_i_particles() {
    int i,j,k;
    
    for(i = 0 ; i < unit->ni; i++) {
        g6_i_particle  * current_i  = unit->i_particles + i;
        
        clear_i_particle(current_i);  
           
        for(j = 0; j < unit->nj ; j++) {
            g6_j_particle  * current_j  = unit->j_particles + j;
            g6_j_predicted_particle * current_jp = unit->jp_particles + j; 
            
            if(current_i->id == current_j->id) {
                continue;
            }
            if(current_j->id ==  -1) {
                continue;
            }
            
            v4df rij = {0.0,0.0,0.0};
            v4df vij = {0.0,0.0,0.0};
            for(k = 0; k < 4; k++) {
                rij[k] = current_jp->x[k] - current_i->x[k];
                vij[k] = current_jp->v[k] - current_i->v[k];
            }
            
            double r_squared = vec_squared(rij);
            
            if( current_i->nearest_j < 0 || r_squared < current_i->nearest_r_squared) {
                current_i->nearest_r_squared = r_squared;
                current_i->nearest_j = current_j->id;
            }
            
            double r_squared_smooth = r_squared + current_i->eps2;
            
            /*
            if(r_squared_smooth == 0.0) {
                fprintf(stderr,"distance between two particles is 0!n");
                continue;
            }
            */
            
            double r_squared_smooth_raised_32 = pow(r_squared_smooth, 3.0 / 2.0);
            double r_squared_smooth_raised_52 = pow(r_squared_smooth, 5.0 / 2.0);
            double r_squared_smooth_raised_12 = pow(r_squared_smooth, 1.0 / 2.0);
                       
            double gmj = current_j->mass;
            
            double r_dot_v = vec_dot(rij, vij);
            
            for(k = 0; k < 3; k++) {
                double acceleration = gmj * rij[k] / r_squared_smooth_raised_32;
                
                double jerk = vij[k] / r_squared_smooth_raised_32;
                jerk -= (3 * r_dot_v * rij[k]) / r_squared_smooth_raised_52;
                jerk *= gmj;
                
                current_i->acc[k] += acceleration;
                current_i->jerk[k] += jerk;
            }
            current_i->pot += - gmj / r_squared_smooth_raised_12;
        }
    }
}




int g6_open_(int *id) {
    int i =0;
    
    if(unit) {
        free(unit);
    }
    unit = (g6_unit *) malloc(sizeof(g6_unit));
    if (!unit) {
        return -1;
    }
    
    g6_reset_(id);
    //fprintf(stderr, "open %d\n", sizeof(g6_unit));
    return 0;
}

int g6_close_(int *id){
    if (unit) {
        free(unit);
        unit = NULL;
    }
    //fprintf(stderr, "close\n");
    return 0;
}

int g6_npipes_(){
    //fprintf(stderr, "npipes\n");
    return 1;
}

int g6_set_tunit_(double* p){
    //fprintf(stderr, "g6_set_tunit_\n");
    return 0;
}
int g6_set_xunit_(double* p){
    //fprintf(stderr, "g6_set_xunit_\n");
    return 0;
}

int g6_set_ti_(int *id, double *ti){
    //fprintf(stderr, "g6_set_ti_ %f\n", ti);
    unit->ti = *ti;
    return 0;
}

int g6_set_j_particle_(int *cluster_id,
         int *address,
         int *index,
         double *tj, double *dtj,
         double *mass,
         double a2by18[3], double a1by6[3],
         double aby2[3], double v[3], double x[3]){
    //fprintf(stderr, "g6_set_j_particle_ %d, %d\n", *address, *index);
    if(*address > MAX_NUMBER_OF_PARTICLES) {
        return -1;
    } 
    
    g6_j_particle * particle = &(unit->j_particles[*address]);
    particle->id = *index;
    particle->tj = *tj;
    particle->dtj = *dtj;
    particle->mass = *mass;

    particle->a2by18[0] = a2by18[0];
    particle->a2by18[1] = a2by18[1];
    particle->a2by18[2] = a2by18[2];
    
    particle->a1by6[0] = a1by6[0];
    particle->a1by6[1] = a1by6[1];
    particle->a1by6[2] = a1by6[2];
    
    particle->aby2[0] = aby2[0];
    particle->aby2[1] = aby2[1];
    particle->aby2[2] = aby2[2];
    
    particle->v[0] = v[0];
    particle->v[1] = v[1];
    particle->v[2] = v[2];
    
    particle->x[0] = x[0];
    particle->x[1] = x[1];
    particle->x[2] = x[2];
    
    //fprintf(stderr,"xj[%d,%d] = %f (%d)\n", *address, 0, particle->x[0], particle );
    return 0;
}

void g6calc_firsthalf_(int *cluster_id,
         int *nj, int *ni,
         int index[], 
         double xi[][3], double vi[][3],
         double aold[][3], double j6old[][3],
         double phiold[3], 
         double *eps2, double h2[]){
    
    //fprintf(stderr, "g6calc_firsthalf_ %d, %d, %d\n", *nj, *ni, index[0]);
    int i = 0;
    unit->ni = *ni;
    unit->nj = *nj;
    for(i=0; i<*ni; i++) {
        g6_i_particle * particle = &(unit->i_particles[i]);
        particle->id = index[i];
        
        particle->v[0] = vi[i][0];
        particle->v[1] = vi[i][1];
        particle->v[2] = vi[i][2];
        
        particle->x[0] = xi[i][0];
        particle->x[1] = xi[i][1];
        particle->x[2] = xi[i][2];
        
        particle->eps2 = *eps2;
        
        particle->h2 = h2[i];
    }
}





int g6calc_lasthalf_(int *cluster_id,
           int *nj, int *ni,
           int index[], 
           double xi[][3], double vi[][3],
           double *eps2, double h2[],
           double acc[][3], double jerk[][3], double pot[]){
    //fprintf(stderr, "g6calc_lasthalf_ %d, %d, %d\n", *nj, *ni, index[0]);
    
    
    predict_positions_and_velocities_for_j_particles();
    calculate_acceleration_jerk_and_potential_for_i_particles();
    
    int i = 0;
    for(i=0; i<*ni; i++) {
        g6_i_particle * particle = &(unit->i_particles[i]);
        index[i] = particle->id;
        
        acc[i][0] = particle->acc[0];
        acc[i][1] = particle->acc[1];
        acc[i][2] = particle->acc[2];
        
        jerk[i][0] = particle->jerk[0];
        jerk[i][1] = particle->jerk[1];
        jerk[i][2] = particle->jerk[2];
        
        pot[i] = particle->pot;
    }
    return 0;
}
int g6calc_lasthalf2_(int *cluster_id,
        int *nj, int *ni,
        int index[], 
        double xi[][3], double vi[][3],
        double *eps2, double h2[],
        double acc[][3], double jerk[][3], double pot[],
        int nnbindex[]){
            
    predict_positions_and_velocities_for_j_particles();
    calculate_acceleration_jerk_and_potential_for_i_particles();
    
    int i = 0;
    for(i=0; i<*ni; i++) {
        g6_i_particle * particle = &(unit->i_particles[i]);
        index[i] = particle->id;
        
        acc[i][0] = particle->acc[0];
        acc[i][1] = particle->acc[1];
        acc[i][2] = particle->acc[2];
        
        jerk[i][0] = particle->jerk[0];
        jerk[i][1] = particle->jerk[1];
        jerk[i][2] = particle->jerk[2];
        
        pot[i] = particle->pot;
        nnbindex[i] = particle->nearest_j;
    }
    return 0;
}

int g6_initialize_jp_buffer_(int* cluster_id, int* buf_size){
    //fprintf(stderr, "g6_initialize_jp_buffer_\n");
    return 0;
}
int g6_flush_jp_buffer_(int* cluster_id){
    //fprintf(stderr, "g6_flush_jp_buffer_\n");
    return 0;
}
int g6_reset_(int* cluster_id){
    //fprintf(stderr, "g6_reset_\n");
    int i;
    
    memset(unit, sizeof(g6_unit), 0);
    for (i = 0;  i < MAX_NUMBER_OF_PARTICLES; i++) {
       unit->j_particles[i].id = -1; 
       unit->i_particles[i].id = -1; 
       unit->jp_particles[i].id = -1;
    }
    
    return 0;
}
int g6_reset_fofpga_(int* cluster_id){
    //fprintf(stderr, "g6_reset_fofpga_\n");
    return 0;
}

int g6_read_neighbour_list_(int* cluster_id){
    //fprintf(stderr, "g6_read_neighbour_list_\n");
    return 0;
}

int g6_get_neighbour_list_(int *cluster_id,
             int *ipipe,
             int *maxlength,
             int *n_neighbours,
             int neighbour_list[]){
    //fprintf(stderr, "g6_get_neighbour_list_\n");
    return 0;
}
