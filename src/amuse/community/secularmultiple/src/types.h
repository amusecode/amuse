#include <math.h>
#include <cstdlib>
#include <map>
#include <vector>

/* constants */
/* units (cf. interface.py): 
 * unit_l = units.AU
 * unit_m = units.MSun
 * unit_t = 1.0e6*units.yr
 */
 
#ifndef __CONSTANTS
#define __CONSTANTS
#define CONST_G			    (double)	3.94852492465e+13
#define CONST_G_P2          (double)    CONST_G*CONST_G
#define CONST_G_P3          (double)    CONST_G*CONST_G_P2
#define CONST_C_LIGHT		(double)	63239726386.8
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.0
#define CONST_R_SUN         (double)    0.00464913034382
#define CONST_L_SUN         (double)    2.71040410975e+14

#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div16            (double)    1.0/16.0
#define c_2div3             (double)    2.0/3.0
#define c_3div2             (double)    3.0/2.0
#define c_3div4             (double)    3.0/4.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_3div32            (double)    3.0/32.0
#define c_3div1024          (double)    3.0/1024.0
#define c_5div2             (double)    5.0/2.0
#define c_5div8             (double)    5.0/8.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_7div8             (double)    7.0/8.0
#define c_8div5             (double)    8.0/5.0
#define c_8div7             (double)    8.0/7.0
#define c_9div2             (double)    9.0/2.0
#define c_9div16            (double)    9.0/16.0
#define c_9div32            (double)    9.0/32.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_16div5            (double)    16.0/5.0
#define c_25div16           (double)    25.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_32div5            (double)    32.0/5.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div5            (double)    64.0/5.0
#define c_73div24           (double)    73.0/24.0
#define c_105div4096        (double)    105.0/4096.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#endif

/*	ODE solver macros	*/
#ifndef __ODE_MACROS
#define __ODE_MACROS
    #define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
    #define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
    #ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
    #endif
    #ifndef min
        #define min(X,Y) ((X) < (Y) ? (X) : (Y))
    #endif
#endif

/* vector operators */
#ifndef __VECTOR_OPERATORS
#define __VECTOR_OPERATORS
inline void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
inline double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
inline double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
#endif

/* classes */
#ifndef __Particle
#define __Particle
class Particle
{
    public:
    /* generic properties */
    int index,child1,child2;
    int parent;
    int sibling;
    std::vector<int> parents;
    std::vector<int> connecting_child_in_parents;
    int level;
    int is_binary;
    double mass,child1_mass,child2_mass;

    /*******************
    /* body properties *
     * ****************/
    /* general */
    double radius;
    double spin_vec_x,spin_vec_y,spin_vec_z;

    /* used in ODE solver only */
    double spin_vec[3],dspin_vec_dt[3]; 
    double spin_vec_norm;
    
    
    void set_ODE_quantities();
    void reset_ODE_quantities();
    
    /*********************
    /* binary properties *
     * ******************/
    /* general */
    double e_vec_x,e_vec_y,e_vec_z;
    double h_vec_x,h_vec_y,h_vec_z;
    
    /* PN terms */
    int include_pairwise_1PN_terms,include_pairwise_25PN_terms;
        
    /* tidal friction */
    int include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms;
    double minimum_eccentricity_for_tidal_precession;
    double tides_Q_prime; /* depricated */
    double tides_apsidal_motion_constant, tides_time_lag, tides_gyration_radius;
    double tides_viscous_time_scale;

    /* root finding */
    int check_for_secular_breakdown,secular_breakdown_has_occurred;
    int check_for_dynamical_instability,dynamical_instability_has_occurred,dynamical_instability_criterion;
    int dynamical_instability_central_particle;
    double dynamical_instability_K_parameter;
    int check_for_physical_collision_or_orbit_crossing,physical_collision_or_orbit_crossing_has_occurred;
    int check_for_minimum_periapse_distance,minimum_periapse_distance_has_occurred;
    double check_for_minimum_periapse_distance_value;
    int check_for_RLOF_at_pericentre,check_for_RLOF_at_pericentre_use_sepinsky_fit,RLOF_at_pericentre_has_occurred;

    /* used in ODE solver only */
    double e_vec[3],h_vec[3];
    double e_vec_unit[3],h_vec_unit[3];    
    double de_vec_dt[3],dh_vec_dt[3];
    double child1_mass_plus_child2_mass,child1_mass_minus_child2_mass,child1_mass_times_child2_mass;
    double e,e_p2;
    double j,j_p2,j_p3,j_p4,j_p5; // j=sqrt(1-e^2)
    double h,a;
    
    Particle(int index, int is_binary) : index(index), is_binary(is_binary) { minimum_eccentricity_for_tidal_precession = 1.0e-3; }
};
#endif

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;

int determine_binary_parents_levels_and_masses(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding);

/* CVODE UserData */
#ifndef __UserData
#define __UserData
typedef struct {
	ParticlesMap *particlesMap;
    double hamiltonian;
    int N_root_finding;
} *UserData;
#endif
