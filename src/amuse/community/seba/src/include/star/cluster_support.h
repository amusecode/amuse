#ifndef    _CLUSTER_SUPPORT
#   define _CLUSTER_SUPPORT

//#include "double_star.h"
#include "stdinc.h"
#include "double_support.h"
#include "stdfunc.h"
#include "star_support.h"
#include "star_state.h"

#define MIN_TYPE_NUMBER  0
#define MAX_TYPE_NUMBER  60

enum profile 	{star_type=0, spectral_type, 
                 spectral_addition, no_of_profile};

/*-----------------------------------------------------------------------------
 *  initial_cluster --	
 *-----------------------------------------------------------------------------
 */

struct initial_cluster
     {

    real start_time;
    real end_time;
    int  field;

    int   start_id;
    int   n_steps;
    int   no_of_binaries;
    int   no_of_singles;
    int   nh_binaries;
    int   nh_singles;

    real  r_core;
    real  r_halfm;
    real  r_tidal;
    real  rho_core;
    real  v_disp;

    real  m_min;
    real  m_max;
    real  m_alpha;

    real  q_min;
    real  q_max;
    real  q_alpha;

    real  a_min;
    real  a_max;
    real  a_alpha;

    real  e_min;
    real  e_max;
    real  e_alpha;
    
    int seed;

    initial_cluster();
    real get_delta_t() {return (end_time - start_time)/n_steps;}

    double_init random_initial_conditions();
}; 

    double_init random_initial_conditions(initial_cluster&);


/*-----------------------------------------------------------------------------
 * cluster_profile --
 *-----------------------------------------------------------------------------
 */
struct cluster_profile
     {

        int init_bins[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
        int fin_bins[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
        int mergers[MAX_TYPE_NUMBER];
        int singles[MAX_TYPE_NUMBER];

        int bins_stt[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
        int runner_stt[MAX_TYPE_NUMBER];
        int merge_stt[MAX_TYPE_NUMBER];
        int single_stt[MAX_TYPE_NUMBER];

        int bins_spt[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
        int runner_spt[MAX_TYPE_NUMBER];
        int merge_spt[MAX_TYPE_NUMBER];
        int single_spt[MAX_TYPE_NUMBER];

        int bins_spa[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER][no_of_spec_type];
        int runner_spa[MAX_TYPE_NUMBER][no_of_spec_type];
        int merge_spa[MAX_TYPE_NUMBER][no_of_spec_type];
        int single_spa[MAX_TYPE_NUMBER][no_of_spec_type];


        int no_of_init_bins;
        int no_of_fin_bins;
        int no_of_runners;
        int no_of_mergers;
        int no_of_singles;

     cluster_profile();
     void enhance_cluster_profile(double_profile&, profile, star_type_spec);
     void enhance_cluster_profile(double_profile&);
     void enhance_cluster_profile(star_state&);
     void star_type_cluster_profile(double_profile&);
     void star_type_cluster_profile(star_state&);
     void spectral_type_cluster_profile(double_profile&);
     void spectral_type_cluster_profile(star_state&);
     void spectral_addition_single_profile(double_profile&, star_type_spec);
     void spectral_addition_cluster_profile(double_profile&, star_type_spec);
     void spectral_addition_cluster_profile(star_state&, star_type_spec);
     void print_profile();
     void print_profile(profile);
     void print_profile(stellar_type);
     void print_profile(spectral_class);
     void print_profile(star_type_spec, profile);
     void print_single_profile(star_type_spec);
};

void make_profile(initial_cluster&, 
                  cluster_profile&, profile, star_type_spec);

real next_output_time(int, int, real, real);

/*
struct ---II--initialization{

    real start_time;
    real end_time;
    int  field;

    int   n_steps;
    int   no_of_binaries;
    int   no_of_singles;

    real  m_min;
    real  m_max;
    real  m_alpha;

    real  q_min;
    real  q_max;
    real  q_alpha;

    real  a_min;
    real  a_max;
    real  a_alpha;

    real  e_min;
    real  e_max;
    real  e_alpha;

    clstr_initialization();
    real get_delta_t() {return (end_time - start_time)/n_steps;}

    double_init random_initial_conditions();
};

    double_init random_initial_conditions(clstr_initialization&);

*/

#endif  // _CLUSTER_SUPPORT


