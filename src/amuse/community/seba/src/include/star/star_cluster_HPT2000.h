/*
 *  star_cluster.h: derived class for element evolution systems.
 *          functions as derived class for the real elements.
 *.............................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class star_cluster
 *
 *.............................................................................
 */
#ifndef     _STAR_CLUSTER
#   define  _STAR_CLUSTER

#include "stdinc.h"
#include "base_element.h"
#include "cluster_support.h"
#include "star.h"
#include "double_star.h"
#include "single_star.h"
#include "scatter3_support.h"

#include "node.h"
//#include "main_sequence.h"
#include "double_support.h"

/*-----------------------------------------------------------------------------
 *  star_cluster  --  class for evolving cluster of stars and double_stars
 *-----------------------------------------------------------------------------
 */
class star_cluster
    {
    protected:

       base_element * cluster;
       base_element * halo;

       bool field;

       real cluster_age;
       real dead_time;

       real core_radius;
       real halfm_radius;
       real tidal_radius;

       real core_density;	// n/pc^3

       real total_mass;
       real v_disp;
       
       int n_singles;
       int n_doubles;
       int nh_singles;
       int nh_doubles;

       int n_single_esc;
       int n_double_esc;
       int n_multiple_esc;
 
       int n_scatter;
       int n_collision;
       int n_capture;
     
       initial_cluster init;

    public:
  
       star_cluster(initial_cluster&);

       base_element* get_cluster()   {return cluster;}
       base_element* get_halo()   {return halo;}
       bool get_field()              {return field;}
       real get_cluster_age()        {return cluster_age;}
       void set_cluster_age(real t)  {cluster_age=t;}
       real get_core_density()       {return core_density;}
       real get_core_radius()        {return core_radius;}
       real get_halfm_radius()        {return halfm_radius;}
       real get_tidal_radius()        {return tidal_radius;}
       real get_core_volume()        {return (4./3.)*PI*pow(core_radius, 3);}
       int  no_of_elements()              {return n_singles+n_doubles;}
       real get_velocity_dispersion()	{return v_disp;}
/*
       void inc_n_singles()		{n_singles++;}
       void inc_n_doubles()		{n_doubles++;}
       void dec_n_singles()		{n_singles--;}
       void dec_n_doubles()		{n_doubles--;}
*/
       void inc_n_scatter()                     {n_scatter++;}
       void inc_n_capture()                     {n_capture++;}
       void inc_n_collision()                   {n_collision++;}
       int get_n_singles()			{return n_singles;}
       int get_n_doubles()			{return n_doubles;}
       int get_n_single_esc()			{return n_single_esc;}
       int get_n_double_esc()			{return n_double_esc;}

       real get_dead_time()		{return dead_time;}
       void inc_dead_time(real dt)	{dead_time += dt;}
       void reset_dead_time()		{dead_time = 0;}

       void set_initial_cluster(initial_cluster& initial) {init=initial;}
       void setup_star_cluster();
       void calculate_equilibrium_velocities();
       real escape_velocity();
       void randomize_maxwellian_velocities();
       real random_maxwellian_velocity(const real);
       void freeze();
       real velocity_dispersion();
       //real get_total_energy();
       //real get_potential_energy();
       //real get_kinetic_energy();
       real get_total_mass();
       real obtain_total_mass();
       real get_scatter_cross_section(base_element*);
       void count_number_of_elements();
       void update_cluster_parameters();
// 		Escaper handling
       void delete_escaper_from_cluster(base_element*);
       void add_new_particle_to_cluster(int);
       void remove_escapers();
       double_init& halo_guest();

       void evolve_element(real);
       void put_cluster() {
                //cluster->put_star(cout, cluster);

       }
       void print_roche();
       void put_state();
       void put_hrd(ostream &);
       void dump(ostream &); //char*);
       void dump(char*); //char*);
       void test_maxwellian_velocities();

//		Profiler
       void make_profile();

 
    };

#define  put_element  put_node

#endif      // _STAR_CLUSTER
