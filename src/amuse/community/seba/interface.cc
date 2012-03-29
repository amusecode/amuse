#include "node.h"
#include "single_star.h"
#include "main_sequence.h"
#include "worker_code.h"

local void evolve_star_until_next_time(node* bi, const real out_time, const int n_steps) {
    ofstream starev("starev.data", ios::app|ios::out);
    bi->get_starbase()->dump(starev, false);  
    real current_time = ((star*)bi->get_starbase())->get_current_time();
    real time_step    =  bi->get_starbase()->get_evolve_timestep();

    while (out_time>current_time+time_step ) {
        bi->get_starbase()->evolve_element(current_time+time_step);
        bi->get_starbase()->dump(starev, false);                
        current_time = ((star*)bi->get_starbase())->get_current_time();
        time_step    =  bi->get_starbase()->get_evolve_timestep();
        
        star_state ss(dynamic_cast(star*, bi->get_starbase()));
    }
    
      
    bi->get_starbase()->evolve_element(out_time);
    bi->get_starbase()->dump(starev, false);
    bi->get_starbase()->dump(cerr, false);
    print_star(bi->get_starbase(), cerr);
    starev.close();
}

int evolve_star(double mass, double endtime, double metal, double * resulttime, double * end_mass, double * end_radius, double * end_luminosity, double * end_temperature){
    stellar_type type = Main_Sequence;
    char * star_type_string;
    int  c;

    bool  t_flag = FALSE;
    bool  S_flag = FALSE;
    bool  c_flag = FALSE;
    bool  M_flag = FALSE;
    bool  n_flag = FALSE;
    bool  R_flag = FALSE;
    bool  I_flag = FALSE;
    real  m_tot;
    real  r_hm = 100;
    real  t_hc = 1;
    real  t_start = 0;           // default value;
    real  t_end = 0;
    int n_steps = 1;
    int n_steps_per_phase = 10;
    int n_init = 0;
    int n =1;
    real z;
    
    char  *comment;
    int input_seed=0, actual_seed;
    *resulttime = 0.0;
    *end_mass = 0.0 ;
    *end_radius = 0.0;
    *end_luminosity = 0.0;
    *end_temperature = 0.0;
    if (metal < 0.00001) {
        return -4;
    }
    
    m_tot = mass;
    t_end = endtime;
    z = metal;
    
    actual_seed = srandinter(input_seed);
    node *root;
    root= mknode(1);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);
    
    addstar(root, t_start, type, z, n_init+0, false);
    root->get_starbase()->set_use_hdyn(false);
    
    real delta_t = t_end/((real)n_steps);
    real out_time; 
    
    for_all_daughters(node, root, bi) {
       out_time = 0;
       do {
            out_time = Starlab::min(out_time+delta_t, t_end);
            evolve_star_until_next_time(bi, out_time, n_steps_per_phase);
       }
       while(out_time < t_end);
    }
    
    for_all_daughters(node, root, bi) {
       *resulttime = bi->get_starbase()->get_current_time();
       *end_mass = bi->get_starbase()->get_total_mass() ;
       *end_radius = bi->get_starbase()->get_effective_radius(); 
       *end_luminosity = bi->get_starbase()->get_luminosity();
       *end_temperature = bi->get_starbase()->temperature();

        //    << "   " << type_string(bi->get_starbase()->get_element_type())
    }
    
    rmtree(root);
    
    return 0;
}

