#include "node.h"
#include "single_star.h"
#include "main_sequence.h"
#include "worker_code.h"

#include <map>

static node * seba_root = 0;
static node * seba_insertion_point = 0;
static int next_seba_id = 1;
static map<int, nodeptr> mapping_from_id_to_node;
static double seba_metallicity = 0.02;
static double seba_time = 0.0;
static stellar_type start_type = Main_Sequence;

local real evolve_star_until_next_time(node* bi, const real out_time, const int n_steps) {
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
    return out_time;
}

local int translate_to_SSE_type(stellar_type stp, const real mass) {

  switch (stp) {
    case Brown_Dwarf:
  case Main_Sequence: 
      if (mass<0.1) {
	return 0;
      }
      return 1;
    case Hertzsprung_Gap:
      return 2;
    case Sub_Giant:
      return 3;
    case Horizontal_Branch:
      return 4;
    case Super_Giant:
      return 5;
    case Hyper_Giant:
      return 6;
    case Carbon_Star:
    case Helium_Star: 
      return 7;
    case Helium_Giant:
      return 9;
    case Helium_Dwarf:
      return 10;
    case Carbon_Dwarf:
      return 11;
    case Oxygen_Dwarf:
      return 12;
    case Xray_Pulsar:
    case Radio_Pulsar:
    case Neutron_Star: 
      return 13;
    case Black_Hole:
      return 14;
    case Disintegrated:
      return 15;
    case Proto_Star:
    case Planet:
    case Static_Star:
    case SPZDCH_Star:
    case NAS: 
    case Thorn_Zytkow:
    case Double:
    case no_of_stellar_type:
      return -1;
  }
}

int evolve_star(double mass, double endtime, double metal, double * resulttime, double * end_mass, double * end_radius, double * end_luminosity, double * end_temperature, double *end_time_step, int *end_stellar_type){
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
       *end_time_step = bi->get_starbase()->get_evolve_timestep();
       *end_stellar_type = translate_to_SSE_type(bi->get_starbase()->get_element_type(), mass);

        //    << "   " << type_string(bi->get_starbase()->get_element_type())
    }
    
    rmtree(root);
    
    return 0;
}


int initialize_code(){
    next_seba_id = 1;
    seba_root = new node();
    seba_root->set_root(seba_root); 
    seba_root->get_starbase()->set_stellar_evolution_scaling(
        1,   // mass units is 1 MSun
        2.255e-8, // radius units is RSun
        1   // time units is Myr
    );
    // XXX AVE
    // in all starlab files this this done after adding stars, 
    // do we need to follow that convention?
    // it looks like this is important for
    // double stars!!
    seba_root->get_starbase()->set_use_hdyn(false); 
    return 0;
}

int cleanup_code(){
    return 0;
}

int commit_parameters(){
    return 0;
}

int recommit_parameters(){
    return 0;
}

int set_metallicity(double metallicity){
    seba_metallicity = metallicity;
    return 0;
}
int get_metallicity(double * metallicity){
    *metallicity = seba_metallicity;
    return 0;
}




int new_particle(int * index_of_the_star, double mass){

    node * new_node = new node();
    new_node->set_label(next_seba_id);
    new_node->set_parent(seba_root);
    new_node->set_mass(mass);
    mapping_from_id_to_node[next_seba_id] = new_node;
    
    if(seba_insertion_point == 0) {
        seba_insertion_point = new_node;
        seba_root->set_oldest_daughter(new_node);
    } else {
        seba_insertion_point->set_younger_sister(new_node);
        new_node->set_elder_sister(seba_insertion_point);
        seba_insertion_point = new_node;
    }
    
    addstar(new_node, seba_time, start_type, seba_metallicity, 0, false);
    *index_of_the_star = next_seba_id;
    
    next_seba_id++;
    
    return 0;
}

int delete_star(int index_of_the_star){
    
    map<int, nodeptr>::iterator i = mapping_from_id_to_node.find(index_of_the_star);
    if(i == mapping_from_id_to_node.end()) {
        return -1;
    } else {
        node * node_to_remove = i->second;
        if (node_to_remove == seba_insertion_point) {
            seba_insertion_point = node_to_remove->get_younger_sister();
        }
        detach_node_from_general_tree(node_to_remove);
        
        mapping_from_id_to_node.erase(i);
        return 0;
    }
    
}

int commit_particles(){
    return 0;
}

int recommit_particles(){
    return 0;
}

node * get_seba_node_from_index(int index_of_the_star, int * errorcode)
{
    map<int, nodeptr>::iterator i = mapping_from_id_to_node.find(index_of_the_star);
    if(i == mapping_from_id_to_node.end()) {
        *errorcode = -1;
        return 0;
    } else {
        *errorcode = 0;
        return i->second;
    }
}

int get_mass(int index_of_the_star, double * mass){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *mass= seba_node->get_starbase()->get_total_mass() ;
    return error_code;
}

int get_temperature(int index_of_the_star, double * temperature){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *temperature= seba_node->get_starbase()->temperature() ;
    return error_code;
}


int get_time_step(int index_of_the_star, double * time_step){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *time_step= seba_node->get_starbase()->get_evolve_timestep() ;
    return error_code;
}

int get_luminosity(int index_of_the_star, double * luminosity){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *luminosity= seba_node->get_starbase()->get_luminosity() ;
    return error_code;
}

int get_age(int index_of_the_star, double * age){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *age= seba_node->get_starbase()->get_current_time() ;
    return error_code;
}

int get_radius(int index_of_the_star, double * radius){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *radius= seba_node->get_starbase()->get_effective_radius() ;
    return error_code;
}


int get_stellar_type(int index_of_the_star, int * stellar_type){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    double mass = seba_node->get_starbase()->get_total_mass();
    *stellar_type = translate_to_SSE_type(seba_node->get_starbase()->get_element_type(), mass);
    return error_code;
}

int evolve_one_step(int index_of_the_star){
    int error_code = 0;
    int n_steps_per_phase = 10;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    
    double out_time = seba_node->get_starbase()->get_current_time();
    out_time += seba_node->get_starbase()->get_evolve_timestep();
    evolve_star_until_next_time(seba_node, out_time, n_steps_per_phase);
    return error_code;
}

int evolve_for(int index_of_the_star, double delta_t){
    int error_code = 0;
    int n_steps_per_phase = 10;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    double out_time = seba_node->get_starbase()->get_current_time();
    out_time += delta_t;
    evolve_star_until_next_time(seba_node, out_time, n_steps_per_phase);
    return error_code;
}

int get_number_of_particles(int * number_of_particles){
    return 0;
}

int evolve_system(double end_time) {
    int n_steps = 1;
    int n_steps_per_phase = 10;
    real delta_t = (end_time - seba_time)/((real)n_steps);
    real out_time; 
    nodeptr bi;
    for_all_daughters(node, seba_root, bi) {
        out_time = seba_time;
        do {
            out_time = Starlab::min(out_time+delta_t, end_time);
            out_time = evolve_star_until_next_time(bi, out_time, n_steps_per_phase);
        }
        while(out_time < end_time);
    }
    return 0;
}
