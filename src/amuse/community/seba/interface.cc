#include "node.h"
#include "single_star.h"
#include "main_sequence.h"
#include "worker_code.h"
#include "double_star.h"

#include <map>

static node * seba_root = 0;
static node * seba_insertion_point = 0;
static int next_seba_id = 1;
static map<int, nodeptr> mapping_from_id_to_node;
static double seba_metallicity = 0.02;
static double seba_time = 0.0;
static stellar_type start_type = Main_Sequence;
static binary_type binary_start_type = Detached;
static bool is_logging_of_evolve_enabled = false;

local void addbinary(
    node *bi, real stellar_time,	// bi is binary CM
    binary_type type,
    real sma, real ecc)
{
    node* od = bi->get_oldest_daughter();
    int id;
    if((od->is_low_level_node() &&
        od->get_younger_sister()->is_low_level_node()) &&
       (od->get_elder_sister() == NULL) &&
       bi->n_leaves()==2 ) {

        if (!has_dstar(od)) {
            //bi->get_parent()->get_starbase()->get_element_type()!=Double) {

            //    if (bi->is_parent() && !bi->is_root())
            story * old_story = bi->get_starbase()->get_star_story();
            bi->get_starbase()->set_star_story(NULL);

            id = bi->get_index();

            // cerr << "Adding binary to "<< id << " at time = "
            //      << stellar_time << endl;

            double_star* new_double
            = new_double_star(bi, sma, ecc, stellar_time, id, type);
            // synchronize_binary_components(dynamic_cast(double_star*,
            // 						  bi->get_starbase()));

            // Give the new binary the old star_story.

            new_double->set_star_story(old_story);

            
        }
        else {
            cerr << "No double_star to node needed in adddouble"<<endl;
        }
    }
    else {
        cerr << "Sorry, no binary node in adddouble"<<endl;
    }
}


local real evolve_star_until_next_time(node* bi, const real out_time, const int n_steps) {
    ofstream starev;
    if(is_logging_of_evolve_enabled) {
        starev.open("starev.data", ios::app|ios::out);
        bi->get_starbase()->dump(starev, false);  
    }
    real current_time = ((star*)bi->get_starbase())->get_current_time();
    real time_step    =  bi->get_starbase()->get_evolve_timestep();

    while (out_time>current_time+time_step ) {
        bi->get_starbase()->evolve_element(current_time+time_step);
        if(is_logging_of_evolve_enabled) {
            bi->get_starbase()->dump(starev, false);     
        }
        current_time = ((star*)bi->get_starbase())->get_current_time();
        time_step    =  bi->get_starbase()->get_evolve_timestep();
        
        star_state ss(dynamic_cast(star*, bi->get_starbase()));
    }
    
      
    bi->get_starbase()->evolve_element(out_time);
    bi->get_starbase()->dump(cerr, false);
    
    if(is_logging_of_evolve_enabled) {
        bi->get_starbase()->dump(starev, false);
        print_star(bi->get_starbase(), cerr);
        starev.close();
    }
    return out_time;
}

local int translate_stellar_type_to_int(stellar_type stp, const real mass) {

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



local int translate_binary_type_to_int(binary_type btp) {
  switch (btp) {
    case Strong_Encounter:
        return -1;
    case Unknown_Binary_Type:
        return 0;
    case Synchronized:
        return 1;
    case Detached:
        return 2;
    case Semi_Detached:
        return 3;
    case Contact:
        return 4;
    case Common_Envelope:
        return 5;
	case Double_Spiral_In:
        return 6;
    case Merged:
        return 7;
    case Disrupted:
        return 8;
    case Spiral_In:
        return 9;
    default:
        return -2;
  }
}

local binary_type translate_int_to_binary_type(int btp) {
  switch (btp) {
    case -1:
        return Strong_Encounter;
    case 0:
        return Unknown_Binary_Type;
    case 1:
        return Synchronized;
    case 2:
        return Detached;
    case 3:
        return Semi_Detached;
    case 4:
        return Contact;
    case 5:
        return Common_Envelope;
	case 6:
        return Common_Envelope;
    case 7:
        return Merged;
    case 8:
        return Disrupted;
    case 9:
        return Spiral_In;
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
       *end_stellar_type = translate_stellar_type_to_int(bi->get_starbase()->get_element_type(), mass);

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
        100, //was 2.255e-8, // radius units is RSun
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

int set_supernova_kick_velocity(double v_disp) {
    cnsts.v_disp = v_disp;
    PRL(v_disp);
    PRL(cnsts.v_disp);
    return 0;
}
int get_supernova_kick_velocity(double * v_disp){
  *v_disp = cnsts.v_disp;
    return 0;
}

int get_is_logging_of_evolve_enabled(int *value){
    *value = is_logging_of_evolve_enabled ? 1 : 0;
    return 0;
}
int set_is_logging_of_evolve_enabled(int value){
    is_logging_of_evolve_enabled = value == 1;
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
    *stellar_type = translate_stellar_type_to_int(seba_node->get_starbase()->get_element_type(), mass);
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



int new_binary(
    int * index_of_the_star, 
    double semi_major_axis, 
    double eccentricity, 
    int index_of_child1, 
    int index_of_child2){
    int error_code = 0;
    node * child1 = get_seba_node_from_index(index_of_child1, &error_code);
    if(error_code < 0) {return error_code;}
    node * child2 = get_seba_node_from_index(index_of_child2, &error_code);
    if(error_code < 0) {return error_code;}
    
    if (child1 == seba_insertion_point) {
        seba_insertion_point = child1->get_younger_sister();
    }
    if (child2 == seba_insertion_point) {
        seba_insertion_point = child2->get_younger_sister();
    }
    if (child1 == seba_insertion_point) {
        seba_insertion_point = child1->get_younger_sister();
    }
    detach_node_from_general_tree(child1);
    detach_node_from_general_tree(child2);   
    
    node * new_node = new node();
    new_node->set_label(next_seba_id);
    new_node->set_parent(seba_root);
    new_node->set_mass(child1->get_mass() + child2->get_mass());
    mapping_from_id_to_node[next_seba_id] = new_node;
    
    if(seba_insertion_point == 0) {
        seba_insertion_point = new_node;
        seba_root->set_oldest_daughter(new_node);
    } else {
        seba_insertion_point->set_younger_sister(new_node);
        new_node->set_elder_sister(seba_insertion_point);
        seba_insertion_point = new_node;
    }
    new_node->set_oldest_daughter(child1);
    child1->set_younger_sister(child2);
    child1->set_elder_sister(0);
    child2->set_younger_sister(0);
    child2->set_elder_sister(child1);
    child1->set_parent(new_node);
    child2->set_parent(new_node);
    
    addbinary(new_node, seba_time, binary_start_type, semi_major_axis, eccentricity);
    *index_of_the_star = next_seba_id;
    
    next_seba_id++;
    
    return 0;
}

int delete_binary(int index_of_the_star){
    
    map<int, nodeptr>::iterator i = mapping_from_id_to_node.find(index_of_the_star);
    if(i == mapping_from_id_to_node.end()) {
        return -1;
    } else {
        node * node_to_remove = i->second;
        node * parent = node_to_remove->get_parent();
        node * younger_sister = node_to_remove->get_younger_sister();
        
        node * child1 = node_to_remove->get_oldest_daughter();
        if(child1 != 0) {
            node * child2 = child1->get_younger_sister();
            
            node * elder_sister = node_to_remove->get_elder_sister();
            
            child1->set_elder_sister(elder_sister);
            child1->set_parent(parent);
            if(elder_sister) {
                elder_sister->set_younger_sister(child1);
            }
            
            if(parent->get_oldest_daughter() == node_to_remove) {
                parent->set_oldest_daughter(child1);
            }
            
            if(child2 != 0) {
                child2->set_younger_sister(younger_sister);
                child2->set_parent(parent);
                
                if(younger_sister) {
                    younger_sister->set_elder_sister(child2);
                }
                
                if (node_to_remove == seba_insertion_point) {
                    seba_insertion_point = child2;
                }
            } else {
                if(younger_sister) {
                    younger_sister->set_elder_sister(child1);
                }
                if (node_to_remove == seba_insertion_point) {
                    seba_insertion_point = child1;
                }
            }
            
        } else {
            detach_node_from_general_tree(node_to_remove);
            
            if (node_to_remove == seba_insertion_point) {
                seba_insertion_point = younger_sister;
            }
        }

        mapping_from_id_to_node.erase(i);
        return 0;
    }
    
}


int get_children_of_binary(
    int index_of_the_star,
    int * child1_index, int * child2_index)
{
    
    *child1_index = -1;
    *child2_index = -1;
    map<int, nodeptr>::iterator i = mapping_from_id_to_node.find(index_of_the_star);
    if(i == mapping_from_id_to_node.end()) {
        return -1;
    } else {
        node * binary = i->second;
        
        node * child1 = binary->get_oldest_daughter();
        if(child1 != 0) {
            node * child2 = child1->get_younger_sister();
            *child1_index = child1->get_index();
            *child2_index = child2->get_index();
        }
        return 0;
    }   
}

int get_eccentricity(int index_of_the_star, double * value){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *value= seba_node->get_starbase()->get_eccentricity() ;
    return error_code;
}

int get_semi_major_axis(int index_of_the_star, double * value){
    int error_code = 0;
    node * seba_node = get_seba_node_from_index(index_of_the_star, &error_code);
    if(error_code < 0) {return error_code;}
    *value= seba_node->get_starbase()->get_semi() ;
    return error_code;
}




