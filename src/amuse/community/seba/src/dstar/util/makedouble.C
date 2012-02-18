//
//  mkdouble.C: construct a linked list of unit-mass nodes.
//
//	       Simon Portegies Zwart, July 1996
//

#include "node.h"
#include "double_star.h"
#include "main_sequence.h"
#include "dstar_to_dyn.h"

#ifndef TOOLBOX

void add_secondary(node* original, real mass_ratio) {

    node* primary = new node;
    node* secondary = new node;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    primary->set_mass(original->get_mass());
    secondary->set_mass(mass_ratio*original->get_mass());
    original->inc_mass(secondary->get_mass());

    // Naming convention:

    if (original->get_name() == NULL)
        if (original->get_index() >= 0) {
            char tmp[64];
            sprintf(tmp, "%d", original->get_index());
            original->set_name(tmp);
        }

    primary->set_name(original->get_name());
    secondary->set_name(original->get_name());
    strcat(primary->get_name(), "a");
    strcat(secondary->get_name(), "b");

   }

void mksecondary(node* b, real binary_fraction, real lower_limit) {

    // For now, use a flat distribution in secondary mass ratio.
    // Assume that the probability of a star being the primary of
    // a binary is independent of mass.

    real sum = 0;
    b->set_mass(0);

    for_all_daughters(node, b, bi) {
        sum += binary_fraction;
        if (sum >= 1) {
            sum -= 1;

            real mass_ratio = randinter(lower_limit, 1);        // Quick fix...
            add_secondary(bi, mass_ratio);

        }
        b->inc_mass(bi->get_mass());
    }
}

#else

local void  evolve_the_stellar_system(node* b, real time) {

      for_all_nodes(node, b, bi) {
         cerr<<bi<< " is_root?: "<<bi->is_root()<<endl;
         cerr<<bi<< " is_parent?: "<<bi->is_parent()<<endl;
         if (!bi->is_root())
            if (bi->get_parent()->is_root()) {
                cerr<<"binary "<<bi<<endl;
                bi->get_starbase()->evolve_element(time);
            }
         
      }
   }

void main(int argc, char ** argv) {

    bool A_flag = false;
    bool a_flag = false;
    bool E_flag = false;
    bool e_flag = false;
    bool reandom_initialization = false;
    int  n = 1;

    int id=1;
    real a_min = 1;
    real a_max = 1.e+6;
    real e_min = 0;
    real e_max = 1;
    real m_tot=1,r_hm=1, t_hc=1;
    stellar_type type = Main_Sequence;
    binary_type bin_type = Detached;
    real binary_fraction = 1.0;
    real lower_limit = 0.0;
    int random_seed = 0;
    char seedlog[64];
    real t_start = 0;
    real t_end   = 100;
    extern char *poptarg;
    int c;
    const char *param_string = "A:a:E:e:f:l:M:n:R:S:s:T:t:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c) {
            case 'A': A_flag = true;
		      a_max = atof(poptarg);
                      break;
            case 'a': a_flag = true;
                      a_min = atof(poptarg);
                      break;
            case 'E': E_flag = true;
                      e_max = atof(poptarg);
                      break;
            case 'e': e_flag = true;
                      a_min = atof(poptarg);
                      break;
            case 'M': m_tot = atof(poptarg);
                      break;
	    case 'n': n = atoi(poptarg);
		      break;
            case 'r': r_hm = atof(poptarg);
                      break;
            case 't': t_hc = atof(poptarg);
                      break;
            case 'T': t_end = atof(poptarg);
                      break;
            case 'S': type = (stellar_type)atoi(poptarg);
                      break;
            case 'f': binary_fraction = atof(poptarg);
                      break;
            case 'l': lower_limit = atof(poptarg);
                      break;
            case 's': random_seed = atoi(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (n <= 0) err_exit("mknodes: N > 0 required!");
    int actual_seed = srandinter(random_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

    if (A_flag || a_flag || E_flag || E_flag)
       reandom_initialization = true;

    // Create flat tree 
//    node *root  = mkstar(n, t_start, type, mf);
    node *root  = mknode(n);
    root->log_history(argc, argv);
    root->log_comment(seedlog);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);
    mksecondary(root, binary_fraction, lower_limit);
    addstar(root, t_start, type);

    cerr<<"calling adddouble...\n";
    adddouble(root, t_start, bin_type, reandom_initialization, 
	      a_min, a_max, e_min, e_max);

    put_node(root);
//	Test pointer structure
    cerr<<"est pointer structure"<<endl;
    node *b = root->get_oldest_daughter();
    starbase *s = b->get_starbase();
    star *st     = (star*)b->get_starbase();
    cerr<<"binary+s:" << b<<" "<<s<<" " <<st<<endl;
    cerr<<"a e m"<<s->get_effective_radius() <<" "
               <<s->get_total_mass()<<endl;
    cerr<<"p, s"<< st->get_primary()<<" "<<st->get_secondary() <<endl;
    cerr<<"M R, m r"<<st->get_primary()->get_total_mass() << " " 
                <<st->get_primary()->get_effective_radius() << " "
                <<st->get_secondary()->get_total_mass() << " " 
                <<st->get_secondary()->get_effective_radius()<< endl;
    cerr<<"companion p, s: "<< 
         st->get_secondary()->get_companion()<<" " 
       <<st->get_primary()->get_companion()<<endl;
   cerr<<"companion p, s: "<<st->get_companion(st->get_secondary())<<" "
                           <<st->get_companion(st->get_primary())<< endl; 
   
    evolve_the_stellar_system(root, t_end);
    put_node(root);
}

#endif

/* end of: mkbinaries.c */
