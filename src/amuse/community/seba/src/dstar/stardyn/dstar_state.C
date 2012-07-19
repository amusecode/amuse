#include "dstar_to_dyn.h"
#include "double_star.h"
#include "dyn.h"

#ifndef TOOLBOX

// Declaration can be found in: dstar_to_dyn.h

void print_binary_dstars(dyn* bi) {	// bi is center-of-mass node

    //  cerr<<"Nl= " << bi->n_leaves()<<endl;

    if (bi->n_leaves() == 2) {
	double_state bst = make_state(bi);
	if(bst.type != Unknown_Binary_Type) {
	    cerr << "                     "; 
	    put_state(bst);
	}
    }
}

double_state make_state(dyn* b) {

    double_state bst;

      //  if (has_dstar(b->get_parent())) {
      //if (has_dstar(b->get_oldest_daughter())) {
      //bst = make_state(dynamic_cast(double_star*,
      //				  b->get_parent()->get_starbase()));
      //  }
    if (has_dstar(b->get_oldest_daughter())) {
	bst = make_state(dynamic_cast(double_star*,
				      b->get_starbase()));
    }
    else if (b->get_oldest_daughter()->get_star_story()!=NULL &&
	     b->get_oldest_daughter()->get_binary_sister()
	                             ->get_star_story()!=NULL) { 

	//cerr<<"                     Has no dstar but is binary node."<<endl;

	dyn *bi = b->get_oldest_daughter();
	dyn *bj = b->get_oldest_daughter()->get_binary_sister();
	real M = bi->get_mass() + bj->get_mass();
	vec dx = bj->get_pos() - bi->get_pos();
	vec dv = bj->get_vel() - bi->get_vel();
	real mu = (bi->get_mass() * bj->get_mass()) / M;
	real E = mu*(0.5*dv*dv - M / abs(dx));

	kepler k;
	k.set_time(0);
	k.set_total_mass(M);
	k.set_rel_pos(dx);
	k.set_rel_vel(dv);

	k.set_circular_binary_limit(MAX_CIRC_ECC);

	k.initialize_from_pos_and_vel(true, false);

	real time = b->get_system_time();
	if (time > -VERY_LARGE_NUMBER)
	    time = b->get_starbase()->conv_t_dyn_to_star(time);
	else
	    time = 0;
	real sma = b->get_starbase()->
	                    conv_r_dyn_to_star(k.get_semi_major_axis());
	real ecc = k.get_eccentricity();
	real sma_c = sma*(1-pow(ecc, 2));

	dyn *pi=bi;
	dyn *si=bj;
	if (si->get_mass() > pi->get_mass()) {
	    pi = bj;
	    si = bi;
	}
	real q   = si->get_mass()/pi->get_mass();
	real Rl1 = 0.46*sma_c/pow(1+q, cnsts.mathematics(one_third));
	real Rl2 = 0.46*sma_c*pow(q/(1+q), cnsts.mathematics(one_third));

	star_state primary = make_star_state(pi);
	star_state secondary= make_star_state(si);
      
	binary_type type = Detached;

	if(ecc==0)  type = Synchronized;

	if (primary.radius>=Rl1) {
	    primary.class_spec[Rl_filling] = 1;
	    type = Semi_Detached;
	}
	if(secondary.radius>=Rl2) {
	    secondary.class_spec[Rl_filling] = 1;
	    if(!primary.class_spec[Rl_filling])
		type = Semi_Detached;
	    else
		type = Contact;
	}

	bst.identity = b->get_index();
	bst.time = time;
	bst.type = type;
	bst.semi = sma;
	bst.ecc  = ecc;
	bst.velocity  = 0;
	bst.total_mass = primary.mass + secondary.mass;

	bst.primary   = primary;
	bst.secondary = secondary;

    }
    else {
	//cerr << "    No stellar information found for: ";
	b->pretty_print_node(cerr);
	return bst;
    }

    return bst;
}

#else


main() {
  cerr <<"Hello"<<endl;
}


#endif

