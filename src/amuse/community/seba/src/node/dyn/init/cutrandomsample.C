
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Create a linked list of equal-mass nodes.
////
//// Usage:  cutrandomsample [OPTIONS]
////
//// Options:
////       -m       specify total mass [1]
////       -n       specify number of nodes [1]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//	      Steve McMillan, July 1996

#include "dyn.h"

#ifndef TOOLBOX
#else

local real get_primary_mass(dyn *b, bool S_flag) {

  real m = -1;
  if(b->is_leaf()) {
    m = b->get_mass();
  }
  else if(!S_flag) {
    for_all_leaves(dyn, b, bi) {
      m = Starlab::max(m, bi->get_mass());
    }
  }

  return m;
}

local dyn* cut_randomsample(dyn * b, int nd, real mmin, bool S_flag) {

  int nb=0;
  for_all_daughters(dyn, b, bi) {
    if(get_primary_mass(bi, S_flag)>mmin) 
      nb++;
  }
  real paccept = 1;
  if(nd>=0)
    paccept = real(nd)/real(nb);
  PRC(nb);PRC(nd);PRL(paccept);

  real mtot=0;
  dyn *d = new dyn();
  dyn *bj=NULL;
  int ndone = 0;
  for_all_daughters(dyn, b, bi) {
    ndone++;
    
    real r = randinter(0, 1);
    //    PRC(paccept);PRL(r);
    if(get_primary_mass(bi, S_flag)>mmin && paccept>=r) {
	bi->pretty_print_node(cerr);

      if (bi != NULL && bi->get_younger_sister()!= NULL) {
	bj = bi->get_younger_sister();
	detach_node_from_general_tree(bi);

	//	bi->pretty_print_node(cerr);
	//	PRI(3);PRL(bi->get_mass());

	add_node(bi, d);
	mtot += bi->get_mass();
	if(bj->get_elder_sister()!=NULL) {
	  bi = bj->get_elder_sister();
	}
	else {
	  cerr << "no elder sister from younger sister in cutrandomsample"
	       <<endl;
	  bi = bj;
	  //  exit(-1);
	}
      }
    }
  }
  cerr << "\n# done: " << ndone << endl;


  //  d->set_mass(b->get_mass());
  //d->set_mass(mtot);
  PRC(mtot);
  PRC(b->get_mass());
  //mtot=0;
  //for_all_leaves(dyn, d, di) {
  //   mtot += di->get_mass();
  //}
  PRL(mtot);
  d->set_mass(mtot);

  d->pretty_print_node(cerr);
  cerr << endl;

  return d;
}

int main(int argc, char ** argv)
{
    bool S_flag = false;
    bool s_flag = false;
    bool leave_unscaled = true;
    bool s = false;
    real mmin = 0;
    int nd = -1;

    check_help();

    int  input_seed, actual_seed;
    extern char *poptarg;
    int c;
    const char *param_string = "m:n:Ss:u";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {

	    case 'm': mmin = atof(poptarg);
		      break;
	    case 'n': nd = atoi(poptarg);
		      break;
   	    case 'u': leave_unscaled = !leave_unscaled;
		      break;
	    case 'S': S_flag = !S_flag;
		      break;
	    case 's': s_flag = !s_flag;
	              input_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    dyn *b = get_dyn();

    b->log_history(argc, argv);

    if (s_flag == false) input_seed = 0;	// default
    actual_seed = srandinter(input_seed);

    char seedlog[64];
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    if(!leave_unscaled) {
      mmin = b->get_starbase()->conv_m_star_to_dyn(mmin);
      PRL(mmin);
    }

    dyn *d = cut_randomsample(b, nd, mmin, S_flag);
    putrq(b->get_log_story(), "initial_mass", b->get_mass());
    if(rmq(b->get_log_story(), "initial_mass")) {
      cerr << "initial_mass" << endl;
    }
    d->set_log_story(b->get_log_story());

    real mtot = d->get_mass();
    PRL(mtot);
      d->get_starbase()->print_stellar_evolution_scaling(cerr);
      real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
      real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
      d->get_starbase()->set_stellar_evolution_scaling(mtot,
						     old_r_vir,
						     old_t_vir);
      d->get_starbase()->print_stellar_evolution_scaling(cerr);


    put_dyn(d);
    return 0;
}

#endif

#if 0
    //    if(s) {
    real mtot = 0;
    real dmtot = 0;
    for_all_daughters(dyn, d, di) {
      mtot += di->get_mass();
      dmtot += di->get_starbase()->conv_m_star_to_dyn(di->get_mass());
    }
    PRC(mtot);PRL(dmtot);

    real eps = 0;
      bool c_flag, e_flag, e, m_flag, q_flag, r_flag;
      c_flag = e_flag = m_flag = q_flag = r_flag = true;
      real m = 1;
      real q = 0.5;
      real r = 1;
      //  scale(b, eps, c_flag, e_flag, e, m_flag, m, q_flag, q, r_flag, r);

      d->get_starbase()->print_stellar_evolution_scaling(cerr);
      //      real mtot = b->get_starbase()->conv_r_star_to_dyn(mtot);
      real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
      real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
      d->get_starbase()->set_stellar_evolution_scaling(mtot/dmtot,
						     old_r_vir,
						     old_t_vir);
      d->get_starbase()->print_stellar_evolution_scaling(cerr);
      //    }
#endif
