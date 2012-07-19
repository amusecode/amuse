
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Double the mass of one or more stars in a snapshot to
//// approximate the presence of a binary.
////
//// Usage:  makeheavystar [OPTIONS]
////
//// Options:
////           -f    specify fraction of stars to double [0.1]
////           -l    specify lower mass limit (unused) [1]
////           -s    specify random seed [random from system clock]
////
//// Written by Kimberly Engle.
////
//// Report bugs to starlab@sns.ias.edu.

// Kimberly Engle, November 1997

#include "node.h"

// All functions are node -- why is this here??  (Steve, 5/03)

#ifdef TOOLBOX

// double_mass: Double the mass of a star to approximate a binary.
//		  Only the mass is modified at this stage.

local void double_mass(node* original, real mass_ratio)
{
    // Set new masses.

    original->set_mass(mass_ratio*original->get_mass());
    //    putiq(original->get_log_story(), "mass_doubled", 1);
}

local void makeheavystar(node* b, real fraction_doubled, real lower_limit,
		         real mass_ratio)
{
 
    real sum = 0;
    b->set_mass(0);

    PRI(6); putrq(b->get_log_story(),
		  "fraction_of_masses_doubled", fraction_doubled);

    for_all_daughters(node, b, bi) {
	sum += fraction_doubled;
	if (sum >= 1) {
	    sum -= 1;

	    double_mass(bi, mass_ratio);
            putrq(bi->get_log_story(), "mass_doubled", 1);
	}
	else putiq(bi->get_log_story(), "mass_doubled", 0);

	b->inc_mass(bi->get_mass());
    }

    // Write essential model information to root dyn story.

    putrq(b->get_log_story(), "initial_mass", b->get_mass());

    // We have probably broken any virialization.  Unset flage here.

    rmq(b->get_log_story(), "initial_rvirial");
    rmq(b->get_log_story(), "initial_rtidal_over_rvirial");
}

int main(int argc, char ** argv)
{
    real fraction_doubled = 0.1;
    real lower_limit = 1;
    int random_seed = 0;
    real mass_ratio = 2;
    //    char seedlog[64];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "f:l:q:s:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.11 $", _SRC_)) != -1)
	switch(c) {

	    case 'f': fraction_doubled = atof(poptarg);
		      break;
	    case 'l': lower_limit = atof(poptarg);
		      break;
	    case 'q': mass_ratio = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}

    if (fraction_doubled < 0 || fraction_doubled > 1)
	err_exit("makeheavystar: Illegal doubling fraction");
    //    if (lower_limit <= 0 || lower_limit > 1)
    //	err_exit("makeheavystar: Illegal mass ratio limit");

    node *b;
    b = get_node();

    b->log_history(argc, argv);

    //    int actual_seed = srandinter(random_seed);

    //    sprintf(seedlog,
    //            "       random number generator seed = %d",actual_seed);
    //    b->log_comment(seedlog);

    makeheavystar(b, fraction_doubled, lower_limit, mass_ratio);

    put_node(b);
    return 0;
}
#endif
