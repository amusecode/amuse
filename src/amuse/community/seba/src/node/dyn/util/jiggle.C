
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add a small, energy-preserving perturbation to all top-level
//// velocities, leaving positions unchanged.
////
//// Usage: jiggle [OPTIONS] < input > output
////
//// Options:     
////              -c    add a comment to the output snapshot [false]
////              -f    specify (fractional) scale of perturbation [0.01]
////              -s    specify random seed [take from system clock]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Jun 2000   Steve McMillan

#include "dyn.h"

#ifdef TOOLBOX

local void jiggle(dyn * b, real f)
{
    if (f >= 1) f = 1;

    for_all_daughters(dyn, b, bi) {

	vec vel = bi->get_vel(), norm;

	while(1) {
	    real costh = randinter(-1, 1);
	    real sinth = sqrt(Starlab::max(0.0, 1-costh*costh));
	    real phi = randinter(0, 2*M_PI);
	    vec rvec = vec(sinth*cos(phi), sinth*cos(phi), costh);
	    norm = rvec ^ vel;
	    if (abs(norm) > 0.01*abs(vel)) break;    // norm not too close to
						     // being parallel to vel
	}

	// Now norm is a random vector perpendicular to vel.

	real v = abs(vel);
	vel = sqrt(1 - f*f) * vel + f * v * norm/abs(norm);
	
	bi->set_vel(vel);
    }
}

main(int argc, char ** argv)
{
    bool c_flag = false, s_flag = false;
    char *comment;
    real f = 0.01;
    int input_seed, actual_seed;
    char seedlog[128];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:f:s:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.10 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'f': f = atof(poptarg);
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if (!s_flag) input_seed = 0;                         	// default
    actual_seed = srandinter(input_seed);

    cerr << "jiggle: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d", actual_seed);

    dyn *b;

    while (b = get_dyn()) {

        b->log_history(argc, argv);
        if (c_flag) b->log_comment(comment);
	b->log_comment(seedlog);

        jiggle(b, f);
	put_dyn(b);
	rmtree(b);
    }
}

#endif
