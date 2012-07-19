
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add Plummer parameters to an input snapshot and write it out again.
//// Do not change the input data. See add_power_law.
////
//// Usage: add_plummer [OPTIONS] < input > output
////
//// Options:      
////		   -c    add comment [none]
////               -C    specify center [(0,0,0)]
////               -f    turn on/off friction [false]
////               -m/M  specify mass [1]
////               -a/R  specify scale [1]
////               -n    force interpretation of parameters in N-body units [no]
////
//// Written by Steve McMillan
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Aug/Sep 2001   Steve McMillan

// This function is just a special case of add_power_law, and is now
// handled as such, except that the optional friction flag applies only
// to Plummer (for now).

#include "dyn.h"

#ifndef TOOLBOX

void toggle_plummer_friction(dyn *b)
{
    // Toggle the current setting of kira_plummer_friction.

    if (find_qmatch(b->get_log_story(), "kira_plummer_friction")) {
	int f = getiq(b->get_log_story(), "kira_plummer_friction");
	f = (f == 0 ? 1 : 0);
	putiq(b->get_log_story(), "kira_plummer_friction", f);
    }
}

void add_plummer(dyn *b,
		 real coeff, real scale,
		 vec center,			// default = (0, 0, 0)
		 bool n_flag,			// default = false
		 bool verbose,			// default = false
		 bool fric_flag)		// default = false
{
    add_power_law(b, coeff, 0, scale, center, n_flag, verbose, false);
    if (fric_flag) toggle_plummer_friction(b);
}

#else

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    real mass = 1, scale = 1;
    vec center = 0;

    bool n_flag = false;
    bool fric_flag = false;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    const char *param_string = "a:c:C:::fm:M:nR:";

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.13 $", _SRC_)) != -1) {
	switch (c) {
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'C':	center = vec(atof(poparr[0]),
					atof(poparr[1]),
					atof(poparr[2]));
			break;
	    case 'f':	fric_flag = true;
			break;
	    case 'm':
	    case 'M':	mass = atof(poptarg);
			break;
	    case 'a':
	    case 'R':	scale = atof(poptarg);
			break;

	    case 'n':	n_flag = true;
	    		break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    dyn *b = get_dyn();
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);
    if (c_flag)	b->log_comment(comment);

    add_plummer(b, mass, scale, center, n_flag, true, fric_flag);
    put_dyn(b);
}

#endif
