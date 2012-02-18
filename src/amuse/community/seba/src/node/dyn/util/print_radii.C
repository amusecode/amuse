
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out indices (if any), masses, and radii (distances from the
//// origin) for an N-body system.
////
//// Usage: print_radii [OPTIONS] < input > output
////
//// Options:     
////		  -p    specify precision of output [6 sig. fig.]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char** argv)
{
    dyn *root, *ni;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "p:";

    int p = STD_PRECISION;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'p': p = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

    cerr.precision(p);
    int i = 0;

    while (root = get_dyn()) {

	for_all_daughters(dyn, root, ni) {

	    cerr << ++i << " "
		 << ni->get_index() << " "
	         << ni->get_mass() << " "
		 << abs(ni->get_pos()) << endl;

	}

	rmtree(root);
    }
}

#endif
