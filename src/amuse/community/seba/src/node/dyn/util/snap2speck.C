
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Dump out an N-body system in a dumb format suitable for digestion
//// by partiview.
////
//// Usage: snap2speck [OPTIONS] < input > output
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
    dyn *root;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "p:";

    int p = STD_PRECISION;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'p': p = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

    cout.precision(p);

    while (root = get_dyn()) {

	real n = 0;
	for_all_daughters(dyn, root, ni) n += 1;

	cout << "datavar 0 colorb_v" << endl;
	cout << "datavar 1 lum" << endl;

	for_all_daughters(dyn, root, ni)

	    cout << ni->get_pos() << " "
		 << ni->get_index()/n << " "	// colors 0 to 1
		 << ni->get_mass()
		 << endl;

	rmtree(root);
    }
}

#endif
