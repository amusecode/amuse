
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out the energy of an N-body system.  Does not include 
//// the effects of any external tidal field.
////
//// Usage: energy [OPTIONS] < input > output
////
//// Options:   
////		  -e    specify softening parameter [0]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char **argv)
{
    dyn *root;
    int  e_flag = 0;
    real eps = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "e:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c) {
	    case 'e': e_flag = 1;
                      eps = atof(poptarg);
	              break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    // Loop over input until no more data remain.

    while ( (root = get_dyn()) != NULL) {

	real kinetic_energy, potential_energy;
	get_top_level_energies(root, eps*eps, potential_energy, kinetic_energy);

	cout << "Top-level energies (T, U, E):  "
	     << kinetic_energy << "  "
	     << potential_energy << "  "
	     << kinetic_energy + potential_energy 
	     << endl;
    }
}

#endif

// endof: energy.C
