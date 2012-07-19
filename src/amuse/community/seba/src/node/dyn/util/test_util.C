
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// -----------------------------------------------------------------
//// Handy "placeholder" program to test new tools.  Always compiled,
//// but the content may vary according to needs.
//// -----------------------------------------------------------------
////
//// Print the half-mass radii corresponding to various mass ranges.
//// Masses are specified on the command line as fractions or percentiles.
//// The default is fractions, but if a mass > 1 is detected we switch
//// to percentiles.  The input masses are scaled if necessary and sorted,
//// and a leading 0 and trailing 100 are added.  No checks for duplicates
//// or percentiles out of range are performed.  Thus, if no masses are
//// specified, we print the half-mass radius for all stars.
////
//// Usage: test_util mass1 mass2 mass3 ... < input > output
////
//// Written by Steve McMillan, July 2005.
////
//// Report bugs to steve@physics.drexel.edu.

#include "dyn.h"
#include <vector>
#include <algorithm>

main(int argc, char *argv[])
{
    // Read and sort the list of masses.

    vector<real> mass;
    bool percent = false;
    mass.push_back(0);
    for (int i = 1; i < argc; i++) {
	real m = atof(argv[i]);
	mass.push_back(m);
	if (m > 1) percent = true;	// default interpretation is fraction,
					// >1 ==> switch to percentile
    }
    if (percent)
	for (int i = 0; i < mass.size(); i++) mass[i] = 0.01*mass[i];
    sort(mass.begin(), mass.end());
    mass.push_back(1.0);

    dyn *b;				// root node
    while (b = get_dyn()) {

	// Compute the quartiles (not actually printed, but stored
	// in the root dyn story) for each mass range.

	cerr << b->get_system_time() << " ";

	for (int i = 1; i < mass.size(); i++) {

	    if (i == 1)
		set_lagr_cutoff_mass(b, mass[i-1], mass[i]);	// sort masses
	    else
		reset_lagr_cutoff_mass(b, mass[i-1], mass[i]);	// no sort

	    real rhalf = compute_lagrangian_radii(b, 0.5,
						  false,	// not verbose
						  4);		// mass filter

	    // Alternatively, if several Lagrangian radii are needed:
	    //
	    // real lagr_array[3] = {0.1, 0.5, 0.9};
	    // compute_lagrangian_radii(b,
	    //			     lagr_array, 3,
	    //			     false,		// not verbose
	    //			     4);		// mass filter
	    // rhalf = lagr_array[1];

	    cerr << rhalf<< " ";
	}
	cerr << endl;

	rmtree(b);
    }
}
