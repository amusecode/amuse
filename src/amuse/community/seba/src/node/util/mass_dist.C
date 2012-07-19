
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute a mass histogram from the input snapshot(s).
////
//// Usage: mass_dist [OPTIONS]
////
//// Options:
////             -b    specify number of bins [25; 0 is OK]
////             -v    toggle verbose mode [off]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "node.h"

#ifdef TOOLBOX

local void get_mass_dist(node* b, int num_of_bins, bool verbose)
{
    int n = b->n_leaves();
    real * mass = new real[n];

    real min_mass = VERY_LARGE_NUMBER;
    real max_mass = -VERY_LARGE_NUMBER;
    real mean_mass = 0;

    int i = 0;
    int j = 0;

    for_all_leaves(node, b, bb) {

        mass[i] = bb->get_mass();
	if ( mass[i] <= min_mass ) min_mass = mass[i];
	if ( mass[i] >= max_mass ) max_mass = mass[i];
	mean_mass += mass[i];

	mass[i] = log10(mass[i]);		      // NB mass --> log10(mass)
	i++;

    }

    if (verbose) {
	PRL(min_mass);
	PRL(max_mass);
	mean_mass /= i;
	PRL(mean_mass);
    }

    if (num_of_bins <= 0) return;

    min_mass = log10(min_mass);
    max_mass = log10(max_mass);

    // Analyze the mass distribution...

    real range = max_mass - min_mass;
    real width = range / (float)num_of_bins;

    real * bin = new real[num_of_bins + 1];
    int * n_in_bin = new int[num_of_bins];

    bin[0] = min_mass;
    for ( j = 1; j <= num_of_bins; j++ ) {
	bin[j] = bin[0] + (real)j * width;
    }

    // Binning:
    //
    //	bin[0]	bin[1]	bin[2]	 . . . . .	bin[nb-1]	bin[nb]
    //	  |	  |	  |			  |		  |
    //	     n[0]    n[1]    n[2]   ...  n[nb-2]	n[nb-1]

    for ( i = 0; i < n; i++ )
	for ( j = 0; j < num_of_bins; j++ )
	    if ( mass[i] >= bin[j] && mass[i] <= bin[j+1] ) {
	        n_in_bin[j]++;
	        break;
	    }

    // Current output is:	log10(m)	log10(dN/dm).
    //
    // Binning is linear in log10(m); m is mass at logarithmic bin center.
    // dN is normalized to N_total = 1.

    for ( j = 0; j < num_of_bins; j++ ) {
	real m = 0.5 * (bin[j] + bin[j+1]);
	real dm = pow(10.0, bin[j+1]) - pow(10.0, bin[j]);
	cerr << m << "  ";
	if (n_in_bin[j] > 0)
	    cerr << log10(n_in_bin[j]/(n*dm));
	else
	    cerr << "---";
	cerr << endl;
    }

    delete mass;
    delete bin;
    delete n_in_bin;
}

int main(int argc, char ** argv)
{
    node *b;
    int num_of_bins = 25;
    bool verbose = false;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "b:v";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
        switch(c) {

	case 'b': num_of_bins = atoi(poptarg);
	          break;
	case 'v': verbose = !verbose;
	          break;
	case '?': params_to_usage(cerr, argv[0], param_string);
		  get_help();
		  exit(1);
	}

    while ((b = get_node())) {
        get_mass_dist(b, num_of_bins, verbose); 
        delete b;
    }
    return 0;
}
#endif   

