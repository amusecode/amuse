
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Split the specified particles in input snapshot into binaries having
//// specified orbital properties.
////
//// Naming convention: particle X is split into components X1 and X2.
//// Splitting may be applied recursively to produce multiple systems,
//// e.g.
////
////        split_particles -l 1 -a 0.1 -e 0 -l 1b -a 0.01 -e 0.5
////
//// creates a triple system (11, (121, 122)).
////
//// Unlike mkscat, structure is never implicit, but must be created
//// from the top down.
////
//// Usage:  split_particles  [-s #] -l l1 -a a1 -e e1 -E E1 -q q1  -l l2 (etc.)
////
//// Options:
////             -s    specify random seed [random from system clock]
////             -i    use "a" and "b" instead of "1" and "2" [false]
////             -l    specify label of next particle to split [no default]
////             -a    specify semi-major axis of current binary [0.1]
////             -e    specify eccentricity of current binary [random]
////             -E    specify energy of current binary [0]
////             -q    specify mass ratio (secondary[b]/primary[a]) of
////                   the current binary [1]
////
//// The "-E" option takes precedence over the "-a" option.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// Only the main routine and check_and_initialize (at end) differ
// from the functions in split_particle.
// The repeated functions will end up in a library someday...

// Steve McMillan, August 1996
//		   June 2002

#include "dyn.h"

#ifdef TOOLBOX

// split_particle: Split the specified node into a binary with the specified
//                 parameters.  All unspecified orbit elements are chosen
//                 randomly.  Newly created leaves have names "n1" and "n2",
//                 where "n" is the name of the leaf being split.

local bool split_particle(dyn* bi, real ecc, real sma,
			  real mass_ratio, bool iname)
{
    if (bi->get_oldest_daughter() != NULL) {
	cerr << "Can't split a binary node!\n";
	return false;
    } else if (sma <= 0) {
	cerr << "Must specify semi-major axis > 0.\n";
	return false;
    } else if (ecc < 0 || ecc >= 1) {
	cerr << "Must specify eccentricity in [0,1).\n";
	return false;
    } else if (mass_ratio <= 0 || mass_ratio > 1) {
	cerr << "Must specify mass ratio in (0,1]!\n";
	return false;
    }

    // Update the binary tree structure:

    dyn* d1 = new dyn;
    dyn* d2 = new dyn;

    bi->set_oldest_daughter(d1);

    d1->set_parent(bi);
    d2->set_parent(bi);

    d1->set_younger_sister(d2);
    d2->set_elder_sister(d1);

    // Set new masses and radii:

    real m_total = bi->get_mass();
    real m1 = m_total / (1 + mass_ratio);
    real m2 = m1 * mass_ratio;

    // By convention, the first component has larger mass.

    if (m1 < m2) {
	real temp = m1;
	m1 = m2;
	m2 = temp;
    }

    d1->set_mass(m1);
    d2->set_mass(m2);

    // Set internal orbital elements:

    kepler k;

    real peri = 1; // Default value (unimportant unless ecc = 1).
    if (ecc == 1) peri = 0;

    // For now, binary phase is random.

    real mean_anomaly = randinter(-PI, PI);

    make_standard_kepler(k, 0, m_total, -0.5 * m_total / sma, ecc,
			 peri, mean_anomaly);

    set_random_orientation(k);

    d1->set_pos(-m2 * k.get_rel_pos() / m_total);
    d1->set_vel(-m2 * k.get_rel_vel() / m_total);

    d2->set_pos(m1 * k.get_rel_pos() / m_total);
    d2->set_vel(m1 * k.get_rel_vel() / m_total);

    // Naming convention:

    if (bi->get_name() == NULL)
	if (bi->get_index() >= 0) {
	    char tmp[64];
	    sprintf(tmp, "%d", bi->get_index());
	    bi->set_name(tmp);
	}

    d1->set_name(bi->get_name());
    d2->set_name(bi->get_name());
    if (iname) {
      strcat(d1->get_name(), "1");
      strcat(d2->get_name(), "2");
    } else {
      strcat(d1->get_name(), "a");
      strcat(d2->get_name(), "b");
    }

    return true;
}

// leaf_with_label: Return a pointer to the leaf with the specified name.

#include <string.h>

local dyn* leaf_with_label(dyn* b, char* label)
{
    for_all_leaves(dyn, b, bi) {
	if (bi->get_name() != NULL) {
	    if (strcmp(bi->get_name(), label) == 0)
		return bi;
	} else if (bi->get_index() >= 0) {
	    char id[64];
	    sprintf(id, "%d", bi->get_index());
	    if (strcmp(id, label) == 0)
		return bi;
	}
    }
    return NULL;
}

local void scale_energy(real energy, int N, real M, real E)
{
    // Scale energy to kT = -(2/3) E/N.

    energy *= (-E / (1.5*N));
}

//----------------------------------------------------------------------

local void check_and_initialize(dyn* b, char* label,
				real eccentricity, real energy,
				real semi_major_axis, real mass_ratio,
				bool iname, bool verbose = true)
{
    if (eccentricity < 0)
	eccentricity = sqrt(randinter(0,1));	// Thermal distribution

    if (energy < 0) {

	// Working with energy, rather than semi-major axis.

	// Look in the dyn story to see if a total system energy is known.
	// If it is, scale 'energy' accordingly.

	// Thus, if the total system energy is known, 'energy' is specified
	// in "kT" units.  If not, 'energy' is the binary kinetic energy.

	const char* energy_string = "initial_total_energy";

	if (find_qmatch(b->get_log_story(), energy_string))
	    scale_energy(energy, b->n_daughters(), b->get_mass(),
			 getrq(b->get_log_story(), energy_string));
	else
	    cerr << "split_particle: Using unscaled energy limits.\n";

    }

    // Locate the particle to be split.

    dyn* bi = leaf_with_label(b, label);

    if (bi == NULL) {

	cerr << "split_particles: warning: particle \"" << label
	     << "\" not found." << endl;

    } else {

	// Determine the semi-major axis if necessary.

	if (energy < 0)
	    semi_major_axis = (1 - mass_ratio) * bi->get_mass()
	    			* mass_ratio * bi->get_mass()
				    	/ (-2 * energy);

	// If a real error occurs, terminate the data flow.

	if (!split_particle(bi, eccentricity, semi_major_axis,
			    mass_ratio, iname))
	    err_exit("split_particles: fatal error");

	if (verbose) {
	    cerr << "split particle \"" << label << "\"" << endl;
	    PRC(semi_major_axis); PRC(eccentricity); PRL(mass_ratio);
	}
    }
}

#define DEFAULT_ENERGY		 0.0
#define DEFAULT_ECC		-1.0
#define DEFAULT_SMA		 0.1
#define DEFAULT_Q		 1.0

int main(int argc, char ** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.11 $", _SRC_);

    real energy = DEFAULT_ENERGY;
    real eccentricity = DEFAULT_ECC;
    real semi_major_axis = DEFAULT_SMA;
    real mass_ratio = DEFAULT_Q;

    bool iname = true;

    char label[64];
    label[0] = '\0';
    bool previous_label = false;

    int random_seed = 0;
    char seedlog[64];

    dyn* b;
    b = get_dyn();

    b->log_history(argc, argv);

    // Parse the command line by hand, modifying the system as we go.

    int i = 0;
    while (++i < argc)
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {

		case 'a': semi_major_axis = atof(argv[++i]);
			  break;

		case 'e': eccentricity = atof(argv[++i]);
			  break;

		case 'E': energy = atof(argv[++i]);
			  break;
			  
		case 'i': iname = !iname;
			  break;

		case 'l': if (previous_label)

		    	      // Initialize the previous binary with the
		    	      // present parameters.

		    	      check_and_initialize(b, label,
						   eccentricity, energy,
						   semi_major_axis, mass_ratio,
						   iname);

			  strcpy(label, argv[++i]);

			  // Options: we can either reinitialize to defaults,
			  // or we can pick up the previously defined values...
			  // Do the former, at least for now.

			  energy = DEFAULT_ENERGY;
			  eccentricity = DEFAULT_ECC;
			  semi_major_axis = DEFAULT_SMA;
			  mass_ratio = DEFAULT_Q;

			  previous_label = true;
			  break;

		case 'q': mass_ratio = atof(argv[++i]);
			  break;

		case 's': {
		    	  random_seed = atoi(argv[++i]);
			  int actual_seed = srandinter(random_seed);
			  sprintf(seedlog,
				  "       random number generator seed = %d",
				  actual_seed);
			  b->log_comment(seedlog);
		          }
			  break;

		default:  get_help();
			  exit(1);
	    }

    if (previous_label)
	check_and_initialize(b, label,
			     eccentricity, energy,
			     semi_major_axis, mass_ratio,
			     iname);

    put_dyn(b);
    return 0;
}

#endif
