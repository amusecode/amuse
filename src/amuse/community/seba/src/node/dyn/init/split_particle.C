
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// split_particle.C: Split a specified particle into a binary.
//
//	Mostly stolen from sdyn/init/mkscat.C
//
//	Steve McMillan, August 1996

#include "dyn.h"

#ifdef TOOLBOX

// split_particle: Split the specified node into a binary with the specified
//                 parameters.  All unspecified orbit elements are chosen
//                 randomly.  Newly created nodes have names "n1" and "n2",
//                 where "n" is the name of the node being split.

local bool split_particle(dyn* bi, real ecc, real sma, real mass_ratio)
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
    strcat(d1->get_name(), "1");
    strcat(d2->get_name(), "2");

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

void main(int argc, char ** argv)
{
    real energy = 0.0;
    real eccentricity = -1;
    real semi_major_axis = 0.1;
    real mass_ratio = 1.0;
    int random_seed = 0;
    char seedlog[64];
    char label[64];

    label[0] = '\0';

    extern char *poptarg;
    int c;
    const char *param_string = "a:e:E:l:s:o:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {

	    case 'a': semi_major_axis = atof(poptarg);
		      break;
	    case 'e': eccentricity = atof(poptarg);
		      break;
	    case 'E': energy = atof(poptarg);
		      break;
	    case 'l': strcpy(label, poptarg);
		      break;
	    case 'q': mass_ratio = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (label[0] == '\0')
	err_exit("must specify a particle label");

    dyn* b;
    b = get_dyn();

    // Bookkeeping.

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog,
	    "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    if (eccentricity < 0)
	eccentricity = sqrt(randinter(0,1));	// Thermal distribution

    if (energy < 0) {

	// Working with energy, rather than semi-major axis.

	// Look in the dyn story to see if a total system energy is known.
	// If it is, scale 'energy' accordingly.

	// Thus, if the total system energy is known, 'energy' is specified
	// in "kT" units.  If not, 'energy' is the binary kinetic energy.

	char* energy_string = "initial_total_energy";

	if (find_qmatch(b->get_log_story(), energy_string))
	    scale_energy(energy, b->n_daughters(), b->get_mass(),
			 getrq(b->get_log_story(), energy_string));
	else
	    cerr << "split_particle: Using unscaled energy limits.\n";

    }

    // Locate the particle to be split.

    dyn* bi = leaf_with_label(b, label);

    if (bi == NULL) {

	// Do nothing in this case, but don't stop the flow of data.

	cerr << "Warning: particle " << label << " not found.\n";
	put_dyn(b);

    } else {

	// Determine the semi-major axis if necessary.

	if (energy < 0)
	    semi_major_axis = (1 - mass_ratio) * bi->get_mass()
	    			* mass_ratio * bi->get_mass()
				    	/ (-2 * energy);

	// If a real error occurs, terminate the data flow.

	if (split_particle(bi, eccentricity, semi_major_axis, mass_ratio))
	    put_dyn(b);
    }
}

#endif
