
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct a Plummer model, with a spatial or mass cut-off to ensure
//// finite radius.  The new model system is written to standard output.
//// The model system is shifted to place its center of mass at rest at
//// the origin of coordinates.  Unscaled systems will be in approximate
//// virial equilibrium, based on the continuum limit.
////
//// Usage:   makeplummer [OPTIONS]
////
//// Options:
////          -c    add a comment to the output snapshot [false]
////          -C    output data in 'col' format [no]
////          -i    number the particles sequentially [don't number]
////          -m    specify mass cutoff (for finite radius) [0.999]
////          -n    specify number of particles [no default]
////          -o    echo value of random seed [don't echo]
////          -r    specify radius cutoff [22.804 for default mass cutoff]
////          -R    toggle reshuffle of particles to remove correlation
////                between index and distance from cluster center [true]
////          -s    specify random seed [random from system clock]
////          -u    leave unscaled [scale to E=-1/4, M = 1, R = 1]
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//............................................................................
//   version 1:  July 1989   Piet Hut               email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//   version 2:  Dec  1992   Piet Hut  --  adapted to the new C++-based Starlab
//   version 3:  Sep  2004   Steve McMillan -- corrected long-standing bug!
//............................................................................
//    litt: S.J. Aarseth, M. Henon and R. Wielen (1974),
//          Astron. and Astrophys. 37, p. 183.
//............................................................................

#include "dyn.h"

#ifdef TOOLBOX

#define  MFRAC_DEFAULT  0.999                 // radial cutoff               
#define  RFRAC_DEFAULT  22.8042468            // corresponding radial cutoff 
                                              //   (not used right now)      
#define  SEED_STRING_LENGTH  60

local void  makeplummer(dyn * b, int n, real mfrac, real rfrac, bool u_flag)
{
    int  i;

    real  partmass;		   // equal mass of each particle	     
    real  radius;		   // absolute value of position vector      
    real  velocity;		   // absolute value of velocity vector      
    real  theta, phi;		   // direction angles of above vectors      
    real  x, y;		           // for use in rejection technique         
    real  scalefactor;             // for converting between different units 
    real  inv_scalefactor;         // inverse scale factor                   
    real  sqrt_scalefactor;        // sqare root of scale factor             
    real  mrfrac;                  // m( rfrac )                             
    real  m_min, m_max;            // mass shell limits for quiet start      
    dyn * bi;

    b->set_mass(1.0);
    partmass = 1.0 / ((real) n);

//  Calculating the coordinates is easiest in STRUCTURAL units;
//  conversion to VIRIAL units will be performed below.
//
//  Recipe for scaling to the proper system of units:
//
//  Since G = M = 1, if we switch from a coordinate system with
//  length unit  r_old  to a coordinate system with length unit  r_new ,
//  the length units simply scale by a factor  C = r_new / r_old .
//  Consequently, the coordinate values of physical quantities
//  such as positions should transform inversely to maintain the same
//  coordinate-invariant meaning. Similarly, the square of the velocities
//  should transform inversely proportional to the positions,
//  since  GM = 1  (cf. a relation such as  v*v = G*M/r ).
//
//  Thus, if
//              r_unit(new) = C * r_unit(old),
//  then
//              pos(new) = (1/C) * pos(old)
//  and
//              vel(new) = sqrt(C) * vel(old).

    scalefactor = 16.0 / (3.0 * PI);
    inv_scalefactor = 1.0 / scalefactor;
    sqrt_scalefactor = sqrt( scalefactor );

    // Now convert  rfrac  into an equivalent mfrac, if necessary:

    if (rfrac < VERY_LARGE_NUMBER) { // Note: the following powers can cause
	                             // problems on some compilers...

	rfrac *= scalefactor;        // Convert from VIRIAL to STRUCTURAL units
	mrfrac = rfrac*rfrac*rfrac / pow(1.0 + rfrac*rfrac, 1.5);
	if (mrfrac < mfrac)
	    mfrac = mrfrac;          // mfrac = min(mfrac, m(rfrac))
    }

    // Now construct the individual particles.

    for (i = 0, bi = b->get_oldest_daughter(); i < n;
         i++, bi = bi->get_younger_sister()) {

	bi->set_mass(partmass);

	// The position coordinates are determined by inverting the cumulative
	// mass-radius relation, with the cumulative mass drawn randomly from
	// [0, mfrac]; cf. Aarseth et al. (1974), eq. (A2).

        m_min = (i * mfrac)/n;
        m_max = ((i+1) * mfrac)/n;
	real rrrr = randinter(m_min, m_max);
	radius = 1 / sqrt( pow (rrrr, -2.0/3.0) - 1);

	// Note that this procedure arranges the particles approximately
	// in order of increasing distance fro the cluster center, which may
	// sometimes be desirable, and sometimes not.  Extra option added
	// by Steve (7/99) to allow radii to be randomized (reshuffled).

	theta = acos(randinter(-1.0, 1.0));
	phi = randinter(0.0, TWO_PI);

        bi->set_pos(vec(radius * sin( theta ) * cos( phi ),
			radius * sin( theta ) * sin( phi ),
			radius * cos( theta )));

	// The velocity coordinates are determined using von Neumann's
	// rejection technique, cf. Aarseth et al. (1974), eq. (A4,5).
	// First we take initial values for x, the ratio of velocity and
	// escape velocity (q in Aarseth et al.), and y, as a trick to
	// enter the body of the while loop.

	x = 0.0;
	y = 0.1;

	// Then we keep spinning the random number generator until we find a
	// pair of values (x,y), so that y < g(x) = x*x*pow( 1.0 - x*x, 3.5).
	// Whenever an y-value lies above the g(x) curve, the (x,y) pair is
	// discarded, and a new pair is selected. The value 0.1 is chosen as
	// a good upper limit for g(x) in [0,1] : 0.1 > max g(x) = 0.092 for
	// 0 < x < 1.

	while (y > x*x*pow( 1.0 - x*x, 3.5)) {
	    x = randinter(0.0,1.0);
	    y = randinter(0.0,0.1);
	}

	// If y < g(x), proceed to calculate the velocity components:

	velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
	theta = acos(randinter(-1.0, 1.0));
	phi = randinter(0.0,TWO_PI);

        bi->set_vel(vec(velocity * sin( theta ) * cos( phi ),
			   velocity * sin( theta ) * sin( phi ),
			   velocity * cos( theta )));
    }

    // The "raw" Plummer model just constructed has a virial radius of
    // approximately 1.695 units and is close to virial equilibrium.
    // Incorporate this factor into the basic model so that the unscaled
    // model (-u on the command line) is already almost in standard
    // units.  (The correct value is 16 / 3 PI = 1.6977, but the default
    // choice of mfrac makes the system slightly too compact.)  This is
    // convenient in situations where exact scaling may be too expensive.

    // Note from Steve (9/04).  The old way of doing this was INCORRECT,
    // as the scaled value of radius was used in computing the velocity.
    // As a result, the velocity distribution was not quite right...
    // Separating the rescaling completely from the generation of the
    // particle data eliminates this error.

    real xfac = 1/1.695;
    real vfac = 1/sqrt(xfac);

    for_all_daughters(dyn, b, bi) {
	bi->set_pos(xfac*bi->get_pos());
	bi->set_vel(vfac*bi->get_vel());
    }

    // Transform to center-of-mass coordinates, and scale to standard
    // parameters, assuming a softening parameter of zero.

    b->to_com();

    putrq(b->get_log_story(), "initial_mass", 1.0);

    if (!u_flag && n > 1) {

        real kinetic, potential;

	// Note: scale_* operates on internal energies.

	get_top_level_energies(b, 0.0, potential, kinetic);
	scale_virial(b, -0.5, potential, kinetic);	// scales kinetic
	real energy = kinetic + potential;
	scale_energy(b, -0.25, energy);			// scales energy

	putrq(b->get_log_story(), "initial_total_energy", -0.25);
	putrq(b->get_log_story(), "initial_rvirial", 1.0);
	putrq(b->get_dyn_story(), "total_energy", -0.25);
    }
}

local void swap(dyn* bi, dyn* bj)
{
    if (bi != bj) {

	real mass  = bi->get_mass();
	vec pos = bi->get_pos();
	vec vel = bi->get_vel();

	bi->set_mass(bj->get_mass());
	bi->set_pos(bj->get_pos());
	bi->set_vel(bj->get_vel());

	bj->set_mass(mass);
	bj->set_pos(pos);
	bj->set_vel(vel);
    }
}

// Original reshuffling code was O(N^2)!  Fixed by Steve, July 2000.

local void reshuffle_all(dyn* b, int n)
{
    // Reorder the tree by swapping mass, pos, and vel of each node
    // with a randomly chosen node.

    // First make a list of nodes.

    int i = 0;
    dynptr *list = new dynptr[n];
    for_all_daughters(dyn, b, bi)
	list[i++] = bi;

    // Then go down the list and randomize.

    for (i = 0; i < n; i++)
	swap(list[i], list[(int)randinter(0, n)]);

    delete [] list;
}

//-----------------------------------------------------------------------------
//  main  --
//      usage:
//                makeplummer -n# [options]  ,
// 
//            where # is the number of bodies in the Nbody system.
//    options:
//            The following options are allowed:
//       seed:
//            -s #   where # stands for an integer either in the range
//                   [1,2147483647], or # = 0 which causes a seed to be
//                   chosen as the value of the UNIX clock in seconds; this
//                   guarantees that no two calls will give the same value for
//                   the seed if they are more than 2.0 seconds apart.
//     cutoff:
//            -m #   mass fraction of the (infinitely extended) Plummer model
//            -r #   radius fraction of the (infinitely extended) Plummer model
//
//            If the mass fraction and radius fraction are both 1.0 then 
//            particles will be sprinkled in all over space. If mfrac < 1.0 or
//            rfrac < 1.0 then each particle is constrained to lie within
//            both the radial and (cumulative) mass bound.
//               For example, if rfrac( mfrac ) > rfrac then rfrac is the 
//            limiting factor, but if rfrac( mfrac ) < rfrac then mfrac limits
//            the extent of the Plummer realization.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    int  i;
    int  n;
    int  input_seed, actual_seed;
    real  mfrac;
    real  rfrac;

    bool  c_flag = false;
    bool  C_flag = false;
    bool  i_flag = false;
    bool  m_flag = false;
    bool  n_flag = false;
    bool  o_flag = false;
    bool  r_flag = false;
    bool  R_flag = true;
    bool  s_flag = false;
    bool  u_flag = false;

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:Cim:n:os:r:Ru";

    mfrac = rfrac = VERY_LARGE_NUMBER;

    while ((c = pgetopt(argc, argv, param_string,
			"$Revision: 1.20 $", _SRC_)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'm': m_flag = true;
		      mfrac = atof(poptarg);
		      break;
	    case 'n': n_flag = true;
		      n = atoi(poptarg);
		      break;
	    case 'o': o_flag = true;
		      break;
	    case 'r': r_flag = true;
		      rfrac = atof(poptarg);
		      break;
	    case 'R': R_flag = !R_flag;
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'u': u_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    if (!n_flag) {
        cerr << "makeplummer: specify the number # of";
	cerr << " particles with -n#\n";
	exit(1);
    }
    
    dyn *b, *by, *bo;
    b = new dyn();
    b->set_root();

    if (C_flag) b->set_col_output(true);

    if (n > 0) {
	bo = new dyn();
	if (i_flag)
	    bo->set_label(1);
	b->set_oldest_daughter(bo);
	bo->set_parent(b);
    }

    for (i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag)
	    by->set_label(i+1);
	by->set_parent(b);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    if (c_flag)
        b->log_comment(comment);
    b->log_history(argc, argv);

    if (!s_flag)
        input_seed = 0;                         	// default
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makeplummer: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    if (m_flag == FALSE && r_flag == FALSE)
        mfrac = MFRAC_DEFAULT;                   	 // default

    if (n != 0) {
        makeplummer(b, n, mfrac, rfrac, u_flag);
	if (R_flag) reshuffle_all(b, n);
    }

    put_dyn(b);
    rmtree(b);
}

#endif

