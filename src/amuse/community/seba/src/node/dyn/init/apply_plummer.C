
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Apply a Plummer-model spatial and velocity distribution to an
//// existing snapshot, without scaling.  Masses are left unchanged.
////
//// Usage:  apply_plummer [OPTIONS] < input > output
////
//// Options:
////         -s    specify random seed [random from system clock]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// For details on model construction, see notes in makeplummer.C
// Steve McMillan, July 1996

#include "dyn.h"

#ifdef TOOLBOX

#define  MFRAC_DEFAULT  0.999                 /* radial cutoff               */
#define  RFRAC_DEFAULT  22.8042468            /* corresponding radial cutoff */
                                              /*   (not used right now)      */
#define  SEED_STRING_LENGTH  60

local void  addplummer(dyn * b, real mfrac, real rfrac)
{
    real  radius;		   /* absolute value of position vector      */
    real  velocity;		   /* absolute value of velocity vector      */
    real  theta, phi;		   /* direction angles of above vectors      */
    real  x, y;		           /* for use in rejection technique         */
    real  scalefactor;             /* for converting between different units */
    real  inv_scalefactor;         /* inverse scale factor                   */
    real  sqrt_scalefactor;        /* sqare root of scale factor             */
    real  m_min, m_max;            /* mass shell limits for quiet start      */

    scalefactor = 16.0 / (3.0 * PI);
    inv_scalefactor = 1.0 / scalefactor;
    sqrt_scalefactor = sqrt(scalefactor);

    //  Set up positions and velocities for top-level nodes.

    m_max = 0;
    for_all_daughters(dyn, b, bi) {

        m_min = m_max;
        m_max += mfrac * bi->get_mass() / b->get_mass();
	real m = randinter(m_min, m_max);
	radius = 1.0 / sqrt( pow (m, -2.0/3.0) - 1.0);

	theta = acos(randinter(-1.0, 1.0));
	phi = randinter(0.0, TWO_PI);

        bi->set_pos(vec(radius * sin( theta ) * cos( phi ),
			   radius * sin( theta ) * sin( phi ),
			   radius * cos( theta )));

	x = 0.0;
	y = 0.1;

	while (y > x*x*pow( 1.0 - x*x, 3.5)) {
	    x = randinter(0.0, 1.0);
	    y = randinter(0.0, 0.1);
	}

	velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
	theta = acos(randinter(-1.0, 1.0));
	phi = randinter(0.0, TWO_PI);

        bi->set_vel(vec(velocity * sin( theta ) * cos( phi ),
			   velocity * sin( theta ) * sin( phi ),
			   velocity * cos( theta )));
    }

    //  Transform to center-of-mass coordinates.

    b->to_com();
}

main(int argc, char ** argv)
{
    int  random_seed;
    real  mfrac = MFRAC_DEFAULT;
    real  rfrac = RFRAC_DEFAULT;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "s:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c) {

	    case 's': random_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
		      exit(1);
	}            
    
    dyn *b;
    b = get_dyn();

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    char seedlog[64];
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    addplummer(b, mfrac, rfrac);

    put_dyn(b);
}

#endif
