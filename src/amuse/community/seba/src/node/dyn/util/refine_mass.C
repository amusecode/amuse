
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// refine_mass.C:  More careful determination of cluster mass.
//		   Also flag escapers.
//
//		   Separate from refine_mass2 because the solution
//		   for the Jacobi radius is analytical, and the
//		   form of the potential is particularly simple.
//
// Externally visible functions:
//
//	void refine_cluster_mass

#include "dyn.h"

#define G_TOL 1.e-6

local real solve3(real a, real b, real c)
{
    // Iteratively solve
    //
    //		a z^3 + b z + c  =  0
    //
    // where we know a > 0, b > 0, c < 0.

    // Proceed by Newton-Raphson iteration.  Since a > 0, b > 0, so
    // there are are no turning points to complicate the logic.

    real z = 2*Starlab::max(pow(abs(c/a), 1.0/3), sqrt(abs(b/a)));
    real g = c + z * (b + a*z*z);
    real gpr = b + 3*a*z*z;		// always > 0

    while (abs(g) > G_TOL) {
	z = z - g/gpr;
	g = c + z * (b + a*z*z);
	gpr = b + 3*a*z*z;
    }

    return z;
}

#define M_TOL 1.e-4
#define M_ITER_MAX 20

void refine_cluster_mass(dyn *b,
			 int verbose)		// default = 0
{
    if (b->get_external_field() == 0) return;

    if (b->get_tidal_field() == 0) {
	refine_cluster_mass2(b, verbose);	// external, not tidal
	return;
    }

    // Self-consistently determine the total mass within the outermost
    // closed zero-velocity surface.  Use a point-mass approximation for
    // the cluster potential and iterate until the actual mass within
    // the surface agrees with the mass used to generate the surface.
    // The total potential is
    //
    //		phi  =  -GM/r + (alpha1 x^2 + alpha3 z^2) / 2
    //
    // where we measure everything relative to the standard center (which
    // is the density center, if known, and the modified center of mass
    // otherwise).  The Jacobi radius for this field is
    //
    //		r_J  =  (-GM/alpha1)^{1/3}
    //
    // and the potential of the last closed surface is
    //
    //		phi_J  =  1.5 alpha1 r_J^2.

    vec center, vcenter;
    int which = get_std_center(b, center, vcenter);

    center -= b->get_pos();			// std_center quantities
    vcenter -= b->get_vel();			// include the root node

    // which = 1 for density center, 2 for mcom.

    // Note from Steve (7/01):  The "center" should really be the
    // center of mass of particles within the Jacobi surface,
    // determined self-consistently.
    //
    // To do...  (See also check_and_remove_escapers().)

    real M_inside = total_mass(b), M = -1, M0 = M_inside;
    real r_J, r_x2, r_y2, r_z2, r_max_inside;
    int  N_inside, iter = 0;

    real phi_J;

    real M_J, M_x, M_y, M_z;
    int  N_J, N_x, N_y, N_z;

    if (verbose)
	cerr << endl << "  refine_cluster_mass: getting M by iteration"
	     << endl << "  initial total system mass = " << M_inside
	     << endl;

    while (iter++ < M_ITER_MAX
	   && M_inside > 0
	   && abs(M_inside/M - 1) > M_TOL) {

	M = M_inside;
	M_inside = 0;
	M_J = M_x = M_y = M_z = 0;
	N_J = N_x = N_y = N_z = 0;

	r_J = pow(-M/b->get_alpha1(), 0.3333333);

	phi_J = 1.5 * b->get_alpha1() * r_J * r_J;
	real r_max = 0;

	r_x2 = square(r_J);		// zero-velocity surface crosses x-axis
	r_y2 = square(-M/phi_J);	// zero-velocity surface crosses y-axis

	// zero-velocity surface crosses z-axis where
	//
	//	-GM/z + 0.5 alpha3 z^2  =  phi_J
	// so
	//	0.5 alpha3 z^3 -phi_J z -GM  =  0

	r_z2 = square(solve3(0.5*b->get_alpha3(), -phi_J, -M));

	N_inside = 0;
	r_max_inside = 0;

	for_all_daughters(dyn, b, bb) {

	    vec dx = bb->get_pos() - center;

	    real x = dx[0];
	    real y = dx[1];
	    real z = dx[2];
	    real r = abs(dx);

	    if (r < r_J) {

		N_J++;
		M_J += bb->get_mass();

		if (r == 0
		    || -M/r + 0.5 * (b->get_alpha1()*x*x + b->get_alpha3()*z*z)
			    < phi_J) {
		    N_inside++;
		    M_inside += bb->get_mass();
		    r_max_inside = Starlab::max(r, r_max_inside);
		}

	    }
	    r_max = Starlab::max(r, r_max);

	    // Count projected masses and numbers.

	    real xy = x*x + y*y;
	    real xz = x*x + z*z;
	    real yz = y*y + z*z;

	    if (xy < r_x2) {		// z projection
		M_x += bb->get_mass();
		N_x++;
	    }
	    if (xz < r_x2) {		// y projection
		M_x += bb->get_mass();
		N_x++;
	    }

	    if (yz < r_y2) {		// x projection
		M_y += bb->get_mass();
		N_y++;
	    }
	    if (xy < r_y2) {		// z projection
		M_y += bb->get_mass();
		N_y++;
	    }

	    if (yz < r_z2) {		// x projection
		M_z += bb->get_mass();
		N_z++;
	    }
	    if (xz < r_z2) {		// y projection
		M_z += bb->get_mass();
		N_z++;
	    }
	}

	M_x /= 2;
	N_x /= 2;
	M_y /= 2;
	N_y /= 2;
	M_z /= 2;
	N_z /= 2;

#if 0
	fprintf(stderr, "    %2d  ", iter);
	PRC(M); PRC(r_J); PRL(r_max);
	PRI(8); PRC(N_inside); PRC(M_inside); PRL(r_max_inside);
#endif

    }

    if (iter >= M_ITER_MAX) {
	if (verbose) PRI(2);
	cerr << "warning: refine_cluster_mass: too many iterations at time "
	     << b->get_system_time() << endl;
    }

    if (verbose)
	cerr << "  within last zero-velocity surface:"
	     << "  M = " << M_inside << "  N = " << N_inside
	     << endl
	     << "                                    "
	     << "  r_J = " << r_J << "  r_max = " << r_max_inside
	     << endl
	     << "  within r_J:                       "
	     << "  M = " << M_J << "  N = " << N_J
	     << endl
	     << "  within r_J (projected):           "
	     << "  M = " << M_x << "  N = " << N_x
	     << endl
	     << "  within r_y (projected):           "
	     << "  M = " << M_y << "  N = " << N_y << "  r_y = " << sqrt(r_y2)
	     << endl
	     << "  within r_z (projected):           "
	     << "  M = " << M_z << "  N = " << N_z << "  r_z = " << sqrt(r_z2)
	     << endl;

    // Repeat the inner loop above and flag stars as escapers or not.
    // Too many iterations probably means a limit cycle of some sort;
    // accept the results in that case.

    bool disrupted = (M_inside < 0.01*M0);

    for_all_daughters(dyn, b, bb) {

	vec dx = bb->get_pos() - center;

	real x = dx[0];
	real z = dx[2];
	real r = abs(dx);

	bool escaper = true;
	if (!disrupted) {
	    if (r < r_J
		&& (r == 0 || -M/r + 0.5 * (b->get_alpha1()*x*x
					     + b->get_alpha3()*z*z) < phi_J))
		escaper = false;
	}
	putiq(bb->get_dyn_story(), "esc", escaper);

	// Print out t_esc if this is a new escaper.

	if (escaper && !find_qmatch(bb->get_dyn_story(), "t_esc"))
	    putrq(bb->get_dyn_story(), "t_esc", (real)b->get_system_time());

	// Delete t_esc if this is not an escaper.

	if (!escaper && find_qmatch(bb->get_dyn_story(), "t_esc"))
	    rmq(bb->get_dyn_story(), "t_esc");
    }
}
