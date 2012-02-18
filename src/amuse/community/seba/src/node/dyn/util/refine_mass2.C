
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// refine_mass2.C:  More careful determination of cluster mass and
//		    center, for use with an external non-tidal field.
//		    Also flag escapers.
//
// Externally visible functions:
//
//	void refine_cluster_mass2

#include "dyn.h"

// In the following "ext" functions, b is any dyn, used simply to access
// the global class data.  The vector pos is relative to the root node.

local inline vec ext_acc(dyn *b, vec pos)
{
    vec v = 0, acc = 0, j = 0;
    real p = 0;
    get_external_acc(b, pos, v, p, acc, j);		// (discard p and j)
    return acc;
}

local inline real ext_pot(dyn *b, vec pos)
{
    vec v = 0, a = 0, j = 0;
    real pot = 0;
    get_external_acc(b, pos, v, pot, a, j, true);	// (discard a and j)
    return pot;
}

// Handy, for now...

static vec ext_center, Rhat, acc_CM;	// ext_center is the absolute position
					// of the external center, for use by
					// the "solve" functions below.
							
static real R;

local real g1(dyn *b, real r)
{
    vec pos = ext_center + r*Rhat - b->get_root()->get_pos();
    return (ext_acc(b, pos) - acc_CM) * Rhat + pow(r-R, -2);
}

local real g2(dyn *b, real r)
{
    vec pos = ext_center + r*Rhat - b->get_root()->get_pos();
    return (ext_acc(b, pos) - acc_CM) * Rhat - pow(r-R, -2);
}

#define EPS_R 1.e-4

local inline real solve(real (*g)(dyn*, real), dyn *b, real r1, real r2)
{
    // Return the root of g(r) = 0 (if any) between r1 and r2.
    // Avoid evaluation of g at r1 or r2.  Start the search at r2,
    // proceed toward r1 (could have r1 > r2), and return the first
    // zero found.

    int fac = (r1 < r2);
    fac = 2*fac - 1;		// = +1 if r1 < r2, -1 otherwise
    r1 *= 1 + fac*EPS_R;
    r2 *= 1 - fac*EPS_R;

    real dr = 0.01 * (r2 - r1);

    // Search for a zero, starting at r2.

    real g2 = g(b, r2);
    real r = r2 - dr;
    while (fac*(r - r1) >= 0 && g(b, r)*g2 > 0) r -= dr;

    if (fac*(r - r1) < 0) return r;

    // Refine by bisection.

    r1 = r;
    r2 = r + dr;
    real g1 = g(b, r1);
    g2 = g(b, r2);

    while (abs(r2/r1 - 1) < EPS_R) {
	r = 0.5*(r1+r2);
	real gr = g(b, r);
	if (gr == 0) return r;
	if (gr*g1 > 0) {
	    r1 = r;
	    g1 = gr;
	} else {
	    r2 = r;
	    g2 = gr;
	}
    }

    // Final refinement by linear interpolation.

    r = r1 + (r2-r1)*(0-g1)/(g2-g1);

    return r;
}

local void get_rL(dyn *b,		// root node
		  real M,		// current mass estimate
		  vec center,		// current center estimate (relative)
		  real& r_L1,
		  real& r_L2)
{
    // Find the inner Lagrangian point in the total external field.
    // Construction of the external field functions is such that
    // it is just as easy to work directly with 3-D vectors...

    // Note that on entry center is relative to the root node.

    acc_CM = ext_acc(b, center);

    center += b->get_pos();		// center now in absolute coordinates
    ext_center = b->get_external_center();

    R = abs(center - ext_center);
    Rhat = (center - ext_center)/(R*M);	// factor of M avoids passing M to g

    real scale
	= sqrt(b->get_external_scale_sq());  // starting point for L1 search

    // Scale is based on the first external field found.  Only used to
    // set an overall scale for use in solve(g1), so probably no need
    // to check for multiple internal fields.  The accelerations used
    // in solve(g?) include all fields.

    r_L1 = solve(g1, b, scale, R);
    r_L2 = solve(g2, b, R, 10*R);
}

local int bitcount(unsigned int i)
{
    // Count nonzero bits.  From K&R.

    int b;
    for (b = 0; i != 0; i >>= 1)
	if (i & 01) b++;
    return b;
}

#define M_TOL 1.e-4
#define M_ITER_MAX 20

#define TTOL 1.e-12				// arbitrary tolerance

void refine_cluster_mass2(dyn *b,		// root node
			  int verbose)		// default = 0
{
    // Self-consistently determine the total mass within the outermost
    // closed zero-velocity surface under the specified external field(s).
    // Use a point-mass approximation for the cluster potential and iterate
    // until the actual mass within the surface agrees with the mass used
    // to generate the surface.
    //
    // This code is most likely to be called from refine_cluster_mass().

    // Quit if internal forces are to be neglected.

    if (b->get_ignore_internal()) return;

    // Do nothing if all we want is to set the dyn story and the current
    // values are up to date.

    if (verbose == 0
	&& twiddles(getrq(b->get_dyn_story(), "bound_center_time"),
		    b->get_system_time(), TTOL))
	return;

    unsigned int ext = b->get_external_field();
    if (!ext || b->get_tidal_field()) return;

    // Method assumes that the external field can be characterized as
    // depending on distance from a single point -- i.e. that it has a
    // single center -- and that it is independent of velocity.  OK to
    // have multiple external fields (e.g. central mass + power law),
    // so long as they have a common center.

    int bits_set = bitcount(ext);
    vec external_center = 0;

    if (bits_set != 1) {

	if (bits_set < 1)
	    
	    return;

	else {

	    // We have multiple external fields.  Check for common center.

	    cerr << "  refine_cluster_mass2: "; PRC(ext); PRL(bits_set);

	    int bit = 0, cbit = -1;
	    for (unsigned int i = ext; i != 0; i >>= 1, bit++) {
		if (i & 01) {
		    if (cbit < 0) {
			cbit = bit;
			external_center = b->get_external_center(cbit);
		    } else {
			if (!twiddles(square(b->get_external_center(bit)
					      - external_center),
				      0)) {
			    cerr << "  refine_cluster_mass2: center " << bit
				 << " inconsistent with center " << cbit
			         << endl;
			    return;
			}
		    }
		}
	    }
	    cerr << "  common center = " << external_center << endl;
	}

    } else

	external_center = b->get_external_center();

    // Function may be called repeatedly for a series of systems.
    // Make sure the root node is correctly set for this one.

    b->set_root(b);
    vec bpos = b->get_pos(), bvel = b->get_vel();

    // Use the standard center as our starting point.  The center will be
    // redetermined self-consistently, along with the mass.

    vec center, vcenter;
    get_std_center(b, center, vcenter);		// std_center quantities
						// include the root node

    // Choose the initial mass to include only the volume between the
    // standard center and the center of the external potential.

    real M_inside = 0;
    R = abs(center - external_center);			// global R

    center  -= bpos;					// remove the root from
    vcenter -= bvel;					// "center" quantities

    for_all_daughters(dyn, b, bb)
	if (abs(bb->get_pos() - center) <= R)
	    M_inside += bb->get_mass();

    // The starting quantities for the loop are not determined in the
    // same way as those calculated within the loop -- possible room
    // for inconsistency?

    real M = -1;					// (to start the loop)
    int N_inside;
    vec cen = center, vcen = vcenter;			// relative to root

    if (verbose) {
	cerr << endl << "  refine_cluster_mass2: getting mass by iteration"
	     << endl << "  initial total system mass = " << M_inside
	     << endl << "  initial center (abs) = " << center + bpos
	     << endl;
    }

    int iter = 0;
    real r_L = 0, phi_lim = 0;

    // Iteration can (should!) run away to zero mass if the cluster
    // density is too low relative to the local tidal field strength.
    // Keep track of the "50% center of mass" during the iteration as
    // a convenient measure of the cluster center in this case.

    real M0 = M_inside;
    vec cen50 = center, vcen50 = vcenter;		// relative to root
    bool set50 = false;

    while (iter++ < M_ITER_MAX
	   && M_inside > 0
	   && abs(M_inside/M - 1) > M_TOL) {

	// Set current mass and center (relative to root):

	M = M_inside;
	center = cen;
	vcenter = vcen;

	// Reinitialize:

	M_inside = 0;
	N_inside = 0;
	cen = vcen = 0;

	// Determine the Lagrangian points of the (point-mass) cluster
	// in the *total* external field, then count stars within the
	// limiting equipotential.

	real r_L1, r_L2;
	get_rL(b, M, center, r_L1, r_L2);		// NB relative center
							// needed here
	if (verbose > 1) {
	    cerr << endl;
	    PRC(R); PRC(r_L1); PRC(r_L2);
	}

	// Limiting potential (ext_pot takes absolute position):

	real phi_L1 = ext_pot(b, ext_center + r_L1*Rhat - bpos) - M/(R-r_L1);
	real phi_L2 = ext_pot(b, ext_center + r_L2*Rhat - bpos) - M/(r_L2-R);
	// PRC(phi_L1); PRC(phi_L2);

	phi_lim = Starlab::max(phi_L1, phi_L2);    // maximize the cluster mass
	r_L = Starlab::max(R-r_L1, r_L2-R);

	if (verbose > 1) PRL(r_L);

	for_all_daughters(dyn, b, bb) {
	    real r = abs(bb->get_pos() - center);
	    if (r < r_L) {
		if (r == 0
		    || -M/r + ext_pot(b, bb->get_pos()) < phi_lim) {
		    N_inside++;
		    M_inside += bb->get_mass();
		    cen += bb->get_mass()*bb->get_pos();
		    vcen += bb->get_mass()*bb->get_vel();
		}
	    }
	}

	if (M_inside > 0) {
	    cen /= M_inside;
	    vcen /= M_inside;
	}

	// Linearly interpolate an estimate of the "50% center."

	if ((M > 0.5*M0 && M_inside <= 0.5*M0)
	    || (M < 0.5*M0 && M_inside >= 0.5*M0)) {
	    cen50 = center + (0.5*M0-M)*(cen-center)/(M_inside-M);
	    vcen50 = vcenter + (0.5*M0-M)*(vcen-vcenter)/(M_inside-M);
	    set50 = true;
	}

	if (verbose > 1) {
	    PRI(2); PRC(iter); PRC(N_inside); PRL(M_inside);
	    PRI(2); PRL(cen + bpos);
	}
    }

    if (iter >= M_ITER_MAX) {
	if (verbose) PRI(2);
	cerr << "warning: refine_cluster_mass: too many iterations at time "
	     << b->get_system_time() << endl;
    }

    if (verbose == 1) {
	PRI(2); PRC(iter); PRC(N_inside); PRL(M_inside);
	PRI(2); PRL(cen + bpos);
	PRI(2); PRL(vcen + bvel);
    }

    // Too many iterations probably means a limit cycle of some sort;
    // accept the results in that case.

    bool disrupted = (M_inside < 0.01*M0);

    if (disrupted) {

	// Looks like the cluster no longer exists.  Use 50% or ext center.

	if (set50) {
	    center = cen50;
	    vcenter = vcen50;
	} else {
	    center = ext_center - bpos;
	    vcenter = -bvel;			// so v = 0...
	}

	if (verbose) {
	    PRI(2); PRL(center);
	    PRI(2); PRL(vcenter);
	}

    } else {

	center = cen;
	vcenter = vcen;
    }

    // Now center and vcenter should be usable, even if M_inside and
    // N_inside aren't very meaningful.

    // Convert the center to absolute coordinates (add root quantities)
    // and write our best estimate to the dyn story.

    vec bc_pos = bpos + center;
    vec bc_vel = bvel + vcenter;

    putrq(b->get_dyn_story(), "bound_center_time", (real)b->get_system_time(),
	  HIGH_PRECISION);
    putvq(b->get_dyn_story(), "bound_center_pos", bc_pos);
    putvq(b->get_dyn_story(), "bound_center_vel", bc_vel);

    // Send the relevant data to the dynamical friction code, and compute
    // the current frictional acceleration.  Could combine these calls, but
    // keep as is for now.

    set_friction_mass(M_inside);
    set_friction_vel(bc_vel);
    set_friction_acc(b, abs(bc_pos - external_center));

    // Repeat the inner loop above and flag stars as escapers or not.

    int n_mem = 0;
    for_all_daughters(dyn, b, bb) {
	bool escaper = true;
	if (!disrupted) {
	    real r = abs(bb->get_pos() - center);
	    if (r < r_L
		&& (r == 0 || -M/r + ext_pot(b, bb->get_pos()) < phi_lim))
		escaper = false;
	}
	putiq(bb->get_dyn_story(), "esc", escaper);
	if (!escaper) n_mem++;

	// Print out t_esc if this is a new escaper.

	if (escaper && !find_qmatch(bb->get_dyn_story(), "t_esc"))
	    putrq(bb->get_dyn_story(), "t_esc", (real)b->get_system_time());

	// Delete t_esc if this is not an escaper.

	if (!escaper && find_qmatch(bb->get_dyn_story(), "t_esc"))
	    rmq(bb->get_dyn_story(), "t_esc");

	// Note that the esc flag is redundant -- all necessary
	// information is contained in the t_esc line.
    }

#if 0
    if (verbose) PRI(2);
    cerr << "refine_mass2: "; PRC(b->get_system_time()); PRL(n_mem);
#endif

}
