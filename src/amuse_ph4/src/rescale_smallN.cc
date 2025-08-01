// Rescale an hdyn system to a specified radius, for reinsertion into
// parallel_hermite_4.  Functions here repeat thse found elsewhere in
// smallN...

#include "hdyn.h"

static void move_to_com(hdyn *b)
{
    // Force the system under b to the CM frame.

    real mtot = 0;
    vec cmpos = 0, cmvel = 0;
    for_all_daughters(hdyn, b, bb) {
	mtot += bb->get_mass();
	cmpos += bb->get_mass()*bb->get_pos();
	cmvel += bb->get_mass()*bb->get_vel();
    }
    cmpos /= mtot;
    cmvel /= mtot;
    for_all_daughters(hdyn, b, bb) {
	bb->inc_pos(-cmpos);
	bb->inc_vel(-cmvel);
    }
}

static void get_energy_and_scale(hdyn *b, real& kin, real& pot, real& rmax)
{
    kin = 0, pot = 0, rmax = 0;
    for_all_daughters(hdyn, b, bi) {
	real mi = bi->get_mass();
	vec xi = bi->get_pos();
	real r2 = square(xi);
	if (r2 > rmax) rmax = r2;
	kin += mi*square(bi->get_vel());
	real ppot = 0;
	for (hdyn *bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	    ppot -= bj->get_mass()/abs(bj->get_pos()-xi);
	pot += mi*ppot;
    }
    kin /= 2;
    rmax = sqrt(rmax);
}

void rescale(hdyn *b, real size)
{
    move_to_com(b);

    // Compute the energies and find the top-level node farthest from
    // the CM.  See top_level_energy() in smallN.cc for redundancy.

    real kin, pot, rmax;
    get_energy_and_scale(b, kin, pot, rmax);
    PRC(rmax); PRC(-kin/pot); PRL(kin+pot);

    real rscale = size/rmax;
    real vscale = sqrt(1 - pot*(1/rscale-1)/kin);	// check this!

    // Rescale all positions by the same factor to place the farthest
    // particle at distance size, and rescale velocities to preserve
    // the total energy.  ** Doesn't preserve angular momentum... **
    // TBD

    for_all_daughters(hdyn, b, bb) {
	bb->set_pos(rscale*bb->get_pos());
	bb->set_vel(vscale*bb->get_vel());
    }

    // Debugging: check energies and scale...

    get_energy_and_scale(b, kin, pot, rmax);
    PRC(rmax); PRC(-kin/pot); PRL(kin+pot);
}
