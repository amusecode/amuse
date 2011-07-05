
// Functions to handle partial or full unperturbed motion in smallN.
// Structure shamelessly cribbed from Starlab::hdyn_unpert.C by the
// author.  Retained Starlab data structures to simplify code
// importation.
//
// Author: Steve McMillan.
//
// Externally visible functions:
//
//	real get_tidal_potential(hdyn *b, hdyn *bi, hdyn *bj, hdyn *cm)
//	bool is_unperturbed_and_approaching(hdyn *b, hdyn *bi, hdyn *bj)
//	bool create_binary(hdyn *bi, hdyn *bj, bool verbose = false)
//	bool extend_or_end_binary(hdyn*& bi, bool verbose = false)

#include "hdyn.h"

real get_tidal_potential(hdyn *b, hdyn *bi, hdyn *bj, hdyn *cm,
			 bool absolute)		// default = true
{
    // Calculate the tidal potential of the system on (bi,bj) with
    // top-level center of mass cm.  If absolute is true (the
    // default), the component positions are assumed to be absolute;
    // if false, they are relative to their parent node.

    real mi = bi->get_mass(), mj = bj->get_mass(), mcm = cm->get_mass();
    vec xi = bi->get_pos(), xj = bj->get_pos(), xcm = cm->get_pos();

    if (!absolute) {
	xi += xcm;
	xj += xcm;
    }

    real phi_t = 0;
    for_all_daughters(hdyn, b, bk)
	if (bk != bi && bk != bj && bk != cm) {
	    vec xk = bk->get_pos();
	    real dphi_k = - mi/abs(xi-xk) - mj/abs(xj-xk) + mcm/abs(xcm-xk);
	    phi_t += bk->get_mass()*dphi_k;
	}
    return phi_t;
}

bool is_unperturbed_and_approaching(hdyn *b, real dt, hdyn *bi, hdyn *bj)
{
    // Check whether the binary (bi,bj) satisfies the requirements for
    // unpeturbed motion.  OK to use pos and vel here because pred_pos
    // and pos are the same at the end of a step, when this function
    // is called.
 
    if (dt > b->get_dt_crit()) return false;
    if (bi->is_parent() || bj->is_parent()) return false;
    if (square(bi->get_pos()-bj->get_pos()) > b->get_r2_crit()) return false;
    if ((bi->get_pos()-bj->get_pos()) * (bi->get_vel()-bj->get_vel()) > 0)
	return false;

    // Particles are within r_crit and approaching.  Check the
    // perturbation due to the next nearest neignbor.  Work with
    // |potential| rather than simple distance to take mass into
    // account.

    vec cm = (bi->get_mass()*bi->get_pos() + bj->get_mass()*bj->get_pos())
		/ (bi->get_mass() + bj->get_mass());
    hdyn *bk_nn = NULL;
    real max_pot2 = 0, nn_dist2 = 0;
    for_all_leaves(hdyn, b, bk)
	if ((bk != bi) && (bk != bj) ) {
	    real dist2 = square(bk->get_pos() - cm);
	    real pot2 = pow(bk->get_mass(), 2)/dist2;
	    if (pot2 > max_pot2) {
		max_pot2 = pot2;
		nn_dist2 = dist2;
		bk_nn = bk;
	    }
	}

    // Particle bk_nn has the largest potential relative to the CM.
    // Compute its perturbation, and return true if it is small.

    real dist2 = square(bi->get_pos()-bj->get_pos());
    real mass = bi->get_mass() + bj->get_mass();
    real gamma2 = 4*pow(bk_nn->get_mass()/mass, 2)*pow(dist2/nn_dist2, 3);

    return (gamma2 < b->get_gamma2_unpert());
}



static inline vec get_cm_acc(hdyn *b, hdyn *bi, hdyn *bj, vec pos)
{
    // Return the acceleration at location pos due to all nodes
    // excluding bi and bj.

    vec acc = 0;
    for_all_daughters(hdyn, b, bk)
	if (bk != bi && bk != bj) {
	    vec dr = bk->get_pos() - pos;
	    acc += bk->get_mass()*dr/pow(dr*dr, 1.5);
	}
    return acc;
}

static inline real integrand(real t, real e)
{
    return cos(2*t)/pow(1+e*cos(t),4);
}

static inline real num_int(kepler *k)
{
    // Numerial determination of the orbital integral using the
    // trapezoid rule (intended as a test of Steve's algebra -- seems
    // to be OK).

    real M = k->get_total_mass();
    real a = k->get_semi_major_axis();
    real e = k->get_eccentricity();
    real J = k->get_angular_momentum();
    real th0 = fabs(k->get_true_anomaly());

    real dth = 0.01*th0, sum = 0.5*(integrand(0,e)+integrand(th0,e));
    for (real t = dth; t < th0-0.5*dth; t += dth) sum += integrand(t,e);

    if (fabs(1-e) < 1.e-6)
	return 2*pow(J,7)*sum*dth/pow(M,4);
    else
	return 2*pow(a*fabs(1-e*e),4)*sum*dth/J;
}

static inline real orb_int(kepler *k)
{
    // Compute the angular integral
    //
    //		\int (x^2 - y^2) dt
    //
    // that appears in the calculation of the change in the angular
    // momentum.  Transform it to
    //
    //		\int r^4 cos(2 theta) d(theta)/J
    //
    // from th0 to -th0, where th0 < 0 is the incoming true anomaly
    // and J is the angular momentum.

    //  Angular integral (for e not equal to 1) is
    //
    //		(A^4/J) \int cos(2 theta) / (1 + e cos theta)^4
    //
    // from th0 to -th0, where th0 < 0 is the true anomaly, A =
    // a(1-e^2) and J is angular momentum.

    real M = k->get_total_mass();
    real a = k->get_semi_major_axis();
    real e = k->get_eccentricity();
    real J = k->get_angular_momentum();
    real r = k->get_separation();
    real th0 = fabs(k->get_true_anomaly());
    real s = sin(th0);
    real c = cos(th0);

    // The following integrals are taken from Mathematica, with
    // additional algebra to simplify some expressions.

    real integral;

    if (fabs(1-e) < 1.e-6) {

	// ~Parabolic (not linear) orbit.  In this case,
	//
	//	r (1 + cos(theta)) = J^2/M
	//
	// and the integral becomes
	//
	//	(J^7/M^4) \int cos(2 theta) / (1 + cos(theta)) d(theta)

	integral = 4*pow(J, 7)*pow(M, -4)
			* (-th0 + sin(th0) + 0.5*tan(th0/2));

    } else {

	// Elliptical or hyperbolic orbit, with
	//
	//	r (1 + e cos(theta)) = A = a |1-e^2|
	//
	// and the integral is
	//
	//	(A^4/J) \int cos(2 theta) / (1 + e cos theta)^4
	//
	// Simplify the Mathematica expressions using the above
	// expression for r and
	//
	//	A^4/J = sqrt(a^7/M) (1-e^2)^(7/2)

	real fac1 = -2*pow(e,5)+9*pow(e,3)+8*e
		    + 3*(3*pow(e,4)+9*e*e-2)*c
		    + e*(8*pow(e,4)+9*e*e-2)*c*c;
	real fac2 = s * J * pow(r,3) / (3*pow(fabs(1-e*e),3)*M);
	real fac3 = 10 * e * e * sqrt(pow(a,7)/M);
	real fac4 = sqrt(fabs((1-e)/(1+e)))*tan(th0/2);
    
	if (e < 1)
	    integral = -fac1*fac2 + fac3*atan(fac4);
	else
	    integral =  fac1*fac2 - fac3*atanh(fac4);
    }

    // cout << "orb_int: " << e << " " << integral << " "
    //	 << num_int(k) << endl << flush;

    return integral;
}

bool create_binary(hdyn *bi, hdyn *bj, bool verbose)
{
    // Combine two nodes, replacing the first by the center of mass of
    // the two, and moving both to lie below the new center of mass
    // node.  On entry we have just completed a step, so {pos,vel} and
    // pred_{pos,vel} should be the same.

    if (!bi || !bj) return false;
    if (bi->get_kepler()
	&& bi->get_kepler()->tidal_potential != 0) return false;

    // Construct a kepler structure describing the relative motion.

    real time = bi->get_system_time();
    real mi = bi->get_mass();
    real mj = bj->get_mass();

    kepler *k = new kepler;
    k->set_time(time);
    k->set_total_mass(mi+mj);
    k->set_rel_pos(bj->get_pos() - bi->get_pos());
    k->set_rel_vel(bj->get_vel() - bi->get_vel());
    k->initialize_from_pos_and_vel();

    // Set a time step.  Choose the time to separation at the same
    // radius.  Could use pred_advance_to_periastron(), but this may
    // be quite inaccurate.  Better to use mean motion.

    real mean_anomaly = k->get_mean_anomaly();
    if (k->get_energy() < 0)
	mean_anomaly = sym_angle(mean_anomaly);	    // place in (-PI, PI]

    real peri_time;
    if (mean_anomaly >= 0) {			    // should never happen
	delete k;
	return false;
    } else
	peri_time = time - mean_anomaly / k->get_mean_motion();

    // Create a new center of mass node.

    hdyn *cm = new hdyn;
    cm->set_mass(mi+mj);
    cm->set_pos((mi*bi->get_pos()+mj*bj->get_pos())/cm->get_mass());
    cm->set_vel((mi*bi->get_vel()+mj*bj->get_vel())/cm->get_mass());

    // Calculate and save data needed to correct the energy and
    // angular momentun after the unperturbed segment.  See
    // partial_unperturbed.pdf for details.

    // The energy correction is just the change in the tidal potential
    // (components - CM).  Store the initial value of the tidal
    // potential in k->tidal_potential.

    k->tidal_potential = get_tidal_potential(bi->get_parent(), bi, bj, cm);

    // The angular momentum correction (neglecting reorientation of
    // the orbit plane) requires an integral over the orbit and
    // computation of a spatial derivative of the acceleration.
    // Specifically, we need to calculate d(a_y)/dx, where x and y are
    // coordinates in the orbit plane, with x along the major axis of
    // the orbit and y perpendicular to it.  We could compute the
    // Jacobi matrix \partial da_i/\partial dx_j analytically and
    // transform it to the orbital frame, but since we need only a
    // single component it is easier to calculate it numerically.

    // Compute the "x" component of the acceration at the CM and at a
    // point slightly offset (by a distance equal to the current
    // particle separation) in the "y" direction.

    vec cmpos = cm->get_pos();
    vec acc = get_cm_acc(bi->get_parent(), bi, bj, cmpos);
    real dx = k->get_separation();
    vec cmpos1 = cmpos + dx*k->get_longitudinal_unit_vector();
    vec acc1 = get_cm_acc(bi->get_parent(), bi, bj, cmpos1);
    real daydx = (acc1-acc)*k->get_transverse_unit_vector()/dx;

    // Compute the angular integral orb_int(theta) to the reflection
    // point and calculate the change in the angular momentum.  Note
    // that we expect dJ/J to be less than dE/E by a factor of
    // sqrt(gamma), and we keep only the first-order term parallel to
    // J, neglecting (second-order) changes in the orbit orientation.

    k->delta_angular_momentum = daydx*orb_int(k);

    if (verbose) {
	cout << endl << "too close: " << bi->get_index();
	cout << " and " << bj->get_index() << endl << flush;
	PRL(total_energy(bi->get_parent()));
    }

    //-----------------------------------------------------------------

    // Offset components to the center of mass frame.

    bi->inc_pos(-cm->get_pos());
    bi->inc_vel(-cm->get_vel());
    bj->inc_pos(-cm->get_pos());
    bj->inc_vel(-cm->get_vel());

    // Restructure the tree to create the binary (bi,bj).  The new CM
    // node is cm.

    create_binary_node(cm, bi, bj);

    // Leaf indices are meainingful and retain their meaning
    // throughout.  CM node indices are generated internally, and are
    // discarded when the node is destroyed.  Desirable to have the CM
    // index reflect the indices of the components, but this is in
    // general not feasible if the leaf indices are large (e.g. if
    // they come from a larger simulation).  Instead of a formula to
    // compute the indices, just use a manager as in ph4.

    // cm->set_index(100000+1000*bi->get_index()+10*bj->get_index());
    int cm_index = bi->get_cm_index();
    cm->set_index(cm_index);
    bi->set_cm_index(cm_index+1);

    // The relevant kepler is attached to (only) the first component.

    bi->set_kepler(k);
    bi->set_t_pred(2*peri_time-time);	// end of unperturbed segment

    // Start with pericenter reflection.  Later, allow possibility of
    // completely unperturbed motion extending over several orbits, as
    // in kira.

    cm->set_fully_unperturbed(false);
    bi->set_fully_unperturbed(false);

    if (verbose) {
	cout << "created new binary " << cm->get_index()
	     << " at time " << bi->get_system_time() << endl;
	k->print_all();

	real dt_unpert = bi->get_t_pred() - time;
	int p = cout.precision(10);
	PRC(dt_unpert); PRL(bi->get_t_pred());
	PRC(peri_time); PRL(peri_time-time);
	PRL(bi->get_t_pred()-time);
	cout.precision(p);
	PRL(k->tidal_potential);
	PRC(daydx); PRL(k->delta_angular_momentum);

	// Special planar case: angle to 3rd star; acc should be in
	// the x-y plane.

	// real phi3 = atan2(acc[1], acc[0]);
	// PRL(phi3);
    }

    // Make sure pos and pred_pos agree.

    for_all_nodes(hdyn, cm, bb) {
	bb->set_pred_pos(bb->get_pos());
	bb->set_pred_vel(bb->get_vel());
    }

    return true;
}



static inline real time_to_radius(real dr,	// distance to cover
				  real vr,	// relative radial velocity
				  real ar)	// relative radial acceleration
{
    // Estimate the time required for two particles, of given relative
    // velocity and acceleration, to decrease their separation by dr.
    // Use the following time scales:
    //
    //		t1  =  dr / |vr|		("crossing" time)
    //		t2  =  sqrt (2 * dr / |ar|)	("free-fall" time)
    //
    // If the particles are in the same clump, then ar is the two-body
    // acceleration and t2 really is the free-fall time.  Otherwise,
    // ar is the relative acceleration in the global cluster field.

    if (dr <= 0) return 0;

    real dt = _INFINITY_;

    if (vr < 0) dt = -dr / vr;
    if (ar < 0) dt = fmin(dt, sqrt (-2 * dr / ar));

    return dt;
}

static inline real dt_perturbers(hdyn *bcm, hdyn*& bmin)
{
    // Return an estimate of the time needed for any particle to
    // exceed the threshhold for unperturbed multiple motion of
    // the multiple CM bcm.

    if (!bcm->get_allow_full_unperturbed())	 // turn off extended motion
	return -1;

    real t_min = _INFINITY_;
    bmin = NULL;
    
    hdyn *od = bcm->get_oldest_daughter();

    if (od && od->get_kepler()) {

	real scale = od->get_kepler()->get_semi_major_axis();
	real mass = od->get_kepler()->get_total_mass();

	for_all_daughters(hdyn, bcm->get_parent(), bb)
	    if (bb != bcm) {

		// Estimate the time needed for node bb to perturb bcm
		// at level gamma.  The critical distance for bb to be
		// such a perturber is rpert, defined by
		//
		//	  2 * m_bb * scale		 mass
		//	 ------------------  =  gamma * -------
		//	       rpert^3			scale^2
		//
		// or
		//
		//    rpert  =  gamma^{-1/3} * scale
		//    			     * (2 * m_bb / mass)^{1/3}.

		real rpert = bcm->get_gamma_inv3() * scale
				        * pow(2*bb->get_mass()/mass, 1./3);

		// Time estimate is a combination of the free-fall and
		// crossing times.

		vec dpos = bb->get_pos() - bcm->get_pos();
		vec dvel = bb->get_vel() - bcm->get_vel();
		vec dacc = bb->get_acc() - bcm->get_acc();

		real dr = abs(dpos);
		real vr = dvel * dpos / dr;
		real ar = dacc * dpos / dr;

		real tr = time_to_radius(dr - rpert, vr, ar);

		if (tr < t_min) {
		    t_min = tr;
		    bmin = bb;
		}
		if (tr <= 0) break;
	    }

    } else

	t_min = -1;

    return t_min;
}

static inline real relative_energy(hdyn *bi, hdyn *bj)
{
    // Return the relative energy (per unit reduced mass) of nodes bi
    // and bj.  Assume they have the same parent.

    return 0.5*square(bi->get_vel()-bj->get_vel())
	      -(bi->get_mass()+bj->get_mass())
		  / abs(bi->get_pos()-bj->get_pos());
}

static void adjust_binary(hdyn *cm, real& tidal_potential, real dphi_tidal)
{
#if 0
    // EXPERIMENTAL CODE: Try to absorb the accumulated tidal error
    // dphi_tidal by rotating the binary without changing its core
    // orbital elements.

    // *** ---> Looks like changing the orientation of the binary ***
    // *** ---> generally doesn't work.                           ***

    hdyn *root = cm->get_root();
    hdyn *od = cm->get_oldest_daughter();
    hdyn *yd = od->get_younger_sister();
    real fac = yd->get_mass()/cm->get_mass();

    // Save the component positions, just in case.

    vec xo = od->get_pos();
    vec xy = yd->get_pos();
    vec dx = xy - xo;
    real sep = abs(dx);

    // Calculate the unit vector in the direction of the acceleration.
    // To the extent that the acceleration is predominantly due to the
    // nearest neighbor, we expect the tidal potential to scale as
    //
    //		sep^2 [3 cos^2(theta) - 1]
    //
    // where theta is the angle between the separation vector dx and
    // the external acceleration.  Specifically, this function is
    // linear in cos^2(theta) and monotonic decreasing with theta for
    // theta between 0 and pi/2.

    // Construct unit vectors orb in the direction of the initial
    // orbital separation, acc in the direction of the acceleration, n
    // perpendicular to the two, and t such that (acc, t, n) form a
    // right-handed triad.

    vec orb = dx/sep;
    vec acc = cm->get_acc() / abs(cm->get_acc());
    vec n = orb ^ acc;
    n /= abs(n);
    vec t = -acc^n;

    real c20 = pow(orb*acc, 2);

    // Starting point: cos(theta) = c0 gives a tidal potential of
    // tidal_potential.  We want to find a new value of theta that
    // gives tidal_potential - dphi_tidal.  Return the new value of
    // the tidal potential (or at least the best we can do).

    real phi_target = tidal_potential - dphi_tidal;
    real dphi0 = dphi_tidal;
    PRC(c20); PRL(tidal_potential);

    PRC(tidal_potential); PRL(phi_target);

    // Try to bracket the solution.  Approximate limits on the
    // possible values of the tidal potential are obtained with theta
    // = 0 (lower) and theta = pi/2 (upper).

    real c21, dphi1;
    // if (dphi_tidal < 0) {

    // Rotate perpendicular to the acceleration.

    od->set_pos(-fac*sep*t);
    yd->set_pos((1-fac)*sep*t);
    c21 = 0;
    dphi1 = get_tidal_potential(root, od, yd, cm, false);
    PRC(c21); PRL(dphi1);

    // } else {

    // Rotate parallel to the acceleration.

    od->set_pos(-fac*sep*acc);
    yd->set_pos((1-fac)*sep*acc);
    c21 = 1;
    dphi1 = get_tidal_potential(root, od, yd, cm, false);
    // }

    PRC(c21); PRL(dphi1);

    // Restore the positions.

    od->set_pos(xo);
    yd->set_pos(xy);
#endif
}

bool extend_or_end_binary(hdyn*& bi, bool verbose)
{
    // See if unperturbed motion can be extended, or terminate it.

    if (!bi) return false;
    hdyn *od = bi->get_oldest_daughter();
    if (!od) return false;
    hdyn *yd = od->get_younger_sister();
    if (!yd) return false;
    kepler *k = od->get_kepler();
    if (!k) err_exit("smallN: binary with no kepler.");

    // Correct energy and angular momentum after partial unperturbed
    // motion before considering whether to extend the motion.
    // Recompute the tidal potential in all cases.

    real old_tidal_potential = k->tidal_potential;
    real new_tidal_potential = get_tidal_potential(bi->get_parent(), od, yd, bi,
						   false);
    real dphi_tidal = new_tidal_potential - old_tidal_potential;

    if (!od->get_fully_unperturbed()) {

	if (verbose) {
	    cout << endl << "correcting binary " << bi->get_index()
		 << " at time " << bi->get_system_time()
		 << endl << flush;
	    PRL(total_energy(bi->get_parent()));
	}

	// Compute the change in the tidal potential (components -
	// CM) due to the rest of the system.

	if (verbose) {
	    PRC(old_tidal_potential); PRL(new_tidal_potential);
	    PRL(dphi_tidal);
	    int p = cout.precision(12);
	    PRL(relative_energy(od, yd));
	    cout.precision(p);
	}

	// Apply the changes to the kepler structure, then transmit
	// it to the components, to keep the kepler consistent in
	// case it is needed later.

	// Energy change is -dphi_tidal, or -dphi_tidal/mu in kepler.

	real mu = od->get_mass()*yd->get_mass()/bi->get_mass();

#if 0
	// Old code simply absorbed dphi_total into the binary
	// energy by rescaling the relative velocities of the
	// components.  This changes the angular momentum by an
	// amount proportional to vfac.

	vec vrel = k->get_rel_vel();
	real vfac = sqrt(1 - 2*dphi_tidal/(mu*square(vrel)));
	k->set_rel_vel(vfac*vrel);
	k->initialize_from_pos_and_vel();
	
	if (verbose) {
	    cout << "corrected component velocities: "; PRL(vfac);
	}
#else
	// New code corrects the energy by -dphi_tidal/mu and the
	// angular momentum by k->delta_angular_momentum, reinitalizes
	// the orbit, and updates the components.

	k->set_energy(k->get_energy() - dphi_tidal/mu);
	k->set_angular_momentum(k->get_angular_momentum()
				+ k->delta_angular_momentum);
	k->initialize_from_integrals_and_separation();

	if (verbose) {
	    vec cmpos = bi->get_pos();
	    vec acc = get_cm_acc(bi->get_parent(), bi, od, cmpos);
	    real dx = k->get_separation();
	    vec cmpos1 = cmpos + dx*k->get_longitudinal_unit_vector();
	    vec acc1 = get_cm_acc(bi->get_parent(), bi, od, cmpos1);
	    real daydx = (acc1-acc)*k->get_transverse_unit_vector()/dx;
	    PRL(daydx);
	    PRL(k->delta_angular_momentum);
	    int p = cout.precision(12);
	    cout << "new energy = " << k->get_energy()
		 << " angular momentum = " << k->get_angular_momentum()
		 << endl << flush;
	    cout.precision(p);
	}
#endif

	// Update the components using standard tools.

	advance_components_to_time(bi, bi->get_system_time());
	update_components_from_pred(bi);

	if (verbose) {
	    k->print_all();
	    int p = cout.precision(12);
	    PRL(relative_energy(od, yd));
	    cout.precision(p);
	    PRL(total_energy(bi->get_parent()));
	}

    } else {

	// Experimental code to try to absorb the change in the tidal
	// potential.  Not necessary if we freeze and resolve fully
	// unperturbed binaries.

	adjust_binary(bi, new_tidal_potential, dphi_tidal);
    }

    // Check whether extension is possible.  Extension will be for an
    // integral number of periods, at fixed binary phase (components
    // should be receding).

    hdyn *bmin = NULL;
    real dtp = 0.5*dt_perturbers(bi, bmin);	// 0.5 is conservative

    real period = k->get_period();
    if (bi->get_allow_full_unperturbed() && dtp >= period) {
	float np = floor(dtp/period);
	od->set_t_pred(od->get_t_pred()+np*period);
	if (verbose) {
	    cout << endl << "extending binary " << bi->get_index()
		 << " by " << np << " periods = " << np*period
		 << " at time " << bi->get_system_time()
		 << endl << flush;
	    if (bmin)
		cout << "NN is " << bmin->get_index()
		     << ": dist = " << abs(bmin->get_pos()-bi->get_pos())
		     << " sma = " << k->get_semi_major_axis()
		     << endl << flush;
	    PRL(total_energy(bi->get_parent()));
	}
	bi->set_fully_unperturbed(true);
	od->set_fully_unperturbed(true);

	// Save the tidal potential for correction at the end of this
	// segment.

	k->tidal_potential = new_tidal_potential;

	return false;
    }

    // Terminate the binary.

    if (verbose) {
	cout << endl << "terminating binary " << bi->get_index()
	     << " at time " << bi->get_system_time()
	     << endl << flush;
	if (bmin)
	    cout << "NN is " << bmin->get_index()
		 << ": dist = " << abs(bmin->get_pos()-bi->get_pos())
		 << " sma = " << k->get_semi_major_axis()
		 << endl << flush;
    }

    // Clear kepler pointers.

    od->set_kepler(NULL);
    yd->set_kepler(NULL);
    delete k;

    // Include the center of mass position and velocity in the
    // component quantities.

    od->inc_pos(bi->get_pos());
    od->inc_vel(bi->get_vel());
    yd->inc_pos(bi->get_pos());
    yd->inc_vel(bi->get_vel());

    // Update the tree.

    add_node(od, bi);
    add_node(yd, bi);
    detach_node(bi);

    // Delete the CM node (careful to detach the daughters to
    // avoid recursive deletion!).

    bi->set_oldest_daughter(NULL);
    delete(bi);

    // Clear flags.

    od->set_fully_unperturbed(false);
    yd->set_fully_unperturbed(false);

    // Make sure pos and pred_pos agree.  Note that init_pred sets
    // t_pred = time as well as pred_pos = pos.

    od->init_pred();
    yd->init_pred();

    // Change bi to yd, so next node is correct in loop below.

    bi = yd;

    if (verbose) PRL(total_energy(bi->get_parent()));
    return true;
}
