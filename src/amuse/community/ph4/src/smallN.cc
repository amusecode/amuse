
// Self-contained few-body integrator, using a fourth-order Hermite
// scheme, and incorporating a modified unperturbed treatment of close
// approaches, with time symmetrization.  Binary trees are used to
// implement unperturbed motion.  The basic data structure is a flat
// tree, with some nodes possibly having substructure representing
// unperturbed binaries.
//
// Shamelessly copied from Starlab::kira_smallN.C by the author.
// Retained Starlab data structures to simplify code importation.
//
// Notes:
//
//	* the system is integrated in isolation -- no external tidal
//	  field,
//
//	* the integration ends after a specified time, or when a size
//	  criterion is exceeded, or when the analysis routines
//	  determine that the interaction is over,
//
//	* partially unperturbed nodes (pericenter reflection) are
//	  advanced in the center of mass approximation -- they are
//	  *not* resolved in computing forces, but *are* resolved when
//	  computing the total energy, both to include the binary
//	  energy and to make the tidal correction continuous,
//
//	* partially unperturbed nodes are predicted using the kepler
//	  package prior to being resolved,
//
//	* both the energy and the angular momentum are corrected at
//	  the end of partial unperturbed motion, using expressions
//	  integrated along the unperturbed (kepler) orbit,
//
//	* fully unperturbed nodes are resolved into components when
//	  computing both the acceleration *and* the potential, in
//	  order to properly incorporate the tidal potential; no
//	  correction is applied at the end of the unperturbed motion,
//
//	* fully unperturbed nodes are *not* predicted prior to being
//	  resolved -- their positions and velocities are held fixed
//	  during the unperturbed motion.
//
// Author: Steve McMillan.
//
// TO DO:
//
//	* implement treatment of unperturbed multiples, based on the
//	  output of check_structure() -- same approach as unperturbed
//	  binaries,
//
//	* include secular changes to orbital eccentricity during fully
//	  unperturbed motion,
//
//	* add an external tidal field.
//
// Externally visible functions:
//
//	void advance_components_to_time(hdyn *bi, real t);
//	void update_components_from_pred(hdyn *bi);
//	int smallN_evolve(hdyn *b,
//			  real t_end = _INFINITY_,
//			  real dt_check = _INFINITY_,
//			  real break_r2 = _INFINITY_,
//			  real dt_log = _INFINITY_,
//			  int  verbose = 0);

#include "hdyn.h"

#ifndef TOOLBOX

// Global pointers to the closest pair (shortest mutual time step).

static hdyn *bi_min = NULL, *bj_min = NULL;



static inline real kepler_step_sq(real distance2, real distance3,
				  real mass, real vel2, hdyn *b)
{
    real dtff2 = 0.5 * distance3 / mass;
    real dtv2 = distance2 / abs(vel2);
    real dt2 = fmin(dtff2, dtv2);

#if 0

    // Experimental!  But really unnecessary to suppress the maximum
    // error during an orbit since the symmetrization means no
    // cumulative effect.

    const real dt2crit = 1.e-2;
    if (b->get_system_time() > 10000 && dt2 < dt2crit) {
	real dt2fac = pow(dt2/dt2crit, 0.1);
	if (dt2fac < 0.25) dt2fac = 0.25;
	dt2 *= dt2fac;
    }
#endif

    return dt2;
}

static inline real kepler_step_sq(real distance, real mass, real vel2, hdyn *b)
{
    real d2 = distance*distance;
    real d3 = distance*d2;
    return kepler_step_sq(d2, d3, mass, vel2, b);
}

typedef struct {real m; vec x; vec v;} body;

static body bbi[128], bbj[128];	// use fixed arrays for efficiency;
				// 128 should be larger than n

static inline real get_pairwise_acc_and_jerk_CPT(hdyn *bi, hdyn *bj,
						 vec &force, vec &jerk)
{
    // Compute the force and jerk on top-level node bi due to
    // top-level node bj by summing all pairwise component forces.
    // Return the unscaled time step appropriate to the minimum
    // distance between any component of bi and any component of bj.

    force = jerk = 0;
    real min_distance2 = _INFINITY_;
    real min_distance3 = _INFINITY_;

#if 1

    // NOTE: assuming here that trees are just a single level deep.
    // MUST be generalized (see below) if we choose to allow
    // unperturbed multiples.

    vec delx = 0, delv = 0;
    if (bi->is_parent()) {
	delx -= bi->get_pred_pos();
	delv -= bi->get_pred_vel();
    }
    if (bj->is_parent()) {
	delx += bj->get_pred_pos();
	delv += bj->get_pred_vel();
    }

    for_all_leaves(hdyn, bi, bbi) {

	// Compute the force and jerk on leaf bbi due to all leaves under bj.

	vec iforce = 0, ijerk = 0;
	for_all_leaves(hdyn, bj, bbj) {
	    vec dx = bbj->get_pred_pos() - bbi->get_pred_pos() + delx;
	    vec dv = bbj->get_pred_vel() - bbi->get_pred_vel() + delv;
	    real distance2 = dx*dx;
	    real distance  = sqrt(distance2);
	    real distance3 = distance2*distance;

	    iforce += bbj->get_mass() * dx / distance3;
	    ijerk  += bbj->get_mass() * (dv / distance3
					 - 3*dx*(dx*dv)/(distance3*distance2));

	    if (distance2 < min_distance2) {
		min_distance2 = distance2;
		min_distance3 = distance3;
	    }
	}

	force += bbi->get_mass() * iforce;
	jerk  += bbi->get_mass() * ijerk;
    }

#else

    // More general: create arrays containing the masses, absolute
    // positions, and velocities of the leaves of bi and bj.
    // *** NOT TESTED!! ***

    int ni = 0;
    for_all_leaves(hdyn, bi, bb) {
	vec xi = bb->get_pred_pos();
	vec vi = bb->get_pred_vel();
	hdyn *pi = bb->get_parent();
	while (pi->get_parent()) {
	    xi += pi->get_pred_pos();
	    vi += pi->get_pred_vel();
	    pi = pi->get_parent();
	}
	bbi[ni].m = bb->get_mass();
	bbi[ni].x = xi;
	bbi[ni++].v = vi;
    }

    int nj = 0;
    for_all_leaves(hdyn, bj, bb) {
	vec xj = bb->get_pred_pos();
	vec vj = bb->get_pred_vel();
	hdyn *pj = bb->get_parent();
	while (pj->get_parent()) {
	    xj += pj->get_pred_pos();
	    vj += pj->get_pred_vel();
	    pj = pj->get_parent();
	}
	bbj[nj].m = bb->get_mass();
	bbj[nj].x = xj;
	bbj[nj++].v = vj;
    }

    for (int i = 0; i < ni; i++) {

	// Compute the force and jerk on leaf bbi[i] due to all leaves
	// in bbj.

	real mi = bbi[i].m;
	vec  xi = bbi[i].x;
	vec  vi = bbi[i].v;
	vec iforce = 0, ijerk = 0;

	for (int j = 0; j < nj; j++) {
	    real mj = bbj[j].m;
	    vec  dx = bbj[j].x - xi;
	    vec  dv = bbj[j].v - vi;
	    real distance2 = dx*dx;
	    real distance  = sqrt(distance2);
	    real distance3 = distance2*distance;

	    iforce += mj * dx / distance3;
	    ijerk  += mj * (dv / distance3
			     - 3*dx*(dx*dv)/(distance3*distance2));

	    if (distance2 < min_distance2) {
		min_distance2 = distance2;
		min_distance3 = distance3;
	    }
	}

	force += mi * iforce;
	jerk  += mi * ijerk;
    }

#endif

    // Time step criterion:

    real timestep2 = kepler_step_sq(min_distance2, min_distance3,
				    bj->get_mass() + bi->get_mass(),
				    square(bj->get_pred_vel()
					   - bi->get_pred_vel()),
				    bi);
    return timestep2;
}

static inline real get_pairwise_acc_and_jerk_CM(hdyn *bi, hdyn *bj,
						vec &force, vec &jerk)
{
    // Compute the force and jerk on top-level node bi due to
    // top-level node bj in the center of mass approximation.  Return
    // the unscaled time step appropriate to the distance between bi
    // and bj.

    real mi = bi->get_mass();
    real mj = bj->get_mass();
    vec dx = bj->get_pred_pos() - bi->get_pred_pos();
    vec dv = bj->get_pred_vel() - bi->get_pred_vel();
    real distance2 = dx*dx;
    real distance  = sqrt(distance2);
    real distance3 = distance2*distance;

    force = mi*mj * dx/distance3;
    jerk  = mi*mj * (dv/distance3 - 3*dx*(dx*dv) / (distance3*distance2));

    // Time step criterion:

    return kepler_step_sq(distance2, distance3, mi+mj, dv*dv, bi);
}

static inline real get_pairwise_acc_and_jerk(hdyn *bi, hdyn *bj,
					     vec &force, vec &jerk)
{
    // Compute the force and jerk on top-level node bi due to
    // top-level node bj.  Return the unscaled time step appropriate
    // to the distance between bi and bj.  Switch between pure
    // center-of-mass approximation and resolution of unperturbed
    // binaries.

#if 1
    return get_pairwise_acc_and_jerk_CPT(bi, bj, force, jerk);
#else
    return get_pairwise_acc_and_jerk_CM(bi, bj, force, jerk);
#endif
}



real calculate_top_level_acc_and_jerk(hdyn *b)
{
    // Compute the acc and jerk on all top-level nodes.  All nodes are
    // resolved into components for purposes of computing the acc and
    // jerk.  Return the minimum time step associated with any
    // top-level pair.

    for_all_daughters(hdyn, b, bi) {
	bi->set_acc(0);
	bi->set_jerk(0);
    }

    real min_timestep2 = _INFINITY_;
    bi_min = bj_min = NULL;

    for_all_daughters(hdyn, b, bi)
	for (hdyn *bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {

	    vec force, jerk;
	    real timestep2 = get_pairwise_acc_and_jerk(bi, bj, force, jerk);

	    real mi = 1 / bi->get_mass();
	    bi->inc_acc(mi*force);
	    bi->inc_jerk(mi*jerk);

	    real mj = 1 / bj->get_mass();
	    bj->inc_acc(-mj*force);
	    bj->inc_jerk(-mj*jerk);

	    if (timestep2 < min_timestep2) {
		min_timestep2 = timestep2;
		bi_min = bi;
		bj_min = bj;
	    }
	}

    real dt = b->get_eta()*sqrt(min_timestep2);		// natural time step
    return dt;
}



//----------------------------------------------------------------------
//
// Advance the components of a binary to the specified time.

void advance_components_to_time(hdyn *bi, real t)	// unperturbed
{							// "predictor"
    hdyn *od = bi->get_oldest_daughter();

    if (od) {

	// Only way this can occur is for unperturbed motion.
	// Check and flag if no kepler found.

	kepler *k = od->get_kepler();

	if (!k) err_exit("smallN_evolve: daughter node with no kepler.");

	else {

	    // Advance the components to time t.  We won't actually
	    // integrate the internal motion, so could just set pos =
	    // pred_pos (etc.) here (no corrector), but for now defer
	    // that until the corrector step.

	    // Note that partially unperturbed binaries are treated in
	    // the center of mass approximation for purposes of force
	    // calculation, so their components are resolved only when
	    // computing the total energy.

	    // In the case of fully unperturbed motion, we resolve the
	    // components but don't predict them -- instead, keep the
	    // components' positions and velocities unchanged.

	    hdyn *yd = od->get_younger_sister();

	    if (!od->get_fully_unperturbed())
		k->transform_to_time(t);

	    real fac = yd->get_mass()/bi->get_mass();
	    od->set_pred_pos(-fac*k->get_rel_pos());
	    od->set_pred_vel(-fac*k->get_rel_vel());
	    yd->set_pred_pos((1-fac)*k->get_rel_pos());
	    yd->set_pred_vel((1-fac)*k->get_rel_vel());
	}
    }
}


void update_components_from_pred(hdyn *bi)		// unperturbed
{							// "corrector"
    // Called after a step is completed.

    hdyn *od = bi->get_oldest_daughter();

    if (od) {
	od->set_pos(od->get_pred_pos());
	od->set_vel(od->get_pred_vel());
	hdyn *yd = od->get_younger_sister();
	yd->set_pos(yd->get_pred_pos());
	yd->set_vel(yd->get_pred_vel());
    }
}



void hdyn::correct_acc_and_jerk(const real new_dt,
				const real prev_new_dt)
{
    // Correct the values of acc and jerk from time prev_new_dt to
    // new_dt.  We simply fit a polynomial to old_acc and old_jerk at
    // time 0 and acc and jerk at time prev_new_dt, then evaluate it
    // at time new_dt.

    if (new_dt == prev_new_dt) return;

    real dt_off = new_dt - 0.5 * prev_new_dt;
                                    // offset from midpoint of prev_new_dt step
    real theta = 0.5 * prev_new_dt;
    real tau = dt_off / theta;	    // equals 1 if new_dt = prev_new_dt

    real inv_theta = 1 / theta;
    real tau2 = tau * tau;
    real tau3 = tau2 * tau;

    vec prev_acc = acc;
    vec prev_jerk = jerk;

    acc = 0.25 * (old_acc * (2 - 3 * tau + tau3)
		   + prev_acc * (2 + 3 * tau - tau3)
		   + old_jerk * theta * (1 - tau - tau2 + tau3)
		   + prev_jerk * theta * (-1 - tau + tau2 + tau3));

    jerk = 0.25 * (old_acc * inv_theta * 3 * (-1 + tau2)
		    + prev_acc * inv_theta * 3 * (1 - tau2)
		    + old_jerk * (-1 - 2*tau + 3*tau2)
		    + prev_jerk * (-1 + 2*tau + 3*tau2));
}

void hdyn::correct_pos_and_vel(const real new_dt)
{
    // Apply a corrector in the form presented by Hut et al. (1995).
    // The "pred" quantities are those at the end of the step.

    real new_dt2 = new_dt * new_dt;

    pred_vel = vel + new_dt * (acc + old_acc)/2
		   - new_dt2 * (jerk - old_jerk)/12;
    pred_pos = pos + new_dt * (pred_vel + vel)/2
		   - new_dt2 * (acc - old_acc)/12;
}



// Take the next step.  Return the actual system time at the end of
// the step, after symmetrization if specified.  Starting time step is
// dt.

static real take_a_step(hdyn *b,	// root node
			real &dt)	// natural step at start/end
{
    real t = b->get_system_time();

    // Impose an absolute limit on the step if unperturbed motion is
    // underway.

    real dt_unpert_limit = _INFINITY_;
    for_all_daughters(hdyn, b, bi) {
	hdyn *od = bi->get_oldest_daughter();
	if (od) {

	    // Kepler termination time is stored in od->t_pred.

	    real dt_term = od->get_t_pred() - t;
	    if (dt_term < dt_unpert_limit) dt_unpert_limit = dt_term;
	}
    }

    // If unperturbed motion is due to end, force t to that time and
    // don't iterate.

    int n_iter = b->get_n_iter();
    if (dt >= dt_unpert_limit) {
	dt = dt_unpert_limit;
	n_iter = 0;
    }

    // Predict all to time t + dt -- set "pred" quantities throughout.

    predict_loworder_all(b, t+dt);	// top-level nodes
    for_all_daughters(hdyn, b, bi) {
	bi->store_old_force();
	if (bi->is_parent())		// components
	    advance_components_to_time(bi, t+dt);
    }

    // Compute forces and correct, iterating if desired.  Note that
    // dt is the natural time step associated with the state of the
    // system on entry, and the natural time step of the updated
    // system on exit.  We use "pred" quantities to represent the
    // current iterate of quantities the end of the time step (even
    // after correction; these values are called "new" in the sdyn and
    // sdyn3 versions of the symmetrization scheme.)  The acc and jerk
    // at the start of the step are "old_acc" etc.; acc and jerk are
    // defined at the predicted time.

    real new_dt = dt;
    real end_point_dt = dt;

    // Ultimately, new_dt will be the actual step taken, while
    // end_point_dt will be the natural time step at the end of the
    // step.

    for (int i = 0; i <= n_iter; i++) {

	real prev_new_dt = new_dt;
	end_point_dt = calculate_top_level_acc_and_jerk(b);

	// Obtain the next iterate of the time step (actually, do two
	// at a time).

	if (i < n_iter) {

	    // First iteration is just Newton's method for solving
	    //
	    //	0.5*(step(t) + step(t+dt))  =  dt

	    new_dt = 0.5 * (dt + end_point_dt);

	    // Second iteration assumes step(t) is differentiable.

	    new_dt = dt + 0.5 * (end_point_dt - dt) * (new_dt/prev_new_dt);

	    // See if we have exceeded the unperturbed limit.

	    if (new_dt >= dt_unpert_limit) {
		new_dt = dt_unpert_limit;
		n_iter = 0;
	    }
	}

	// Extrapolate acc and jerk to the end of the new step, and
	// apply the corrector for this iteration.

	for_all_daughters(hdyn, b, bi) {
	    if (new_dt != prev_new_dt)
		bi->correct_acc_and_jerk(new_dt, prev_new_dt);
	    bi->correct_pos_and_vel(new_dt);	// sets pred_xxx
	    if (bi->is_parent()) advance_components_to_time(bi, t+new_dt);
	}
     }

    // Complete the step.

    for_all_daughters(hdyn, b, bi) {
	bi->set_pos(bi->get_pred_pos());
	bi->set_vel(bi->get_pred_vel());
	if (bi->is_parent()) update_components_from_pred(bi);
    }

    dt = end_point_dt;
    return t + new_dt;
}



static void print_most_bound(hdyn *b)
{
    real Emin = _INFINITY_, mu = 0;
    int imin = 0, jmin = 0;
    for_all_leaves(hdyn, b, bi) {
	real mi = bi->get_mass();
	vec xi = bi->get_pos();
	vec vi = bi->get_vel();
	if (!bi->is_top_level_node()) {
	    xi += bi->get_parent()->get_pos();
	    vi += bi->get_parent()->get_vel();
	}
	for (hdyn *bj = bi->next_node(b); bj != NULL;
	     bj = bj->next_node(b))
	    if (bj->is_leaf()) {
		real mj = bj->get_mass();
		vec xj = bj->get_pos();
		vec vj = bj->get_vel();
		if (!bj->is_top_level_node()) {
		    xj += bj->get_parent()->get_pos();
		    vj += bj->get_parent()->get_vel();
		}
		real Eij = 0.5*square(vj-vi) - (mi+mj)/abs(xj-xi);
		if (Eij < Emin) {
		    Emin = Eij;
		    imin = bi->get_index();
		    jmin = bj->get_index();
		    mu = mi*mj/(mi+mj);
		}
	    }
    }
    cout << imin << " " << jmin << " " << mu*Emin << " ";
}

static real top_level_energy(hdyn *b)
{
    // Return the total energy of the system, treating all top-level
    // nodes in the senter of mass approximation.

    real kin = 0, pot = 0;
    for_all_daughters(hdyn, b, bi) {
	real mi = bi->get_mass();
	vec xi = bi->get_pos();
	kin += mi*square(bi->get_vel());
	real ppot = 0;
	for (hdyn *bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	    ppot -= bj->get_mass()/abs(bj->get_pos()-xi);
	pot += mi*ppot;
    }
    return kin/2 + pot;
}

real get_energies(hdyn *b, real& kin, real& pot)
{
    // Return the total energy of the system, resolving all binary
    // nodes.

    kin = pot = 0;
    for_all_leaves(hdyn, b, bi) {
	real mi = bi->get_mass();
	vec xi = bi->get_pos();
	vec vi = bi->get_vel();

	// Make all quantities absolute (not very efficient).

	hdyn *pi = bi->get_parent();
	while (pi->get_parent()) {
	    xi += pi->get_pos();
	    vi += pi->get_vel();
	    pi = pi->get_parent();
	}
	kin += mi*square(vi);
	real ppot = 0;
	for (hdyn *bj = bi->next_node(b); bj != NULL;
	     bj = bj->next_node(b))
	    if (bj->is_leaf()) {
		vec dx = bj->get_pos() - xi;
		hdyn *pj = bj->get_parent();
		while (pj->get_parent()) {
		    dx += pj->get_pos();
		    pj = pj->get_parent();
		}
		ppot -= bj->get_mass()/abs(dx);
	    }
	pot += mi*ppot;
    }
    kin /= 2;
    return kin + pot;
}

real total_energy(hdyn *b)
{
    real kin, pot;
    return get_energies(b, kin, pot);
}

static void log_output(hdyn *b, int n_steps)
{
    // One-line essential output.

    static real E0 = 0;
    real E = total_energy(b);
    if (E0 == 0) E0 = E;
    real Etop = top_level_energy(b);
    int n_leaves = 0;
    for_all_leaves(hdyn, b, bb) n_leaves++;
    int n_unp = 0;
    for_all_daughters(hdyn, b, bb) if (bb->is_parent()) n_unp++;

    // real phi_tidal = 0;
    // for_all_daughters(hdyn, b, bb) {
    // 	if (bb->is_parent()) {
    // 	    hdyn *od = bb->get_oldest_daughter();
    // 	    hdyn *yd = od->get_younger_sister();
    // 	    phi_tidal += get_tidal_potential(b, od, yd, bb, false);
    // 	    phi_tidal -= od->get_kepler()->tidal_potential;
    // 	}
    // }

    int p = cout.precision(10);
    cout << "smallN%% " << b->get_system_time() << " ";
    cout.precision(p);

    cout << E << " " << E-E0
      	 << " " << Etop << " ";
      //	 << " " << E-E0-phi_tidal << " ";
    print_most_bound(b);
    real user, sys;
    get_cpu_time(user, sys);
    cout << n_leaves << " " << n_unp << " " << user << " " << n_steps;
    cout << endl << flush;
}

// Evolve the system to time t_end, using the input data and settings
// without modification.  Assume that we start with a flat tree.  Stop
// and return if
//
//	(1) t >= t_end (checked at start of each iteration),
//	(2) any particle gets too far (break_r) from the origin
//	    (checked after every NCHECK steps)
//
// The return value is 1 or 2 for these two cases.  A return value
// of 0 means that the interaction is over, as determined by
// check_structure().

#define NCHECK 100

static void two_body(hdyn *b, real time, real radius)
{
    // Follow the motion of two bodies using kepler.  Advance to time
    // or (outgoing) radius, whichever is sooner.

    hdyn *od = b->get_oldest_daughter();
    if (!od) return;
    hdyn *yd = od->get_younger_sister();
    if (!yd) return;

    kepler *k = hdyn_to_kepler(b);
    // k->print_all(cout);

    if (k->get_energy() >= 0) {
	k->transform_to_time(time);
	if (k->get_separation() > radius
	    && k->get_rel_pos()*k->get_rel_vel() > 0)
	    k->return_to_radius(radius);
    } else {
	real apo = k->get_apastron();
	if (apo <= radius)
	    k->transform_to_time(time);
	else {
	    k->advance_to_radius(radius);
	    if (k->get_time() > time)
		k->transform_to_time(time);
	    else if (k->get_rel_pos()*k->get_rel_vel() < 0) {
		real peri = k->get_apastron();
		k->advance_to_radius(peri + 0.999*(radius-peri));
		k->advance_to_radius(radius);
	    }
	}
    }

    b->set_system_time(k->get_time());

    // Update the daughters with the new orbital data.

    real total_mass = od->get_mass() + yd->get_mass();
    vec cmpos = (od->get_mass()*od->get_pos()
		 + yd->get_mass()*yd->get_pos()) / total_mass;
    vec cmvel = (od->get_mass()*od->get_vel()
		 + yd->get_mass()*yd->get_vel()) / total_mass;
    real fac = yd->get_mass() / total_mass;

    od->set_pos(cmpos-fac*k->get_rel_pos());
    od->set_vel(cmvel-fac*k->get_rel_vel());
    yd->set_pos(cmpos+(1-fac)*k->get_rel_pos());
    yd->set_vel(cmvel+(1-fac)*k->get_rel_vel());

    delete k;
}

void spaces(int n) {for (int i = 0; i < n; i++) cout << " ";}

void print(hdyn *b, int level = 0)
{
    spaces(4*level);
    cout << b->get_index() << "  " << "mass = " << b->get_mass();
    if (b->get_kepler()) cout << "  kepler";
    cout << endl;
    spaces(4*level);
    cout << "    pos = " << b->get_pos() << endl;
    spaces(4*level);
    cout << "    vel = " << b->get_vel() << endl;
    for_all_daughters(hdyn, b, bb)
	print(bb, level+1);
}

int smallN_evolve(hdyn *b,
		  real t_end,		// default = _INFINITY_
		  real break_r2,	// default = _INFINITY_
		  real dt_check,	// default = _INFINITY_
		  real dt_log,		// default = _INFINITY_
		  int verbose)		// default = 0
{
    set_kepler_tolerance(2);	// energy corrections may force orbital
				// separations outside allowed limits

    // print(b);

    // Treat special cases (that may come from AMUSE).

    int n_leaves = 0;
    for_all_leaves(hdyn, b, bi) n_leaves++;

    // cout << "In smallN_evolve: "; PRL(n_leaves);

    if (n_leaves == 1)
	return 0;
    else if (n_leaves == 2) {
	// cout << "smallN: two-body encounter" << endl << flush;
	two_body(b, t_end, sqrt(break_r2));
	return 0;
    }

    // cout << "smallN: direct integration" << endl << flush;
    int n_steps = 0;
    for_all_daughters(hdyn, b, bi)
	bi->init_pred();

    if (b->get_cm_index() <= 0) {	// default is -1

	// Make up a value that is well removed from existing indices.

	cout << "auto definition of cm_index: ";
	int cm_index = 0;
	for_all_daughters(hdyn, b, bi) {
	    int i = bi->get_index();
	    if (i > cm_index) cm_index = i;
	}
	cm_index = pow(10, (int)log10((real)cm_index+1) + 2.) + 1;
	PRL(cm_index);
	b->set_cm_index(cm_index);
    }

    real dt = calculate_top_level_acc_and_jerk(b);

    real tmp;
    get_cpu_time(tmp, tmp);

    real t_log = b->get_system_time();
    if (dt_log < _INFINITY_) log_output(b, n_steps);
    t_log += dt_log;

    real t_check = b->get_system_time() + dt_check;

    while (b->get_system_time() < t_end) {

	// Take a step.  Don't set the end time in advance, as the
	// symmetrization process will determine the actual time step.
	// During the entire step, all times refer to the time at the
	// *start* of the step.  The time step dt on entry is the
	// natural time step for the system at the current time.  On
	// return it is the new natural time step for the system at
	// the end of the step.  The return value of the function is
	// the new system time.

	b->set_system_time(take_a_step(b, dt));
	n_steps++;

	// The time step dt was set by bi_min and bj_min during the
	// last acc and jerk calculation.

	// Check second (size) termination criterion.

	if (n_steps%NCHECK == 0)
	    for_all_daughters(hdyn, b, bi) {
		real r2 = square(bi->get_pos());
		if (r2 > break_r2) {
		    cout << "smallN: returning with rmax > break_r"
			 << endl << flush;
		    return 2;
		}
	    }

	// Check for the start of unperturbed motion.  Use various
	// thresholds to avoid this check at the end of every step:
	//
	//	time step dt < dt_crit
	//	bi_min and bj_min are leaves
	//	bi_min and bj_min are within distance r_crit
	//	bi_min and bj_min are approaching
	//	others...
	//
	// Once the check fails, a more clever search would defer
	// further checks until the components have approached each
	// other significantly, but this could entail significant
	// bookeeping.

	bool tree_changed = false;

	// Check for new binary creation.  Only allow one new binary
	// at a time.  The first step is always to the reflection
	// point.

	if (is_unperturbed_and_approaching(b, dt, bi_min, bj_min))
	    tree_changed = create_binary(bi_min, bj_min, verbose > 1);

	// Check for extension or termination of unperturbed motion.
	// By construction, if unperturbed motion is due to end, we
	// should have just taken a step to that time.

	for_all_daughters(hdyn, b, bi)
	    if (bi->get_oldest_daughter()
		&& bi->get_oldest_daughter()->get_t_pred()
		    <= b->get_system_time())
	      tree_changed |= extend_or_end_binary(bi, verbose > 1);

	// Recompute accs, jerks, and the time step, if necessary.

	if (tree_changed) dt = calculate_top_level_acc_and_jerk(b);

	// Basic diagnostic output.

	if (dt_log == 0
	    || (dt_log > 0 && b->get_system_time() >= t_log)) {
	    log_output(b, n_steps);
	    if (dt_log > 0)
		while (b->get_system_time() >= t_log) t_log += dt_log;
	}

	// Structure analysis:

	if (dt_check > 0 && b->get_system_time() >= t_check) {
	    bool over = check_structure(b, break_r2, verbose);
	    if (over) return 0;
	    while (b->get_system_time() >= t_check) t_check += dt_check;
	}
    }

    real rmax2 = 0;
    for_all_daughters(hdyn, b, bi) {
	real r2 = square(bi->get_pos());
	if (r2 > rmax2) rmax2 = r2;
    }
    cout << "smallN: returning with "; PRC(rmax2); PRL(break_r2);

    return 1;
}



#else

#include <string>
#include <unistd.h>

template <class Q>
bool get_quantity(string s, const char *q, Q& value)
{
    size_t i = s.find(q);
    if (i != string::npos) {
	i = s.find("=", i);
	if (i != string::npos) {
	    s = s.substr(i+1);
	    istringstream ss(s);
	    ss >> value;
	    return true;
	}
    }
    return false;
}

void initialize_special(hdyn *b)
{
    // Set up a hierarchical 3-body system.

    real m1 = 0.75, m2 = 0.25, m3 = 1;
    real sma_in = 1, ecc_in = 0.95;
    real sma_out = 10, ecc_out = 0;

    real m12 = m1 + m2;
    real m123 = m12 + m3;

    // Inner binary lies in the x-y plane and has its long axis along
    // the x-axis.  It starts at apastron.

    kepler k_in;
    k_in.set_time(0);
    k_in.set_total_mass(m12);
    k_in.set_semi_major_axis(sma_in);
    k_in.set_eccentricity(ecc_in);
    k_in.set_mean_anomaly(M_PI);	// apastron
    k_in.align_with_axes(1);
    k_in.initialize_from_shape_and_phase();

    // Outer binary has random orientation and phase.

    kepler k_out;
    k_out.set_time(0);
    k_out.set_total_mass(m123);
    k_out.set_semi_major_axis(sma_out);
    k_out.set_eccentricity(ecc_out);
    k_out.set_mean_anomaly(randinter(0, 2*M_PI));
    set_random_orientation(k_out);
    k_out.initialize_from_shape_and_phase();

    hdyn *b1 = new hdyn;
    hdyn *b2 = new hdyn;
    hdyn *b3 = new hdyn;
    b1->set_parent(b);
    b2->set_parent(b);
    b3->set_parent(b);
    b->set_oldest_daughter(b1);
    b1->set_younger_sister(b2);
    b2->set_younger_sister(b3);
    b2->set_older_sister(b1);
    b3->set_older_sister(b2);

    b1->set_index(1);
    b2->set_index(2);
    b3->set_index(3);
    b1->set_mass(m1);
    b2->set_mass(m2);
    b3->set_mass(m3);

    b3->set_pos(m12 * k_out.get_rel_pos() / m123);
    b3->set_vel(m12 * k_out.get_rel_vel() / m123);
    b1->set_pos(-m3 * k_out.get_rel_pos() / m123);
    b1->set_vel(-m3 * k_out.get_rel_vel() / m123);
    b2->set_pos(b1->get_pos());
    b2->set_vel(b1->get_vel());

    // Then set up the inner binary.

    b1->inc_pos(-m2 * k_in.get_rel_pos() / m12);    // rel_pos is from 1 to 2
    b1->inc_vel(-m2 * k_in.get_rel_vel() / m12);
    b2->inc_pos( m1 * k_in.get_rel_pos() / m12);
    b2->inc_vel( m1 * k_in.get_rel_vel() / m12);

    b->set_system_time(0);
    PRL(k_in.pred_advance_to_periastron());
}

void initialize_from_stdin(hdyn *b)
{
    int n = 0;
    hdyn *bp = NULL, *bb = NULL;
    string s;

    while (getline(cin, s)) {
	if (s[0] != ';') {
	    istringstream ss(s);
	    int i = -1;
	    real m;
	    vec x, v;
	    ss >> i >> m >> x >> v;
	    if (i < 0) break;
	    bb = new hdyn();
	    n++;
	    if (bp) {
		bp->set_younger_sister(bb);
		bb->set_older_sister(bp);
	    } else
		b->set_oldest_daughter(bb);
	    bb->set_parent(b);
	    bb->set_index(i);
	    bb->set_mass(m);
	    bb->set_pos(x);
	    bb->set_vel(v);
	    bp = bb;
	} else {
	    real t;
	    if (get_quantity(s, "system_time", t)) b->set_system_time(t);
	    int seed;
	    if (get_quantity(s, "seed", seed)) {PRL(seed); b->set_seed(seed);}
	}
    }
    PRC(n); PRL(b->get_system_time());
}

int main(int argc, char *argv[])
{
    const char *param_string = "a:c:d:fg:r:s:t:Tv:";

    real eta = 0.03;
    real gamma = 1.e-6;
    real break_r2 = _INFINITY_;
    real t_end = 100;
    real dt_log = _INFINITY_;
    real dt_check = 50;
    bool allow_full = true;
    int verbose = 0;
    bool test = false;

    int c;
    extern char *optarg;
    while ((c = getopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'a': eta = atof(optarg);
		      break;
	    case 'c': dt_check = atof(optarg);
		      break;
	    case 'd': dt_log = atof(optarg);
		      break;
	    case 'f': allow_full = !allow_full;
		      break;
	    case 'g': gamma = atof(optarg);
		      break;
	    case 'r': break_r2 = pow(atof(optarg), 2);
		      break;
	    case 's': srandom(atoi(optarg));
		      break;
	    case 't': t_end = atof(optarg);
		      break;
	    case 'T': test = true;
		      break;
	    case 'v': verbose = atoi(optarg);
		      break;
        }            

    hdyn *b = new hdyn;
    if (test)
	initialize_special(b);
    else
	initialize_from_stdin(b);
    b->set_eta(eta);
    b->set_gamma(gamma);
    b->set_allow_full_unperturbed(allow_full);
    int ret = smallN_evolve(b, t_end, break_r2, dt_check, dt_log, verbose);
    cout << "end at time " << b->get_system_time()
	 << "; return status = " << ret << endl;

    return 0;
}

#endif
