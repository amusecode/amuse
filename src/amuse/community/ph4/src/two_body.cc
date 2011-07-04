
// Manage true two-body encounters and multiple encounters.
//
// Global function:
//
//	void jdata::two_body()
//	void jdata::multiple()

#include "hdyn.h"
#include "jdata.h"
#include "scheduler.h"
#include "idata.h"
#include "debug.h"
#include <vector>
#include <algorithm>
#include <unistd.h>

// Simple dynamical structure for internal use.

typedef struct {
    int jindex;
    real mass;
    vec pos;
    vec vel;
} dyn;

// Note the uncomfortable mix of arrays and vectors below.  The jdata
// functions all use arrays of ints/reals, while we use vectors
// internally for programming convenience.

static real partial_potential(vector<dyn> list1,
			      vector<int> list2, jdata& jd)
{
    // Return the potential energy of list1 (vector of dyns, separate
    // from the jdata arrays) relative to list2 (vector of jdata
    // indices).

    real pot = 0;
    for (unsigned int l1 = 0; l1 < list1.size(); l1++) {
	real pot1 = 0;
	int j1 = list1[l1].jindex;
	vec pos = list1[l1].pos;
	for (unsigned int l2 = 0; l2 < list2.size(); l2++) {
	    int j2 = list2[l2];
	    if (j2 != j1) {
		real r2 = jd.eps2;
		for (int k = 0; k < 3; k++)
		    r2 += pow(pos[k]-jd.pos[j2][k], 2);
		pot1 += jd.mass[j2]/sqrt(r2);
	    }
	}
	pot -= list1[l1].mass*pot1;
    }
    return pot;
}



static bool reflect_or_merge_orbit(real total_mass,
				   vec& rel_pos, vec& rel_vel,
				   real& energy, real& semi_major_axis,
				   real& eccentricity,
				   bool allow_mergers,
				   real rmin = _INFINITY_,
				   bool verbose = false)
{
    // Advance a two-body orbit past pericenter out to the same
    // separation.  We only need the unit vectors for the orbit in
    // order to perform the reflection, but these entail solving for
    // most of the orbital elements...

    // Code stolen from Starlab/kepler, with some fine points omitted.
    // *** Should use the standalone kepler code now. ***

    // Dynamics and geometry.

    real separation = abs(rel_pos);
    energy = 0.5 * rel_vel * rel_vel - total_mass / separation;

    vec normal_unit_vector = rel_pos ^ rel_vel;
    real angular_momentum = abs(normal_unit_vector);

    real periastron;

    if (energy != 0) {

        semi_major_axis = -0.5 * total_mass / energy;

        real rdotv = rel_pos * rel_vel;
        real temp = 1 - separation / semi_major_axis;
        eccentricity = sqrt(rdotv * rdotv / (total_mass * semi_major_axis)
                            + temp * temp);

	// Deal with rounding error, and avoid problems with nearly
	// circular or nearly linear orbits:

	if (eccentricity < 1.e-6) eccentricity = 0;

	if (energy > 0 && fabs(eccentricity-1) < 1.e-6) {

	    // Force a linear orbit.

	    eccentricity = 1;
	    angular_momentum = 0;
	}

        // Convention: semi_major_axis is always > 0.

        semi_major_axis = fabs(semi_major_axis);
	periastron = semi_major_axis * fabs(1 - eccentricity);

    } else {

        eccentricity = 1;
        semi_major_axis = _INFINITY_;
        periastron = 0.5 * angular_momentum * angular_momentum / total_mass;
    }

    if (verbose) {PRC(energy); PRC(semi_major_axis); PRL(eccentricity);}

    if (allow_mergers
	&& eccentricity < 1 && semi_major_axis <= rmin) return true;

    vec r_unit = rel_pos / separation;

    if (angular_momentum != 0) 
        normal_unit_vector /= angular_momentum;

    else {

        eccentricity = 1;
	periastron = 0;
        vec temp = vec(1,0,0);  	// construct an arbitrary normal vector
        if (fabs(r_unit[0]) > .5) temp = vec(0,1,0);
        normal_unit_vector = r_unit ^ temp;
        normal_unit_vector /= abs(normal_unit_vector);
    }

    vec t_unit = normal_unit_vector ^ r_unit;

    // Determine cosine and sine of the true anomaly.

    real cos_true_an = 1, sin_true_an = 0;

    if (eccentricity > 0) {

	if (eccentricity != 1) {

	    cos_true_an = ((periastron/separation) * (1 + eccentricity) - 1)
				/ eccentricity;
	    if (cos_true_an >  1 - 1.e-6) {
		cos_true_an = 1;
		sin_true_an = 0;
	    } else if (cos_true_an < -1 + 1.e-6) {
		cos_true_an = -1;
		sin_true_an = 0;
	    } else {
		sin_true_an = sqrt(1 - pow(cos_true_an,2));
		if (rel_pos * rel_vel < 0) sin_true_an = -sin_true_an;
	    }

	} else {

	    if (angular_momentum > 0) {

		// Special case: see McCuskey, p. 54.

		real true_anomaly = 2*acos(sqrt(periastron/separation));
		if (rel_pos * rel_vel < 0) true_anomaly = -true_anomaly;
		cos_true_an = cos(true_anomaly);
		sin_true_an = sin(true_anomaly);

	    } else {

		// Linear orbit.

		cos_true_an = -1;    // to get the unit vectors right below
		sin_true_an = 0;
	    }
	}
    }

    // if (verbose) {PRC(cos_true_an); PRL(sin_true_an);}

    vec longitudinal_unit_vector = cos_true_an * r_unit - sin_true_an * t_unit;
    vec transverse_unit_vector = sin_true_an * r_unit + cos_true_an * t_unit;

    // Reflecting the orbit simply entails reversing the transverse
    // component of rel_pos and the longitudinal component of dv.

    real dr_trans = rel_pos*transverse_unit_vector;
    real dv_long = rel_vel*longitudinal_unit_vector;

    // if (verbose) {PRC(abs(rel_pos)); PRL(abs(rel_vel));}
    rel_pos -= 2*dr_trans*transverse_unit_vector;
    rel_vel -= 2*dv_long*longitudinal_unit_vector;
    // if (verbose) {PRC(abs(rel_pos)); PRL(abs(rel_vel));}

    return false;
}



#define REVERSE 0

void jdata::two_body(int j1, int j2, vector<int> nbrlist)
{
    // Evolve the two-body orbit of j1 and j2 past periastron, and
    // correct all bookkeeping.

    const char *in_function = "jdata::two_body";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    int icomp1 = id[j1], icomp2 = id[j2];
    real mass1 = mass[j1], mass2 = mass[j2];
    real total_mass = mass1+mass2;
    real m2 = mass2/total_mass;
    real reduced_mass = mass1*m2;

    //-----------------------------------------------------------------
    // Calculate the center of mass and relative coordinates of
    // particles j1 and j2.

    vec cmpos, cmvel;
    for (int k = 0; k < 3; k++) {
	cmpos[k] = (mass[j1]*pos[j1][k]
		     + mass[j2]*pos[j2][k]) / total_mass;
	cmvel[k] = (mass[j1]*vel[j1][k]
		     + mass[j2]*vel[j2][k]) / total_mass;
    }

    vec dr = vec(pos[j1][0]-pos[j2][0],
		 pos[j1][1]-pos[j2][1],
		 pos[j1][2]-pos[j2][2]);
    vec dv = vec(vel[j1][0]-vel[j2][0],
		 vel[j1][1]-vel[j2][1],
		 vel[j1][2]-vel[j2][2]);

    // Note: by choice of sign, pos[j1] = cmpos + m2*dr,
    //                          pos[j2] = cmpos - (1-m2)*dr

    //-----------------------------------------------------------------
    // Calculate the potential energy of the (j1,j2) pair relative to
    // the neighbors.

    int jcm_init[2] = {j1, j2};
    int ncm_init = 2;

    // Express the cm data as a dyn array.

    vector<dyn> cm_init;
    for (int icm = 0; icm < ncm_init; icm++) {
	int j = jcm_init[icm];
	dyn cm;
	cm.jindex = j;
	cm.mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm.pos[k] = pos[j][k];
	    cm.vel[k] = vel[j][k];
	}
	cm_init.push_back(cm);
    }

    real pot_init = partial_potential(cm_init, nbrlist, *this);
    // PRL(pot_init);

    //-----------------------------------------------------------------
    // Advance the relative orbit past pericenter and out to the same
    // separation, or collapse the pair into a single particle.  The
    // factor of 2 in the merger condition is TBD: TODO.

    vec dr_old = dr, dv_old = dv;
    real energy, semi_major_axis, eccentricity;

    bool merge = reflect_or_merge_orbit(total_mass, dr, dv, energy,
					semi_major_axis, eccentricity,
					manage_encounters > 1,
					2*rmin, mpi_rank == 0);
    // PRC(merge); PRL(nbrlist.size());
    // PRC(dr); PRL(abs(dr));
    // PRC(dv); PRL(abs(dv));

    if (merge) {

	// Suppress merger if NN is too close.  Factor TBD.

	int jnn = nbrlist[0];		// list is sorted
	real dr2 = 0;
	for (int k = 0; k < 3; k++)
	    dr2 += pow(pos[jnn][k]-cmpos[k], 2);
	if (dr2 < rmin*rmin) {
	    if (mpi_rank == 0)
		cout << "suppressing merger because " << jnn
		     << " (ID =" << id[jnn] << ") is too close"
		     << endl << flush;
	    merge = false;
	}
    }

    //-----------------------------------------------------------------
    // Update the jdata arrays (mass, pos, and vel) to their final
    // values, and list the particles left after the encounter.

    int *jcm_final;
    int ncm_final;

    if (merge) {

	// Define CM quantities.

	real newstep = fmin(timestep[j1], timestep[j2]);
	real newrad = radius[j1]+radius[j2];
	int newid = binary_base + binary_count++;

	// Remove both components from the jdata arrays, add the CM,
	// and correct nbrlist.  We need to keep track of the changes
	// in order to recompute the potential energy, update the GPU
	// with new j-data, and recompute the forces on the
	// components/CM and neighbors.  The scheduler is updated by
	// remove_particle().

	// Temporarily convert nbrlist indices into IDs, to be
	// converted back at the end, so we don't have to know how
	// particles are added or removed, only that IDs are
	// preserved.

	for (unsigned int i = 0; i < nbrlist.size(); i++)
	    nbrlist[i] = id[nbrlist[i]];

	if (mpi_rank == 0)
	    cout << "removing " << j1 << " (ID = " << id[j1] << ")"
		 << endl << flush;
	remove_particle(j1);	

	j2 = inverse_id[icomp2];

	if (mpi_rank == 0)
	    cout << "removing " << j2 << " (ID = " << id[j2] << ")"
		 << endl << flush;
	remove_particle(j2);

	add_particle(total_mass, newrad, cmpos, cmvel, newid, newstep);
	int jcm = inverse_id[newid];
	if (mpi_rank == 0)
	    cout << "added " << jcm << " (ID = " << newid << ")"
		 << endl << flush;

	// Convert nbrlist back to indices.

	for (unsigned int i = 0; i < nbrlist.size(); i++)
	    nbrlist[i] = inverse_id[nbrlist[i]];

	// Save the binary properties for later use.

	binary_list.push_back(binary(newid, icomp1, icomp2, mass1, mass2,
				     semi_major_axis, eccentricity));

	// Put the CM on the final list.

	ncm_final = 1;
	jcm_final = (int*)malloc(ncm_final*sizeof(int));
	jcm_final[0] = jcm;

    } else {

	for (int k = 0; k < 3; k++) {
#if REVERSE == 0

	    // Use the reflected orbit just computed.

	    pos[j1][k] = cmpos[k] + m2*dr[k];
	    pos[j2][k] = cmpos[k] - (1-m2)*dr[k];
	    vel[j1][k] = cmvel[k] + m2*dv[k];
	    vel[j2][k] = cmvel[k] - (1-m2)*dv[k];
#else

	    // *** EXPERIMENT: Simply reverse the velocities in the CM
	    // *** frame.  No energy error, and statistically OK, even if
	    // *** it is dynamically completely wrong.

	    vel[j1][k] = cmvel[k] - m2*dv_old[k];	  // note sign change
	    vel[j2][k] = cmvel[k] + (1-m2)*dv_old[k];
#endif

	    // Affected j-data locations are j1, j2 only.
	}

	// Put the components on the final list.

	ncm_final = 2;
	jcm_final = (int*)malloc(ncm_final*sizeof(int));
	jcm_final[0] = j1;
	jcm_final[1] = j2;
    }

    //-----------------------------------------------------------------
    // Calculate the new potential energy and the energy error.

    vector<dyn> cm_final;
    for (int icm = 0; icm < ncm_final; icm++) {
	int j = jcm_final[icm];
	dyn cm;
	cm.jindex = j;
	cm.mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm.pos[k] = pos[j][k];
	    cm.vel[k] = vel[j][k];
	}
	cm_final.push_back(cm);
    }

    real pot_final = partial_potential(cm_final, nbrlist, *this);
    // PRL(pot_final);

    real de = pot_final - pot_init;
    if (mpi_rank == 0) {PRC(pot_init); PRC(pot_final); PRL(de);}

    //-----------------------------------------------------------------
    // Absorb the energy error into the relative motion of the
    // components.  For mergers, simply log the error.

    if (merge) {

	// Report the merger and the error, and save the energy change
	// (NB dE currently includes both internal and tidal
	// components).

	de -= reduced_mass*energy;
	if (mpi_rank == 0)
	    cout << "merged "
		 << j1 << " (ID = " << icomp1 << ") and "
		 << j2 << " (ID = " << icomp2
		 << ") at time " << system_time
		 << "  dE = " << de
		 << endl << flush;
	update_merger_energy(-de);

    } else {

	// Correct the error.

	if (de != 0) {	// de should be zero in the REVERSE case

	    // Redistribution is rather ad hoc.  Simplest approach is
	    // to modify the relative velocities of the interacting
	    // particles.

	    real kin = 0.5*reduced_mass*dv*dv;
	    real vfac2 = 1 - de/kin;

	    if (vfac2 < 0.25) {

		// We'll need to be cleverer in this case.  Let's see how
		// often it occurs...

		if (mpi_rank == 0)
		    cout << "warning: can't correct component velocities."
			 << endl;

	    } else {
		real v_correction_fac = sqrt(vfac2);
		if (mpi_rank == 0) PRL(v_correction_fac);
		dv *= v_correction_fac;
		for (int k = 0; k < 3; k++) {
		    vel[j1][k] = cmvel[k] + m2*dv[k];
		    vel[j2][k] = cmvel[k] - (1-m2)*dv[k];
		}
	    }
	}
    }

    //-----------------------------------------------------------------
    // Send new data on all affected j-particles to GPU.

    if (use_gpu) {

	// First make sure all pending data are properly flushed.

#ifndef NOSYNC		// set NOSYNC to omit this call;
	sync_gpu();	// hardly a true fix for the GPU update
			// problem, since we don't understand why it
			// occurs, but this does seem to work...
#endif

	if (merge && mpi_size > 1) {

	    // We may have changed the distribution of particles
	    // across domains.  Simplest to reinitialize all worker
	    // processes.

	    initialize_gpu(true);

	} else {
	    
	    // No change to the distribution of particles among
	    // domains.  Update the GPU data only for the affected
	    // locations.  Function update_gpu() will ignore any
	    // references to j >= nj.

	    // In all cases, the original j1 and j2 locations have
	    // changed.

	    update_gpu(jcm_init, ncm_init);

	    // If a merger has occurred, the CM location has changed
	    // too.

	    if (merge) update_gpu(jcm_final, ncm_final);
	}
    }

    //-----------------------------------------------------------------
    // Recompute forces on pair/CM and neigbors.  The following code
    // is taken from jdata::synchronize_list() and idata::advance().
    // Retain current time steps and scheduling.  Note that in the
    // REVERSE case, the accelerations should not change.

    int nnbr = nbrlist.size();
    int ni = nnbr + ncm_final;
    int ilist[ni];
    for (int i = 0; i < nnbr; i++) ilist[i] = nbrlist[i];
    for (int icm = 0; icm < ncm_final; icm++) ilist[nnbr+icm] = jcm_final[icm];

    if (!use_gpu) predict_all(system_time);
    idat->set_list(ilist, ni);
    idat->gather();
    idat->predict(system_time);
    idat->get_acc_and_jerk();		// compute iacc, ijerk
    idat->scatter();			// j acc, jerk <-- iacc, ijerk

    if (use_gpu) idat->update_gpu();

    // Could cautiously reduce neighbor steps here (and reschedule),
    // but that seems not to be necessary.

    delete [] jcm_final;
}



// *** Multiple code under development. ***

int jdata::find_binary(int id)
{
    for (unsigned int ib = 0; ib < binary_list.size(); ib++)
	if (binary_list[ib].binary_id == id) return (int)ib;
    return -1;
}

void jdata::expand_binary(hdyn *bb,
			  real *time_scale,	// default = NULL
			  real *length_scale2)	// default = NULL
{
    // Recursively expand binary components.  Note the implicit loop
    // over the entire binary list in find_binary -- to be improved.

    if (!bb) return;
    int ibin = find_binary(bb->get_index());
    if (ibin < 0) return;
    binary bin = binary_list[ibin];

    // Create a binary with bb as center of mass.

    real total_mass = bin.mass1 + bin.mass2;
    real m2 = bin.mass2/total_mass;

    kepler k;
    k.set_time(0);
    k.set_total_mass(total_mass);
    k.set_semi_major_axis(bin.semi);
    k.set_eccentricity(bin.ecc);
    set_random_orientation(k);
    k.set_mean_anomaly(randinter(0, 2*M_PI));
    k.initialize_from_shape_and_phase();

    if (time_scale) *time_scale = fmax(*time_scale, k.get_period());
    if (length_scale2) *length_scale2 = fmax(*length_scale2,
					     pow(k.get_semi_major_axis(),2));

    // Recursively create the component nodes.

    hdyn *b1 = new hdyn;
    hdyn *b2 = new hdyn;
    bb->set_oldest_daughter(b1);

    b1->set_parent(bb);
    b1->set_younger_sister(b2);
    b1->set_index(bin.comp1);
    b1->set_mass(bin.mass1);
    b1->set_pos(-m2*k.get_rel_pos());
    b1->set_vel(-m2*k.get_rel_vel());
    expand_binary(b1);

    b2->set_parent(bb);
    b2->set_older_sister(b1);
    b2->set_index(bin.comp2);
    b2->set_mass(bin.mass2);
    b2->set_pos((1-m2)*k.get_rel_pos());
    b2->set_vel((1-m2)*k.get_rel_vel());
    expand_binary(b2);
}

void jdata::remove_binary(int id)
{
    // Recursively remove a binary and substructure from the binary list.

    int ibin = find_binary(id);
    if (ibin >= 0) {
	binary bin = binary_list[ibin];
	remove_binary(bin.comp1);
	remove_binary(bin.comp2);
	binary_list.erase(binary_list.begin()+ibin);
    }
}

void jdata::add_binary(hdyn *b1)
{
    // Recursively add a binary and substructure to the binary list.

    if (!b1) return;
    hdyn *b2 = b1->get_younger_sister();
    if (!b2) return;

    add_binary(b1->get_oldest_daughter());
    add_binary(b2->get_oldest_daughter());

    int newid = binary_base + binary_count++;
    int icomp1 = b1->get_index();
    int icomp2 = b2->get_index();
    real mass1 = b1->get_mass();
    real mass2 = b2->get_mass();

    kepler k;
    k.set_time(0);
    k.set_total_mass(mass1+mass2);
    k.set_rel_pos(b2->get_pos()-b1->get_pos());
    k.set_rel_vel(b2->get_vel()-b1->get_vel());
    k.initialize_from_pos_and_vel();
    
    b1->get_parent()->set_index(newid);
    binary_list.push_back(binary(newid,
				 icomp1, icomp2, mass1, mass2,
				 k.get_semi_major_axis(),
				 k.get_eccentricity()));
}

void flatten(hdyn *b)
{
    // Flatten the hdyn tree under b.  Replace each center of mass by
    // the list of leaves under it.

    hdyn *bb = b->get_oldest_daughter();
    while (bb) {
	hdyn *bn = bb->get_younger_sister();
	hdyn *b1 = bb->get_oldest_daughter();
	if (b1) {

	    // Promote the daughter nodes.

	    for_all_daughters(hdyn, bb, d) {
		d->inc_pos(bb->get_pos());
		d->inc_vel(bb->get_vel());
	    }

	    // Replace bb by b1 in the tree.

	    if (b->get_oldest_daughter() == bb) b->set_oldest_daughter(b1);
	    b1->set_parent(b);
	    b1->set_older_sister(bb->get_older_sister());
	    hdyn *d = b1;
	    while (d->get_younger_sister()) d = d->get_younger_sister();
	    d->set_younger_sister(bn);

	    // Delete bb (setting all pointers NULL is overkill).

	    bb->set_parent(NULL);
	    bb->set_oldest_daughter(NULL);
	    bb->set_older_sister(NULL);
	    bb->set_younger_sister(NULL);
	    delete bb;

	    bb = b1;

	} else
	    bb = bn;
    }
}

nbody_descriptor jdata::integrate_multiple(vector<int> jcomp,
					   real dt_fac, real r2_fac)
{
    // Create and integrate to completion the multiple system
    // consisting of the listed components.  Return a structure
    // describing the system.  Do not rescale the top-level nodes.

    // Simple structure representing the multiple system (returned by
    // this function).

    nbody_descriptor s = {-1, 0, 0, 0, NULL};

    if (jcomp.size() < 2) return s;
    const char *in_function = "jdata::integrate_multiple";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Determine local IDs, masses, center of mass, and relative
    // components.

    int ncomp = (int)jcomp.size();
    int icomp[ncomp], mcomp[ncomp], total_mass = 0;
    vec rel_pos[ncomp], rel_vel[ncomp], cmpos = 0, cmvel = 0;
    for (int i = 0; i < ncomp; i++) {
	int j = jcomp[i];
	real mj = mass[j];
	icomp[i] = id[j];
	mcomp[i] = mj;
	total_mass += mj;
	cmpos += mj*newvec(pos[j]);
	cmvel += mj*newvec(vel[j]);
    }

    cmpos /= total_mass;
    cmvel /= total_mass;

    for (int i = 0; i < ncomp; i++) {
	int j = jcomp[i];
	rel_pos[i] = newvec(pos[j]) - cmpos;
	rel_vel[i] = newvec(vel[j]) - cmvel;
    }

    // Base initial length and time scales on only the first two
    // elements on the list (the interacting pair -- any others are
    // neighbors).

    real m12 = mass[0] + mass[1];
    vec dr12 = rel_pos[1] - rel_pos[0];
    vec dv12 = rel_vel[1] - rel_vel[0];
    real length_scale2 = square(dr12);
    real time_scale = sqrt(fmax(length_scale2/square(dv12),
				pow(length_scale2, 1.5)/m12));

    // Build a flat hdyn tree from the listed components.
    
    hdyn *b = new hdyn;
    b->set_mass(total_mass);
    b->set_pos(0);
    b->set_vel(0);
    hdyn *bp = NULL;
    for (int i = 0; i < ncomp; i++) {
	hdyn *bb = new hdyn;
	bb->set_parent(b);
	if (!bp)
	    b->set_oldest_daughter(bb);
	else {
	    bb->set_older_sister(bp);
	    bp->set_younger_sister(bb);
	}
	bb->set_index(icomp[i]);	// index is the jdata index, note
	bb->set_mass(mcomp[i]);
	bb->set_pos(rel_pos[i]);
	bb->set_vel(rel_vel[i]);
	bp = bb;
    }

    // Recursively expand binary components.

    for_all_daughters(hdyn, b, bb)
	expand_binary(bb, &time_scale, &length_scale2);

    s.cmpos = cmpos;
    s.cmvel = cmvel;
    s.b = b;
    s.scale2 = length_scale2;

    // Flatten the tree and send it to smallN_evolve().

    real t_end = _INFINITY_;	// use radius or structure for termination
    real dt_log = _INFINITY_;	// suppress log output
    bool verbose = false;
    b->set_eta(0.03);
    b->set_gamma(1.e-6);
    b->set_allow_full_unperturbed(true);
    flatten(b);

    s.status = smallN_evolve(b, t_end, dt_log,
			     dt_fac*time_scale, r2_fac*length_scale2,
			     verbose);

    // The state of the input tree after smallN_evolve isn't quite
    // what we want.  Construct a tree reflecting the hierarchical
    // structure of the system and delete b.

    s.b = get_tree(b);
    rmtree(b);
    return s;
}

bool jdata::check_add_neighbors(nbody_descriptor &s,
				vector<int>mult, vector<int> nbrlist)
{
    real scale2 = 9*s.scale2;	// conservative, but not foolproof
    bool changed = false;

    for (unsigned int i = 0; i < nbrlist.size(); i++) {
	int j = nbrlist[i];
	vec rrel, vrel;
	for (int k = 0; k < 3; k++) {
	    rrel[k] = pos[j][k] - s.cmpos[k];
	    vrel[k] = vel[j][k] - s.cmvel[k];
	}
	real r2 = square(rrel);
	real vr = rrel*vrel;
	if (vr < 0 && acos(sqrt(-vr/(r2*square(vrel)))) < atan(scale2/r2)) {

	    // Do some more work to refine the disance and time of
	    // closest approach.

	    if (0) {

		// TODO

	    } else {

		// Move neighbor i into the multiple.  Note that
		// removing this element means that i now refers to
		// the next entry, so we must reduce the value of i
		// here.

		mult.push_back(j);
		nbrlist.erase(find(nbrlist.begin(), nbrlist.end(), j));
		i--;
		changed = true;
	    }
	}
     }
    return changed;
}

real jdata::rescale(hdyn *b, real sep2)
{
    // Rescale the top-level nodes under b to squared maximum
    // separation sep2, preserving energy, angular momentum, and the
    // center of mass position and velocity (both 0).  Return a time
    // step to use in the jdata arrays.

    // The quick way to get this working *** FOR CODE INTEGRATION
    // PURPOSES ONLY *** is just to rescale the positions to maximum
    // separation sep and velocities to preserve the total energy.
    // *** THIS WILL GET THE ANGUAR MOMENTUM WRONG ***.  In the
    // two-body case it is trivial to use kepler to rescale to sep or
    // perastron, whichever is larger, to preserve angular momentum.
    // In the 3+-body case, more care is needed, and it simply may not
    // be possible to preserve everything.  Quite debatable how
    // important angular momentum conservation really is, since it is
    // small in any case. *** TODO ***

    return 0;
}

void jdata::multiple(int j1, int j2, vector<int> nbrlist_in)
{
    // Realize j1 and j2 as multiple particles and integrate the
    // internal motion to completion.  Check for inclusion of
    // neighbors, then rescale the final system, compute and correct
    // the tidal error, update all bookkeeping, and return.

    const char *in_function = "jdata::multiple";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Make a local copy of the neighbor list.

    vector<int> nbrlist;
    for (unsigned int i = 0; i < nbrlist.size(); i++)
	nbrlist.push_back(nbrlist_in[i]);

    // Integrate the multiple system in isolation.  Start by
    // establishing time and length scalings for the integration.

    real sep2 = square(newvec(pos[j2])-newvec(pos[j1]));
    real dt_fac = 5;			// ~arbitrary
    real r2_fac = _INFINITY_;		// 25

    vector<int> mult;
    mult.push_back(j1);
    mult.push_back(j2);
    nbody_descriptor s = integrate_multiple(mult, dt_fac, r2_fac);

    // See if any neighbors could have come close during the
    // interaction.  If so, move them into the multiple and repeat
    // the calculation.

    if (check_add_neighbors(s, mult, nbrlist)) {
	rmtree(s.b);
	s = integrate_multiple(mult, dt_fac, r2_fac);
    }

    // There are two possibilities on return: (1) the interaction is
    // really over, with all top-level components unbound and
    // receding; (2) some component is outside the limiting distance
    // from the center of mass.  In case (1) we are done, and can
    // rescale and exit.  In case (2) we have to be more cautious in
    // handling substructure in the multiple system... TODO...

    // *** For now, suppress possibility (2) by setting r2_fac large
    // above. ***

    // The hdyn pointer s.b returned by integrate_multiple() comes
    // directly from get_tree() in analyze.cc, which uses names
    // internally to refer to CM nodes and never changes indices.
    // Thus all CM nodes in s.b have index -1, while leaf indices
    // still refer directly to the jdata structures.  TODO: CHECK
    // THIS!

    //-----------------------------------------------------------------
    // Calculate the initial potential energy of the mult list
    // relative to the neighbors.

    vector<dyn> cm_init;
    for (unsigned int icm = 0; icm < mult.size(); icm++) {
	int j = mult[icm];
	dyn cm;
	cm.jindex = j;
	cm.mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm.pos[k] = pos[j][k];
	    cm.vel[k] = vel[j][k];
	}
	cm.pos += s.cmpos;
	cm.vel += s.cmvel;
	cm_init.push_back(cm);
    }

    real pot_init = partial_potential(cm_init, nbrlist, *this);
    // PRL(pot_init);

    //-----------------------------------------------------------------
    // Delete all top-level initial CMs from the jdata arrays.
    // Temporarily convert mult and nbrlist indices into IDs, to be
    // converted back later, so we don't have to know how particles
    // are added or removed, only that IDs are preserved.

    for (unsigned int i = 0; i < mult.size(); i++)
	mult[i] = id[mult[i]];
    for (unsigned int i = 0; i < nbrlist.size(); i++)
	nbrlist[i] = id[nbrlist[i]];

    // Recursively remove binaries in mult from jdata::binary_list.

    for (unsigned int i = 0; i < mult.size(); i++)
	remove_binary(mult[i]);

    // Remove *all* top-level particles in mult from the jdata arrays.

    for (unsigned int i = 0; i < mult.size(); i++) {
	if (mpi_rank == 0)
	    cout << "removing " << mult[i] << " (ID = " << id[mult[i]] << ")"
		 << endl << flush;
	remove_particle(inverse_id[mult[i]]);
    }

    // Convert nbrlist IDs back to indices.

    for (unsigned int i = 0; i < nbrlist.size(); i++)
	nbrlist[i] = inverse_id[nbrlist[i]];

    //-----------------------------------------------------------------
    // Define new CM particles and update jdata::binary_list.
    // Henceforth neglect low-level structure in s.b.

    for_all_daughters(hdyn, s.b, bb)
	add_binary(bb->get_oldest_daughter());

    //-----------------------------------------------------------------
    // Scale the new CMs (i.e. s.b top level) back to the initial sphere.

    real newstep = rescale(s.b, sep2);

    //-----------------------------------------------------------------
    // 
    // Put the top-level particles in a dyn vector and calculate the
    // tidal error.

    vector<dyn> cm_final;
    real kin = 0;
    for_all_daughters(hdyn, s.b, bb) {
	dyn cm;
	cm.jindex = bb->get_index();	// CM index will be -1
	cm.mass = bb->get_mass();
	cm.pos = bb->get_pos() + s.cmpos;
	vec vel = bb->get_vel();
	cm.vel = vel + s.cmvel;
	kin += cm.mass*square(vel);
	cm_final.push_back(cm);
    }
    kin /= 2;

    real pot_final = partial_potential(cm_final, nbrlist, *this);
    real de = pot_final - pot_init;
    if (mpi_rank == 0) {PRC(pot_init); PRC(pot_final); PRL(de);}

    rmtree(s.b);

    //-----------------------------------------------------------------
    // Correct the tidal error by rescaling the top-level velocities.

    real vfac = sqrt(1-de/kin);
    for (vector<dyn>::iterator i = cm_final.begin(); i != cm_final.end(); i++)
	i->vel = s.cmvel + vfac*(i->vel-s.cmvel);

    //-----------------------------------------------------------------
    // Add the new CMs to the jdata arrays.  The neighbor indices
    // shouldn't be affected by this (with the current "add" scheme),
    // but don't rely on this, so convert indices to IDs and back, as
    // above.

    for (unsigned int i = 0; i < nbrlist.size(); i++)
	nbrlist[i] = id[nbrlist[i]];

    // ***** NOTE that ncm_final could be 1.  CHECK - TODO *****

    int ncm_final = cm_final.size();
    int jcm_final[ncm_final];
    for (int icm = 0; icm < ncm_final; icm++) {
	dyn cm = cm_final[icm];
	int newid = cm.jindex;
	if (newid < 0)
	    newid = binary_base + binary_count++;
	else
	    newid = id[newid];
	real newrad = 0;	// true radius; close encounters use rmin - TODO
	add_particle(cm.mass, newrad, cm.pos, cm.vel, newid, newstep);
	jcm_final[icm] = inverse_id[newid];
	if (mpi_rank == 0)
	    cout << "added " << jcm_final[icm] << " (ID = " << newid << ")"
		 << endl << flush;
    }

    // Convert nbrlist IDs back to indices.

    for (unsigned int i = 0; i < nbrlist.size(); i++)
	nbrlist[i] = inverse_id[nbrlist[i]];

    //-----------------------------------------------------------------
    // Update the GPU(s).  Maybe too complicated to try to keep track
    // of which indices are affected, in general.  For now, just
    // reinitialize the GPU data in all worker processes (even if
    // there is only one worker, and the domains are unchanged).

    if (use_gpu) {

	// Make sure all pending data are properly flushed.

#ifndef NOSYNC		// set NOSYNC to omit this call;
	sync_gpu();	// not a true fix for the GPU update
			// problem, since we don't understand why it
			// occurs, but this does seem to work...
#endif

	initialize_gpu(true);
    }

    //-----------------------------------------------------------------
    // Recompute forces on the CM system and neigbors.  The following
    // code follows that in two_body().  Retain current time steps and
    // scheduling.

    // Make a list of all potentially affected particles -- the
    // neighbor list and the final CM particles -- and recompute their
    // accs and jerks.

    int nnbr = nbrlist.size();
    int ni = nnbr + ncm_final;
    int ilist[ni];
    for (int i = 0; i < nnbr; i++) ilist[i] = nbrlist[i];
    for (int icm = 0; icm < ncm_final; icm++) ilist[nnbr+icm] = jcm_final[icm];

    if (!use_gpu) predict_all(system_time);
    idat->set_list(ilist, ni);
    idat->gather();
    idat->predict(system_time);
    idat->get_acc_and_jerk();		// compute iacc, ijerk
    idat->scatter();			// j acc, jerk <-- iacc, ijerk

    if (use_gpu) idat->update_gpu();

    // Could cautiously reduce neighbor steps here (and reschedule),
    // but that seems not to be necessary.
}
