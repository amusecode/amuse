
// Manage true two-body encounters.
//
// Global function:
//
//	void jdata::two_body()

#include "jdata.h"
#include "scheduler.h"
#include "idata.h"
#include "debug.h"
#include <vector>
#include <algorithm>
#include <unistd.h>

typedef struct {
    int jindex;
    real mass;
    vec pos;
    vec vel;
} dyn;

static real partial_potential(dyn list1[], int n1,
			      int list2[], int n2,
			      jdata& jd)
{
    // Return the potential energy of list1 (array of dyns) relative
    // to list2 (array of jdata indices).

    real pot = 0;
    for (int l1 = 0; l1 < n1; l1++) {
	real pot1 = 0;
	int j1 = list1[l1].jindex;
	vec pos = list1[l1].pos;
	for (int l2 = 0; l2 < n2; l2++) {
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

void jdata::two_body(int j1, int j2, int nbrlist[], int nnbr)
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
    dyn cm_init[ncm_init];

    for (int icm = 0; icm < ncm_init; icm++) {
	int j = jcm_init[icm];
	cm_init[icm].jindex = j;
	cm_init[icm].mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm_init[icm].pos[k] = pos[j][k];
	    cm_init[icm].vel[k] = vel[j][k];
	}
    }

    real pot_init = 0;
    if (nnbr > 2)
	pot_init = partial_potential(cm_init, ncm_init, 
				     nbrlist, nnbr, *this);
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
    // PRC(merge); PRL(nnbr);
    // PRC(dr); PRL(abs(dr));
    // PRC(dv); PRL(abs(dv));

    if (merge) {

	// Suppress merger if NN is too close.  Factor TBD.

	int jnn = nbrlist[0];
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
	int newid = binary_base + binary_list.size();

	// Remove both components from the jdata arrays, add the CM,
	// and correct nbrlist.  We need to keep track of the changes
	// in order to recompute the potential energy, update the GPU
	// with new j-data, and recompute the forces on the
	// components/CM and neighbors.  The scheduler is updated by
	// remove_particle().

	// Temporarily convert nbrlist indices into IDs, to be
	// converted back at the end, so we don't have to know how
	// particles are added or removed, only that IDs will be
	// preserved.

	for (int i = 0; i < nnbr; i++)
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

	for (int i = 0; i < nnbr; i++)
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

    dyn cm_final[ncm_final];

    for (int icm = 0; icm < ncm_final; icm++) {
	int j = jcm_final[icm];
	cm_final[icm].jindex = j;
	cm_final[icm].mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm_final[icm].pos[k] = pos[j][k];
	    cm_final[icm].vel[k] = vel[j][k];
	}
    }
    real pot_final = 0;
    if (nnbr > 2)
	pot_final = partial_potential(cm_final, ncm_final,
				      nbrlist, nnbr, *this);
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

    int ni = nnbr+2-merge;
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



// Multiple code under development.  Template is two_body() above.

void jdata::multiple(int j1, int j2, int nbrlist[], int nnbr)
{
    // Realize j1 and j2 as multiple particles and integrate the
    // internal motion to completion.  Update bookkeeping, store all
    // internal information, scale the CMs back to the initial radius,
    // and correct the dynamics.

    const char *in_function = "jdata::multiple";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

#if 0

    int icomp1 = id[j1], icomp2 = id[j2];
    real mass1 = mass[j1], mass2 = mass[j2];
    real total_mass = mass1+mass2;
    real m2 = mass2/total_mass;

    //-----------------------------------------------------------------
    // Calculate the center of mass and relative coordinates of
    // particles j1 and j2.  OK.  Same as above.

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
    // the neighbors.  OK.  Same as above.

    int jcm_init[2] = {j1, j2};

    int ncm_init = 2;
    dyn cm_init[ncm_init];

    for (int icm = 0; icm < ncm_init; icm++) {
	int j = jcm_init[icm];
	cm_init[icm].jindex = j;
	cm_init[icm].mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm_init[icm].pos[k] = pos[j][k];
	    cm_init[icm].vel[k] = vel[j][k];
	}
    }

    real pot_init = 0;
    if (nnbr > 2)
	pot_init = partial_potential(cm_init, ncm_init, 
				     nbrlist, nnbr, *this);
    // PRL(pot_init);

    //-----------------------------------------------------------------
    // Realize and integrate the smallN system.  Build an hdyn tree
    // representing the multiple, flatten it, and send it to
    // smallN_evolve.  TODO.  Recursively expand binaries and remove
    // them from binary_list as we use them.  Also update
    // UpdatedParticles as needed.  Correct the jdata arrays
    // (add/remove/update) directly from the final tree at the end of
    // the calculation.

    //-----------------------------------------------------------------
    // Scale the surviving tree CMs back to the initial sphere.  TODO.

    //-----------------------------------------------------------------
    // Update jd.pos and jd.vel to their final values.  TODO.
    // Procedure: calculate the tidal error, correct the tree, then
    // repopulate the jdata arrays from the tree.

    int *jcm_final;
    int ncm_final;

    // Make a list of final CMs.  Add and remove jdata particles
    // accordingly from the top level of the tree.  Correct binary
    // bookkeeping.  TODO.


// *** OLD:

	// Define CM quantities.

	real newstep = fmin(timestep[j1], timestep[j2]);
	real newrad = radius[j1]+radius[j2];
	int newid = binary_base + binary_list.size();

	// Remove both components from the jdata arrays, add the CM,
	// and correct nbrlist.  We need to keep track of the changes
	// in order to recompute the potential energy, update the GPU
	// with new j-data, and recompute the forces on the
	// components/CM and neighbors.  The scheduler is updated by
	// remove_particle().

	// Temporarily convert nbrlist indices into IDs, to be
	// converted back at the end, so we don't have to know how
	// particles are added or removed, only that IDs will be
	// preserved.

	for (int i = 0; i < nnbr; i++)
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

	for (int i = 0; i < nnbr; i++)
	    nbrlist[i] = inverse_id[nbrlist[i]];

	// Save the binary properties for later use.

	binary_list.push_back(binary(newid, icomp1, icomp2, mass1, mass2,
				     semi_major_axis, eccentricity));

	// Put the CM on the final list.

	ncm_final = 1;
	jcm_final = (int*)malloc(ncm_final*sizeof(int));
	jcm_final[0] = jcm;


    //-----------------------------------------------------------------
    // Calculate the new potential energy and the energy error.  OK.
    // Same as above.

    dyn cm_final[ncm_final];

    for (int icm = 0; icm < ncm_final; icm++) {
	int j = jcm_final[icm];
	cm_final[icm].jindex = j;
	cm_final[icm].mass = mass[j];
	for (int k = 0; k < 3; k++) {
	    cm_final[icm].pos[k] = pos[j][k];
	    cm_final[icm].vel[k] = vel[j][k];
	}
    }
    real pot_final = 0;
    if (nnbr > 2)
	pot_final = partial_potential(cm_final, ncm_final,
				      nbrlist, nnbr, *this);
    // PRL(pot_final);

    real de = pot_final - pot_init;
    if (mpi_rank == 0) {PRC(pot_init); PRC(pot_final); PRL(de);}

    //-----------------------------------------------------------------
    // Redistribute the energy error among the final top-level CMs.
    // TODO.

    //-----------------------------------------------------------------
    // Send new data on all affected j-particles to the GPU(s).
    // Almost OK.  Mostly same as above.

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
	    update_gpu(jcm_final, ncm_final);	// check final != init TODO
	}
    }

    //-----------------------------------------------------------------
    // Recompute forces on pair/CM and neigbors.  The following code
    // is taken from jdata::synchronize_list() and idata::advance().
    // Retain current time steps and scheduling.  OK.  Same as above.

    int ni = nnbr+2-merge;
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

#endif
}
