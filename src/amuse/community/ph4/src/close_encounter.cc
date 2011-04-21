
// Special approximate treatment of close encounters.  Suitable for
// runs in which we need to resolve relaxation but don't want/care
// about binaries.
//
// The code in this file is currently not used by AMUSE, but it does
// contain one jdata member function, only referenced in the
// standalone program.
//
// Global function:
//
//	void jdata::resolve_encounter()

#include "jdata.h"
#include "scheduler.h"
#include "idata.h"
#include "debug.h"
#include <vector>
#include <algorithm>
#include <unistd.h>

class rdata {
  public:
    int jindex;
    real r_sq;
};

inline bool operator < (const rdata& x, const rdata& y)
{
    return x.r_sq < y.r_sq;
}

static vector<rdata> rlist;

static void sort_neighbors(jdata& jd, vec center)
{
    // Direct single-process neighbor computation.  Speed up!  TODO.

    rlist.clear();
    rdata element;
    for (int j = 0; j < jd.nj; j++) {
	element.jindex = j;
	element.r_sq = 0;
	for (int k = 0; k < 3; k++)
	    element.r_sq += pow(jd.pos[j][k]-center[k],2);
	rlist.push_back(element);
    }
    sort(rlist.begin(), rlist.end());	// NB precise ordering of identical
					//    elements is unpredictable
}

static inline void swap(int list[], int j1, int j2)
{
    if (j1 != j2) {
	int l = list[j1];
	list[j1] = list[j2];
	list[j2] = l;
    }
}

static real partial_potential(int list1[], int n1,
			      int list2[], int n2,
			      jdata& jd)
{
    // Return the potential energy of list1 relative to list2.

    real pot = 0;
    for (int l1 = 0; l1 < n1; l1++) {
	real pot1 = 0;
	int j1 = list1[l1];
	for (int l2 = 0; l2 < n2; l2++) {
	    int j2 = list2[l2];
	    if (j2 != j1) {
		real r2 = jd.eps2;
		for (int k = 0; k < 3; k++)
		    r2 += pow(jd.pos[j1][k]-jd.pos[j2][k], 2);
		pot1 += jd.mass[j2]/sqrt(r2);
	    }
	}
	pot -= jd.mass[j1]*pot1;
    }
    return pot;
}



static bool reflect_or_merge_orbit(real total_mass,
				   vec& rel_pos, vec& rel_vel,
				   real& energy, real& semi_major_axis,
				   real& eccentricity,
				   real rmin = _INFINITY_,
				   bool verbose = false)
{
    // Advance a two-body orbit past pericenter out to the same
    // separation.  We only need the unit vectors for the orbit in
    // order to perform the reflection, but these entail solving for
    // most of the orbital elements...

    // Code stolen from Starlab/kepler, with some fine points omitted.
    // Should use the standalone kepler code now...

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

    if (eccentricity < 1 && semi_major_axis <= rmin) return true;

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

bool jdata::resolve_encounter()
{
    const char *in_function = "jdata::resolve_encounter";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    bool status = false;
    //PRRL(1);
    if (close1 < 0 || close2 < 0 || eps2 > 0) return status;

    //-----------------------------------------------------------------
    // We will treat this encounter as an unperturbed two-body event
    // and absorb the tidal errors into the nearby motion if:
    //
    // (1) particles close1 and close2 are approaching, and
    // (2) the next nearest particle is more than twice as far away
    //     as the separation between close1 and close2.
    //
    // We will improve on these criteria (e.g. to handle a mass
    // spectrum) soon.  TODO.

    int comp1 = close1, comp2 = close2;
    int j1 = inverse_id[comp1];
    int j2 = inverse_id[comp2];

    //PRRC(2); PRC(comp1); PRC(j1); PRC(comp2); PRL(j2);
    if (j1 < 0 || j2 < 0) return status;

    // Make j1 < j2, but note we may have to repeat this process with
    // the lists as constructed below.

    if (j1 > j2) {
	int temp = j1; j1 = j2; j2 = temp;
	temp = comp1; comp1 = comp2; comp2 = temp;
    }

    //PRRC(21); PRC(comp1); PRC(j1); PRC(comp2); PRL(j2);
    int pair[2] = {j1, j2};
    //PRRC(22); PRC(pos[j1][0]); PRC(pos[j1][1]); PRL(pos[j1][2]);
    //PRRC(22); PRC(pos[j2][0]); PRC(pos[j2][1]); PRL(pos[j2][2]);
    synchronize_list(pair, 2);
    //PRRC(23); PRC(pred_pos[j1][0]); PRC(pred_pos[j1][1]); PRL(pred_pos[j1][2]);
    //PRRC(23); PRC(pred_pos[j2][0]); PRC(pred_pos[j2][1]); PRL(pred_pos[j2][2]);

    predict(j1, system_time);
    predict(j2, system_time);

    vec dr = vec(pred_pos[j1][0]-pred_pos[j2][0],
		 pred_pos[j1][1]-pred_pos[j2][1],
		 pred_pos[j1][2]-pred_pos[j2][2]);
    vec dv = vec(pred_vel[j1][0]-pred_vel[j2][0],
		 pred_vel[j1][1]-pred_vel[j2][1],
		 pred_vel[j1][2]-pred_vel[j2][2]);

    //PRRC(3); PRL(dr*dv);
    //PRL(dr);
    //PRL(dv);
    if (dr*dv >= 0) return status;		// criterion (1)
    //PRRL(dr*dv);

    //-----------------------------------------------------------------
    // We will probably need to list neighbors soon anyway, so just
    // predict all particles and make a list of indices sorted by
    // distance away.  O(N) front-end operations -- could be
    // parallelized and sped up using the GPU.  TODO.

    predict_all(system_time, true);

    real mass1 = mass[j1], mass2 = mass[j2], total_mass = mass1 + mass2;
    real reduced_mass = mass1*mass2/total_mass;
    vec cmpos;
    for (int k = 0; k < 3; k++)
	cmpos[k] = (mass[j1]*pred_pos[j1][k]
		     + mass[j2]*pred_pos[j2][k]) / total_mass;

    //PRRC(31); PRL(cmpos);
    sort_neighbors(*this, cmpos);		// sorted list is rlist

    real dr2 = dr*dr;
    //PRRC(4); PRC(rlist[2].r_sq); PRL(9*dr2);
    if (rlist[2].r_sq < 9*dr2) return status;	// criterion (2): factor TBD

    if (mpi_rank == 0)
	cout << endl << "managing two-body encounter of "
	     << j1 << " (ID = " << comp1 << ") and "
	     << j2 << " (ID = " << comp2 
	     << ") at time " << system_time
	     << endl << flush;

    // Ordering of j1 and j2 elements in rlist is not clear.  Force
    // element 0 to be j1, since later lists depend on this ordering.

    if (rlist[0].jindex != j1) {
	rlist[0].jindex = j1;
	real tmp = rlist[0].r_sq;
	rlist[0].r_sq = rlist[1].r_sq;
	rlist[1].jindex = j2;
	rlist[1].r_sq = tmp;
    }

#if 0
    if (mpi_rank == 0) {
	cout << "neighbor distances (rmin = " << rmin << "):" << endl;
	int nl = rlist.size();
	if (nl > 5) nl = 5;
	for (int il = 0; il < nl; il++)
	    cout << "    " << rlist[il].jindex << " " << sqrt(rlist[il].r_sq)
		 << endl << flush;
    }
#endif

    status = true;

    //-----------------------------------------------------------------
    // Prepare for two-body motion by synchronizing the neighbors.
    // Make a list of all particles within 100 |dr|, and a sub-list of
    // particles that need to be synchronized.  Optimal factor TBD.
    // TODO.

    int *nbrlist0 = new int[nj+1];	// nbrlist0 leaves room for CM
    int *nbrlist = nbrlist0 + 1;	// nbrlist contains pair and neighbors
    int *synclist = new int[nj];
    int nnbr = 0, nsync = 0;

    for (int jl = 0; jl < nj; jl++) {
	int j = rlist[jl].jindex;
	if (rlist[jl].r_sq <= 1.e4*dr2) {
	    nbrlist[nnbr++] = j;
	    if (time[j] < system_time) synclist[nsync++] = j;
	} else
	    break;
    }

    synchronize_list(synclist, nsync);
    delete [] synclist;

    // Note that the lists are based on pred_pos, but now the
    // positions have been advanced and corrected.  Still use the
    // indices based on rlist (but don't trust the distances).

    //-----------------------------------------------------------------
    // Recalculate the center of mass and relative coordinates of
    // particles j1 and j2.

    vec cmvel;
    for (int k = 0; k < 3; k++) {
	cmpos[k] = (mass[j1]*pos[j1][k]
		     + mass[j2]*pos[j2][k]) / total_mass;
	cmvel[k] = (mass[j1]*vel[j1][k]
		     + mass[j2]*vel[j2][k]) / total_mass;
    }

    dr = vec(pos[j1][0]-pos[j2][0],
	     pos[j1][1]-pos[j2][1],
	     pos[j1][2]-pos[j2][2]);
    dv = vec(vel[j1][0]-vel[j2][0],
	     vel[j1][1]-vel[j2][1],
	     vel[j1][2]-vel[j2][2]);
    dr2 = dr*dr;
    //PRRL(dr2);

    // PRL(cmpos);
    // PRL(cmvel);
    // PRC(dr); PRL(abs(dr));
    // PRC(dv); PRL(abs(dv));

    real m2 = mass[j2]/total_mass;

    // Note: by choice of sign, pos[j1] = cmpos + m2*dr,
    //                          pos[j2] = cmpos - (1-m2)*dr

    //-----------------------------------------------------------------
    // Make sure j1 and j2 are at the start of nbrlist (shouldn't be
    // necessary).

    int loc = 0;
    for (int jl = 0; jl < nnbr; jl++) {
	if (nbrlist[jl] == j1 || nbrlist[jl] == j2)
	    swap(nbrlist, loc++, jl);
	if (loc > 1) break;
    }
    if (loc < 2 && mpi_rank == 0) cout << "nbrlist: huh?" << endl << flush;

    //-----------------------------------------------------------------
    // Calculate the potential energy of the (j1,j2) pair relative to
    // the neighbors (uses pos).

    real pot_init = 0;
    if (nnbr > 2)
	pot_init = partial_potential(nbrlist, 2, nbrlist+2, nnbr-2, *this);
    // PRL(pot_init);
    // real total_init = total_energy(nbrlist, nnbr, *this);

    //-----------------------------------------------------------------
    // Advance the relative orbit past pericenter and out to the same
    // separation, or collapse the pair into a single particle.  The
    // factor of 2 in the merger condition is TBD: TODO.

    //cout << "))) " << mpi_rank << " " << j1 << " ";
    //for (int k = 0; k < 3; k++) cout << " " << pos[j1][k];
    //cout << endl << flush;
    //cout << "))) " << mpi_rank << " " << j2 << " ";
    //for (int k = 0; k < 3; k++) cout << " " << pos[j2][k];
    //cout << endl << flush;

    vec dr_old = dr, dv_old = dv;
    real energy, semi_major_axis, eccentricity;
    bool merge = reflect_or_merge_orbit(total_mass, dr, dv, energy,
					semi_major_axis, eccentricity,
					2*rmin, mpi_rank == 0);
    // PRC(merge); PRL(nnbr);
    // PRC(dr); PRL(abs(dr));
    // PRC(dv); PRL(abs(dv));

    //PRRC(dr_old); PRL(dr);

    if (merge) {

	// Suppress merger if next NN is too close.  Factor TBD.

	real dr2 = 0;
	for (int k = 0; k < 3; k++)
	    dr2 += pow(pos[nbrlist[2]][k]-cmpos[k], 2);
	if (dr2 < rmin*rmin) {
	    if (mpi_rank == 0)
		cout << "suppressing merger because " << nbrlist[2]
		     << " (ID =" << id[nbrlist[2]] << ") is too close"
		     << endl << flush;
	    merge = false;
	}
    }

    //-----------------------------------------------------------------
    // Update jd.pos and jd.vel to their final values.

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
	// remove_particle().  Note that, as written, the code below
	// needs to know explicitly how removal and addition affect
	// the internal j-data -- probably not good.

	// Removal of particle j swaps j with the last particle and
	// reduces nj.

	if (mpi_rank == 0)
	    cout << "removing " << j1 << " (ID = " << id[j1] << ")"
		 << endl << flush;
	remove_particle(j1);	
	for (int jl = 1; jl < nnbr; jl++)	// recall 0, 1 are j1, j2
	    if (nbrlist[jl] == nj) nbrlist[jl] = j1;

	j2 = nbrlist[1];

	if (mpi_rank == 0)
	    cout << "removing " << j2 << " (ID = " << id[j2] << ")"
		 << endl << flush;
	remove_particle(j2);
	for (int jl = 2; jl < nnbr; jl++)
	    if (nbrlist[jl] == nj) nbrlist[jl] = j2;

	add_particle(total_mass, newrad, cmpos, cmvel, newid, newstep);
	if (mpi_rank == 0)
	    cout << "added " << nj-1 << " (ID = " << id[nj-1] << ")"
		 << endl << flush;

	// Strange new storage order preserves contiguous lists
	// below.  Affected j-data locations are j1, j2, nj-1.

	nbrlist[-1] = j1;
	nbrlist[0]  = j2;
	nbrlist[1]  = nj-1;

	// Save the binary properties for later use.

	binary_list.push_back(binary(newid, comp1, comp2, mass1, mass2,
				     semi_major_axis, eccentricity));

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
    }

    //-----------------------------------------------------------------
    // Calculate the new potential energy and the energy error.  Use
    // j1 and j2 (0 and 1) for a flyby, jcm (1) for a merger.

    real pot_final = 0;
    if (nnbr > 2)
	pot_final = partial_potential(nbrlist+merge, 2-merge,
				      nbrlist+2, nnbr-2, *this);
    real de = pot_final - pot_init;
    if (mpi_rank == 0) {PRC(pot_init); PRC(pot_final); PRL(de);}
    // real dtotal = total_energy(nbrlist+merge, nnbr-merge, *this)-total_init;
    // PRL(dtotal);

    //-----------------------------------------------------------------
    // Redistribute the energy error among the components and the
    // neighbors.  For mergers, simply report and live with the error,
    // for now.

    if (merge) {

	// Simply report the merger and the error (NB dE currently
	// includes both internal and tidal components).

	PRC(de);
	de -= reduced_mass*energy;
	PRL(de);
	update_merger_energy(-de);
	if (mpi_rank == 0)
	    cout << "merged "
		 << j1 << " (" << comp1 << ") and "
		 << j2 << " (" << comp2
		 << ") at time " << system_time
		 << "  dE = " << de
		 << endl << flush;

    } else {

	// Correct the error.

	if (de != 0) {	// de should be zero in the REVERSE case

	    // Redistribution is rather ad hoc.  Simplest approach is
	    // to modify the relative velocities of the interacting
	    // particles.

	    real kin = 0.5*mass[j1]*mass[j2]*dv*dv/total_mass;
	    real vfac2 = 1 - de/kin;
	    PRC(mpi_rank); PRL(vfac2);
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

	    update_gpu(nbrlist-merge, 2+merge);  // nbrlist-1 is nbrlist0
	}
    }

    //-----------------------------------------------------------------
    // Recompute forces on pair and neigbors.  The following code is
    // taken from jdata::synchronize_list() and idata::advance().
    // Retain current time steps and scheduling.  Note that in the
    // REVERSE case, the accelerations should not change.

    if (!use_gpu) predict_all(system_time);
    idat->set_list(nbrlist+merge, nnbr-merge);
    idat->gather();
    idat->predict(system_time);
    idat->get_acc_and_jerk();		// compute iacc, ijerk
    idat->scatter();			// j acc, jerk <-- iacc, ijerk

    if (use_gpu) idat->update_gpu();

    // Could cautiously reduce neighbor steps here (and reschedule),
    // but that seems not to be necessary.

    delete [] nbrlist0;
    return status;
}
