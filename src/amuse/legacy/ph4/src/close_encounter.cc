
// Special approximate treatment of close encounters.  Suitable for
// runs in which we need to resolve relaxation but don't care about
// binaries.
//
// The code in this file is currently not used by AMUSE, but it does
// contain 1 jdata member function.
//
// Global function:
//
//	void jdata::resolve_encounter(idata& id, scheduler& sched)

#include "jdata.h"
#include "scheduler.h"
#include "idata.h"
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
    for (int j = 0; j < jd.get_nj(); j++) {
	element.jindex = j;
	element.r_sq = 0;
	for (int k = 0; k < 3; k++)
	    element.r_sq += pow(jd.pos[j][k]-center[k],2);
	rlist.push_back(element);
    }
    sort(rlist.begin(), rlist.end());
}

static inline void swap(int list[], int j1, int j2)
{
    if (j1 != j2) {
	int l = list[j1];
	list[j1] = list[j2];
	list[j2] = l;
    }
}

static real partial_potential(int list[], int nj, jdata& jd)
{
    // Return the potential energy of the first two elements of the
    // list relative to the rest.

    int j1 = list[0], j2 = list[1];
    real pot1 = 0, pot2 = 0;
    for (int jl = 2; jl < nj; jl++) {
	int j = list[jl];
	real r1 = 0, r2 = 0;
	for (int k = 0; k < 3; k++) {
	    r1 += pow(jd.pos[j][k]-jd.pos[j1][k], 2);
	    r2 += pow(jd.pos[j][k]-jd.pos[j2][k], 2);
	}
	pot1 += jd.mass[j]/sqrt(r1);
	pot2 += jd.mass[j]/sqrt(r2);
    }
    return -jd.mass[j1]*pot1 - jd.mass[j2]*pot2;
}

static void reflect_orbit(real total_mass, vec& rel_pos, vec& rel_vel)
{
    // Advance a two-body orbit past pericenter out to the same
    // separation.  We only need the unit vectors for the orbit in
    // order to perform the reflection, but these entail solving for
    // most of the orbital elements...

    // Code stolen from Starlab/kepler, with some fine points omitted.

    // Dynamics and geometry.

    real separation = abs(rel_pos);
    real energy = 0.5 * rel_vel * rel_vel - total_mass / separation;

    vec normal_unit_vector = rel_pos ^ rel_vel;
    real angular_momentum = abs(normal_unit_vector);

    real semi_major_axis, eccentricity, periastron;

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

        semi_major_axis = abs(semi_major_axis);
	periastron = semi_major_axis * abs(1 - eccentricity);

    } else {

        eccentricity = 1;
        semi_major_axis = huge;
        periastron = 0.5 * angular_momentum * angular_momentum / total_mass;
    }

    PRC(semi_major_axis); PRL(eccentricity);

    vec r_unit = rel_pos / separation;

    if (angular_momentum != 0) 
        normal_unit_vector /= angular_momentum;

    else {

        eccentricity = 1;
	periastron = 0;
        vec temp = vec(1,0,0);  	// construct an arbitrary normal vector
        if (abs(r_unit[0]) > .5) temp = vec(0,1,0);
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

    // PRC(cos_true_an); PRL(sin_true_an);

    vec longitudinal_unit_vector = cos_true_an * r_unit - sin_true_an * t_unit;
    vec transverse_unit_vector = sin_true_an * r_unit + cos_true_an * t_unit;

    // Reflecting the orbit simply entails reversing the transverse
    // component of rel_pos and the longitudinal component of dv.

    real dr_trans = rel_pos*transverse_unit_vector;
    real dv_long = rel_vel*longitudinal_unit_vector;

    // PRC(abs(rel_pos)); PRL(abs(rel_vel));
    rel_pos -= 2*dr_trans*transverse_unit_vector;
    rel_vel -= 2*dv_long*longitudinal_unit_vector;
    // PRC(abs(rel_pos)); PRL(abs(rel_vel));
}

#define REVERSE 0

bool jdata::resolve_encounter(idata& id, scheduler& sched)
{
    const char *infunction = "jdata::resolve_encounter";
    if (DEBUG > 2) PRL(in_function);

    bool status = false;
    if (close1 < 0 || close2 < 0 || eps2 > 0) return status;

    // We will treat this encounter as an unperturbed two-body event
    // and absorb the tidal errors into the nearby motion if:
    //
    // (1) particles close1 and close2 are approaching, and
    // (2) the next nearest particle is more than twice as far away
    //     as the separation between close1 and close2.
    //
    // We will improve on these criteria (e.g. to handle a mass
    // spectrum) soon.  TODO.

    int j1 = inverse_id[close1];
    int j2 = inverse_id[close2];

    if (j1 < 0 || j2 < 0) return status;	// caution!

    predict(j1, system_time);
    predict(j2, system_time);

    vec dr = vec(pred_pos[j1][0]-pred_pos[j2][0],
		 pred_pos[j1][1]-pred_pos[j2][1],
		 pred_pos[j1][2]-pred_pos[j2][2]);
    vec dv = vec(pred_vel[j1][0]-pred_vel[j2][0],
		 pred_vel[j1][1]-pred_vel[j2][1],
		 pred_vel[j1][2]-pred_vel[j2][2]);

    if (dr*dv >= 0) return status;		// criterion (1)

    // We will probably need to list neighbors soon anyway, so just
    // predict all particles and make a list of indices sorted by
    // distance away.  O(N) front-end operations -- could be
    // parallelized and sped up using the GPU.  TODO.

    predict_all(system_time, true);

    real mtot = mass[j1] + mass[j2];
    vec cmpos;
    for (int k = 0; k < 3; k++)
	cmpos[k] = (mass[j1]*pred_pos[j1][k]
		     + mass[j2]*pred_pos[j2][k]) / mtot;

    sort_neighbors(*this, cmpos);		// sorted list is rlist

    real dr2 = dr*dr;
    if (rlist[2].r_sq < 4*dr2) return status;	// criterion (2)

    cout << endl << "managing two-body encounter of "
	 << j1 << " (" << close1 << ") and "
	 << j2 << " (" << close2 << ") at time " << system_time
	 << endl << flush;

    status = true;

    // Prepare for two-body motion by synchronizing the neighbors.
    // Make a list of all particles within 100 |dr|, and a sub-list of
    // particles that need to be synchronized.  Factor TBD.  TODO.

    int *nbrlist = new int[nj];
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

    synchronize_list(synclist, nsync, id, &sched);
    delete [] synclist;

    // Note that the lists are based on pred_pos, but now the
    // positions have been advanced and corrected.  Still use the
    // indices based on rlist (but don't trust the distances).

    vec cmvel;
    for (int k = 0; k < 3; k++) {
	cmpos[k] = (mass[j1]*pos[j1][k]
		     + mass[j2]*pos[j2][k]) / mtot;
	cmvel[k] = (mass[j1]*vel[j1][k]
		     + mass[j2]*vel[j2][k]) / mtot;
    }

    dr = vec(pos[j1][0]-pos[j2][0],
	     pos[j1][1]-pos[j2][1],
	     pos[j1][2]-pos[j2][2]);
    dv = vec(vel[j1][0]-vel[j2][0],
	     vel[j1][1]-vel[j2][1],
	     vel[j1][2]-vel[j2][2]);
    dr2 = dr*dr;

    real mu = mass[j2]/mtot;

    // Note: by choice of sign, pos[j1] = cmpos + mu*dr,
    // pos[j2] = cmpos - (1-mu)*dr

    // PRL(cmpos+mu*dr);
    // PRC(pos[j1][0]); PRC(pos[j1][1]); PRL(pos[j1][2]); 
    // PRL(cmvel+mu*dv);
    // PRC(vel[j1][0]); PRC(vel[j1][1]); PRL(vel[j1][2]); 

    // Make sure j1 and j2 are at the start of nbrlist (shouldn't be
    // necessary).

    int loc = 0;
    for (int jl = 0; jl < nnbr; jl++) {
	if (nbrlist[jl] == j1 || nbrlist[jl] == j2)
	    swap(nbrlist, loc++, jl);
	if (loc > 1) break;
    }
    if (loc < 2) cout << "nbrlist: huh?" << endl << flush;

    // Calculate the potential energy of the (j1,j2) pair relative to
    // the neighbors (uses pos).

    real pot_init = partial_potential(nbrlist, nnbr, *this);

    // Advance the relative orbit past pericenter out to the same
    // separation.

    vec dr_old = dr, dv_old = dv;
    reflect_orbit(mtot, dr, dv);

    // Update jd.pos and jd.vel.

    for (int k = 0; k < 3; k++) {
#if REVERSE == 0

	// Use the reflected orbit just computed.

 	pos[j1][k] = cmpos[k] + mu*dr[k];
 	pos[j2][k] = cmpos[k] - (1-mu)*dr[k];
 	vel[j1][k] = cmvel[k] + mu*dv[k];
 	vel[j2][k] = cmvel[k] - (1-mu)*dv[k];
#else

	// *** EXPERIMENT: Simply reverse the velocities in the CM
	// *** frame.  No energy error, and statistically OK, even if
	// *** it is dynamically completely wrong.

 	vel[j1][k] = cmvel[k] - mu*dv_old[k];	  // note sign change
 	vel[j2][k] = cmvel[k] + (1-mu)*dv_old[k];
#endif
    }

    if (use_gpu) update_gpu(nbrlist, 2);	// modify pos and vel

    // Calculate the new potential energy and the energy error.

    real pot_final = partial_potential(nbrlist, nnbr, *this);
    real de = pot_final - pot_init;
    PRC(nnbr); PRC(pot_init); PRC(pot_final); PRL(de);

    // Redistribute the energy error among the components and the
    // neighbors.

    if (de != 0) {	// de should be zero in the REVERSE case

	// Redistribution is rather ad hoc.  Simplest approach is
	// to modify the relative velocities of the interacting
	// particles.

	real kin = 0.5*mass[j1]*mass[j2]*dv*dv/mtot;
	real vfac2 = 1 - de/kin;
	if (vfac2 < 0.25)

	    // We'll need to be cleverer in this case.  Let's see how
	    // often it occurs...

	    cout << "warning: can't correct component velocities." << endl;

	else {
	    real v_correction_fac = sqrt(vfac2);
	    PRL(v_correction_fac);
	    dv *= v_correction_fac;
	    for (int k = 0; k < 3; k++) {
		vel[j1][k] = cmvel[k] + mu*dv[k];
		vel[j2][k] = cmvel[k] - (1-mu)*dv[k];
	    }
	}
    }

    // Recompute all forces.  The following code is taken from
    // jdata::synchronize_list() and idata::advance().  Retain current
    // time steps and scheduling.  Note that in the REVERSE case, the
    // accelerations should not change.

    if (!use_gpu) predict_all(system_time);
    id.set_list(nbrlist, nnbr);
    id.gather(*this);
    id.predict(system_time, *this);
    id.get_acc_and_jerk(*this);		// compute iacc, ijerk
    id.scatter(*this);			// j acc, jerk <-- iacc, ijerk
    if (use_gpu) id.update_gpu(*this);

    // Could cautiously reduce neighbor steps here (and reschedule),
    // but that seems not to be necessary.

    delete [] nbrlist;
    return status;
}
