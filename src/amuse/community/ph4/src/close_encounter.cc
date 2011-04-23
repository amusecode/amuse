
// Special treatment of close encounters.  Current options, selected
// by jdata.manage_encounters, are:
//
//	0	no special treatment
//	1	analytic pericenter reflection
//	2	(1) and merge close binaries
//	3	(2) and full integration of multiple (TODO)
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



#define REVERSE 0

bool jdata::resolve_encounter()
{
    const char *in_function = "jdata::resolve_encounter";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    bool status = false;
    if (!manage_encounters || eps2 > 0
	|| close1 < 0 || close2 < 0) return status;

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

    if (j1 < 0 || j2 < 0) return status;

    // Make j1 < j2, but note we may have to repeat this process with
    // the lists as constructed below.

    if (j1 > j2) {
	int temp = j1; j1 = j2; j2 = temp;
	temp = comp1; comp1 = comp2; comp2 = temp;
    }

    int pair[2] = {j1, j2};
    synchronize_list(pair, 2);

    predict(j1, system_time);
    predict(j2, system_time);

    vec dr = vec(pred_pos[j1][0]-pred_pos[j2][0],
		 pred_pos[j1][1]-pred_pos[j2][1],
		 pred_pos[j1][2]-pred_pos[j2][2]);
    vec dv = vec(pred_vel[j1][0]-pred_vel[j2][0],
		 pred_vel[j1][1]-pred_vel[j2][1],
		 pred_vel[j1][2]-pred_vel[j2][2]);

    if (dr*dv >= 0) return status;		// criterion (1)

    //-----------------------------------------------------------------
    // We will probably need to list neighbors soon anyway, so just
    // predict all particles and make a list of indices sorted by
    // distance away.  O(N) front-end operations -- could be
    // parallelized and sped up using the GPU.  TODO.

    predict_all(system_time, true);

    real mass1 = mass[j1], mass2 = mass[j2], total_mass = mass1 + mass2;
    vec cmpos;
    for (int k = 0; k < 3; k++)
	cmpos[k] = (mass[j1]*pred_pos[j1][k]
		     + mass[j2]*pred_pos[j2][k]) / total_mass;

    sort_neighbors(*this, cmpos);		// sorted list is rlist

    real dr2 = dr*dr;
    if (rlist[2].r_sq < 9*dr2) return status;	// criterion (2): factor TBD

    if (mpi_rank == 0)
	cout << endl << "managing two-body encounter of "
	     << j1 << " (ID = " << comp1 << ") and "
	     << j2 << " (ID = " << comp2 
	     << ") at time " << system_time
	     << endl << flush;

    // Ordering of j1 and j2 elements in rlist is not clear.  Force
    // element 0 to be j1, since later lists (may) depend on this
    // ordering.

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
    // Prepare for two-body/multiple motion by synchronizing the
    // neighbors.  Make a list of all particles within 100 |dr|, and a
    // sub-list of particles that need to be synchronized.  Optimal
    // factor TBD.  TODO.

    int *nbrlist = new int[nj];
    int *synclist = new int[nj];
    int nnbr = 0, nsync = 0;

    for (int jl = 0; jl < nj; jl++) {
	int j = rlist[jl].jindex;
	if (rlist[jl].r_sq <= 1.e4*dr2) {
	    if (j != j1 && j != j2) nbrlist[nnbr++] = j;
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
    // Hand off the rest of the calculation to the two_body()
    // function.  Later, we will add multiple() functionality, and
    // reinstate duplicated code at the end of this function.

    if (0 && (is_multiple(id[j1]) || is_multiple(id[j2])))
	multiple(j1, j2, nbrlist, nnbr);
    else
	two_body(j1, j2, nbrlist, nnbr);

    delete [] nbrlist;
    return status;
}
