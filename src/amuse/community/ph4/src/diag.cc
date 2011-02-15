
// ***************
// * diagnostics *
// ***************
//
// Global functions:
//
//	real get_elapsed_time()
//	real get_cpu_time(real& user_time, real& system_time)
//	vec jdata::get_center()
//	void jdata::print_percentiles(vec center)

#include "jdata.h"
#include <sys/time.h>
#include <sys/resource.h>

#include <vector>
#include <algorithm>

// Note: the first call to a timer function sets its zero point.

static bool etset = false;
static timeval et0;
real get_elapsed_time()
{
    if (!etset) {
	gettimeofday(&et0, NULL);
	etset = true;
	cout << "initialized elapsed time counter" << endl << flush;
	return 0;
    } else {
	timeval et;			// et.tv_set = time(NULL), note
	gettimeofday(&et, NULL);
	return et.tv_sec - et0.tv_sec + 1.e-6*(et.tv_usec - et0.tv_usec);
    }
}

static bool ctset = false;
static real user0, sys0;
void get_cpu_time(real& user_time, real& system_time)
{
    struct rusage tim;
    getrusage(RUSAGE_SELF, &tim);

    if (!ctset) {
	user0 = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec;
	sys0 = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec;
	user_time = system_time = 0;
	ctset = true;
	cout << "initialized CPU time counter" << endl << flush;
    } else {
	user_time = tim.ru_utime.tv_sec + 1.e-6*tim.ru_utime.tv_usec - user0;
	system_time = tim.ru_stime.tv_sec + 1.e-6*tim.ru_stime.tv_usec - sys0;
    }
}

vec jdata::get_center()
{
    return vec(0,0,0);		// TODO
}

class mrdata {
  public:
    int index;
    real mass;
    real r_sq;
    mrdata(int ii, real mm, real rr) {index = ii; mass = mm; r_sq = rr;}
};

bool operator < (const mrdata& x, const mrdata& y) {return x.r_sq < y.r_sq;}

void jdata::get_lagrangian_radii(vector<real>& mlist,
				 vector<real>& rlist,
				 vec center)	// default = (0,0,0)
{
    // Compute lagrangian radii corresponding to the input lagrangian
    // masses.  Probably should parallelize this...

    vector<mrdata> mrlist;
    real total_mass = 0;
    for (int j = 0; j < nj; j++) {
	total_mass += mass[j];
	real r_sq = 0;
	for (int k = 0; k < 3; k++)
	    r_sq += pow(pos[j][k] - center[k], 2);
	mrlist.push_back(mrdata(id[j], mass[j], r_sq));
    }

    sort(mrlist.begin(), mrlist.end());

    // Extract the desired radii.

    vector<real>::const_iterator miter = mlist.begin();
    real mcurr = 0;
    rlist.clear();

    for (vector<mrdata>::const_iterator mriter = mrlist.begin();
	 mriter != mrlist.end(); mriter++) {
	mcurr += mriter->mass;

	if (mcurr >= *miter) {
	    rlist.push_back(sqrt(mriter->r_sq));
	    miter++;
	    if (miter == mlist.end()) break;
	}
    }
    mpi_comm.Barrier();
}

void jdata::print_percentiles(vec center)	// default = (0,0,0)
{
    // Print selected percentiles.

    vector<real> mlist, rlist;

    mlist.push_back(0.01);
    mlist.push_back(0.02);
    mlist.push_back(0.05);
    mlist.push_back(0.10);
    mlist.push_back(0.25);
    mlist.push_back(0.50);
    mlist.push_back(0.75);
    mlist.push_back(0.90);
    rlist.clear();

    get_lagrangian_radii(mlist, rlist, center);

    cout << "lagrangian radii (";
    for (unsigned int ip = 0; ip < mlist.size(); ip++) {
	if (ip > 0) cout << " ";
	printf("%.2f", mlist[ip]);
    }
    cout << "):" << endl << "   " << flush;
    for (unsigned int ip = 0; ip < mlist.size(); ip++)
	printf(" %.4f", rlist[ip]);
    cout << endl << flush;
}
