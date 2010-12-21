#ifndef IDATA_H
#define IDATA_H

// Define the idata class: data and methods operating on the
// current interaction list.

#include "stdinc.h"

class jdata;
class scheduler;

// i-data:

class idata {

  private:
    int ni;

  public:
    real ti;				// most recent update time

    int *iid, *ilist, *inn, *pnn;
    real *imass, *iradius, *itime, *itimestep, *ipot, *ppot, *idnn, *pdnn;
    real2 old_pos, old_vel, ipos, ivel,
	  old_acc, old_jerk, iacc, ijerk, pacc, pjerk;

    int get_ni() {return ni;}

    idata(jdata* jd = NULL) {
	iid = ilist = inn = pnn = NULL;
	imass = iradius = itime = itimestep
	    = ipot = ppot = idnn = pdnn = NULL;
	old_pos = old_vel = ipos = ivel
	    = old_acc = old_jerk = iacc = ijerk
	    = pacc = pjerk = NULL;
	if (jd) setup(*jd);
    }

    void cleanup();			// (in idata.cc)
    ~idata() {cleanup();}

    // In idata.cc:

    void setup(jdata& jd);
    void set_ni(int n = 0);
    void gather(jdata& jd);
    void scatter(jdata& jd);
    void set_list_all(jdata& jd);
    void set_list_sync(jdata& jd);
    real set_list(scheduler& sched);
    void set_list(int list[], int n);
    void get_partial_acc_and_jerk(jdata& jd);
    void get_acc_and_jerk(jdata& jd);
    real get_pot(jdata& jd);
    real get_kin(jdata& jd);
    void predict(real t, jdata& jd);
    void correct(real tnext, real eta);
    void advance(jdata& jd, real tnext, real eta);
    void check_encounters(jdata& jd);

    // In gpu.cc:

    void update_gpu(jdata& jd);
    void get_partial_acc_and_jerk_on_gpu(jdata& jd, bool pot = false);
};

#endif
