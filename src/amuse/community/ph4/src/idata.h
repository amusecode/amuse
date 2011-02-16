#ifndef IDATA_H
#define IDATA_H

// Define the idata class: data and methods operating on the
// current interaction list.

#include "stdinc.h"

class jdata;
class scheduler;

// Note: after proper initialization:
//
//     jdata and scheduler have pointers to each other
//     jdata and idata have pointers to each other
//     idata has a pointer to scheduler
//
// Order of initialization: jdata, idata(jdata), scheduler(jdata).

class idata {

  public:
    int ni;				// public, but note set_ni() below
    real ti;				// most recent update time

    int *iid, *ilist, *inn, *pnn, *lnn;
    real *imass, *iradius, *itime, *itimestep, *ipot, *ppot, *idnn, *pdnn;
    real2 old_pos, old_vel, ipos, ivel,
	  old_acc, old_jerk, iacc, ijerk, pacc, pjerk;

    jdata *jdat;
    scheduler *sched;

    idata(jdata* jd = NULL) {
	iid = ilist = inn = pnn = lnn = NULL;
	imass = iradius = itime = itimestep
	    = ipot = ppot = idnn = pdnn = NULL;
	old_pos = old_vel = ipos = ivel
	    = old_acc = old_jerk = iacc = ijerk
	    = pacc = pjerk = NULL;
	jdat = jd;
	setup();
    }

    void cleanup();			// (in idata.cc)
    ~idata() {cleanup();}

    // In idata.cc:

    void setup();
    void set_ni(int n = 0);
    void gather();
    void scatter();
    void set_list_all();
    void set_list_sync();
    real set_list();
    void set_list(int jlist[], int njlist);
    void get_partial_acc_and_jerk();
    void get_acc_and_jerk();
    real get_pot();
    real get_kin();
    void predict(real t);
    void correct(real tnext);
    void advance(real tnext);
    void check_encounters();

    // In gpu.cc:

    void update_gpu();
    void get_partial_acc_and_jerk_on_gpu(bool pot = false);
};

#endif
