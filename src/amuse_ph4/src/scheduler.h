#ifndef SCHEDULER_H
#define SCHEDULER_H

// Define the scheduler class.

#include "stdinc.h"
#include "jdata.h"

// Note: after proper initialization:
//
//     jdata and scheduler have pointers to each other
//     jdata and idata have pointers to each other
//     idata has a pointer to scheduler
//
// Order of initialization: jdata, idata(jdata), scheduler(jdata).

#include <list>
#include <deque>
#include <algorithm>

struct jtime {		// really only want to compare t_next, but the list
    int  jindex;	// compare function argument apparently can't be a
    real t_next;	// member function of class scheduler...
    jtime(){}
    jtime(int j, real t): jindex(j), t_next(t) {}
    ~jtime(){}
};

inline bool compare(const jtime& j1, const jtime& j2) {
    return j1.t_next < j2.t_next;
}
typedef std::list<jtime>::iterator li;

class scheduler {
public:

    jdata *jdat;
    list<jtime> blist;
    deque<li> bp;

    // jdat is a pointer to the jdata structure served by the scheduler
    //
    // blist contains indices of particles in the jdata array, sorted
    // by t_next
    //
    // bp contains pointers (iterators) to the first element of each
    // blist block; bp[0] is blist.begin(), for uniformity of
    // expression
    //
    // Thus, bp[k] points to the list element at the start of block k,
    // *bp[k] is the corresponding jtime structure, and bp[k]->jindex
    // is the actual index into the jdata.  A drawback of this
    // approach, apparently forced on us by the sort() argument, is
    // that we have to keep jlist.t_next up to date -- done in
    // update().

    void clear() {
	blist.clear();
	bp.clear();
    }

    real t_next(const int j) const {return jdat->time[j] + jdat->timestep[j];}

    void initialize(jdata *jd = NULL);
    scheduler(jdata *jd = NULL) {initialize(jd);}
    ~scheduler() {cleanup();}

    // In scheduler.cc:

    real get_list(int *ilist, int& n) const;
    int find_block(real t, int i1 = 0) const;
    void update();

    void add_particle(int j);
    bool remove_particle(int j);

    void print_blist();
    void print_bp(bool verbose = false);
    void print(bool verbose = false);
    void check(bool repair = true, bool verbose = false);

    void cleanup();
};

#endif
