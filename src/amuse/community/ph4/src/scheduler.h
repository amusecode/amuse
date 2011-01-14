#ifndef SCHEDULER_H
#define SCHEDULER_H

// Define the scheduler class.

#include "stdinc.h"
#include "jdata.h"

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

    jdata *jd;
    list<jtime> blist;
    deque<li> bp;

    // jd is a pointer to the jdata structure served by the scheduler
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

    real t_next(const int j) const {return jd->time[j] + jd->timestep[j];}

    void initialize();
    void set_jd(jdata *j = NULL) {jd = j; initialize();}
    scheduler(jdata *j = NULL) {set_jd(j);}
    ~scheduler() {clear();}

    real get_list(int *ilist, int& n) const;
    int find_block(real t, int i1 = 0) const;
    void update();

    void add_particle(int j);
    void remove_particle(int j);

    void print_blist();
    void print_bp(bool verbose = false);
    void print(bool verbose = false);
    void check(bool repair = true, bool verbose = false);
};

#endif
