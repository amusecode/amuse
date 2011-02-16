
// *****************************************
// * Implementation of the scheduler class *
// *****************************************
//
// Operations:
//
//	set/change the jdata pointer
//	(re)initialize from a jdata structure
//	return a list of "next" particle indices in the jdata array
//	update after a step -- partially ordered lists
//	print the state of the scheduler
//	check the state of the scheduler
//
// Global functions:
//
//	void scheduler::initialize()
//	real scheduler::get_list(int *ilist, int& n) const
//	int scheduler::find_block(real t, int i1) const
//	void scheduler::update()
//	void scheduler::add_particle(int j)
//	void scheduler::remove_particle(int j)
//	void scheduler::print_blist()
//	void scheduler::print_bp(bool verbose)
//	void scheduler::print(bool verbose)
//	void check(bool repair, bool verbose)   (NOT IMPLEMENTED)

#include "idata.h"
#include "scheduler.h"
#include <unistd.h>

void scheduler::initialize(jdata *jd)	// (re)initialize based on the
					// time step data in jd
{
    const char *in_function = "scheduler::initialize";
    if (DEBUG > 2 && jd->mpi_rank == 0) PRL(in_function);

    if (jd) {
	jdat = jd;
	jdat->sched = this;
	if (jdat->idat)
	    jdat->idat->sched = this;
	else
	    cout << "scheduler::initialize: jdata has no idata pointer"
		 << endl << flush;
    }

    if (jdat) {
	clear();		// scheduler member function
	for (int j = 0; j < jdat->nj; j++)
	    blist.push_back(jtime(j, t_next(j)));

	blist.sort(compare);	// --> the main reason to use stl::list

	li ib = blist.begin(), jb = ib;
	bp.push_back(ib);
	while (++ib != blist.end()) {
	    if (t_next(ib->jindex) > t_next(jb->jindex)) bp.push_back(ib);
	    jb = ib;
	}
	bp.push_back(ib);	// list runs from blist.begin() to blist.end()
    }
}

real scheduler::get_list(int *ilist, int& n) const
{
    const char *in_function = "scheduler::get_list";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    li ib;
    for (n = 0, ib = bp[0]; ib != bp[1]; ib++) ilist[n++] = ib->jindex;
    return bp[0]->t_next;
}

int scheduler::find_block(real t,
			  int i1) const		// default = 0
{
    const char *in_function = "scheduler::find_block";
    if (DEBUG > 3 && jdat->mpi_rank == 0) PRL(in_function);

    // Return the starting index of the first block with t_next >= t.
    // Start the search at block i1.  Return i1 if t is less than any
    // t_next, and nbp-1 if t is greater than any t_next.
    //
    //            current/next
    //            active list          rest of the system
    //           t_next = time            t_next > time
    //          ________________________________________________
    // blist:  |_______________|________________________________|
    //
    //          ^               ^       ^      ...      ^        ^
    // bp:      0               1       2      ...    nbp-2    nbp-1
    //
    // Bisection below assumes that there are enough blocks to start,
    // so i1 <= i2, or bp.size() >= 3.  If the active particle list is
    // the entire j system, then bp.size() = 2 and we have to treat
    // this separately.

    int nbp = bp.size();		// should *always* be >= 2
    if (nbp == 2) return i1;

    real t1 = t_next(bp[i1]->jindex);
    if (t <= t1) return i1;

    int i2 = nbp - 2;
    real t2 = t_next(bp[i2]->jindex);
    if (t == t2) return i2;
    if (t > t2) return i2+1;	// force a new entry at the end of bp

    // Now t1 < t < t2 strictly.  Use bisection to locate the
    // desired block.

    while (i2 > i1+1) {		// only way for i2 = i1 + 1 without breaking
				// from the loop is if t1 < t < t2
	int im = (i1+i2)/2;
	real tm = t_next(bp[im]->jindex);
	if (t == tm)
	    return im;
	else if (t < tm)
	    i2 = im;
	else
	    i1 = im;
    }
    return i2;			// should check if t < t2 on return
}

void scheduler::update()
{
    const char *in_function = "scheduler::update";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // The portion of blist from bp[0] to bp[1]-1 is assumed to
    // have just been advanced, and is now out of order.  Merge it
    // back into the rest of the list.
    //
    //------------------------------------------------------------
    // The most efficient strategy depends on the length n of the
    // disordered part relative to the total length N of the list.
    //
    // 1. We could reinsert the elements one by one.  Cost ~ n log
    // log N, assuming len(bp) ~ log N.
    //
    // 2. If the disordered part is long, n ~ N, it might be
    // faster just to sort the entire list and reestablish the
    // pointers (reinitialize).  Cost ~ N log N + N.
    //
    // 3. Or it might be better to sort the disordered part, then
    // merge the range that overlaps the >bp[1] ordered region, as
    // in kira.  Cost ~ n log n + <n.
    //
    // Seems clear that (2) is never the best strategy (but beware
    // factors of order unity...).  Looks like (3) < (2), but (3)
    // --> (2) as n --> N, not surprisingly.  Also seems that (1)
    // is always faster...  Go with (1) for now.
    // ------------------------------------------------------------

    // For each element from bp[0] to bp[1]-1, update its t_next,
    // determine its new block, and insert it as the first element
    // of that block.

#if 0
    cout << "in update..." << endl << flush;
    for (li ib = bp[0]; ib != bp[1]; ib++) {
	cout << ib->jindex << "  ib->t_next = " << ib->t_next
	     << "  jdat->t_next = " << t_next(ib->jindex) << endl << flush;
    }
    print_blist();
    print_bp();
#endif

    for (li ib = bp[0]; ib != bp[1]; ib++) {
	ib->t_next = t_next(ib->jindex);
#if 0
	cout << ib->jindex << ":  set t_next = " << ib->t_next
	     << endl << flush;
#endif
	// Return the index in array bp of the first block with t_next
	// >= t_next(ib->jindex) after a time step.  Search only from
	// bp[1] to the end, since bp[0] to bp[1]-1 is the current
	// time step list, and will be removed.

	int ibp = find_block(ib->t_next, 1);
#if 0
	cout << "    found block = " << ibp << endl << flush;
#endif
	// Insert *ib before bp[ibp], which may run from start to end.

	li newbp = blist.insert(bp[ibp], *ib);

	// Update the pointers.

	if (bp[ibp] == blist.end() || ib->t_next < bp[ibp]->t_next) {
	    bp.insert(bp.begin()+ibp, newbp);
	} else
	    bp[ibp] = newbp;
    }

    // Clean up blist and bp.

    blist.erase(bp[0], bp[1]);
    bp.erase(bp.begin());
#if 0
    print_blist();
    print_bp();
#endif
}

void scheduler::add_particle(int j)
{
    const char *in_function = "scheduler::add_particle";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Add particle j to the scheduler.  Code is essentially the same
    // as in update().

    // Identify the block to which particle belongs.

    int ibp = find_block(t_next(j));

    // Insert the particle as the first member of that block.  Don't
    // check whether j is already in blist.

    li newbp = blist.insert(bp[ibp], jtime(j, t_next(j)));

    // Update the block pointers.

    if (bp[ibp] == blist.end() || t_next(j) < bp[ibp]->t_next) {
	bp.insert(bp.begin()+ibp, newbp);
    } else
	bp[ibp] = newbp;
}

bool scheduler::remove_particle(int j)
{
    const char *in_function = "scheduler::remove_particle";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Remove particle j from the scheduler.

    // Locate the particle in blist.  First, identify the block, then
    // search linearly for j.

    int ibp = find_block(t_next(j));
    
    for (li ib = bp[ibp]; ib != blist.end() && ib != bp[ibp+1]; ib++) {

	//if (jdat->system_time >= 2.19141) {
	//    PRC(j); PRC(ibp); PRL(ib->jindex);
	//}

	if (ib->jindex == j) {

	    // Correct pointers.

	    if (ib == bp[ibp]) {
		li jb = ib;
		bp[ibp] = ++jb;		// can't say ib+1
		if (bp[ibp] == bp[ibp+1]) bp.erase(bp.begin()+ibp);
	    }

	    // Remove ib.

	    blist.erase(ib);
	    return true;
	}
    }

    // Didn't find j.  Flag it.

    cout << "scheduler::remove_particle(): " << j << " not found, ";
    PRL(ibp);
    return false;
}

void scheduler::print_blist()
{
    const char *in_function = "scheduler::print_blist";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    if (jdat->mpi_rank == 0) {
	cout << "blist (size = " << blist.size() << "):" << endl << flush;
	for (li ib = blist.begin(); ib != blist.end(); ib++)
	    cout << "    " << ib->jindex << "  " << ib->t_next
		 << "  " << jdat->time[ib->jindex]
		 << "  " << jdat->timestep[ib->jindex]
		 << endl << flush;
    }
}

void scheduler::print_bp(bool verbose)		// default = false
{
    const char *in_function = "scheduler::print_bp";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    if (jdat->mpi_rank == 0) {

	int total = 0;
	for (int kb = 0; kb < (int)bp.size()-1; kb++)
	    for (li ib = bp[kb]; ib != bp[kb+1]; ib++) total++;
	cout << "block pointers (size = "		   // NB last element
	     << bp.size() << ", total = " << total << "):" //    is blist.end()
	     << endl << flush;	

	for (int kb = 0; kb < (int)bp.size()-1; kb++) {
	    cout << "    block " << kb << ", t_next = " << bp[kb]->t_next;
	    if (verbose)
		cout << endl << "     ";
	    else
		cout << "  ";
	    cout << flush;
	    int nb = 0;
	    for (li ib = bp[kb]; ib != bp[kb+1]; ib++) {
		nb++;
		if (verbose) cout << " " << ib->jindex;
	    }
	    if (verbose) cout << "  ";
	    cout << "[nb = " << nb << "]" << endl << flush;
	}
    }
}

void scheduler::print(bool verbose)		// default = false
{
    const char *in_function = "scheduler::print";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    if (jdat->mpi_rank == 0 && verbose) {
	cout << endl << "----------"
	     << endl << "system time = " << jdat->system_time << ", " << flush;
	print_blist();
    }
    print_bp(verbose);
    if (jdat->mpi_rank == 0 && verbose) cout << "----------"
					     << endl << flush;
}

void scheduler::check(bool repair,		// default = true
		      bool verbose)		// default = false
{
    const char *in_function = "scheduler::check";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Checks that should be carried out:
    //	 every particle is on the list
    //	 everything on the list is a particle
    //	 blist is ordered
    //	 blist jindex pointers are within the jdata range
    //	 blist elements are ordered in t_next
    //	 blist elements are consistent with jdata t_next
    //	 bp pointers are correct
    //
    // TODO
}

void scheduler::cleanup()
{
    clear();
}
