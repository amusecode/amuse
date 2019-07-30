
// *********************
// * i-data operations *
// *********************
//
// Global functions:
//
//	void idata::get_partial_acc_and_jerk()
//	void idata::get_acc_and_jerk()				* MPI *
//	void idata::predict(real t)
//	void idata::correct(real tnext, real eta)
//	real idata::get_pot()
//	void void idata::gather()
//	void idata::scatter()
//	void idata::check_encounters()
//	void idata::set_list_all()
//	void idata::set_list_sync()
//	real idata::set_list()
//	void set_list(int jlist[], int njlist)
//	void idata::set_list(int jlist[], int njlist)
//	void idata::advance(real tnext, real eta)
//
// OpenMP acceleration in get_partial_acc_and_jerk is for the purpose
// of illustration/experimentation.  Preferably use GPU acceleration
// for prediction and j-force calculation.

#include "jdata.h"
#include "idata.h"
#include "scheduler.h"

void idata::set_ni(int n)
{
    const char *in_function = "idata::set_ni";
    if (DEBUG > 2) PRL(in_function);

    if (n > 0) {
	cleanup();
	iid = new int[n];
	ilist = new int[n];
	inn = new int[n]();
	pnn = new int[n];
	lnn = NULL;		// force recomputation in get_partial_acc_and_jerk()
	imass = new real[n];
	iradius = new real[n];
	itime = new real[n];
	itimestep = new real[n];
	ipot = new real[n];
	ppot = new real[n];
	idnn = new real[n];
	pdnn = new real[n];
	old_pos = new real[n][3];
	old_vel = new real[n][3];
	ipos = new real[n][3];
	ivel = new real[n][3];
	old_acc = new real[n][3];
	old_jerk = new real[n][3];
	iacc = new real[n][3];
	ijerk = new real[n][3];
	pacc = new real[n][3];
	pjerk = new real[n][3];
    }
}

void idata::setup()
{
    const char *in_function = "idata::setup";
    if (DEBUG > 2 && jdat && jdat->mpi_rank == 0) PRL(in_function);
    if (jdat) {
	ti = jdat->system_time;
	set_ni(jdat->nj);
	set_list_all();
	gather();
	get_acc_and_jerk();			// set i pot, acc, jerk
	scatter();				// set j pot, acc, jerk
	jdat->idat = this;
    } else
	cout << "idata::setup: jdata pointer is NULL" << endl << flush;
}

void idata::cleanup()
{
    const char *in_function = "idata::cleanup";
    if (DEBUG > 2) PRL(in_function);

    if (iid) delete [] iid;
    if (ilist) delete [] ilist;
    if (inn) delete [] inn;
    if (pnn) delete [] pnn;
    if (imass) delete [] imass;
    if (iradius) delete [] iradius;
    if (itime) delete [] itime;
    if (itimestep) delete [] itimestep;
    if (ipot) delete [] ipot;
    if (ppot) delete [] ppot;
    if (idnn) delete [] idnn;
    if (pdnn) delete [] pdnn;
    if (old_pos) delete [] old_pos;
    if (old_vel) delete [] old_vel;
    if (ipos) delete [] ipos;
    if (ivel) delete [] ivel;
    if (old_acc) delete [] old_acc;
    if (old_jerk) delete [] old_jerk;
    if (iacc) delete [] iacc;
    if (ijerk) delete [] ijerk;
    if (pacc) delete [] pacc;
    if (pjerk) delete [] pjerk;
}

// The i-data temporarily mirror a subset of the j-data during a block
// time step.  Between idata::gather() and idata::scatter(), the
// following data are currently used:
//
//	j-data			i-data
//
//	id	     -->	iid
//	nn	    <--		inn
//	pot	    <-->	ipot
//	dnn	    <--		idnn
//	mass	     -->	imass
//	radius	     -->	iradius
//	time	    <-->	itime
//	timestep    <-->	itimestep
//	acc	     -->	old_acc
//		    <--		iacc
//	jerk	     -->	old_jerk
//		    <--		ijerk
//	pos	     -->	old_pos
//		    <--		ipos
//	vel	     -->	old_vel
//		    <--		ivel
//  no GPU:
//
//	pred_pos     -->	ipos
//	pred_vel     -->	ivel
//
// The j-data are predicted prior to gather if GPU is not used.
// Otherwise, the j-prediction is done on the GPU, but the i-data are
// predicted on the CPU.  Note that pred_pos and pred_vel are needed
// ONLY if the force calculation is done on the front end.  If a GPU
// is used, they are never used.

void idata::get_partial_acc_and_jerk()
{
    const char *in_function = "idata::get_partial_acc_and_jerk";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Compute the partial pots, accs, and jerks on the i-particles.
    //
    // Without a GPU, this function calculate all quantities with
    // respect to the pred_pos and pred_vel quantities in the
    // j-domain.
    //
    // With a GPU, function get_partial_acc_and_jerk_on_gpu() uses
    // sapporo/GRAPE library calls to accomplish the same task.
    //
    // Note that all "p" lists are contiguous.  Sets pnn, pdnn, ppot,
    // pacc, and pjerk.

    // Avoid unnecessary reductions in the 1-process case.  These
    // pointers only need be recomputed after setup() is called, which
    // sets lnn = NULL.

    real *lpot, *ldnn;
    real2 lacc, ljerk;

    if (jdat->mpi_size == 1) {
        lnn = inn;
        ldnn = idnn;
        lpot = ipot;
        lacc = iacc;
        ljerk = ijerk;
    } else {
        lnn = pnn;
        ldnn = pdnn;
        lpot = ppot;
        lacc = pacc;
        ljerk = pjerk;
    }
    
    // Define the j-domains.  These only need be recomputed if nj
    // changes.

    static int curr_nj = 0;
    static int j_start, j_end;

    if (jdat->nj != curr_nj) {
	jdat->define_domain(j_start, j_end);
	curr_nj = jdat->nj;
    }

    // Calculate the gravitational forces on the i-particles.

    real dx[3], dv[3], r2, xv;
    real r2i, ri, mri, mr3i, a3;
    real eps2 = jdat->eps2;

    for (int i = 0; i < ni; i++) {
	lpot[i] = 0;
	ldnn[i] = _INFINITY_;
	for (int k = 0; k < 3; k++) lacc[i][k] = ljerk[i][k] = 0;
	for (int j = j_start; j < j_end; j++) {
	    r2 = xv = 0;
	    for (int k = 0; k < 3; k++) {
		dx[k] = jdat->pred_pos[j][k] - ipos[i][k];
		dv[k] = jdat->pred_vel[j][k] - ivel[i][k];
		r2 += dx[k]*dx[k];
		xv += dx[k]*dv[k];
	    }
	    r2i = 1/(r2+eps2+_TINY_);
	    ri = sqrt(r2i);
	    mri = jdat->mass[j]*ri;
	    mr3i = mri*r2i;
	    a3 = -3*xv*r2i;
	    // PRC(jdat->mpi_rank); PRC(ri); PRL(mri);
	    if (r2 > _TINY_) {
		lpot[i] -= mri;
		if (r2 < ldnn[i]) {
		    ldnn[i] = r2;
		    lnn[i] = j;
		}
	    }
	    for (int k = 0; k < 3; k++) {
		lacc[i][k] += mr3i*dx[k];
		ljerk[i][k] += mr3i*(dv[k]+a3*dx[k]);
	    }
	}
	ldnn[i] = sqrt(ldnn[i]);
    }
}

void idata::get_acc_and_jerk()
{
    const char *in_function = "idata::get_acc_and_jerk";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Compute the accs and jerks on the i-particles with respect to
    // the predicted quantities in the j system.  Every node has a
    // complete copy of the j-data, but computes only its part of the
    // total.  Combine in process 0 to get the final result.  Sets
    // inn, idnn, ipot, iacc, and ijerk.

    if (!jdat->use_gpu) {

	// Calculate partial accs and jerks on the front end.

	get_partial_acc_and_jerk();

    } else {

	// Use the GPU to compute the forces.  Note that "pred"
	// quantities are never used.

	if (DEBUG > 1) cout << "getting acc and jerk on GPU" << endl << flush;
	get_partial_acc_and_jerk_on_gpu();
	if (DEBUG > 1) cout << "acc and jerk on GPU done" << endl << flush;

    }

#ifndef NOMPI
    //cout << "idata Barrier 1 for " << jdat->mpi_rank << endl << flush;
    //jdat->mpi_comm.Barrier();
    MPI_Barrier(jdat->mpi_comm);
    //cout << "idata Barrier 1a for " << jdat->mpi_rank << endl << flush;
#endif

    if (jdat->mpi_size > 1) {

	// Combine partial potentials, accs, and jerks to get the totals.
	// If size = 1, the data have already been saved in ipot, etc.

#ifndef NOMPI
	//jdat->mpi_comm.Allreduce(ppot, ipot, ni, MPI_DOUBLE, MPI_SUM);
	//jdat->mpi_comm.Allreduce(pacc, iacc, 3*ni, MPI_DOUBLE, MPI_SUM);
	//jdat->mpi_comm.Allreduce(pjerk, ijerk, 3*ni, MPI_DOUBLE, MPI_SUM);

	MPI_Allreduce(ppot, ipot, ni, MPI_DOUBLE, MPI_SUM, jdat->mpi_comm);
	MPI_Allreduce(pacc, iacc, 3*ni, MPI_DOUBLE, MPI_SUM, jdat->mpi_comm);
	MPI_Allreduce(pjerk, ijerk, 3*ni, MPI_DOUBLE, MPI_SUM, jdat->mpi_comm);
	// Update and distribute the nn pointers.  Start by determining
	// the global nearest neighbor distances.

	MPI_Allreduce(pdnn, idnn, ni, MPI_DOUBLE, MPI_MIN, jdat->mpi_comm);
#else
    for(int i = 0; i < ni; i++)
    {
        ipot[i] = ppot[i];
        for(int j = 0; j < 3; j++)
        {
            iacc[i][j] = pacc[i][j];
            ijerk[i][j] = pjerk[i][j];
        }
    }
#endif
	// Currently, pnn[i] is the nearest neighbor of i in the local j
	// domain.  Use the idnn array to decide if the local nearest
	// neighbor is the global one.  (We are assuming here that the
	// distances are correctly preserved by the above reduction, so
	// only one pnn will be nonzero...)

	for (int i = 0; i < ni; i++)
	    if (pdnn[i] > idnn[i]) pnn[i] = 0;

#ifndef NOMPI
	//jdat->mpi_comm.Allreduce(pnn, inn, ni, MPI_INT, MPI_SUM);
	MPI_Allreduce(pnn, inn, ni, MPI_INT, MPI_SUM, jdat->mpi_comm);
#else
    for(int i = 0; i < ni; i++)
    {
        inn[i] = pnn[i];
    }
#endif
    }

#ifndef NOMPI
    //cout << "idata Barrier 2 for " << jdat->mpi_rank << endl << flush;
    //jdat->mpi_comm->Barrier();
    MPI_Barrier(jdat->mpi_comm);
    //cout << "idata Barrier 2a for " << jdat->mpi_rank << endl << flush;
#endif
}

void idata::predict(real t)
{
    const char *in_function = "idata::predict";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Predict the i system to time t (old_pos, oldvel --> ipos,
    // ivel).  Note: should skip ipos, ivel in gather in this case.
    // We must do this for all i-particles in the GPU case (no
    // j-prediction on front end), and for all i-particles not in the
    // current j-domain in the non-GPU case.  In the latter case, ipos
    // is already jpred_pos (etc.) for particles in the current
    // j-domain.

    // Don't clutter the GPU code with extra checks used only in the
    // non-GPU case.  Just write the loop twice.

    if (jdat->use_gpu) {
	for (int i = 0; i < ni; i++) {
	    real dt = t - itime[i];
	    if (dt == 0) {
		for (int k = 0; k < 3; k++) {
		    ipos[i][k] = old_pos[i][k];
		    ivel[i][k] = old_vel[i][k];
		}
	    } else {
		for (int k = 0; k < 3; k++) {
		    ipos[i][k] = old_pos[i][k]
				    + dt*(old_vel[i][k]
					+ 0.5*dt*(old_acc[i][k]
					    + dt*old_jerk[i][k]/3));
		    ivel[i][k] = old_vel[i][k]
				    + dt*(old_acc[i][k]
					+ 0.5*dt*old_jerk[i][k]);
		}
	    }
	}
    } else {

	// Define the j-domains.  Logic would be easier if ilist was
	// sorted.

	int j_start, j_end;
	jdat->define_domain(j_start, j_end);

	for (int i = 0; i < ni; i++) {
	    int j = ilist[i];
	    if (j < j_start || j >= j_end) {
		real dt = t - itime[i];
		if (dt == 0) {
		    for (int k = 0; k < 3; k++) {
			ipos[i][k] = old_pos[i][k];
			ivel[i][k] = old_vel[i][k];
		    }
		} else {
		    for (int k = 0; k < 3; k++) {
			ipos[i][k] = old_pos[i][k]
					+ dt*(old_vel[i][k]
					      + 0.5*dt*(old_acc[i][k]
							+ dt*old_jerk[i][k]/3));
			ivel[i][k] = old_vel[i][k]
					+ dt*(old_acc[i][k]
					      + 0.5*dt*old_jerk[i][k]);
		    }
		}
	    }
	}
    }
}

void idata::correct(real tnext)
{
    const char *in_function = "idata::correct";
    if (DEBUG > 10) PRL(in_function);

    // Compute and apply the ilist corrector.  Return a new time step.
    //
    //		old_pos and old_vel are base quantities
    //		ipos and ivel are predicted quantities
    //		old_acc and old_jerk are at the base time
    //		iacc and ijerk are at the predicted time
    //
    // Note that tnext-time may not equal timestep if synchronization
    // has been forced for some reason (e.g. a collision).  However,
    // in ALL cases, the new time step will be commensurate with tnext
    // on exit.  (See the warnings in starlab/hdyn_ec.C, but note that
    // the problems should be less severe here with effective
    // softening).

#if 0

    // This elegant form of the corrector comes from Makino & Hut
    // (ACS)...

    for (int i = 0; i < ni; i++) {
	real dt = tnext - itime[i];
	real dt2 = dt*dt;
	for (int k = 0; k < 3; k++) {
	    ivel[i][k] = old_vel[i][k] + (old_acc[i][k] + iacc[i][k])*dt/2
				+ (old_jerk[i][k] - ijerk[i][k])*dt2/12;
	    ipos[i][k] = old_pos[i][k] + (old_vel[i][k] + ivel[i][k])*dt/2
				+ (old_acc[i][k] - iacc[i][k])*dt2/12;
	}
	timestep[i] = ???;
    }

#else

    // ...but this ugly version is actually more useful for purposes
    //    of computing the Aarseth time step.

    real eta = jdat->eta;

    for (int i = 0; i < ni; i++) {

	real dt = tnext - itime[i];
	real dt2 = dt*dt;
	real a2 = 0, j2 = 0, k2 = 0, l2 = 0;

	for (int k = 0; k < 3; k++) {
	    real alpha = -3*(old_acc[i][k] - iacc[i][k])
			 - dt * (2*old_jerk[i][k] + ijerk[i][k]);  // a'' dt^2/2
	    real beta =  2*(old_acc[i][k] - iacc[i][k])
			 + dt * (old_jerk[i][k] + ijerk[i][k]);	   // a'''dt^3/6
	    ipos[i][k] += (alpha/12 + beta/20)*dt2;
	    ivel[i][k] += (alpha/3 + beta/4)*dt;
	    a2 += pow(iacc[i][k], 2);				// (a)^2
	    j2 += pow(dt*ijerk[i][k], 2);			// (a')^2 dt^2
	    k2 += pow(2*alpha, 2);				// (a'')^2 dt^4
	    l2 += pow(6*beta, 2);				// (a''')^2 dt^6
	}

	// Aarseth step:

	real newstep = eta * dt * sqrt((sqrt(a2 * k2) + j2)
				       / (sqrt(j2 * l2) + k2));
	// Force the time step down to a power of 2 commensurate with
	// tnext.  Start by finding the first power of 2 below the old
	// timestep.  This is redundant if the old step is already a
	// power of 2, but will become necessary if non-power-of-two
	// time steps are allowed, or forced.

	int exponent;
	real oldstep2 = itimestep[i]/(2*frexp(itimestep[i], &exponent));

#if 0
	cout << "newstep:" << newstep << " , oldstep2:"<< oldstep2 <<endl;
	cout << "eta, dt:"<< eta <<", "<< dt<<endl;
	cout << "eta, dt:"<< eta <<", "<< dt<<endl;
#endif

	// At this stage it is possible that oldstep2 is not
	// commensurate with the new time tnext (e.g. if this step was
	// forced by system synchronization, rather than "occurring
	// naturally").  Force it to be commensurate.  (This is
	// unlikely to drop the step too low for softened motion, but
	// see the starlab warning mentioned above.)

	while (fmod(tnext, oldstep2) != 0) oldstep2 /= 2;

	if (newstep < oldstep2) {

	    // Halving the timestep is always OK.  Note that we don't
	    // necessarily drop the step below newstep -- we assume
	    // that a single halving is sufficient.

	    newstep = oldstep2 / 2;

	} else {

	    // To preserve the synchronization, doubling is OK only if
	    // the current time could be reached by an integral number
	    // of doubled time steps (use of fmod below).

	    real t2 = 2 * oldstep2;
	    if (newstep >= t2 && fmod(tnext, t2) == 0)
		newstep = t2;
	    else
		newstep = oldstep2;
	}

	itime[i] = tnext;
	itimestep[i] = newstep;
    }
#endif
}

real idata::get_pot()
{
    const char *in_function = "idata::get_pot";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Return the total potential energy of the i system by summing
    // the potentials of the individual particles.  Note that the pots
    // refer to the last acc and jerk calculation and hence will omit
    // any corrections.  This is fast but not necessarily accurate.
    // For it to return useful results, we should just have predicted
    // and called get_acc_and_jerk() for the entire system.

    real pot2 = 0;
    for (int i = 0; i < ni; i++)
	pot2 += jdat->mass[ilist[i]]*ipot[i];
    return pot2/2;
}

real idata::get_kin()
{
    const char *in_function = "idata::get_kin";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Return the total kinetic energy of the i system by summing over
    // individual particles.  See note in get_pot().

    real kin2 = 0;
    for (int i = 0; i < ni; i++) {
	real v2 = 0;
	for (int k = 0; k < 3; k++) v2 += pow(ivel[i][k],2);
	kin2 += jdat->mass[ilist[i]]*v2;
    }
    return kin2/2;
}

void idata::gather()
{
    const char *in_function = "idata::gather";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    for (int i = 0; i < ni; i++) {
	int j = ilist[i];
	itime[i] = jdat->time[j];
	itimestep[i] = jdat->timestep[j];
	ipot[i] = jdat->pot[j];
	for (int k = 0; k < 3; k++) {
	    old_pos[i][k] = jdat->pos[j][k];
	    old_vel[i][k] = jdat->vel[j][k];

	    // Need id and radius for encounter checking.  Don't
	    // scatter back.

	    iid[i] = jdat->id[j];
	    iradius[i] = jdat->radius[j];

	    if (!jdat->use_gpu) {

		ipos[i][k] = jdat->pred_pos[j][k];  // NB i??? = jdat->pred_???
		ivel[i][k] = jdat->pred_vel[j][k];

	    } else {

		// Don't use pred_pos/vel in the GPU case.  Also, need
		// mass for GPU update.  Don't scatter back.

		imass[i] = jdat->mass[j];
		ipos[i][k] = jdat->pos[j][k];
		ivel[i][k] = jdat->vel[j][k];
	    }

	    old_acc[i][k] = jdat->acc[j][k];
	    old_jerk[i][k] = jdat->jerk[j][k];
	}
    }
}

void idata::scatter()
{
    const char *in_function = "idata::scatter";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Scatter the i-particles back to the j-arrays.

    for (int i = 0; i < ni; i++) {
	int j = ilist[i];
	jdat->nn[j] = inn[i];
	jdat->pot[j] = ipot[i];
	jdat->dnn[j] = idnn[i];
	jdat->time[j] = itime[i];
	jdat->timestep[j] = itimestep[i];
	for (int k = 0; k < 3; k++) {
//	    jdat->pred_pos[j][k] = 
	    jdat->pos[j][k] = ipos[i][k];
//	    jdat->pred_vel[j][k] = 
	    jdat->vel[j][k] = ivel[i][k];
	    jdat->acc[j][k] = iacc[i][k];
	    jdat->jerk[j][k] = ijerk[i][k];
	}
    }
}

void idata::check_encounters()
{
    const char *in_function = "idata::check_encounters";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Search the current i-list to find the IDs of the particles (if
    // any) involved in the a close encounter or collision.  This
    // function is routinely called from idata::advance() after every
    // step.

    jdat->close1 = jdat->close2 = -1;
    jdat->coll1 = jdat->coll2 = -1;
    real rmin = jdat->rmin;

    if (!jdat->use_gpu || NN > 0) {

	// All relevant particles have nearest neighbor pointers.
	// Search those to identify encounters and collisons.

	real rmax_close = 0, rmax_coll = 0;
	int imax_close = -1, imax_coll = -1;

	for (int i = 0; i < ni; i++) {

	    int jnn = inn[i];			// j index, note
	    if (jnn >= 0) {			// valid neighbor

		// Note no dtmin test.

		real r = rmin/idnn[i];

		if (r > rmax_close) {
		    rmax_close = r;
		    imax_close = i;
		}

		// Hmmm.  Still have to go into jdata for radius[jnn].
		// Not clear if this is expensive.  Should monitor,
		// and suppress this check if not needed.  TODO.

		r = (iradius[i]+jdat->radius[jnn])/idnn[i];

		if (r > rmax_coll) {
		    rmax_coll = r;
		    imax_coll = i;
		}
	    }
	}

	if (rmax_close >= 1) {			// close criterion

	    int close1 = iid[imax_close];
	    int close2 = jdat->id[inn[imax_close]];

	    // Pair (close1,close2) is closest and inside rmin.
	    // Refine the test here before managing the close
	    // encounter.

	    // Default rmin is <m>/<v^2>.

	    int jclose1 = jdat->inverse_id[close1];
	    int jclose2 = jdat->inverse_id[close2];
	    //real mtot = jdat->mass[jclose1] + jdat->mass[jclose2];
	    //real r = idnn[imax_close];
	    real v2 = 0;
	    for (int k = 0; k < 3; k++)
		v2 += pow(jdat->vel[jclose2][k]-jdat->vel[jclose1][k], 2);

	    // PRC(r); PRC(rmin); PRC(mtot/v2); PRL(r*v2/mtot);

	    jdat->close1 = close1;
	    jdat->close2 = close2;
	}

	if (rmax_coll >= 1) {			// collision criterion

// 	    cout << endl;
// 	    int jnn = -1;
// 	    real dnn = 1.e300;
// 	    for (int j = 0; j < jdat->nj; j++)
// 		if (jdat->id[j] != iid[imax_coll]) {
// 		    real dij2 = 0;
// 		    for (int k = 0; k < 3; k++)
// 			dij2 += pow(jdat->pos[j][k]-ipos[imax_coll][k], 2);
// 		    if (dij2 < dnn) {
// 			dnn = dij2;
// 			jnn = j;
// 		    }
// 		}
// 	    dnn = sqrt(dnn);
// 	    PRC(ipos[imax_coll][0]); PRC(ipos[imax_coll][1]); PRL(ipos[imax_coll][2]);
// 	    PRC(jnn); PRL(dnn);

	    jdat->coll1 = iid[imax_coll];
	    jdat->coll2 = jdat->id[inn[imax_coll]];

// 	    int jcoll1 = jdat->inverse_id[jdat->coll1];
// 	    int jcoll2 = jdat->inverse_id[jdat->coll2];
// 	    real jdnn = 0;
// 	    for (int k = 0; k < 3; k++)
// 		jdnn += pow(jdat->pos[jcoll2][k]-jdat->pos[jcoll1][k], 2);
// 	    jdnn = sqrt(jdnn);
// 	    PRC(imax_coll); PRL(inn[imax_coll]);
// 	    cout << "j indices: "; PRC(jcoll1); PRL(jcoll2);
// 	    cout << "j IDs:     "; PRC(jdat->coll1); PRL(jdat->coll2);
// 	    cout << "j IDs:     "; PRC(jdat->id[jcoll1]); PRL(jdat->id[jcoll2]);
// 	    cout << "j radii:   "; PRC(jdat->radius[jcoll1]); PRL(jdat->radius[jcoll2]);
// 	    PRC(jdat->pos[jcoll1][0]); PRC(jdat->pos[jcoll1][1]); PRL(jdat->pos[jcoll1][2]);
// 	    PRC(jdnn); PRL(idnn[imax_coll]);

// 	    real dnnmax = 1.e300;
// 	    int jnnmax = -1;
// 	    for (int j = 0; j < jdat->nj; j++) {
// 		dnn = 0;
// 		for (int k = 0; k < 3; k++)
// 		    dnn += pow(jdat->pos[j][k]-ipos[imax_coll][k], 2);
// 		if (dnn < dnnmax) {
// 		    dnnmax = dnn;
// 		    jnnmax = j;
// 		}
// 	    }
// 	    dnnmax = sqrt(dnnmax);
// 	    PRC(jnnmax); PRL(dnnmax);
// 	    PRC(jdat->pos[jnnmax][0]); PRC(jdat->pos[jnnmax][1]); PRL(jdat->pos[jnnmax][2]);

// 	    for (int i = 0; i < ni; i++) PRC(iid[i]);
// 	    PRL(jdat->id[jnnmax]);
//	    //for (int j = 0; j < jdat->nj; j++) PRC(jdat->id[j]); cout << endl;

// 	    for (int j = 0; j < jdat->nj; j++) {
// 		int jd = jdat->id[j];
// 		if (jdat->inverse_id[jd] != j) {
// 		    PRC(j); PRC(jd); PRL(jdat->inverse_id[jd]);
// 		}
// 	    }
	}
    }
}

void idata::set_list_all()
{
    const char *in_function = "idata::set_list_all";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Make a list of the entire system.

    ni = jdat->nj;
    for (int j = 0; j < ni; j++) ilist[j] = j;
    ti = jdat->system_time;
}

void idata::set_list_sync()
{
    const char *in_function = "idata::set_list_sync";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    // Make a list of particles not already at system_time.

    ni = 0;
    for (int j = 0; j < jdat->nj; j++)
	if (jdat->time[j] < jdat->system_time)
	    ilist[ni++] = j;
    ti = jdat->system_time;
}

real idata::set_list()
{
    const char *in_function = "idata::set_list";
    if (DEBUG > 10) PRL(in_function);

    // System scheduler.  Determine the set of particles with the
    // minimum time + timestep, and return that time.

    // Note that sched already contains a pointer to the jdata, so no
    // need to pass that information twice.

    real tnext = _INFINITY_;

#if 0

    // Do an O(N) search (don't do this!).

    for (int j = 0; j < jdat->nj; j++)
	if (jdat->time[j]+jdat->timestep[j] < tnext)
	    tnext = jdat->time[j]+jdat->timestep[j];
    ni = 0;
    for (int j = 0; j < jdat->nj; j++)
	if (jdat->time[j]+jdat->timestep[j] == tnext) ilist[ni++] = j;

#else

    // Use the scheduler.

    tnext = sched->get_list(ilist, ni);

#endif

    return tnext;
}

void idata::set_list(int jlist[], int njlist)
{
    const char *in_function = "idata::set_list";
    if (DEBUG > 10) PRL(in_function);

    ni = 0;
    for (int jj = 0; jj < njlist; jj++)
	if (jlist[jj] < jdat->nj) ilist[ni++] = jlist[jj];
}

void idata::advance(real tnext,
		    bool zero_step_mode)	// default = false
{
    const char *in_function = "idata::advance";
    if (DEBUG > 2 && jdat->mpi_rank == 0) PRL(in_function);

    //cout << endl; PRC(jdat->system_time); PRL(ni);
    //for (int i = 0; i < ni; i++) PRC(ilist[i]);
    //PRL(tnext);

    gather();				// j pos,vel --> old_ipos, old_ivel
					// (j pred_pos, pred_vel --> ipos, ivel)
					// j acc, jerk --> i old_acc, old_jerk
    ti = tnext;
    //cout << "dt = " << ti - jdat->system_time << endl << flush;

    // Predict i-particles only (if using GPU; the GPU will handle the
    // j-particles), or i-particles not in the current j-domain (no
    // GPU).

    predict(tnext);
    get_acc_and_jerk();			// compute iacc, ijerk
    if (!zero_step_mode)
	correct(tnext);			// --> new ipos, ivel
    scatter();				// j pos, vel <-- ipos, ivel
 					// j acc, jerk <-- iacc, ijerk

    if (jdat->use_gpu) {

	// Send updated information to the GPU.

	update_gpu();
    }

    // Check for close encounters and collisions.

    check_encounters();

    // jdat->coll1 >= 0 ==> take some action in the calling function, or
    // return to AMUSE.

    if (DEBUG > 1) {
	if (jdat->close1 >= 0) {
	    PRC(jdat->system_time); PRC(jdat->close1); PRL(jdat->close2);
	}
	if (jdat->coll1 >= 0) {
	    PRC(jdat->system_time); PRC(jdat->coll1); PRL(jdat->coll2);
	}
    }
}
