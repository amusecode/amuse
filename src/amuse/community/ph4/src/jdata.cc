
// *********************
// * j-data operations *
// *********************
//
// Global functions:
//
//	void jdata::setup_mpi(MPI::Intracomm comm)
//	int  jdata::add_particle(real pmass, real pradius,
//				 vec ppos, vec pvel, int pid)
//	void jdata::remove_particle(int j)
//	void jdata::initialize_arrays()
//	int  jdata::get_inverse_id(int i)
//	void jdata::check_inverse_id()
//	void jdata::set_initial_timestep()
//	real jdata::get_pot(bool reeval)			* MPI *
//	real jdata::get_kin(bool reeval)
//	real jdata::get_energy(bool reeval)
//	void jdata::predict(int j, real t)
//	void jdata::predict_all(real t, bool full)
//	void jdata::advance()
//	void jdata::synchronize_all()
//	void jdata::synchronize_list(int jlist[], int njlist, idata& id)
//	void jdata::print(idata *id)
//	void jdata::cleanup()

#include "jdata.h"
#include "scheduler.h"
#include "idata.h"

#ifdef GPU
#include "grape.h"
#endif

void jdata::setup_mpi(MPI::Intracomm comm)
{
    const char *in_function = "jdata::setup_mpi";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    mpi_comm = comm;
    mpi_size = mpi_comm.Get_size();
    mpi_rank = mpi_comm.Get_rank();
}

static int init_id = 1, next_id = init_id;
int jdata::get_particle_id(int offset)		// default - 0
{
    // Return a unique ID for a new particle.  Just number
    // sequentially, for now.  Don't reuse IDs, and allow the
    // possibility of an additive offset.

    int pid = next_id++ + offset;
    
    if (user_specified_id.size() > 0)
	while (find(user_specified_id.begin(), user_specified_id.end(), pid)
	       != user_specified_id.end()) pid = next_id++ + offset;

    return pid;
}

template <class T>
static string ToString(const T &arg)
{
	ostringstream out;
	out << arg;
	return out.str();
}

int jdata::add_particle(real pmass, real pradius,
			vec ppos, vec pvel,
			int pid,		// default = -1
			real dt)		// default = -1
{
    const char *in_function = "jdata::add_particle";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    if (nj >= njbuf) {

	// Extend the work arrays.

	int dnjbuf = JBUF_INC;
	if (dnjbuf < njbuf/4) dnjbuf = njbuf/4;
	njbuf += dnjbuf;

	// Must preserve name, id, mass, radius, pos, vel already set.
	// All other arrays will be created at the end.

	string *name0 = name;
	int *id0 = id;
	real *mass0 = mass, *radius0 = radius;
	real *time0 = time, *timestep0 = timestep;
	real2 pos0 = pos, vel0 = vel;

	name = new string[njbuf];
	id = new int[njbuf];			// gather
	mass = new real[njbuf];			// gather
	radius = new real[njbuf];
	time = new real[njbuf];			// gather, scatter
	timestep = new real[njbuf];		// gather, scatter
	pos = new real[njbuf][3];		// gather, scatter
	vel = new real[njbuf][3];		// gather, scatter

	for (int j = 0; j < nj; j++) {
	    name[j] = name0[j];
	    id[j] = id0[j];
	    mass[j] = mass0[j];
	    radius[j] = radius0[j];
	    time[j] = time0[j];
	    timestep[j] = timestep0[j];
	    for (int k = 0; k < 3; k++) {
		pos[j][k] = pos0[j][k];
		vel[j][k] = vel0[j][k];
	    }
	}

	if (name0) delete [] name0;
	if (id0) delete [] id0;
	if (mass0) delete [] mass0;
	if (radius0) delete [] radius0;
	if (time0) delete [] time0;
	if (timestep0) delete [] timestep0;
	if (pos0) delete [] pos0;
	if (vel0) delete [] vel0;
    }

    // If a particle ID is specified, check that it isn't already in
    // use.

    if (pid >= 0) {
	if (pid >= init_id && pid < next_id) {
	    cout << "user specified ID " << pid << " already in use (1)."
		 << endl;
	    pid = -1;
	} else if (user_specified_id.size() > 0) {
	    if (find(user_specified_id.begin(), user_specified_id.end(), pid)
		!= user_specified_id.end()) {
		cout << "user specified ID " << pid << " already in use (2)."
		     << endl;
		pid = -1;
	    }
	}
    }

    if (pid < 0) pid = get_particle_id();

    // Place the new particle at the end of the list.

    id[nj] = pid;
    time[nj] = system_time;
    name[nj] = ToString(pid);
    mass[nj] = pmass;
    radius[nj] = pradius;
    for (int k = 0; k < 3; k++) {
	pos[nj][k] = ppos[k];
	vel[nj][k] = pvel[k];
    }

//     if (mpi_rank == 0) {
// 	cout << "add_particle: "; PRC(system_time);
// 	PRC(pmass); PRC(nj); PRL(pid);
//     }

    // Update the inverse ID structure, such that inverse[id[j]] = j.
    // Use an STL map, so we don't have to worry about the range or
    // sparseness of the IDs in use.

    inverse_id[pid] = nj;

    // Set the timestep if required, and update the scheduler, if it
    // exists.

    if (dt > 0) {
	timestep[nj] = dt;
	if (sched) {
	    // cout << "sched.add " << nj << " (" << pid << ")"
	    //	 << endl << flush;
	    sched->add_particle(nj);
	}
    }

    nj++;

    // Return the particle id.

    return pid;    
}

void jdata::remove_particle(int j)
{
    const char *in_function = "jdata::remove_particle";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Remove particle j from the system by swapping it with the last
    // particle and reducing nj.

    // Don't update the GPU(s), if any, yet, as repeated removals will
    // probably change the j-domains, so this must be done later,
    // e.g. in recommit_particles(), which also recomputes forces and
    // time steps.  Do update the scheduler, however: remove j and
    // remove/replace nj, if necessary.

    inverse_id.erase(id[j]);
    // cout << "sched.remove " << j << " (" << id[j] << ")"
    //	 << endl << flush;
    sched->remove_particle(j);

    nj--;
    if (j < nj) {
	// cout << "sched.remove " << nj << " (" << id[nj] << ")"
	//     << endl << flush;
	sched->remove_particle(nj);
	id[j] = id[nj];
	name[j] = name[nj];
	time[j] = time[nj];
	timestep[j] = timestep[nj];
	mass[j] = mass[nj];
	radius[j] = radius[nj];
	pot[j] = pot[nj];
	nn[j] = nn[nj];
	dnn[j] = dnn[nj];
	for (int k = 0; k < 3; k++) {
	    pos[j][k] = pos[nj][k];
	    vel[j][k] = vel[nj][k];
	    acc[j][k] = acc[nj][k];
	    jerk[j][k] = jerk[nj][k];
	    pred_pos[j][k] = pred_pos[nj][k];
	    pred_vel[j][k] = pred_vel[nj][k];
	}
	inverse_id[id[j]] = j;
	// cout << "sched.add " << j << endl << flush;
	sched->add_particle(j);
    } 
}

void jdata::initialize_arrays()
{
    const char *in_function = "jdata::initialize_arrays";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Complete the initialization of all arrays used in the calculation.

    nn = new int[njbuf];			// scatter
    pot = new real[njbuf];			// scatter
    dnn = new real[njbuf];			// scatter
    acc = new real[njbuf][3];			// gather, scatter
    jerk = new real[njbuf][3];			// gather, scatter
    pred_pos = new real[njbuf][3];		// gather, scatter
    pred_vel = new real[njbuf][3];		// gather, scatter

    for (int j = 0; j < nj; j++) {
	nn[j] = 0;
	pot[j] = dnn[j] = timestep[j] = 0;
	time[j] = system_time;
	for (int k = 0; k < 3; k++)
	    acc[j][k] = jerk[j][k] = pred_pos[j][k] = pred_vel[j][k] = 0;
    }

    // Check the integrity of the inverse ID structure (should have
    // inverse[id[j]] = j for all j).

    check_inverse_id();

    // *** To be refined... ***

    rmin = 2./nj;		// 90 degree turnaround in standard units
    dtmin = eta*(2./nj) * 4;	// for nn check: final 4 is a fudge factor

    if (mpi_rank == 0) {
	PRC(rmin); PRL(dtmin);
    }

    if (!use_gpu) {

	if (have_gpu && mpi_rank == 0)
	    cout << endl << "GPU use disabled" << endl << flush;

	predict_all(system_time, true);	 // set all predict_time, pred_xxx

    } else {

	// Copy all data into the GPU.

	if (mpi_rank == 0) cout << endl << "initializing GPU" << endl << flush;
	initialize_gpu();
	if (DEBUG > 1 && mpi_rank == 0) cout << "GPU initialization done"
					     << endl << flush;
    }
}

int jdata::get_inverse_id(int i)
{
    const char *in_function = "jdata::get_inverse_id";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    map<int,int>::iterator ii = inverse_id.find(i);
    if (ii == inverse_id.end()) return -1;
    else return ii->second;
}

void jdata::check_inverse_id()
{
    const char *in_function = "jdata::check_inverse_id";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    if (mpi_rank == 0) {

	cout << "checking inverse IDs... " << flush;
	int id_min = 10000000, id_max = -1;
	bool ok = true;

	for (int j = 0; j < nj; j++) {
	    int i = id[j];
	    if (i < id_min) id_min = i;
	    if (i > id_max) id_max = i;
	    int jj = get_inverse_id(i);
	    if (jj != j) {
		PRC(j); PRC(id[j]); PRL(get_inverse_id(id[j]));
		ok = false;
	    }
	}
	if (ok) cout << "OK  ";
	PRC(id_min); PRL(id_max);
    }
}

void jdata::set_initial_timestep()
{
    const char *in_function = "jdata::set_initial_timestep";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Assume acc and jerk have already been set.

    for (int j = 0; j < nj; j++) {
	real a2 = 0, j2 = 0;
	for (int k = 0; k < 3; k++) {
	    a2 += pow(acc[j][k], 2);
	    j2 += pow(jerk[j][k], 2);
	}

	real firststep = 0.0625 * eta * sqrt(a2/j2);	// conservative
    
    if(a2 == 0.0 || j2 == 0.0)
    {
        firststep = 0.0625 * eta;
    }
    if(eta == 0.0)
    {
        firststep = 0.0625;
    }
    
	// Force the time step to a power of 2 commensurate with
	// system_time.

	int exponent;
	firststep /= 2*frexp(firststep, &exponent);
	while (fmod(system_time, firststep) != 0) firststep /= 2;

	timestep[j] = firststep;
    }
}

real jdata::get_pot(bool reeval)		// default = false
{
    const char *in_function = "jdata::get_pot";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Compute the total potential energy of the entire j-particle
    // set.  If no idata structure is specified, then just do a
    // (parallel) calculation using the predicted j-data.  Otherwise,
    // use the i-data (and GPU) to predict and compute the potential.

    real total_pot = 0;

    if (reeval) {

	// Direct parallel calculation on the j-data, with no GPU.
	// Note: always predict first, so use pred_pos in this case.

	predict_all(system_time, true);	 // all nodes update all j-particles

	// Compute the partial potentials by j-domain, then combine to
	// get the total and distribute to all processes.

	// Define the j-domains.

	int n = nj/mpi_size;
	if (n*mpi_size < nj) n++;
	int j_start = mpi_rank*n;
	int j_end = j_start + n;
	if (mpi_rank == mpi_size-1) j_end = nj;

	real mypot = 0;

	for (int jd = j_start; jd < j_end; jd++) {
	    real dpot = 0;
	    for (int j = 0; j < nj; j++) {
		if (j == jd) continue;
		real r2 = 0;
		for (int k = 0; k < 3; k++)
		    r2 += pow(pred_pos[j][k] - pred_pos[jd][k], 2);
		if (r2 > tiny)
		    dpot -= mass[j]/sqrt(r2+eps2);
	    }
	    mypot += mass[jd]*dpot;
	}

	mpi_comm.Allreduce(&mypot, &total_pot, 1, MPI_DOUBLE, MPI_SUM);

	total_pot /= 2;				// double counting

    } else {

	// Use pot data.  This approach uses the GPU if available, and
	// is always preferable to the above.

	// Predict and recompute individual potentials.

	if (!use_gpu) predict_all(system_time);
	idat->set_list_all();			// select the entire j system
	idat->gather();				// copy into the i system
	idat->predict(system_time);
	idat->get_acc_and_jerk();		// parallel, may use GPU
						// computes acc, jerk, pot
	total_pot = idat->get_pot();
    }

    return total_pot;
}

real jdata::get_kin(bool reeval)		// default = false
{
    const char *in_function = "jdata::get_kin";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Return the total kinetic energy of the j system.

    // *** ASSUME that get_pot() has just been called, with the same
    // argument, so use pred_vel if id == NULL, and idata::get_kin()
    // otherwise.

    if (reeval) {
	real kin2 = 0;
	for (int j = 0; j < nj; j++) {
	    real v2 = 0;
	    for (int k = 0; k < 3; k++) v2 += pow(pred_vel[j][k],2);
	    kin2 += mass[j]*v2;
	}
	return kin2/2;
    } else
	return idat->get_kin();
}

real jdata::get_energy(bool reeval)		// default = false
{
    const char *in_function = "jdata::get_energy";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    real pot = get_pot(reeval);
    real kin = get_kin(reeval);			// order is important!
    return pot + kin;
}

void jdata::predict(int j, real t)
{
    const char *in_function = "jdata::predict";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Predict particle j to time t.  Do so even if dt = 0, since we
    // will want to use pred_ quantities.

    real dt = t - time[j];
    for (int k = 0; k < 3; k++) {
	pred_pos[j][k] = pos[j][k]
			     + dt*(vel[j][k]
				   + 0.5*dt*(acc[j][k]
					+ dt*jerk[j][k]/3));
	pred_vel[j][k] = vel[j][k]
			     + dt*(acc[j][k]
				   + 0.5*dt*jerk[j][k]);
    }
}

void jdata::predict_all(real t,
			bool full_range)	// default = false
{
    const char *in_function = "jdata::predict_all";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Predict the entire j-system to time t (pos, vel --> pred_pos,
    // pred_vel).  This function is only called in the non-GPU case.

    // Define the j-domains.  Predict only the j-domain of the current
    // process unless full_range = true (only on initialization).

    int j_start = 0, j_end = nj;

    if (!full_range) {
	int n = nj/mpi_size;
	if (n*mpi_size < nj) n++;
	j_start = mpi_rank*n;
	j_end = j_start + n;
	if (mpi_rank == mpi_size-1) j_end = nj;
    }

    if (predict_time == t) {
	for (int j = j_start; j < j_end; j++) {
	    for (int k = 0; k < 3; k++) {
		pred_pos[j][k] = pos[j][k];
		pred_vel[j][k] = vel[j][k];
	    }
	}
	return;
    }

    for (int j = j_start; j < j_end; j++) {
	real dt = t - time[j];
	for (int k = 0; k < 3; k++) {
	    pred_pos[j][k] = pos[j][k]
				+ dt*(vel[j][k]
				      + 0.5*dt*(acc[j][k]
						+ dt*jerk[j][k]/3));
	    pred_vel[j][k] = vel[j][k]
				+ dt*(acc[j][k]
				      + 0.5*dt*jerk[j][k]);
	}
    }

    predict_time = t;
}

void jdata::advance()
{
    const char *in_function = "jdata::advance";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Advance the system by a single step.
    // Procedure:
    //		determine the new i-list and next time tnext
    //		predict all particles to time tnext
    //		gather the i-list
    //		calculate accelerations and jerks on the i-list
    //		correct the i-list
    //		scatter the i-list

    real tnext = idat->set_list();	// determine next time, make ilist

    if (!use_gpu) {

	// GPU will do this, if present.  Note that this only predicts
	// the range of particles associated with the current process.

	predict_all(tnext);		// j pos, vel --> j pred_pos, pred_vel
    }

    // All of the actual work is done by idata::advance.

    idat->advance(tnext);
    block_steps += 1;
    total_steps += idat->ni;

    // Note that system_time remains unchanged until the END of the step.

    system_time = tnext;
    sched->update();
    // sched->print(true);
}

void jdata::synchronize_all()
{
    const char *in_function = "jdata::synchronize_all";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Synchronize all particles at the current system time.

    // Make an ilist of all particles not already at system_time.
    // Don't worry about disrupting the scheduler, as we will
    // recompute it after the synchronization is done.

    idat->set_list_sync();

    if (!use_gpu) {

	// GPU will do this, if present.  Note that this only predicts
	// the range of particles associated with the current process.

	predict_all(system_time);	// j pos, vel --> j pred_pos, pred_vel
    }

    // All of the actual work is done by idata::advance.

    idat->advance(system_time);
    block_steps += 1;
    total_steps += idat->ni;

    if (sched) sched->initialize();	// default is to reinitialize later
}

void jdata::synchronize_list(int jlist[], int njlist)
{
    const char *in_function = "jdata::synchronize_list";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Advance the particles on the list to the current system time.
    // They will be out of the normal schedule, so to avoid corrupting
    // the scheduler, remove the particles, advance them, then
    // reinsert them into the scheduler arrays.

    // Check that the particles all need to be advanced, and shrink
    // the list accordingly.

    int joffset = 0;
    for (int jj = 0; jj < njlist; jj++) {
	if (joffset > 0) jlist[jj-joffset] = jlist[jj];
	if (time[jlist[jj]] == system_time) joffset++;
    }
    njlist -= joffset;

    // Transmit the jlist to the idata routines.

    idat->set_list(jlist, njlist);

    if (!use_gpu) {

	// GPU will do this, if present.  Note that this only predicts
	// the range of particles associated with the current process.

	predict_all(system_time);	// j pos, vel --> j pred_pos, pred_vel
    }

    // Remove the particles from the scheduler.

    bool ok = true;
    for (int jj = 0; jj < njlist; jj++)
	ok &= sched->remove_particle(jlist[jj]);
    if (!ok) sched->print(true);

    // Advance the particles as usual.  Count this as one block step.

    idat->advance(system_time);
    block_steps += 1;
    total_steps += idat->ni;

    // Put the particles back into the scheduler.

    for (int jj = 0; jj < njlist; jj++) sched->add_particle(jlist[jj]);
}

static real E0 = 0;
static real Emerge = 0;
void update_merger_energy(real dEmerge) {Emerge += dEmerge;}

void jdata::print()
{
    const char *in_function = "jdata::print";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    real E = get_energy();
    if (E0 == 0) E0 = E;
    real pe = get_pot();
    real total_mass = 0;
    for (int j = 0; j < nj; j++) total_mass += mass[j];
    int n_async = 0;
    for (int j = 0; j < nj; j++) if (time[j] < system_time) n_async++;

    // Note: get_energy() predicts all particles and recomputes the
    // total potential and kinetic energies.  For the most accurate
    // results, synchronize all particles first, but note that this
    // will also tie the results of the N-body calculation to the
    // output interval chosen.

    real dnnmin = huge;
    int jmin = 0;
    for (int j = 0; j < nj; j++)
 	if (dnn[j] < dnnmin) {
 	    jmin = j;
 	    dnnmin = dnn[j];
 	}

    if (mpi_rank == 0) {

	cout << endl;
	PRC(system_time); PRC(nj); PRC(total_mass); PRL(n_async);
	real rvir = -0.5*total_mass*total_mass/pe;
	PRC(pe); PRL(rvir);
	int p = cout.precision(12);
	PRC(E); PRL(E-E0);
	PRC(Emerge); PRL(E-E0-Emerge);
	cout.precision(p);
	print_percentiles();			// OK: runs on process 0 only

	if (NN) {
	    int jnn = nn[jmin];
	    PRC(jmin); PRC(jnn); PRL(dnnmin);	// GPU implementation
	}					// of nn is incomplete

	real elapsed_time = get_elapsed_time();
	real Gflops_elapsed = 6.e-8*(nj-1)*total_steps/(elapsed_time+tiny);
	real user_time, stime;
	get_cpu_time(user_time, stime);
	real Gflops_user = 6.e-8*(nj-1)*total_steps/(user_time+tiny);
	PRC(block_steps); PRL(total_steps); PRL(total_steps/block_steps);
	PRC(elapsed_time); PRL(user_time);
	PRC(Gflops_elapsed); PRL(Gflops_user);

#ifdef GPU		// needed because of explicit g6_npipes() below
	if (use_gpu) {
	    PRC(gpu_calls); PRL(gpu_total);
	    int npipes = g6_npipes();
	    real Gflops_elapsed_max = 0;
	    if (elapsed_time > 0)
		Gflops_elapsed_max
		    = 6.e-8*(nj-1)*gpu_calls*npipes/(elapsed_time+tiny);
	    real Gflops_user_max = 0;
	    if (user_time > 0)
		Gflops_user_max = 6.e-8*(nj-1)*gpu_calls*npipes
							/(user_time+tiny);
	    PRC(Gflops_elapsed_max); PRL(Gflops_user_max);
	}
#endif
    }

    if (use_gpu) {
	// get_densities_on_gpu();		// not tested yet
    }
}

void jdata::log_output()
{
    // Rudimentary log output...

    if (mpi_rank == 0) {
	vector<real> mlist, rlist;
	mlist.push_back(0.1);
	mlist.push_back(0.5);
	mlist.push_back(0.9);
	rlist.clear();
	get_lagrangian_radii(mlist, rlist);
	cout << "%%% " << system_time << " "
	     << rlist[0] << " " << rlist[1] << " " << rlist[2]
	     << endl << flush;
	    }
}

void jdata::cleanup()
{
    const char *in_function = "jdata::cleanup";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    if (name) delete [] name;
    if (id) delete [] id;
    if (nn) delete [] nn;
    if (mass) delete [] mass;
    if (radius) delete [] radius;
    if (pot) delete [] pot;
    if (dnn) delete [] dnn;
    if (time) delete [] time;
    if (timestep) delete [] timestep;
    if (pos) delete [] pos;
    if (vel) delete [] vel;
    if (acc) delete [] acc;
    if (jerk) delete [] jerk;
    if (pred_pos) delete [] pred_pos;
    if (pred_vel) delete [] pred_vel;
    id = nn = NULL;
    name = NULL;
    mass = radius = pot = dnn = time = timestep = NULL;
    pos = vel = acc = jerk = pred_pos = pred_vel = NULL;
    nj = 0;
    njbuf = 0;
    inverse_id.clear();
    user_specified_id.clear();
}
