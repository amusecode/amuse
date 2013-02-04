
// *********************
// * j-data operations *
// *********************
//
// Global functions:
//
//	void jdata::setup_mpi(MPI::Intracomm comm)
//	void jdata::setup_gpu()
//	void jdata::set_manage_encounters(int m)
//	int  jdata::add_particle(real pmass, real pradius,
//				 vec ppos, vec pvel, int pid)
//	void jdata::remove_particle(int j)
//	void jdata::initialize_arrays()
//	int  jdata::get_inverse_id(int i)
//	void jdata::check_inverse_id()
//	void jdata::set_initial_timestep()
//	real jdata::get_pot(bool reeval)			* MPI *
//	real jdata::get_kin()
//	real jdata::get_energy(bool reeval)
//	real jdata::get_total_mass()
//	void jdata::predict(int j, real t)
//	void jdata::predict_all(real t, bool full)
//	void jdata::advance()
//	void jdata::synchronize_all()
//	void jdata::synchronize_list(int jlist[], int njlist)
//	bool jdata::is_multiple(int i)
//	void jdata::print(idata *id)
//	void jdata::cleanup()

#include "jdata.h"
#include "scheduler.h"
#include "idata.h"

#ifdef GPU
#include "grape.h"
#endif

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

#ifndef NOMPI
void jdata::setup_mpi(MPI::Intracomm comm)
{
    const char *in_function = "jdata::setup_mpi";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    mpi_comm = comm;
    mpi_size = mpi_comm.Get_size();
    mpi_rank = mpi_comm.Get_rank();
}
#endif

void jdata::setup_gpu()
{
    const char *in_function = "jdata::setup_gpu";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Determine GPU availability and use, based on the option
    // actually compiled into the module.

#ifdef GPU
    have_gpu = use_gpu = true;
#else
    have_gpu = use_gpu = false;
#endif
}

void jdata::set_manage_encounters(int m)
{
    // Setter for public data, but also checks the AMUSE stopping
    // conditions.  Use this instead of jd.manage_emcounters = ...

    // AMUSE STOPPING CONDITIONS SUPPORT
    int is_collision_detection_enabled;
    is_stopping_condition_enabled(COLLISION_DETECTION,
				  &is_collision_detection_enabled);
				  
    if (is_collision_detection_enabled) {
	manage_encounters = 4;	// unnecessary
    }

    if (manage_encounters == 4 && m != 4){
	if (mpi_rank == 0) {
	    cout << "warning: setting manage_encounters = " << m
		 << " overrides default stopping condition"
		 << endl << flush;
	}
    }

    manage_encounters = m;
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

    int input_pid = pid;
    if (pid < 0) pid = get_particle_id();

    if (mpi_rank == 0 && input_pid >= 0 && pid == 0)
	cout << "jdata:add_particle: adding pid = 0" << endl << flush;

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

    // Update the inverse ID structure, such that inverse[id[j]] = j.
    // Use an STL map, so we don't have to worry about the range or
    // sparseness of the IDs in use.

    inverse_id[pid] = nj;

    // Set the timestep if required, and update the scheduler, if it
    // exists.

    timestep[nj] = dt;
    if (dt > 0) {
	if (sched) {
	    // cout << "sched.add " << nj << " (" << pid << ")"
	    //	 << endl << flush;
	    sched->add_particle(nj);
	}
    }

    nj++;

    if (0 && system_time > 0 && mpi_rank == 0) {
 	cout << "add_particle: "; PRC(system_time);
 	PRC(pmass); PRC(pid); PRL(nj);
	cout << "    pos:";
	for (int k = 0; k < 3; k++) cout << " " << pos[nj-1][k];
	cout << endl << "    vel:";
	for (int k = 0; k < 3; k++) cout << " " << vel[nj-1][k];
	cout << endl << flush;
    }

    // AMUSE bookkeeping:

    UpdatedParticles.push_back(UpdatedParticle(pid, 2));

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

    // AMUSE bookkeeping.  Place the particle on the UpdatedParticles
    // list for later removal, or delete it from the list if it has
    // recently been added.

    bool add_to_list = true;
    for (unsigned int iup = 0; iup < UpdatedParticles.size(); iup++)
	if (UpdatedParticles[iup].index_of_particle == id[j]) {
	    UpdatedParticles.erase(UpdatedParticles.begin()+iup);
	    add_to_list = false;
	    break;
	}
    if (add_to_list)
	UpdatedParticles.push_back(UpdatedParticle(id[j], 1));

    nj--;

    if (0 && system_time > 0 && mpi_rank == 0) {
 	cout << "remove_particle: "; PRC(system_time);
	cout << "id = " << id[j] << ",  mass = " << mass[j] << ",  ";
 	PRL(nj);
	cout << "    pos:";
	for (int k = 0; k < 3; k++) cout << " " << pos[j][k];
	cout << endl << "    vel:";
	for (int k = 0; k < 3; k++) cout << " " << vel[j][k];
	cout << endl << flush;
    }

    if (j < nj) {
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

    // *** To be refined.  Probably want to specify a scaling. ***

    rmin = 2./nj;		// 90 degree turnaround in standard units
    dtmin = eta*(2./nj) * 4;	// for nn check: final 4 is a fudge factor
				// -- not used??
    if (mpi_rank == 0) {
	PRC(rmin); PRL(dtmin);
    }

    // Starting point for binary IDs.

    int bb = log10((real)nj) + 1;
    binary_base = pow(10., bb);

    // Final consistency check on GPU availability/use.

    if (use_gpu && !have_gpu) {
	if (mpi_rank == 0) cout << endl << "GPU unavailable" << endl << flush;
	use_gpu = false;
    }

    if (!use_gpu) {

	if (have_gpu && mpi_rank == 0)
	    cout << endl << "GPU use disabled" << endl << flush;

	predict_all(system_time, true);	 // set all predict_time, pred_xxx

    } else {

	// Copy all data into the GPU.

	cout << endl << "initializing GPU(s) for " << mpi_rank 
	     << endl << flush;
	initialize_gpu();
	if (DEBUG > 1) cout << "GPU initialization done for " << mpi_rank
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

    for (int j = 0; j < nj; j++)
	if (timestep[j] <= 0) {
	    real a2 = 0, j2 = 0;
	    for (int k = 0; k < 3; k++) {
		a2 += pow(acc[j][k], 2);
		j2 += pow(jerk[j][k], 2);
	    }

	    real firststep;
	    if (eta == 0.0)
		firststep = 0.0625;
	    else if (a2 == 0.0 || j2 == 0.0)
		firststep = 0.0625 * eta;
	    else
		firststep = 0.0625 * eta * sqrt(a2/j2);	// conservative
    
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
		if (r2 > _TINY_)
		    dpot -= mass[j]/sqrt(r2+eps2);
	    }
	    mypot += mass[jd]*dpot;
	}

#ifndef NOMPI
	mpi_comm.Allreduce(&mypot, &total_pot, 1, MPI_DOUBLE, MPI_SUM);
#else
    total_pot = mypot;
#endif

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

real jdata::get_kin()
{
    const char *in_function = "jdata::get_kin";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Return the total kinetic energy of the (predicted) j system.

    predict_all(system_time, true);		// make this function parallel?
    real kin2 = 0;
    for (int j = 0; j < nj; j++) {
	real v2 = 0;
	for (int k = 0; k < 3; k++) v2 += pow(pred_vel[j][k],2);
	kin2 += mass[j]*v2;
    }
    return kin2/2;
}

real jdata::get_energy(bool reeval)		// default = false
{
    const char *in_function = "jdata::get_energy";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Note: energy includes Emerge.

    return get_pot(reeval) + get_kin() + Emerge;
    // return get_pot(reeval) + get_kin() + get_binary_energy()
}

real jdata::get_total_mass()
{
    real mtot = 0;
    for (int j = 0; j < nj; j++) mtot += mass[j];
    return mtot;
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

    for (int j = j_start; j < j_end; j++) {
	real dt = t - time[j];
	if (dt == 0) {

	    // Just set pred_??? from ???.

	    for (int k = 0; k < 3; k++) {
		pred_pos[j][k] = pos[j][k];
		pred_vel[j][k] = vel[j][k];
	    }
	} else {
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
	// the particle range associated with the current process.

	predict_all(tnext);		// j pos, vel --> j pred_pos, pred_vel
    }

    // All of the actual work is done by idata::advance.

    idat->advance(tnext);
    block_steps += 1;
    total_steps += idat->ni;

    // Note that system_time remains unchanged until the END of the step.

    system_time = tnext;
    sched->update();
}

#define EPS 0.001	// see couple/multiples.py

bool jdata::advance_and_check_encounter()
{
    bool status = false;
    int collision_detection_enabled;
    advance();

    // Optionally manage close encounters.  AMUSE stopping conditions
    // are enabled with manage_encounters = 4.  Return true if an
    // encounter has been detected (handled by the top level loop or
    // by AMUSE).
    
    // AMUSE STOPPING CONDITIONS SUPPORT

    is_stopping_condition_enabled(COLLISION_DETECTION,
				  &collision_detection_enabled);
    if (collision_detection_enabled) {
        if (coll1 >= 0) {

	    // Duplicate the check made in multiples.py, to
	    // avoid returning too many times.

	    int j1 = get_inverse_id(coll1);
	    int j2 = get_inverse_id(coll2);
	    real r = 0, v = 0, vr = 0;
	    for (int k = 0; k < 3; k++) {
		real dx = pos[j1][k]-pos[j2][k];
		real dv = vel[j1][k]-vel[j2][k];
		r += dx*dx;
		v += dv*dv;
		vr += dx*dv;
	    }

	    if (vr < EPS*sqrt(r*v)) {
		int stopping_index = next_index_for_stopping_condition();
		set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
		set_stopping_condition_particle_index(stopping_index, 0, coll1);
		set_stopping_condition_particle_index(stopping_index, 1, coll2);
		status = true;
	    }
        }
        return status;
    }
    
    
    if (!manage_encounters || eps2 > 1.0e-99 || manage_encounters == 4)
        return false;

    if (close1 >= 0)
        status = resolve_encounter();

    if (status)
        if (0 && mpi_rank == 0) {
            cout << "after resolve_encounter" << endl << flush;
            PRL(get_energy());
        }

    return status;
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

    if (njlist == 0) return;		// nothing to do

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

// Maintain two measures of the energy in merged objects.  The total
// merger energy Emerge includes binary binding energy and the tidal
// error incurred when the merger occurred.  The total binary energy
// omits the tidal effects.

void jdata::update_merger_energy(real dEmerge)
{
    Emerge += dEmerge;
}

real jdata::get_binary_energy()
{
    // Return the total binary binding energy.

    real Eb = 0;
    for (unsigned int ib = 0; ib < binary_list.size(); ib++) {
	real m1m2 = binary_list[ib].mass1 * binary_list[ib].mass2;
	real a = binary_list[ib].semi;
	Eb -= 0.5*m1m2/a;
    }
    return Eb;
}

bool jdata::is_multiple(int i)		// note: i is ID, not j index
{
    for (unsigned int ib = 0; ib < binary_list.size(); ib++)
	if (binary_list[ib].binary_id == i) return true;
    return false;
}

void jdata::print()
{
    const char *in_function = "jdata::print";
    if (DEBUG > 2 && mpi_rank == 0) PRL(in_function);

    // Print standard diagnostic information on the j-data system.

    real E = get_energy();	// note: energy includes Emerge
    if (E0 == 0) E0 = E;
    real pe = get_pot();
    real total_mass = get_total_mass();
    int n_async = 0;
    for (int j = 0; j < nj; j++) if (time[j] < system_time) n_async++;

    // Note: get_energy() predicts all particles and recomputes the
    // total potential and kinetic energies.  For the most accurate
    // results, synchronize all particles first, but note that this
    // will also tie the results of the N-body calculation to the
    // output interval chosen.

    real dnnmin = _INFINITY_;
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
	PRC(Emerge); PRL(get_binary_energy());
	PRC(E); PRL(E-E0);
	cout.precision(p);
	vec cmpos, cmvel;
	get_com(cmpos, cmvel);
	cout << "com:  " << cmpos << endl << flush;
	get_mcom(cmpos, cmvel);
	cout << "mcom: " << cmpos << endl << flush;
	print_percentiles(get_center());	// OK: runs on process 0 only

	if (NN) {
	    int jnn = nn[jmin];
	    PRC(jmin); PRC(jnn); PRL(dnnmin);	// GPU implementation
	}					// of nn is incomplete

	real elapsed_time = get_elapsed_time();
	real Gflops_elapsed = 6.e-8*(nj-1)*total_steps/(elapsed_time+_TINY_);
	real user_time, stime;
	get_cpu_time(user_time, stime);
	real Gflops_user = 6.e-8*(nj-1)*total_steps/(user_time+_TINY_);
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
		    = 6.e-8*(nj-1)*gpu_calls*npipes/(elapsed_time+_TINY_);
	    real Gflops_user_max = 0;
	    if (user_time > 0)
		Gflops_user_max = 6.e-8*(nj-1)*gpu_calls*npipes
							/(user_time+_TINY_);
	    PRC(Gflops_elapsed_max); PRL(Gflops_user_max);
	}
#endif
    }

    if (use_gpu) {
	// get_densities_on_gpu();		// not tested yet
    }
}

void jdata::spec_output(const char *s)		// default = NULL
{
    // Problem-specific output.  Careful with parallel functions!

    real pot = get_pot();

    if (mpi_rank == 0) {
	vector<real> mlist, rlist;
	mlist.push_back(0.005);
	mlist.push_back(0.01);
	mlist.push_back(0.02);
	mlist.push_back(0.05);
	mlist.push_back(0.10);
	mlist.push_back(0.25);
	mlist.push_back(0.50);
	mlist.push_back(0.75);
	mlist.push_back(0.90);
	rlist.clear();
	get_lagrangian_radii(mlist, rlist, get_center());

	real mtot = get_total_mass();
	real rvir = -0.5*mtot*mtot/pot;
	real kin = get_kin();
	real etot = kin + pot;
	real qvir = -kin/pot;

	if (s) cout << s << " ";
	cout << system_time << " " << mtot << " "
	     << kin << " " << pot << " " << etot << " "
	     << rvir << " " << qvir;

	// for (unsigned int i = 0; i < mlist.size(); i++)
	//     cout << " " << rlist[i];
	cout << endl << flush;
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
    binary_list.clear();
}

void jdata::to_com()
{
    // Place the system in the center of mass frame.

    vec cmpos, cmvel;
    get_com(cmpos, cmvel);
    for (int j = 0; j < nj; j++) {
	for (int k = 0; k < 3; k++) {
	    pos[j][k] -= cmpos[k];
	    vel[j][k] -= cmvel[k];
	}
    }
}
