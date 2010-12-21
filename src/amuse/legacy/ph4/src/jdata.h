#ifndef JDATA_H
#define JDATA_H

// Define the jdata class: data and methods operating on the entire
// N-body system.

#include "stdinc.h"
#include <map>

class scheduler;
class idata;

// j-data:

#define JBUF_INC	8192

class jdata {

  private:
    int nj;
    int njbuf;

  public:

    MPI::Intracomm mpi_comm;		// communicator for the N-body system
    int mpi_size;
    int mpi_rank;

    bool have_gpu;
    bool use_gpu;

    real eps2, eta;
    real rmin;				// 90 degree turnaround distance
    real dtmin;				// time step for enabling nn check

    real block_steps, total_steps, gpu_calls, gpu_total;
    real system_time, predict_time;

    int close1, close2;			// close particles (within rmin)
    int coll1, coll2;			// colliding particles

    // NOTE: name and id are unique identifiers for a particle in the
    // j system; id is called "index" in AMUSE, but it doesn't directly
    // index the jdata arrays, as particles may migrate.

    string *name;
    int *id, *nn;
    map<int,int> inverse_id;

    real *mass, *radius, *pot, *dnn, *time, *timestep;
    real2 pos, vel, acc, jerk;
    real2 pred_pos, pred_vel;

    int get_nj() {return nj;}

    jdata() {
	nj = 0;
	njbuf = 0;
	mpi_comm = NULL;
	mpi_size = 0;
	mpi_rank = -1;
#ifdef GPU
	have_gpu = true;
	use_gpu = true;
#else
	have_gpu = false;
	use_gpu = false;
#endif
	eps2 = eta = rmin = dtmin = 0;
	block_steps = total_steps = gpu_calls = gpu_total = 0;
	system_time = predict_time = -1;
	coll1 = coll2 = -1;
	id = nn = NULL;
	inverse_id.clear();
	name = NULL;
	mass = radius = pot = dnn = time = timestep = NULL;
	pos = vel = acc = jerk = pred_pos = pred_vel = NULL;
    }

    void cleanup();		// (in jdata.cc)
    ~jdata() {cleanup();}

    // In jdata.cc:

    void setup_mpi(MPI::Intracomm comm);
    int add_particle(real pmass, real pradius, vec ppos, vec pvel,
		     int pid = -1);
    void remove_particle(int j);
    void initialize_arrays();
    int get_inverse_id(int i);
    void check_inverse_id();
    void set_initial_timestep();
    real get_pot(idata *id = NULL);
    real get_kin(idata *id = NULL);
    real get_energy(idata *id = NULL);
    void predict(int j, real t);
    void predict_all(real t, bool full_range = false);
    void advance(idata& id, scheduler& sched);
    void synchronize_all(idata& id, scheduler *sched = NULL);
    void synchronize_list(int jlist[], int njlist,
			  idata& id, scheduler *sched);
    void print(idata *id = NULL);

    // In gpu.cc:

    void initialize_gpu(bool reinitialize = false);
    void update_gpu(int jlist[], int njlist);
    void get_densities_on_gpu();

    // In diag.cc:

    vec get_center();
    void print_percentiles(vec center = 0);

    // In close_encounter.cc:

    bool resolve_encounter(idata& id, scheduler& sched);
};

#endif
