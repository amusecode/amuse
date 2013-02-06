#ifndef JDATA_H
#define JDATA_H

// Define the jdata class: data and methods operating on the entire
// N-body system.

#include "stdinc.h"
#include <vector>
#include <map>
#include <algorithm>

class scheduler;
class idata;
class hdyn;

typedef struct {
    int status;
    vec cmpos;
    vec cmvel;
    real scale2;
    hdyn *b;
} nbody_descriptor;

class binary {

  // The basic stable multiple structure is a binary tree, where each
  // node is a stable object (star, binary, hierarchical triple,
  // etc.).  Simply store component IDs, masses, and binary
  // properties.  In a hierarchy, one or both IDs refer to another
  // binary.  Orbital phases and orientations are not stored, and are
  // randomized on instantiation.

  public:

    int binary_id;			// CM ID
    int comp1, comp2;			// component IDs
    real mass1, mass2;			// component masses
    real semi, ecc;			// relative orbital elements

    binary(int bid, int id1, int id2, real m1, real m2, real a, real e)
    {binary_id = bid; comp1 = id1; comp2 = id2;
     mass1 = m1; mass2 = m2; semi = a; ecc = e;}
    ~binary() {}
};

class UpdatedParticle {

  // AMUSE bookkeeping.  A particle on the UpdatedParticle list has
  // been added or removed internally.  Use this class to communicate
  // this information back to the python level.
    
  public:

    int index_of_particle;
    int status;		// 1 ==> deleted, 2 ==> added
    
    UpdatedParticle():index_of_particle(-1),status(0) {}
    UpdatedParticle(int i, int s):index_of_particle(i), status(s) {}
    UpdatedParticle(const UpdatedParticle& src)
	:index_of_particle(src.index_of_particle), status(src.status) {}
};

// Note: after proper initialization:
//
//     jdata and scheduler have pointers to each other
//     jdata and idata have pointers to each other
//     idata has a pointer to scheduler
//
// Order of initialization: jdata, idata(jdata), scheduler(jdata).

#define JBUF_INC	8192

class jdata {

  public:

    int nj;
    int njbuf;

#ifndef NOMPI
    MPI::Intracomm mpi_comm;		// communicator for the N-body system
#endif

    int mpi_size;
    int mpi_rank;

    bool have_gpu;			// will be true if -DGPU is compiled in
    bool use_gpu;			// true if actually using GPU

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
    vector<int> user_specified_id;

    real *mass, *radius, *pot, *dnn, *time, *timestep;
    real2 pos, vel, acc, jerk;
    real2 pred_pos, pred_vel;

    // Pointers to other data structures.

    idata *idat;
    scheduler *sched;

    // Initial system energy.

    real E0;

    // Multiple structure management (see also stopping_conditions code):

    real Emerge;
    int manage_encounters;
    int binary_base;
    int binary_count;
    vector<binary> binary_list;

    // Manage internal removal/creation of particles.

    vector<UpdatedParticle> UpdatedParticles;

    jdata() {
	nj = 0;
	njbuf = 0;
#ifndef NOMPI
	mpi_comm = NULL;
	mpi_size = 0;
	mpi_rank = -1;
#else
	mpi_size = 1;
	mpi_rank = 0;
#endif
	have_gpu = false;		// correct values will be set at
	use_gpu = false;		// run time, in setup_gpu()
	eps2 = eta = rmin = dtmin = 0;
	block_steps = total_steps = gpu_calls = gpu_total = 0;
	system_time = predict_time = -1;
	coll1 = coll2 = -1;
	id = nn = NULL;
	inverse_id.clear();
	user_specified_id.clear();
	name = NULL;
	mass = radius = pot = dnn = time = timestep = NULL;
	pos = vel = acc = jerk = pred_pos = pred_vel = NULL;

	idat = NULL;
	sched = NULL;

	E0 = 0;
	Emerge = 0;
	manage_encounters = 1;
	binary_count = 0;
	binary_list.clear();
	UpdatedParticles.clear();
    }

    void cleanup();		// (in jdata.cc)
    ~jdata() {cleanup();}

    // In jdata.cc:

    void set_manage_encounters(int m);

#ifndef NOMPI
    void setup_mpi(MPI::Intracomm comm);
#endif

    void setup_gpu();
    int get_particle_id(int offset = 0);
    int add_particle(real pmass, real pradius, vec ppos, vec pvel,
		     int pid = -1, real dt = -1);
    void remove_particle(int j);
    void initialize_arrays();
    int get_inverse_id(int i);
    void check_inverse_id(const char *s = NULL);
    void set_initial_timestep();
    real get_pot(bool reeval = false);
    real get_kin();
    real get_energy(bool reeval = false);
    real get_total_mass();
    void predict(int j, real t);
    void predict_all(real t, bool full_range = false);
    void advance();
    bool advance_and_check_encounter();
    void synchronize_all();
    void synchronize_list(int jlist[], int njlist);
    void update_merger_energy(real dEmerge);
    real get_binary_energy();
    bool is_multiple(int i);
    void print();
    void spec_output(const char *s = NULL);
    void to_com();

    // In gpu.cc:

    void initialize_gpu(bool reinitialize = false);
    void update_gpu(int jlist[], int njlist);
    void sync_gpu();
    void get_densities_on_gpu();

    // In diag.cc:

    void get_com(vec& pos, vec& vel);
    void get_mcom(vec& pos, vec& vel,
		  real cutoff = 0.9,
		  int n_iter = 2);
    vec get_center();
    void get_lagrangian_radii(vector<real>& mlist,
			      vector<real>& rlist,
			      vec center = 0);
    void print_percentiles(vec center = 0);

    // In close_encounter.cc:

    bool resolve_encounter();

    // In two_body.cc:

    void two_body(int j1, int j2, vector<int> nbrlist);

    int find_binary(int id);
    void expand_binary(hdyn *bb, real *time_scale = NULL,
		       real *length_scale2 = NULL);
    void remove_binary(int id);
    void add_binary(hdyn *b);
    nbody_descriptor integrate_multiple(vector<int> jcomp,
					real dt_fac, real r2_fac);
    bool check_add_neighbors(nbody_descriptor &s,
			     vector<int>mult, vector<int> nbrlist);
    real rescale(hdyn *b, real r2);
    void multiple(int j1, int j2, vector<int> nbrlist);
};

#define PRRC(x) cout << "rank = " << mpi_rank << " " << #x << " = " << x << ",  " << flush
#define PRRL(x) cout << "rank = " << mpi_rank << " " << #x << " = " << x << endl << flush

#endif
