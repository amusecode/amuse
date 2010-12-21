
//----------------------------------------------------------------------
//
// Parallel fourth-order Hermite integrator with GRAPE/GPU
// acceleration.  Run time command-line options are:
//
//	-a eta		set accuracy parameter			[0.14]
//	-d dt_out	set time interval for diagnostic output	[0.5]
//	-e eps		set softening parameter			[0.01]
//	-f infile	specify input file			[none]
//	-g		suppress GPU use			[use GPU]
//	-n N		set number of stars (if no infile)	[1000]
//	-r		recompute the potential in get_energy()	[false]
//	-s seed		set initial random seed (if no infile)	[12345]
//	-t t_max	set integration time			[0.5]
//
// Currently no snapshot output is produced for restart or other
// purposes.  To be fixed...
//
// Version 0, September 2010:	Steve McMillan (steve@physics.drexel.edu)
//
//----------------------------------------------------------------------

#include "stdinc.h"
#include "jdata.h"
#include "idata.h"
#include "scheduler.h"

// Main program sets up initial conditions and parameters, then passes
// them to the worker code.

#include <sys/stat.h> 

static real scale = 1.0/(((unsigned long)1<<31)-1.0);
static real randinter(real a=0, real b=1)
{
    unsigned long r = random();
    real x = scale*r;
    return a + (b-a)*x;
}

void initialize_particles(jdata &jd, int nj, int seed,
			  char *infile, real &system_time)
{
    // Initialize basic j data: mass, radius, position, and velocity,
    // id, and name.

    // Note: every node sees the entire dataset, but operates on only
    // a portion of it.

    // Add particles to the jd structure one by one, managing memory
    // internally as needed.  This approach is adopted in anticipation
    // of the AMUSE new_particle() interface function.  The current
    // AMUSE specification has new_particle() set the id and name via
    // the function add_particle(), but we extend that here to allow
    // the id to be specified by the user.

    if (!infile) {

	// Create a simple homogeneous sphere of unit radius.  The
	// virial mean square velocity is 0.6 (G=1).

	system_time = 0;

	srandom(seed);
	real v2scale = 0.6;
	for (int j = 0; j < nj; j++) {

	    real mass = 1.0/nj;				// total mass = 1
	    real radius = 0;

	    real r = pow(randinter(), 1.0/3);		// uniform sphere
	    real costh = randinter(-1,1);
	    real sinth = sqrt(fmax(0,1-costh*costh));
	    real phi = randinter(0,2*M_PI);
	    vec pos = r*vec(sinth*cos(phi), sinth*sin(phi), costh);
		
	    real v = sqrt(2*v2scale*randinter());	// nonthermal, out of
							// virial equilibrium
	    costh = randinter(-1,1);
	    sinth = sqrt(fmax(0,1-costh*costh));
	    phi = randinter(0,2*M_PI);
	    vec vel = v*vec(sinth*cos(phi), sinth*sin(phi), costh);

	    // Let the system choose the id and name.

	    jd.add_particle(mass, radius, pos, vel);	// (return value is id)
	}

    } else {

	// Read data from a file (parallel read, as here, or serial
	// read from node 0 followed by distribution).  It is
	// convenient to save files in gzipped format, but this may
	// not work well in a parallel environment.  Check and unzip
	// here.  Allow infile to end in .gz or to be the unzipped
	// name.  A pfstream would be handy...

	ifstream s;

	bool zip = false;
	char infile1[1024], infile2[1024];
	strcpy(infile1, infile);		// the unzipped file
	strcpy(infile2, infile);		// the zipped file
	if (strstr(infile1+strlen(infile1)-3, ".gz")) {
	    infile1[strlen(infile1)-3] = '\0';
	    zip = true;
	} else {
	    strcat(infile2, ".gz");
	    s.open(infile2, ifstream::in);
	    if (s) {
		s.close();
		zip = true;
	    }
	}

	if (zip && jd.mpi_rank == 0) {
	    char command[1024];
	    sprintf(command, "gunzip %s", infile1);
	    // cout << command << endl << flush;

	    // System call may be undesirable...

	    system(command);
	}
	jd.mpi_comm.Barrier();

	s.open(infile1, ifstream::in);
	if (!s) {
	    if (jd.mpi_rank == 0)
		cout << "No input file found." << endl << flush;
	    exit(1);
	}

	// Note that all worker processes read the same input file.

	int isnap;
	int id;
	real mass, radius = 0;
	vec pos, vel;
	s >> isnap >> nj >> system_time;
	for (int j = 0; j < nj; j++) {
	    s >> id >> mass;
	    s >> pos[0] >> pos[1] >> pos[2];
	    s >> vel[0] >> vel[1] >> vel[2];

	    // Force the id to be the one in the input file.

	    jd.add_particle(mass, radius, pos, vel, id);
	}

	s.close();
	if (zip && jd.mpi_rank == 0) {
	    char command[1024];
	    sprintf(command, "gzip %s", infile1);
	    // cout << command << endl << flush;

	    // System call may be undesirable...

	    system(command);
	}
	jd.mpi_comm.Barrier();
    }
}



void run_hermite4(int ntotal, int seed, char *file, bool use_gpu,
		  real eps2, real eta, real t_max, real dt_out, bool repot)
{
    real system_time = 0;

    // Set up the jdata parameters and data structures.

    jdata jd;
    jd.setup_mpi(MPI::COMM_WORLD);
    jd.use_gpu = use_gpu;
    jd.eps2 = eps2;
    jd.eta = eta;

    initialize_particles(jd, ntotal, seed, file, system_time);
    jd.system_time = system_time;
    jd.initialize_arrays();
    idata id(&jd);	  // set up idata data structures (sets acc and jerk)
    jd.set_initial_timestep();		// set timesteps (needs acc and jerk)

    if (DEBUG > 1 && jd.mpi_rank == 0)
	cout << "idata and jdata initialization done" << endl << flush;

    // Initialize the scheduler.

    scheduler sched(&jd);
    if (DEBUG > 1 && jd.mpi_rank == 0)
	cout << "scheduler setup done" << endl << flush;

    // Loop to advance the system to time t_max.

    jd.print(&id);	// no repot needed because the system is synchronized
    sched.print();
    real t_out = dt_out;

    int step = 0;
    while (jd.system_time < t_max) {

	jd.advance(id, sched);
#if 0
	cout << endl << "after normal step" << endl;
	int p = cout.precision(10);
	PRC(jd.system_time); PRL(jd.get_energy(&id));
	cout.precision(p);
#endif

	// Special treatment of close encounters (for non-AMUSE code).

	if (jd.close1 >= 0 && jd.eps2 == 0) {
	    bool status = jd.resolve_encounter(id, sched);
	    if (status) {
		int p = cout.precision(10);
		PRC(jd.system_time); PRL(jd.get_energy(&id));
		cout.precision(p);
	    }
	}

	step++;
	if (jd.system_time >= t_out || jd.system_time >= t_max) {
	    jd.print(&id);
	    if (jd.system_time >= t_max) sched.print();
	    if (jd.system_time >= t_out) t_out += dt_out;
	}
    }
}



int main(int argc, char *argv[])
{
    // Defaults:

    real eta = 0.14;
    real dt_out = 0.5;
    real eps2 = 1.e-2;
    char *infile = NULL;
    int ntotal = 1000;
    bool repot = false;
    int seed = 12345;
    real t_max = 0.5;
#ifdef GPU
    bool use_gpu = true;
#else
    bool use_gpu = false;
#endif

    // Parse the command line.

    for (int i = 1; i < argc; i++)
        if (argv[i][0] == '-')
            switch (argv[i][1]) {
                case 'a':	eta = atof(argv[++i]);
				break;
                case 'd':	dt_out = atof(argv[++i]);
				break;
                case 'e':	eps2 = pow(atof(argv[++i]),2);
				break;
		case 'f':	infile = argv[++i];
				break;
		case 'g':	use_gpu = false;
				break;
                case 'n':	ntotal = atoi(argv[++i]);
				break;
                case 'r':	repot = true;
				break;
                case 's':	seed = atoi(argv[++i]);
				break;
                case 't':	t_max = atof(argv[++i]);
				break;
            }

    // Run the integrator.

    MPI::Init(argc, argv);
    run_hermite4(ntotal, seed, infile, use_gpu,
		 eps2, eta, t_max, dt_out, repot);
    MPI_Finalize();
}
