//----------------------------------------------------------------------
//
// Parallel fourth-order Hermite integrator with GRAPE/GPU
// acceleration.  Run time command-line options are:
//
//	-a eta		set accuracy parameter			[0.14]
//	-c close	manage close encounters [0-2]		[1]
//	-d dt_out	set time interval for diagnostic output	[0.5]
//	-e eps		set softening parameter			[0.01]
//	-f infile	specify input file			[none]
//	-g		suppress GPU use			[use GPU]
//	-n N		set number of stars (if no infile)	[1000]
//	-s seed		set initial random seed (if no infile)	[12345]
//	-t t_max	set integration time			[0.5]
//
// Currently no snapshot output is produced for restart or other
// purposes.  To be fixed...
//
// version 0, September 2010:	Steve McMillan (steve@physics.drexel.edu)
// version 1, December 2010:	added treatment of close encounters
// version 2. April 2011:	exposed extra functionality to AMUSE
//
//----------------------------------------------------------------------

#include "stdinc.h"
#include "jdata.h"
#include "idata.h"
#include "scheduler.h"

// Main program sets up initial conditions and parameters, then passes
// them to the worker code.

#include <sys/stat.h> 

void initialize_particles(jdata &jd, int nj, int seed, real vfac,
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

	jd.system_time = system_time = 0;

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
	    v *= vfac;

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
	jd.system_time = system_time;
	for (int j = 0; j < nj; j++) {
	    s >> id >> mass;
	    s >> pos[0] >> pos[1] >> pos[2];
	    s >> vel[0] >> vel[1] >> vel[2];

	    // Force the id to be the one in the input file.

	    jd.add_particle(mass, radius, pos, vfac*vel, id);
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

    // Set the system center of mass position and velocity to zero.

    jd.to_com();
}



void run_hermite4(int ntotal, int seed, char *file, bool use_gpu,
		  real eps2, real eta, real t_max, real dt_out,
		  real vfac, int manage_encounters)
{
    real system_time = 0;

    // Set up the jdata parameters and data structures.

    jdata jd;
    jd.setup_mpi(MPI::COMM_WORLD);
    jd.setup_gpu();		// set have_gpu and default use_gpu
    if (use_gpu && !jd.have_gpu) use_gpu = false;
    jd.use_gpu = use_gpu;
    jd.eps2 = eps2;
    jd.eta = eta;
    jd.set_manage_encounters(manage_encounters);

    initialize_particles(jd, ntotal, seed, vfac, file, system_time);
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

    jd.print();
    sched.print();
    real t_out = dt_out;

    real dt_spec = 1., t_spec = dt_spec;
    jd.spec_output("%%%");

    bool remove = false;	// don't check GPU bug in remove_particle()

    int step = 0;
    while (jd.system_time < t_max) {

	jd.advance_and_check_encounter();
	step++;

	// Problem-specific output.

	if (jd.system_time >= t_spec) {
	    jd.synchronize_all();
	    jd.spec_output("%%%");
	    t_spec += dt_spec;
	}

	if (jd.system_time >= t_out || jd.system_time >= t_max) {

	    // Routine diagnostic output.

	    jd.synchronize_all();
	    jd.print();

	    // End of run output.

	    if (jd.system_time >= t_max) {
		sched.print();
		if (jd.mpi_rank == 0 && jd.binary_list.size() > 0) {
		    cout << "binaries:" << endl;
		    for (unsigned int i = 0; i < jd.binary_list.size(); i++) {
			real semi = jd.binary_list[i].semi;
			real energy = -0.5*jd.binary_list[i].mass1
					  *jd.binary_list[i].mass2/semi;
			cout << "    " << jd.binary_list[i].binary_id
			     << " "    << jd.binary_list[i].comp1
			     << " "    << jd.binary_list[i].comp2
			     << " "    << semi
			     << " "    << jd.binary_list[i].ecc
			     << " "    << energy
			     << endl << flush;
		    }
		}
	    }
	    if (jd.system_time >= t_out) t_out += dt_out;
	}

	// Test code for remove_particle() bug.

	if (remove && jd.system_time > 1) {

	    jd.synchronize_all();
	    jd.print();

	    // Remove some particles.

	    cout << endl;
	    for (int j = 5; j <= 50; j += 5) {
		PRC(jd.system_time); cout << "removing "; PRL(j);
		jd.remove_particle(j);
	    }
	    
	    if (!jd.use_gpu)
		jd.predict_all(jd.system_time, true);	// set pred quantities
	    else
		jd.initialize_gpu(true);		// reload the GPU

	    id.setup();					// compute acc and jerk
	    sched.initialize();				// reconstruct scheduler

	    remove = false;
	    jd.print();
	}
    }
    jd.mpi_comm.Barrier();
}



int main(int argc, char *argv[])
{
    // Defaults:

    real eta = 0.14;
    real dt_out = 0.5;
    real eps = 0.01;
    real eps2 = eps*eps;
    char *infile = NULL;
    int ntotal = 1000;
    int seed = 12345;
    real t_max = 10.0;
    real vfac = 1.0;
#ifdef GPU
    bool use_gpu = true;
#else
    bool use_gpu = false;
#endif
    int manage_encounters = 1;

    // Parse the command line.

    for (int i = 1; i < argc; i++)
        if (argv[i][0] == '-')
            switch (argv[i][1]) {
                case 'a':	eta = atof(argv[++i]);
				break;
                case 'c':	manage_encounters = atoi(argv[++i]);
				break;
                case 'd':	dt_out = atof(argv[++i]);
				break;
                case 'e':	eps = atof(argv[++i]);
				eps2 = pow(eps,2);
				break;
		case 'f':	infile = argv[++i];
				break;
		case 'g':	use_gpu = !use_gpu;
				break;
                case 'n':	ntotal = atoi(argv[++i]);
				break;
                case 's':	seed = atoi(argv[++i]);
				break;
                case 't':	t_max = atof(argv[++i]);
				break;
                case 'v':	vfac = atof(argv[++i]);
				break;
            }

    MPI::Init(argc, argv);

    // Echo the command-line arguments.

    if (MPI::COMM_WORLD.Get_rank() == 0) {
	cout << "Parameters:" << endl << flush;
	cout << "    -a "; PRL(eta);
	cout << "    -c "; PRL(manage_encounters);
	cout << "    -d "; PRL(dt_out);
	cout << "    -e "; PRL(eps);
	if (infile) {
	    cout << "    -f "; PRL(infile);
	} else
	    cout << "    -f (no file)" << endl << flush;
	cout << "    -g "; PRL(use_gpu);
	if (!infile) {
	    cout << "    -n "; PRL(ntotal);
	}
	cout << "    -s "; PRL(seed);
	cout << "    -t "; PRL(t_max);
	cout << "    -v "; PRL(vfac);
    }

    // Run the integrator.

    run_hermite4(ntotal, seed, infile, use_gpu,
		 eps2, eta, t_max, dt_out, vfac,
		 manage_encounters);
    MPI_Finalize();
}
