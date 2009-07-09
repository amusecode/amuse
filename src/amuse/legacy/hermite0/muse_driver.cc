#include <iostream>
#include <cstdlib>
#include "muse_dynamics.h"
#include "local.h"
using namespace std;

// C++ test code for the hermite0 module.
// Run with ../../muse/in.20 for comparison with test.py.
// Default output channel is cout.  In test mode we send
// graphical data to cout, the rest to cerr.

double E0 = 0;

void print_diag(ostream& s)
{
    if (E0 == 0) E0 = get_kinetic_energy() + get_potential_energy();
    s << get_time() << " " << get_time_step() << " "
      << get_kinetic_energy() + get_potential_energy() - E0 << " "
      << get_n_steps() << endl;
    
    // (Note: get_n_steps() is not a standard Gravity function.)
}

main(int argc, char *argv[])
{
    double dt_next = 1.0;
    double t_max = 100;
    double radius = 0.025;
    double eps = 1.e-3;
    double eta = 1./64;
    bool test_mode = false;
    bool reeval = false;

    ostream* sout = &cout;

    for (int i = 1; i < argc; i++)
        if (argv[i][0] == '-')
            switch (argv[i][1]) {
                case 'a':	eta = atof(argv[++i]);
				break;
                case 'd':	dt_next = atof(argv[++i]);
				break;
                case 'e':	eps = atof(argv[++i]);
				break;
                case 'E':	reeval = true;
				break;
                case 'r':	radius = atof(argv[++i]);
				break;
                case 't':	t_max = atof(argv[++i]);
				break;
                case 'T':	test_mode = true;
				sout = &cerr;
				break;
            }

    // Initialize non-default parameters:

    set_dt_dia(1.e9);
    set_dt_param(eta);
    set_eps(eps);
    setup_module(reeval, test_mode);

    int n;
    
    // Read data from cin:

    n = 0;
    while (!cin.eof()) {
	dynamics_state s; 
	cin >> s.id;
	cin >> s.mass;
	s.radius = radius;
	cin >> s.x >> s.y >> s.z;
	cin >> s.vx >> s.vy >> s.vz;
	if (!cin.eof()) n = add_particle(s);
    }

    double t_next = 0;
    initialize_particles(t_next);
    *sout << "n = " << get_number()
	  << ", E = " << get_kinetic_energy()+get_potential_energy()
	  << ", tdyn = " << get_dynamical_time_scale() << endl;
    print_diag(*sout);

    // Run the code:

    while (t_next < t_max) {

	t_next += dt_next;

	// Evolve is supposed to loop to time t_next, but this
	// secondary loop is necessary to handle collisions, which
	// cause an immediate return with the id of the primary.

	bool coll = false;
	while (get_time() + get_time_step() <= t_next) {

	    int id1 = evolve(t_next);

	    // System is at time t_next unless a collision occurred.

	    if (id1 >= 0) {

		evolve(get_time(), 1);	// (sync not actually necessary here)

		int id2 = get_colliding_secondary(id1);
		if (id2 >= 0) {

		    // Only flag one collision per output interval.

		    if (!coll)
			printf(
			"  detected collision of %d and %d at time %.5f\n",
			id1, id2, get_time());
		    coll = true;
		}
	    }
	}

	print_diag(*sout);
    }

    cleanup_module();
}
