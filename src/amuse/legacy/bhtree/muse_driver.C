#include <iostream>
#include "muse_dynamics.h"
#include "local.h"
using namespace std;

// C++ test code for the BHTree module.
// Run with examples/plummer1024.dat for comparison with test.py.

void print_diag()
{
    printf("time = %.5f, dt = %.4e, E = %.6f\n",
	   get_time(), get_time_step(),
	   get_kinetic_energy() + get_potential_energy());
}

int main()
{
    extern BHTC_SYSTEM bhtcs;

    real eps = 0.05;

    set_timestep(0.015625);
    set_eps2_for_gravity(eps*eps);
    set_theta_for_tree(0.75);
    set_use_self_gravity(1);
    set_ncrit_for_tree(1024);
    set_dt_dia(10);

    setup_module();

    printf("eps = %.5f, theta = %.5f, ncrit = %d\n",
	   eps, bhtcs.theta_for_tree, bhtcs.ncrit_for_tree);

    real dt_next = 0.25;
    real t_max = 10;
    real radius = 0.0025;

    // Read data from cin:

    int n = 0;
    while (!cin.eof()) {
	dynamics_state s; 
	cin >> s.id;
	cin >> s.mass;
	s.radius = radius;
	cin >> s.x >> s.y >> s.z;
	cin >> s.vx >> s.vy >> s.vz;
	if (!cin.eof()) n = add_particle(s);
    }

    real t_next = 0;
    initialize_particles(t_next);
    printf("n = %d, tdyn = %.5f\n", get_number(), get_dynamical_time_scale());
    print_diag();

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

		int id2 = find_colliding_secondary(id1);
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

	print_diag();
    }

    cleanup_module();
}
