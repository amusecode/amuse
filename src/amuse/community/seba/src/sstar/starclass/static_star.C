//
// SPZDCH19_star.C
//

#include "static_star.h"

void static_star::instantaneous_element() {

  luminosity = 1;
  effective_radius = radius = 1;
  core_radius = 1;

  envelope_mass = get_total_mass();
  core_mass = 0;

  radius = core_radius;

}

void static_star::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

	next_update_age = relative_age + cnsts.safety(maximum_timestep);

        update();
}


void static_star::update() {

// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

}
