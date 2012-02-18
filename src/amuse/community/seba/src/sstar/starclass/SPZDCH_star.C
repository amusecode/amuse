//
// SPZDCH19_star.C
//

#include "SPZDCH_star.h"

void SPZDCH_star::instantaneous_element() {

  luminosity = 1;
  effective_radius = radius = 1;
  core_radius = 1;

  envelope_mass = get_total_mass();
  core_mass = 0;

  radius = core_radius;

}

void SPZDCH_star::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

	next_update_age = relative_age + cnsts.safety(maximum_timestep);

        update();
	stellar_wind(dt);
}


void SPZDCH_star::update() {

// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

}

void SPZDCH_star::stellar_wind(const real dt) {

    if (!get_use_hdyn()) 
	cerr << " No stellar dynamical information present in SPZDCH_star" 
	     << endl;

    int N = the_node->get_root()->n_leaves();

    real end_time = 1;
    real t_relax = 1;
    if (find_qmatch(the_node->get_root()->get_dyn_story(), "t_relax"))
	t_relax = getrq(the_node->get_root()->get_dyn_story(), "t_relax");
    else if (N>1) {
	real total_mass = the_node->get_root()->get_mass();
	real potential_energy=1;
	if (find_qmatch(the_node->get_root()->get_dyn_story(), 
			"potential_energy"))
	    potential_energy = getrq(the_node->get_root()->get_dyn_story(), 
				     "potential_energy");
	real r_virial = -0.5 * total_mass * total_mass / potential_energy;
	t_relax = 9.62e-2 * sqrt(pow(r_virial, 3) / total_mass)
	    * N / log10(0.4 * N);
    }
    else 
	cerr << "\nwarning: no relaxation time available in SPZDCH_star"
	     << endl;

    real dm = cnsts.parameters(relaxation_driven_mass_loss_constant) 
            * envelope_mass / t_relax; 

  real wind_mass = dm * dt;

  if (wind_mass >= envelope_mass) 
    wind_mass = envelope_mass;

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_mass(wind_mass);

  return;
}

real SPZDCH_star::accretion_limit(const real mdot, const real dt) {

        return 0;
     }

star* SPZDCH_star::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = 0;
      return this;
}

void SPZDCH_star::adjust_accretor_age(const real mdot,
				      const bool rejuvenate=true) {

  return;

}

real SPZDCH_star::zeta_thermal() {
        return 0;
     }

star* SPZDCH_star::reduce_mass(const real mdot) {

  if (envelope_mass>=mdot) {
    envelope_mass -= mdot;

  } 
  else {
    cerr << "WARNING: star* SPZDCH_star::reduce_mass(const real mdot=" 
	 << mdot<< ")"<< endl; 
    cerr << "SPZDCH_star has negative mass..." << endl;
    envelope_mass = 0;
  }

  // Note that the core_mass = 0 at initialization.
  if (envelope_mass<=0) {
      cerr << "SPZDCH_star reduced to zero mass at t="
	   << get_current_time() << endl;
      cerr << "            last mdot was: " << mdot << endl;
  }
           
  return this;
}

void SPZDCH_star::adjust_next_update_age() {
     
// (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;
  next_update_age = relative_age + cnsts.safety(maximum_timestep);
}

real SPZDCH_star::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}


void SPZDCH_star::update_wind_constant(const real tscale) {

  wind_constant = tscale;

}

