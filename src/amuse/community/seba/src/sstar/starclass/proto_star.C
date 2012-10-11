//
// proto_star.C
//
// see also src/star/sstar/stardyn/proto_sat_dyn.C

#include "proto_star.h"


// required only for making double star from proto_star
//#include "double_star.h"

void proto_star::instantaneous_element() {


//    real m = get_total_mass()*cnsts.parameters(solar_mass);
//    real r = effective_radius*cnsts.parameters(solar_radius);
//    rotation_period = 2*cnsts.mathematics(pi)
//                    * gyration_radius_sq()*m*pow(r, 2)/angular_momentum;

    rotation_period = 1000;

    real m_tot = core_mass + envelope_mass;
    core_mass = m_tot * cnsts.parameters(star_formation_efficiency);
    envelope_mass = m_tot - core_mass;

    real alpha, beta, gamma, delta;
 
    real log_mass = log10(core_mass);
    
    if (relative_mass > 1.334) {

	alpha = 0.1509 + 0.1709*log_mass;
	beta  = 0.06656 - 0.4805*log_mass;
	gamma = 0.007395 + 0.5083*log_mass;
	delta = (0.7388*pow(relative_mass, 1.679)
			 - 1.968*pow(relative_mass, 2.887))
	    	     / (1.0 - 1.821*pow(relative_mass, 2.337));

    } else {
      
	alpha = 0.08353 + 0.0565*log_mass;
	beta  = 0.01291 + 0.2226*log_mass;
	gamma = 0.1151 + 0.06267*log_mass;
	delta = pow(relative_mass, 1.25)
	    		* (0.1148 + 0.8604*relative_mass*relative_mass)
		    / (0.04651 + relative_mass*relative_mass);
    }
    cerr<<"check if this is still correct!!"<<endl;
    cerr<<"as base_main_sequence_luminosity changed from EFT to Tout 1996"<<endl;
    //luminosity = base_main_sequence_luminosity(core_mass);
    core_radius = delta;

    // proto stars are real big!
    radius = 10000*core_radius;
    effective_radius = max(effective_radius, radius);

    dump(cerr, false);

}

// evolve a proto_star star upto time argument according to
// the model discribed by Eggleton et al. 1989.
void proto_star::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {

        // proto_stars are static,
	// but lose mass in wind.

      }
      else {
	create_zero_age_object();
	return;
      }
  
      update();
      stellar_wind(dt);
}

void proto_star::create_zero_age_object() {

//  if (core_mass<=cnsts.parameters(minimum_brown_dwarf_mass)) {
//
//    star_transformation_story(Planet);
//    new planet(*this);
//  }

  // Make sure that effective_radius is resetted.
  effective_radius = core_radius;

  if (core_mass>cnsts.parameters(minimum_main_sequence)) {

    // see src/star/sstar/stardyn/proto_sat_dyn.C
    // no binaries yet, ST 25 feb 2009
    // create_binary_from_proto_star();
    cout << "ERROR: binary needs to be made, according to proto_star.C \n";
     return;
  }
  else {

    star_transformation_story(Brown_Dwarf);
    new brown_dwarf(*this);
    return;
  }
}


void proto_star::stellar_wind(const real dt) {

  real wind_mass = wind_constant 
                 * (pow(relative_age/next_update_age,
			cnsts.parameters(massive_star_mass_loss_law))
	         -  pow((relative_age-dt)/next_update_age,
			cnsts.parameters(massive_star_mass_loss_law)));

      if (is_binary_component())
         get_binary()->adjust_binary_after_wind_loss(this, 
                     wind_mass, dt);
      else
         reduce_mass(wind_mass);
    }


real proto_star::helium_core_mass() {

  real m_core = envelope_mass * cnsts.parameters(star_formation_efficiency);
      
  m_core = min(m_core, get_total_mass());

  return m_core;
}

star* proto_star::reduce_mass(const real mdot) {

      if (envelope_mass<=mdot) {
         envelope_mass = 0;
	 create_zero_age_object();
      }

      envelope_mass -= mdot;
      return this;
   }

star* proto_star::subtrac_mass_from_donor(const real dt, real& mdot) {

      //real mdot_temp = relative_mass*dt/get_binary()->get_donor_timescale();
      real mdot_temp = get_total_mass()*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_temp);

      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;
	 
	 create_zero_age_object();
      }

      envelope_mass -= mdot;
      return this;
   }

void proto_star::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {

  return;
}

void proto_star::adjust_next_update_age() {

     if (cnsts.parameters(star_formation_timescale)<0) 
	 next_update_age = abs(cnsts.parameters(star_formation_timescale))
                         * (radius/core_radius) * kelvin_helmholds_timescale();
     else
	 next_update_age = randinter(0., cnsts.parameters(star_formation_timescale));

    }

void proto_star::update_wind_constant() {

      wind_constant = envelope_mass;
    }

real proto_star::zeta_thermal() {

      real z = -10;

      return z;
   }

real proto_star::gyration_radius_sq() {

  return cnsts.parameters(convective_star_gyration_radius_sq); 
}
