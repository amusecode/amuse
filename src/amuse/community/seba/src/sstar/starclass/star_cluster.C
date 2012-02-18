//
// star_cluster.C
//
// Rudimentary class for simulating the evolution of a cluster of stars.
//

#include "star_cluster.h"

void star_cluster::instantaneous_element() {

  luminosity = 1;
  effective_radius = radius = 1;
  core_radius = 1;

  envelope_mass = get_total_mass();
  core_mass = 0;

  radius = core_radius*Rtidal_over_Rvir_KingW09;

  initialize_star_cluster();

}

void star_cluster::print(ostream& s = cerr) {

  PRC(relative_mass);PRC(get_total_mass());PRL(nstar);
  PRC(relative_age);PRC(next_update_age);PRC(relative_mass);
  PRC(black_hole_mass);PRC(trlx_initial);
  PRC(rvir);PRC(nstar);PRC(x_imf);PRL(minimum_mass);

}

void star_cluster::evolve_element(const real end_time) {

  //        cerr << "XXX  Evolve Star cluster" << endl;
	//	print();

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

	adjust_next_update_age();

        update();

	stellar_wind(dt);
	//	cerr << "Leave evolve_element()" << endl;
	//	print();
}


void star_cluster::update() {

  // do nothing for now.
  // Just override single_star::update();
}

real star_cluster::turn_off_mass(const real t) {
// ST 16 feb 2009 simply added solar metallicity to make
//function compatible with new code 
// is not based on anything yet! 
// 
 real metalicity = cnsts.parameters(solar_metalicity);
  real m = 50;  // first guess;
  real mdot = m;
  real to_time = main_sequence_time(m, metalicity);

  while(abs(t-to_time)>1.e-2) {
    mdot *= 0.5;
    if(to_time>t)
      m += mdot;
    else
      m -= mdot;
    to_time = main_sequence_time(m, metalicity);

    if(m>=100 || m<0.1)
      to_time = t;
  }

  return m;
}

// From McMillan & Portegies Zwart Eq. 30
real star_cluster::mean_stellar_mass(const real t, const real mmin) {

  real mto = turn_off_mass(t);
  real mm = ((1-x_imf)/(2-x_imf))
    * (pow(mto, 2-x_imf) - pow(mmin, 2-x_imf))
    / (pow(mto, 1-x_imf) - pow(mmin, 1-x_imf));

  //  PRC(mto);PRC(x_imf);PRC(mmin);PRL(mm);

  return mm;
}

real star_cluster::mean_stellar_mass() {
  return mean_stellar_mass(relative_age, minimum_mass);
}

real star_cluster::relaxation_time(real n, real m, real r) {

  real t_hc = 57 * sqrt(pow(r, 3)/m);
  real t_rlx = 2.05 * sqrt(pow(r, 3)/m) * n/log(0.11*n);

  return t_rlx;
}

void star_cluster::initialize_star_cluster() {

  //  cerr << "void star_cluster::initialize_star_cluster()" << endl;

  black_hole_mass = 0;
  rvir = 1.0;
  x_imf = 2.35;  // Salpeter
  minimum_mass = 1.0;
  //minimum_mass = 0.63;
  //minimum_mass = 0.50;
  //minimum_mass = 0.1;
  nstar = get_total_mass()/mean_stellar_mass();
  trlx_initial = relaxation_time(nstar, get_total_mass(), 
				     tidal_radius());

  evolve_element(relative_age);
}

real star_cluster::mass_loss_by_evolution(const real dt) {

  /// correcttion factor for emulating turn off mass to nuclear lifetime
  //  real age_correction_factor = 0.85;
  real age_correction_factor = 0.70;

  real mm_old = mean_stellar_mass(age_correction_factor*relative_age, 
				  minimum_mass);
  real mm_new = mean_stellar_mass(age_correction_factor*(relative_age + dt), 
				  minimum_mass);
  real mass_lost = nstar * (mm_old-mm_new);

  //  PRC(dt);PRC(mm_old);PRC(mm_new);PRL(mass_lost);

  return mass_lost;
}


void star_cluster::stellar_wind(const real dt) {

    if (!get_use_hdyn()) 
	cerr << " No stellar dynamical information present in star_cluster" 
	     << endl;

    real wind_mass = mass_loss_by_evolution(dt);
    //    PRC(wind_mass);

  if (wind_mass >= envelope_mass) 
    wind_mass = envelope_mass;

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_mass(wind_mass);

  return;
}

real star_cluster::tidal_radius() {
  return rvir*Rtidal_over_Rvir_KingW09;
}

real star_cluster::accretion_limit(const real mdot, const real dt) {

        return 0;
     }

star* star_cluster::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = 0;
      return this;
}

void star_cluster::adjust_accretor_age(const real mdot,
				      const bool rejuvenate=true) {

  return;

}

real star_cluster::zeta_thermal() {
        return 0;
     }

star* star_cluster::reduce_mass(const real mdot) {

  if (envelope_mass>=mdot) {
    envelope_mass -= mdot;

  } 
  else {
    cerr << "WARNING: star* star_cluster::reduce_mass(const real mdot=" 
	 << mdot<< ")"<< endl; 
    cerr << "star_cluster has negative mass..." << endl;
    envelope_mass = 0;
  }

  // Note that the core_mass = 0 at initialization.
  if (envelope_mass<=0) {
      cerr << "star_cluster reduced to zero mass at t="
	   << get_current_time() << endl;
      cerr << "            last mdot was: " << mdot << endl;
  }
           
  return this;
}

void star_cluster::adjust_next_update_age() {

  //  PRC(relative_age);PRC(last_update_age);PRL(next_update_age);

  last_update_age = next_update_age;
  next_update_age = relative_age + 1;
}

real star_cluster::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}


void star_cluster::update_wind_constant(const real tscale = 0) {

  wind_constant = tscale;

}

