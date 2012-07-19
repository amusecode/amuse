//
// black_hole.C
//

#include "black_hole.h"
#include "hyper_giant.h"
#include "super_giant.h"
#include "thorne_zytkow.h"
#include "helium_giant.h"
#include "neutron_star.h"

// ANSI C++ first creates the base class before the dreived classes are
// created. If this is done in another order we are in deep shit.

black_hole::black_hole(hyper_giant & w) : single_star(w) {

      delete &w;

      suddenly_lost_mass = 0;
      real m_tot = get_total_mass();

      core_mass = black_hole_mass();
      envelope_mass = m_tot - core_mass;
      
// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      relative_age = 0;

      bool hit_companion = super_nova();
      post_supernova_story(); 
      
      instantaneous_element();
      update();
    
    if (hit_companion)
         direct_hit();

      post_constructor();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }
}

black_hole::black_hole(super_giant & g) : single_star(g) {

      delete &g;

      suddenly_lost_mass = 0;
      real m_tot = get_total_mass();
      core_mass = black_hole_mass();
      envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      relative_age = 0;

    bool hit_companion = false;
 //     bool hit_companion = super_nova();
 //     post_supernova_story(); 

      instantaneous_element();
      update();
    
      if (hit_companion)
          direct_hit();
 
      post_constructor();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }
}

black_hole::black_hole(thorne_zytkow & t) : single_star(t) {

      delete &t;

      magnetic_field  = 0;
      rotation_period = 0;

      //  TZO core is the black-hole itself!
      suddenly_lost_mass = 0;
      relative_age = 0;

      bool hit_companion = super_nova();
      post_supernova_story(); 

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      instantaneous_element();
      update();

      if (hit_companion)
         direct_hit();

      post_constructor();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }

}

black_hole::black_hole(helium_giant & h) : single_star(h) {
      delete &h;

      suddenly_lost_mass = 0;
      real m_tot = get_total_mass();
      core_mass = black_hole_mass();
      envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      relative_age = 0;

      bool hit_companion = super_nova();
      post_supernova_story(); 

      instantaneous_element();
      update();

      if (hit_companion)
         direct_hit();

      post_constructor();
      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }
}

black_hole::black_hole(neutron_star & n) : single_star(n) {

      delete &n;

      magnetic_field  = 0;
      rotation_period = 0;

      suddenly_lost_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      relative_age = 0;

      bool hit_companion = super_nova();
      post_supernova_story(); 

      instantaneous_element();

      if (hit_companion)
         direct_hit();

      post_constructor();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }

}

#if 0
// Initial black_hole should be properly initialized.
void black_hole::adjust_initial_star() {

  if (core_mass<=0) {
    core_mass = get_total_mass();
    envelope_mass = 0;
  }
  
  if (relative_age<=0) 
    relative_age = max(current_time, 0.0);
  effective_radius = core_radius = radius
                   = cnsts.parameters(kanonical_neutron_star_radius);

  if (relative_mass<cnsts.parameters(super_giant2black_hole)) {
    cerr<<"Black hole with low mass!"
	<<"\n should be turned into a neutron star"
	<<"\nadopt minimum BH relative mass"<<endl;
    relative_mass = cnsts.parameters(super_giant2black_hole);
    return;
  }
}
#endif

real black_hole::get_radius() {
    cerr<<"enter black_hole::get_radius"<<endl;
    real r_eff = radius;
    if (is_binary_component() && get_companion()->get_element_type() != Black_Hole) {
        real m_sec = get_companion()->get_total_mass();
        real r_sec = get_companion()->get_radius();
        real r_tide = r_sec * (1 + pow(get_total_mass()/m_sec,
                                       cnsts.mathematics(one_third)));
        r_eff = r_tide;

    }

    return r_eff;
}

void black_hole::instantaneous_element() {
      next_update_age = relative_age + cnsts.safety(maximum_timestep);

      magnetic_field  = 0;
      rotation_period = 0;
      
      luminosity = 0.01;

      core_radius = effective_radius = radius = 4.25e-6*core_mass;

      update();

   }

// Not much to evolve about a black hole.
// Fixed radius, no wind, no nothing really low budget.
void black_hole::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      next_update_age = relative_age + cnsts.safety(maximum_timestep);
      accrete_from_envelope(dt);

      luminosity = 0.01;

      core_radius = radius = 4.25e-6*core_mass;

      update();

   }

void black_hole::update() {

      detect_spectral_features();
      
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
      effective_radius = radius;
   }

// Binding energy for a black hole.
// When no mass loss during a super-nova (collapse from a neutron star to
// a blck hole) the binading energy is thought to be released in the 
// explosion.
real black_hole::aic_binding_energy() {

      real GM2_RC2 = 3*cnsts.physics(G)*pow(cnsts.parameters(solar_mass)
	           / cnsts.physics(C), 2)/(5*cnsts.parameters(solar_radius));
      return GM2_RC2*core_mass/(4.25e-6*cnsts.parameters(solar_mass));
   }

real black_hole::black_hole_mass() {

  // Relation used if helium2black_hole=15 and
  //                  super_giant2black_hole=40
  // (SPZ+GN:24 Sep 1998)
  // real m = (0.35*relative_mass - 12);

// (GN+SPZ May 12 1999) was nice try, but seemt to be a bit unrealistic
// Try simple MBH = frac * MHe with frac = 0.5

  real m = 0.5 * get_total_mass();

  return m;
#if 0
  real m = 2;
  if (core_mass >= cnsts.parameters(helium2black_hole))
    m = 0.6 * (core_mass-cnsts.parameters(helium2black_hole)) + 2;
  
    //  if (core_mass>=cnsts.parameters(helium2black_hole))
    //    m = 0.6 * (core_mass-cnsts.parameters(helium2black_hole)) + 2;
    //  else if (core_mass>cnsts.parameters(super_giant2neutron_star))
    //    m *= (cnsts.parameters(super_giant2neutron_star)/core_mass);
    //  else if (core_mass>cnsts.parameters(helium2neutron_star))
    //    m = (cnsts.parameters(helium2neutron_star)/core_mass);
    //  else m = 0.043*pow(core_mass,1.67);

  return min(m, core_mass);
#endif

}

// Allow black_hole star to accrete from a disc.
// Mass accreted by the hole is first put into a storage disk
// and on a slower rate accreted onto the hole.
// The accretion effeciency BH_CORE_ACCRETION is defined in
// constants.h and has a magnitude of a few percent.
void black_hole::accrete_from_envelope(const real dt) {

      real mdot = cnsts.parameters(black_hole_accretion_limit)
	        * accretion_limit(envelope_mass, dt);

      if (mdot>0) {

	set_spec_type(Accreting);
	core_mass += mdot;
	envelope_mass -= mdot;
      }
      else
	set_spec_type(Accreting, false);
}

real black_hole::add_mass_to_accretor(const real mdot) {

//		Increase envelope_mass of black hole.
      envelope_mass += mdot;
      relative_mass = max(relative_mass, get_total_mass());

      return mdot;
   }

real black_hole::add_mass_to_accretor(real mdot, const real dt) {

      mdot = accretion_limit(mdot, dt);

      envelope_mass += mdot;
      relative_mass = max(relative_mass, get_total_mass());

      return mdot;
   }

// Accretion is limited by the Eddington luminosity
real black_hole::accretion_limit(const real mdot, const real dt) {

      real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;
      if (cnsts.parameters(hyper_critical))
	eddington *= 1.e8;
      
      return min(eddington, mdot);
}

// A black hole is generally not a donor.
star* black_hole::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = accretion_limit(envelope_mass, dt);
      mdot = mass_ratio_mdot_limit(mdot);

      envelope_mass -= mdot;
      return this;
}


star* black_hole::merge_elements(star* str) {

      envelope_mass = 0;
      real merger_core = str->get_core_mass();
      real merger_env = str->get_envelope_mass();

      if(get_total_mass()<100) {
	  add_mass_to_accretor(merger_env);
	  core_mass += merger_core;
      }
      else {
	  core_mass += merger_core + merger_env;
      }
      
      relative_mass = max(relative_mass, get_total_mass());

      spec_type[Merger]=Merger;
      instantaneous_element();

      return this;
   }

// Tricky, especially in higher order system.
// One of the few routines that requires the random number generator.
// random numbers:
//      kick_velocity,
//      theta (kick-angel in orbital plane),
//      phi   (kick-angel perpendicular to orbital plane).
// separation is determined by solving Keplers equation for the
// Eccentric anomaly taking the mean anomaly randomly.
//
// Note that a black hole does not receive a kick!!
//
// Return value is boolian whather or not the remnant colescence
// with its companion.
bool black_hole::super_nova() {

      suddenly_lost_mass = envelope_mass;

      bool hit_companion = FALSE;
      //real mass_correction =
      //sqrt(cnsts.parameters(kanonical_neutron_star_mass)/core_mass);
      
      real v_kick  = cnsts.super_nova_kick(no_velocity_kick);
                     // mass_correction*random_kick();
      real theta_kick = acos(1-2*random_angle(0, 1));
      real phi_kick   = random_angle(0, cnsts.mathematics(two_pi));

//      cerr << "Supernova kick v = " << v_kick << " [km/s]" << endl;

     if (is_binary_component()) {
       if(v_kick>0) 
	 get_seba_counters()->snwk_in_dstar++;
       else
	 get_seba_counters()->sn_in_dstar++;
     }
     else if(v_kick>0) 
       get_seba_counters()->snwk_in_sstar++;
     else
       get_seba_counters()->sn_in_sstar++;

    if (get_use_hdyn()) {

      // Transform 1D kick velocity to 3D space velocity for dynamics.
      real x_kick = v_kick*sin(theta_kick)*cos(phi_kick);
      real y_kick = v_kick*sin(theta_kick)*sin(phi_kick);
      real z_kick = v_kick*cos(theta_kick);
      vec kick_velocity(x_kick, y_kick, z_kick);
//      PRL(kick_velocity);
      anomal_velocity = kick_velocity;
// cerr << "Kick   : "<< v_kick <<" " << theta_kick <<" " << phi_kick<< endl;
// cerr << "kick v : " << anomal_velocity << endl;

      velocity = sqrt(pow(v_kick, 2) + pow(velocity, 2)
	       + 2*v_kick*velocity*sin(theta_kick)*cos(phi_kick));
      envelope_mass = 0;
//      cerr << "return hit_companion;" <<hit_companion << flush << endl;
      return hit_companion;
    }
    
    else if(get_binary()->get_bin_type() == Disrupted ||
	    get_binary()->get_bin_type() == Merged) {

      get_companion()->set_spec_type(Accreting, false);

      velocity = sqrt(pow(v_kick, 2) + pow(velocity, 2)
	       + 2*v_kick*velocity*sin(theta_kick)*cos(phi_kick));
      envelope_mass = 0;
      return hit_companion;
    }
    else if (get_binary()->get_bin_type() != Merged &&
	     get_binary()->get_bin_type() != Disrupted) {

      real a_init = get_binary()->get_semi();
      real e_init = get_binary()->get_eccentricity();
      real mtot_0 = get_binary()->get_total_mass();
      real msn_0 = get_total_mass();
      real m_psn = core_mass;
      real m_comp_0 = mtot_0 - msn_0;
      real m_comp = m_comp_0;

      real separation = random_separation(a_init, e_init);
      real a_new = post_sn_semi_major_axis(a_init, e_init, separation,
					   msn_0, m_comp_0, m_psn, m_comp,
					   v_kick, theta_kick, phi_kick);
      real e_new = post_sn_eccentricity(a_init, e_init, separation,
					msn_0, m_comp_0, m_psn, m_comp,
					v_kick, theta_kick, phi_kick);
      real v_cm  = post_sn_cm_velocity(a_init, e_init, separation,
				       msn_0, m_comp_0, m_psn, m_comp,
				       v_kick, theta_kick, phi_kick);

//              System bound after kick?
      if (e_new>=0 && e_new<1.) {
	get_binary()->set_eccentricity(e_new);
	get_binary()->set_semi(a_new);

	get_binary()->set_velocity(v_cm);
	get_binary()->calculate_velocities();
	    
//              Does it hit the companion?
	real pericenter = a_new*(1-e_new);
	if (pericenter<=get_companion()->get_radius())
	  hit_companion = TRUE;
      }
      else {
	spec_type[Runaway]=Runaway;
	
	get_binary()->set_eccentricity(1);
	get_companion()->set_spec_type(Runaway);
	get_binary()->set_bin_type(Disrupted);
	get_binary()->set_semi(0);
	real e2_init = e_init*e_init;
	real vr_mean_0 = sqrt(((cnsts.physics(G)*cnsts.parameters(solar_mass)
		       / cnsts.parameters(solar_radius))
		       * (msn_0+m_comp_0)/a_init)
                       * (1-e2_init)/pow(1+0.5*e2_init, 2));
	vr_mean_0 /= cnsts.physics(km_per_s);
	real mu_0 = get_total_mass()/mtot_0;
	real v_comp = mu_0*vr_mean_0;
//      v_comp = velocity_at_infinity(v_comp, separation, m_comp, m_psn);

	real v_sn_0 = (1-mu_0)*vr_mean_0;
	real v_sn   = sqrt(v_sn_0*v_sn_0 + v_kick*v_kick
		     + 2*v_sn_0*v_kick*sin(theta_kick)*cos(phi_kick));
//      v_sn = velocity_at_infinity(v_sn, separation, m_psn, m_comp);

//              Now corect for binary CM velocity:
	real v_cm = get_binary()->get_velocity();
	v_comp = sqrt(pow(v_comp, 2) + pow(v_cm, 2)
		      + 2*v_comp*v_cm*cos(theta_kick));
	get_companion()->set_velocity(v_comp);
	v_sn = sqrt(pow(v_sn, 2) + pow(v_cm, 2)
		    + 2*v_sn*v_cm*cos(theta_kick));
	velocity = v_sn;
         }
      }
    envelope_mass = 0;

    return hit_companion;
}

// Reduce the mass of the donor when it's filling its Roche-lobe.
// For a black_hole this is, of cource, not appropriate.
star* black_hole::reduce_mass(const real mdot) {

      if (mdot>envelope_mass) 
         envelope_mass = 0;
      else 
         envelope_mass -= mdot;
      return this;
   }

// Colide and merge just created neutron star with its companion.
// Call is rare but could happen is high velocity kicks are included.
void black_hole::direct_hit() {
//      Binary and triple stuff....

      real theta = random_angle(0, cnsts.mathematics(two_pi));
      real v_bin = get_binary()->get_velocity();
      real ek_ns = get_total_mass()*velocity*velocity;
      real ek_comp = get_companion()->get_total_mass()
                   * pow(get_companion()->get_velocity(), 2);
      real v_merger = sqrt((ek_ns+ek_comp)/get_binary()->get_total_mass());
      real v_new = sqrt(pow(v_merger, 2) + pow(v_bin, 2)
                 + 2*v_merger*v_bin*cos(theta));

      get_binary()->set_semi(0);
      if (get_total_mass() >= get_companion()->get_total_mass())
         get_binary()->merge_elements(this, get_companion());
      else
         get_binary()->merge_elements(get_companion(), this);

      get_binary()->set_semi(0);
      get_binary()->set_velocity(v_new);
      get_binary()->get_primary()->set_velocity(v_new);

      get_binary()->dump("hit.data", false);
   }

real black_hole::sudden_mass_loss() {

    real mass_lost = suddenly_lost_mass;
    suddenly_lost_mass = 0;

    return mass_lost;

   }


real black_hole::gyration_radius_sq() {

  real a = 1;
  return a/(cnsts.physics(C)*pow(radius, 2));
}

// Angular momentum of homogeneous sphere.
real black_hole::angular_momentum() {
  
  cerr << "black_hole::angular_momentum()"<<endl;
       
  real a = 1;   // Kerr black hole: maximum rotation
  real m = get_total_mass()*cnsts.parameters(solar_mass);

  return a*m;
	
}


