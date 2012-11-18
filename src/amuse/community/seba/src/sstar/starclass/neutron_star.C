#define DEBUG false
//
// neutron_star.C
//

#include "neutron_star.h"
#include "hyper_giant.h"
#include "super_giant.h"
#include "thorne_zytkow.h"
#include "helium_giant.h"


// ANSI C++ first creates the base class before the dreived classes are
// created. If this is done in another order we are in deep shit.

neutron_star::neutron_star(super_giant & g) : single_star(g) {

      delete &g;

      magnetic_field  = cnsts.parameters(pulsar_magnetic_field);
      // log Gauss
      rotation_period = cnsts.parameters(pulsar_pulse_period);
      // seconde


      suddenly_lost_mass     = 0;
      real m_tot             = get_total_mass();
      core_mass = birth_mass = neutron_star_mass(Super_Giant);
      envelope_mass          = m_tot - core_mass;
      relative_age           = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      bool hit_companion = super_nova();
      post_supernova_story();

      refresh_memory();
      instantaneous_element();
      update();

      post_constructor();

      if (hit_companion)            // Kick may cause coalescence.
         direct_hit();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }
}

neutron_star::neutron_star(hyper_giant & w) : single_star(w) {

      delete &w;

      magnetic_field  = cnsts.parameters(pulsar_magnetic_field);
      // log Gauss
      rotation_period = cnsts.parameters(pulsar_pulse_period);
      // seconde

      suddenly_lost_mass     = 0;
      real m_tot             = get_total_mass();
      core_mass = birth_mass = neutron_star_mass(Super_Giant);
      envelope_mass          = m_tot - core_mass;
      relative_age           = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      bool hit_companion = super_nova();
      post_supernova_story();

      refresh_memory();
      instantaneous_element();
      update();

      post_constructor();

      if (hit_companion) 
         direct_hit();

      post_constructor();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	dump("binev.data", false);
	get_companion()->set_spec_type(Accreting, false);
      }
      else {
	dump("binev.data", false);
      }
}

neutron_star::neutron_star(thorne_zytkow & t) : single_star(t) {

      delete &t;

      suddenly_lost_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      relative_age = 0;
      lose_envelope_decent();

      refresh_memory();
      instantaneous_element();

      post_constructor();

      if (is_binary_component()) 
	get_binary()->dump("binev.data", false);
      else 
	dump("binev.data", false);
}

neutron_star::neutron_star(helium_giant & h) : single_star(h) {

      delete &h;


      magnetic_field  = cnsts.parameters(pulsar_magnetic_field);
      // log Gauss
      rotation_period = cnsts.parameters(pulsar_pulse_period);
      // seconde


      suddenly_lost_mass     = 0;
      real m_tot             = get_total_mass();
      core_mass = birth_mass = neutron_star_mass(Helium_Star);
      envelope_mass          = m_tot - core_mass;
      relative_age           = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

    bool hit_companion = super_nova();
    post_supernova_story();

      refresh_memory();
      instantaneous_element();
      update();

      post_constructor();

      if (hit_companion)            // Kick may cause coalescence.
         direct_hit();

      if (is_binary_component()) {
	get_binary()->set_first_contact(false);
	get_companion()->set_spec_type(Accreting, false);
	get_binary()->dump("binev.data", false);
      }
      else {
	dump("binev.data", false);
      }
}


neutron_star::neutron_star(white_dwarf & w) : single_star(w) {

      delete &w;
      
      magnetic_field  = cnsts.parameters(pulsar_magnetic_field);
      // log Gauss
      rotation_period = cnsts.parameters(pulsar_pulse_period);
      // seconde
      
      suddenly_lost_mass = 0;

      real m_tot = get_total_mass();
      core_mass = birth_mass = max(0.0, m_tot - aic_binding_energy());
      envelope_mass = m_tot - core_mass;
      relative_age = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

        bool hit_companion = super_nova();
        post_supernova_story();
    

      refresh_memory();
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

void neutron_star::instantaneous_element() {
  
  next_update_age = relative_age + cnsts.safety(maximum_timestep);

  core_radius = effective_radius = radius = neutron_star_radius();
  luminosity = spindown_luminosity();

  rotation_period = pulsar_spin_down(0);
}

void neutron_star::evolve_element(const real end_time) {
      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      core_radius = radius = neutron_star_radius();
      luminosity = spindown_luminosity();

      accrete_from_envelope(dt);
      rotation_period = pulsar_spin_down(dt);

      if (core_mass>cnsts.parameters(maximum_neutron_star_mass)) {

	if (is_binary_component()) 
	  get_binary()->dump("binev.data", false);
	else
	  dump("binev.data", false);
	
         star_transformation_story(Black_Hole);
         new black_hole(*this);
         return;
      }

      next_update_age = relative_age + cnsts.safety(maximum_timestep);

      update();
   }

void neutron_star::update() {

      detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
      effective_radius = radius;

//      if(get_element_type() != previous.star_type) {
//
//	if (is_binary_component()) 
//	  get_binary()->dump("SeBa.data", true);
//	else
//	  dump("SeBa.data", true);
//      }
    }

void neutron_star::accrete_from_envelope(const real dt) {

  if (DEBUG)
    cerr<<"neutron_star::accrete_from_envelope "<<dt<<endl;
  
  // real dm = NS_CORE_ACCRETION*accretion_limit(envelope_mass, dt);
  // Generally accretion onto neutron star core proceeds at 5% Eddington
  // For the moment this is 100% Eddington.

  real dm = accretion_limit(envelope_mass, dt);

  if (dm>0) {
    if (!propeller(dm, dt)) {

      magnetic_field  = magnetic_field_decay(dm, magnetic_field, dt);
      rotation_period = pulsar_spin_up(dm, dt);
      luminosity = accretion_luminosity(dm, dt);

      set_spec_type(Accreting);
    }
    else {

      rotation_period = pulsar_propeller_torque(dm, dt);
      real tau = 100;
      magnetic_field  = magnetic_field_decay(magnetic_field, dt/tau);

      set_spec_type(Accreting, false);
      dm = 0;
    }
  }
  else {
    real tau = 100;
    magnetic_field = magnetic_field_decay(magnetic_field, dt/tau);

    set_spec_type(Accreting, false);
  }

  core_mass += dm;
  envelope_mass = 0;
}

// Neutron star's binading energy. Energy fraction is lost during a
// Accretion Induced colaescence of a white dwarf into a neutron star.
// Since whole coalescence is not certain this function only used in
// cases of allowing AIC to occur.
// Typically one of those things that makes programming interesting.
real neutron_star::aic_binding_energy() {

      real GM2_RC2 = 3*cnsts.physics(G)
	           * pow(cnsts.parameters(solar_mass)/cnsts.physics(C), 2)
	           / (5*cnsts.parameters(solar_radius));
      real binding_energy = (GM2_RC2*pow(core_mass,2)
			  / cnsts.parameters(kanonical_neutron_star_radius))
	                  / cnsts.parameters(solar_mass);

      return binding_energy;
   }

// One of the few routines that requires the random number generator.
// random numbers:
//	kick_velocity, 
//	theta (kick-angel in orbital plane),
//	phi   (kick-angel perpendicular to orbital plane).
// separation is determined by solving Keplers equation for the 
// Eccentric anomaly taking the mean anomaly randomly.
//
// Return value is boolian whather or not the remnant colescence 
// with its companion.
bool neutron_star::super_nova() {

     suddenly_lost_mass = envelope_mass;

     bool hit_companion = FALSE;
      
     real v_kick     = cnsts.super_nova_kick(); //random_kick();
     real theta_kick = acos(1-2*random_angle(0, 1));
     real phi_kick   = random_angle(0, 2*PI);

     cerr << "Supernova kick v = " << v_kick << " [km/s]" << endl;

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

	real x_kick = v_kick*sin(theta_kick)*cos(phi_kick);
	real y_kick = v_kick*sin(theta_kick)*sin(phi_kick);
	real z_kick = v_kick*cos(theta_kick);
	vec kick_velocity(x_kick, y_kick, z_kick);
	anomal_velocity = kick_velocity;

	// This is correct if velocity is along the z-axis.
	velocity = sqrt(pow(v_kick, 2) + pow(velocity, 2)
		 + 2*v_kick*velocity*sin(theta_kick)*cos(phi_kick));

	envelope_mass = 0;
	return hit_companion;
      }
      else if(is_binary_component()) { // +++

	// Supernova is performed by the binary evolution

	if(get_binary()->get_bin_type() == Disrupted ||
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

            if (pericenter <= get_companion()->get_radius())
	      hit_companion = TRUE;
	    
	  }
	  else {
            set_spec_type(Runaway);

            get_binary()->set_eccentricity(1);
            get_companion()->set_spec_type(Runaway);
            get_binary()->set_bin_type(Disrupted);
            get_binary()->set_semi(0);

            real e2_init = e_init*e_init;
            real vr_mean_0 = sqrt(((cnsts.physics(G)
			   *        cnsts.parameters(solar_mass)
	                   /        cnsts.parameters(solar_radius))
			   *        (msn_0+m_comp_0)/a_init)
                           *        (1-e2_init)/pow(1+0.5*e2_init, 2));
	    vr_mean_0      /= cnsts.physics(km_per_s);
	    
            real mu_0 = get_total_mass()/mtot_0;
            real v_comp = mu_0*vr_mean_0;

	    //              velocity at infinity.
            real v_sn_0 = (1-mu_0)*vr_mean_0;
            real v_sn   = sqrt(v_sn_0*v_sn_0 + v_kick*v_kick
                        + 2*v_sn_0*v_kick*sin(theta_kick)*cos(phi_kick));

	    //              binary CoM velocity:
            real v_cm = get_binary()->get_velocity();
            v_comp = sqrt(pow(v_comp, 2) + pow(v_cm, 2)
                   + 2*v_comp*v_cm*cos(theta_kick));
            get_companion()->set_velocity(v_comp);
	    
            v_sn = sqrt(pow(v_sn, 2) + pow(v_cm, 2)
                 + 2*v_sn*v_cm*cos(theta_kick));
            velocity = v_sn;
         }
      }
      
      } // +++
      envelope_mass = 0;

      
     return hit_companion;
}

//used by subtrac_mass_from_donor and double_star::perform_mass_transfer

real neutron_star::mdot_limit(const real dt){
 
    return mass_ratio_mdot_limit(accretion_limit(envelope_mass, dt));
}



star* neutron_star::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = mdot_limit(dt);

      //		The disc may lose mass.
      envelope_mass -= mdot;

      return this;
}


// Accrete on neutron star.
// Eddington limited. The accreted mass first
// enters a stallar envelope (disc) and is finally accreted on
// the compact object.
// So no magnetic field decay here.
real neutron_star::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {
    //For neutron stars no difference currently between hydrogen/helium/.. accretion

     if (DEBUG) cerr<<"neutron_star::add_mass_to_accretor "<<dt<<endl;

     mdot = accretion_limit(mdot, dt);
     envelope_mass += mdot;
     relative_mass = max(relative_mass, get_total_mass());

     return mdot;
}

// Limit the accretion onto the nuetron star by the Eddingtopn limit.
real neutron_star::accretion_limit(const real mdot, const real dt) {

  if (dt < 0) return mdot;

     real eddington = eddington_limit(radius, dt); 

   if (cnsts.parameters(hyper_critical))
      return min(mdot, (1.e+8)*eddington); 

   return min(mdot, eddington);

}


star* neutron_star::merge_elements(star* str) {

      envelope_mass = 0;
      real merger_core = str->get_core_mass();

      add_mass_to_accretor(str->get_envelope_mass(), str->hydrogen_envelope_star());
      relative_mass = max(relative_mass, get_total_mass() + merger_core);
      core_mass += merger_core;

      set_spec_type(Merger);
      instantaneous_element();

      return this;

   }

star* neutron_star::reduce_mass(const real mdot) {

      if (mdot>envelope_mass) {
	real dm = mdot -envelope_mass;
	envelope_mass = 0;
	core_mass -= dm;
      }
      else envelope_mass -= mdot;
      
      return this;
}

void neutron_star::direct_hit() {

      real theta = random_angle(0, 2*PI);
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

      get_binary()->set_bin_type(Merged);
      get_binary()->set_semi(0);
      get_binary()->set_velocity(v_new);
      get_binary()->get_primary()->set_velocity(v_new);

      get_binary()->dump("hit.data", false);
   }

real neutron_star::sudden_mass_loss() {

    real mass_lost = suddenly_lost_mass;
    suddenly_lost_mass = 0;

    return mass_lost;

   }


/*----------------------------------------------------------------------
 * Relation concerning the radio pulsar's evolution
 *----------------------------------------------------------------------
 */
//      See Timmes, Woosley, Weaver, 1997, ApJ 457, 834
//      Bimodal distribution for supernova type II progenitors and
//      a single Gaussian for type I's
//      For the moment a random function determines the mass of the remnant.
//      If the mass appears to be above 2.0 Msun, evolve element
//      transformes it into a black hole.
//      Note that a kick is already applied in these cases.
//      Which might as well be true: see Sigurdsson & Brandt (~1995).

real neutron_star::neutron_star_mass(stellar_type stp) {

  real m = get_total_mass();
  if (stp==Helium_Star || stp==Hyper_Giant)
        m = 10 + get_total_mass();
  else
      m = relative_mass;

    real a, b, c, d, e;
    a =  2.25928E+00;
    b = -2.68264E-01;
    c =  2.28206E-02;
    d = -7.06583E-04;
    e =  7.48965E-06;
    real mass = a + m*(b + m*(c + m*(d + m*e)));

  return min(mass, get_total_mass());
}

//Assuming a Newtonian Polytrope with n= 3/2
real neutron_star::neutron_star_radius() {

     return 15.12 * (cnsts.physics(kilometer)/cnsts.parameters(solar_radius))
                  / pow(core_mass, cnsts.mathematics(one_third));   // [rsun]
}

//Is this neutron star spinning so fast that it cannot accrete any matter?
bool neutron_star::propeller(const real mdot, const real dt) {
     if (DEBUG)
       cerr<<"neutron_star::propeller"<<endl;

     real propellor = TRUE;
     real eddington = eddington_limit(radius, dt); 

     //Take the propeller effect into account, by allowing accretion
     //only to occur if mdot>propeller_limit
     if (is_binary_component() && eddington>0 && mdot>0 && dt>0) {
         real f = 1;

         real v_wind = get_companion()->wind_velocity();
         real v_orb  = velocity + get_companion()->get_velocity();
//       real v2_rel  = sqrt(0.5 * (pow(v_wind, 2) + pow(v_orb, 2)));
         real v2_rel  = sqrt(pow(v_wind, 2) + pow(v_orb, 2));
         real propeller_limit = 2.06e-8*f
                              * pow(pow(10., magnetic_moment()-30), 2)
            / (v2_rel*pow(rotation_period, 4));               // [Msun/Myr]

         // no accterion, neutron star spins too rapidly.
         if (mdot>propeller_limit*dt)
	              propellor=FALSE;
     }

     if (DEBUG)
       cerr<<propellor<<endl;

    return propellor;
}

real neutron_star::moment_of_inertia() {
     if (DEBUG)
       cerr<<"neutron_star::moment_of_inertia = "<<endl;

     real MR2 = cnsts.parameters(solar_mass)
              * pow(cnsts.parameters(solar_radius), 2);
     real inertia = MR2*gyration_radius_sq()*core_mass*pow(radius, 2);

     if (DEBUG)
       cerr << inertia << endl;

  return inertia;
}

real neutron_star::period_derivative() {
    if (DEBUG)
      cerr<<"neutron_star::period_derivative"<<endl;

  real pi2R6_3c3 = cnsts.mathematics(one_third)*pow(cnsts.mathematics(pi), 2)
                 * pow( pow(cnsts.parameters(solar_radius), 2)
	         /          cnsts.physics(C), 3);

  real Bfield = pow(10., magnetic_field);
  real inertia = moment_of_inertia();

  real Pdot = 8*pi2R6_3c3*pow(radius, 6)*pow(Bfield, 2)
                / (inertia*rotation_period);


    if (DEBUG)
      cerr << Pdot << endl;

  return Pdot;
}

real neutron_star::spindown_luminosity() {
     if (DEBUG)
       cerr<<"neutron_star::spindown_luminosity"<<endl;

     real L400 = 1;         // Luminosity at 400 Mhz.
     if (rotation_period>0) {
       real Pdot = period_derivative();

       real argument = Pdot/pow(rotation_period, 3);

       if(argument>0)
	 L400 = pow(10., 6.63 + ONE_THIRD*log10(argument));

     }
     
     if (DEBUG)
       cerr<<L400<<endl;

     return L400;
}

real neutron_star::accretion_luminosity(const real dm,
					const real dt) {
    if (DEBUG)
      cerr<<"neutron_star::accretion_luminosity"<<endl;

    real L = spindown_luminosity();

    if (dm>0 &&  dt>0) {
      real dmdt = dm/dt;
      real GM2_RMyr = cnsts.physics(G)*pow(cnsts.parameters(solar_mass), 2)
                    / (cnsts.parameters(solar_radius)*cnsts.physics(Myear));
      L = GM2_RMyr*core_mass*dmdt/radius;
      L /= cnsts.parameters(solar_luminosity);
    }
    
    if (DEBUG)
      cerr<<L<<endl;

    return L;
}

real neutron_star::fastness_parameter(const real dmdt) {
    if (DEBUG)
      cerr<<"neutron_star::fastness_parameter"<<endl;

    real mu30 = pow(10., magnetic_moment()-30);
    real omega_s = 1.19*pow(mu30, 6./7.)
                 / (rotation_period*pow(dmdt, 3./7.)*pow(core_mass, 5./7.));

    if (DEBUG)
      cerr<< min(1-0.01, max(0.1, omega_s))<<endl;

    return min(1-0.01, max(0.1, omega_s));
}

real neutron_star::dimension_less_accretion_torque(const real dmdt) {
    if (DEBUG)
      cerr<<"neutron_star::dimension_less_accretion_torque"<<endl;
    
    real omega_s = fastness_parameter(dmdt);
    real n_omega = ((7./6.) - (4./3.)*omega_s + (1./9)*pow(omega_s, 2))
                 / (1-omega_s);
    
    if (DEBUG)
      cerr<<n_omega<<endl;
    
    return n_omega;
}


real neutron_star::pulsar_propeller_torque(const real dm,
                                           const real dt) {
     if (DEBUG)
       cerr<<"neutron_star::pulsar_propeller_torque"<<endl;

//Pulsar torque due to drag should be taken into account.
// the prescription below is for accreting pulsars, not propellering ones.
    
    real dmEdd = 1;
    real dmdt=1;
    if (dm>0 &&  dt>0) {
      dmdt = dm/dt;
      real Eddington = eddington_limit(radius, 1.0); 
      //real Eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius; // [Msun/Myr]
      dmEdd = dmdt/Eddington;
    }
    
    real n_omega = dimension_less_accretion_torque(dmdt);
    real mu30 = pow(10., magnetic_moment()-30);

    real Pdot = -7.2 * n_omega * pow(mu30, 2./7.)
              * pow(dmEdd, 6./7.) * pow(rotation_period, 2); // [s/Myr]
    Pdot = max(0., Pdot);

    real delta_P = rotation_period + Pdot*dt;
    
    if (DEBUG)
      cerr<<delta_P<<endl;

    return delta_P;
}

real neutron_star::magnetic_moment() {
  
     real r3 = pow(radius*cnsts.parameters(solar_radius), 3);
               // log r^3 in [cm]

     return magnetic_field + log10(r3);
}

// Spontanious magnetic field decay.
// the constant for field decay is tau=100 Myr.
// The parameter dt_tau is actually dt/tau. Taken as a single parameter
// to make function overloading workable.
real neutron_star::magnetic_field_decay(const real log_B,
					const real dt_tau=0) {
  
     real B0 = pow(10., log_B);
     real B1 = max(1.,  B0 / exp(dt_tau));

     return log10(B1);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// local neutron star functions.
// Currently only for an accretion rate of 10^-9 Msun/yr.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real neutron_star::magnetic_field_decay(const real delta_m,
                                        const real log_B,
					const real dt=0) {
     real dmEdd = 1.0;
     if (delta_m>0 &&  dt>0) {
       real dmdt = delta_m/dt;
       real Eddington = eddington_limit(radius, 1.0); 

       dmEdd = dmdt/Eddington;
     }
     
     real B1 = magnetic_field_strength(core_mass, dmEdd);
     real B2 = magnetic_field_strength(core_mass+delta_m, dmEdd);
     real delta_B = min(0.0, B2-B1);

     return log_B + delta_B;
}

//The lineair interpolation is from :
// Konar, S., 1997, PhD thesis, Indian Institute of Science,
// Bangalorem India.
// The lineair interpolated curve presented is for Mdor = 10^-9 Msun/yr.
// For lower mass loss rates a linair shift is applied of
// log dM = 0.5 and log dB = 0.3 per order of magnitude decrease
// in accretion rate.
real neutron_star::magnetic_field_strength(const real ns_mass,
					   const real dmEdd = 1) {

      real delta_BBo = -0.2 * log10(dmEdd);
      real log_dm    = log10(ns_mass - birth_mass + 1.0e-10);
      log_dm        -= max(0.0, -0.5*(1+log10(dmEdd)));

      real log_dm_min  = -8.0;
      real log_dm_max  = -2;
      real log_BBo_min = -5.0;
      real log_BBo_max =  0.0;

      if (log_dm>-3)
	        return max(0.0, log_BBo_min + delta_BBo);
      if (log_dm<-7.5)
	        return log_BBo_max;

      if(log_dm>=-3.5) {
          log_dm_min  = -3.5;
          log_dm_max  = -3.0;
          log_BBo_min = -4.5;
          log_BBo_max = -4.9;
      }else if(log_dm>=-4.0) {
          log_dm_min  = -4.0;
          log_dm_max  = -3.5;
          log_BBo_min = -3.1;
          log_BBo_max = -4.5;
      }else if(log_dm>=-4.5) {
          log_dm_min  = -4.5;
          log_dm_max  = -4.0;
          log_BBo_min = -2.4;
          log_BBo_max = -3.1;
      }else if(log_dm>=-5.0) {
          log_dm_min  = -5.0;
          log_dm_max  = -4.5;
          log_BBo_min = -2.0;
          log_BBo_max = -2.4;
      }else if(log_dm>=-5.5) {
          log_dm_min  = -5.5;
          log_dm_max  = -5.0;
          log_BBo_min = -1.5;
          log_BBo_max = -2.0;
      }else if(log_dm>=-6.0) {
          log_dm_min  = -6.0;
          log_dm_max  = -5.5;
          log_BBo_min = -1.1;
          log_BBo_max = -1.5;
      }else if(log_dm>=-6.5) {
          log_dm_min  = -6.5;
          log_dm_max  = -6.0;
          log_BBo_min = -0.7;
          log_BBo_max = -1.1;
      }else if(log_dm>=-7.0) {
          log_dm_min  = -7.0;
          log_dm_max  = -6.5;
          log_BBo_min = -0.4;
          log_BBo_max = -0.7;
      }else if(log_dm>=-7.5) {
          log_dm_min  = -7.5;
          log_dm_max  = -7.0;
          log_BBo_min =  0.0;
          log_BBo_max = -0.4;
}


      real log_B = lineair_interpolation(log_dm, log_dm_min, log_dm_max,
                                               log_BBo_min, log_BBo_max);
      log_B += delta_BBo;

      return min(0.0, log_B);
}

real neutron_star::pulsar_spin_up(const real dm, const real dt) {

     if (!is_binary_component())
       return rotation_period;

     real dmEdd = 1.0;
     if (dm>0 &&  dt>0) {
       real dmdt = dm/dt;
       real Eddington = eddington_limit(radius, 1.0); 
       dmEdd = dmdt/Eddington;
     }
     
     real B = pow(10., magnetic_field);
     real P_eq = 0.0019*pow(B/1.e+9, 6./7.)
               / (pow(core_mass/1.4, 5./7.)*pow(dmEdd, 3./7.));  // seconds

     return min(rotation_period, P_eq);
}

real neutron_star::pulsar_spin_down(const real dt) {
   if (DEBUG)
     cerr << "enter pulsar_spin_down( dt="<<dt<<")"<<endl;

   real B = pow(10., magnetic_field);
   real P0 = rotation_period;

   real Pdot = pow(B/3.2e+19, 2)/P0;
   real delta_P = Pdot * dt * cnsts.physics(Myear);

   return rotation_period + delta_P;
}

bool neutron_star::dead_pulsar() {

   bool dead = FALSE;

   real B = pow(10., magnetic_field);
   real Pdeath = sqrt(B/0.17e+12);

   if (rotation_period >=Pdeath)
          dead = TRUE;

   return dead;
}

real neutron_star::gyration_radius_sq() {

  // Gunn & Ostriker 1969, Nat. 221, 454
  // Momentum of inertia of 12km neutron star = 10^45 gr cm^2
  // this results in:
  
  return 0.25; 
  
}

stellar_type neutron_star::get_element_type() {

  // (SPZ+GN:  1 Aug 2000)
  // Currently not stable due to SeBa.data which would produce too much output.
  // Fix: implement next_update_age in neutron stars.C
#if 0
  if (spec_type[Accreting])
    return Xray_Pulsar;
  else if (!dead_pulsar())
    return Radio_Pulsar;
  else
    return Neutron_Star;
#endif

    return Neutron_Star;

}




real neutron_star::get_evolve_timestep() {
    
    return max(next_update_age - relative_age, 0.0001);
}

