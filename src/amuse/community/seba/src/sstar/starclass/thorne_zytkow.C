//
// thorne_zytkow.C

#include "thorne_zytkow.h"
#include "main_sequence.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.

thorne_zytkow::thorne_zytkow(main_sequence & m) : single_star(m) {

      delete &m;

      // (SPZ: Removed 23March2000)
      //      birth_mass = 0;
      if (is_binary_component()) {

	birth_mass      = get_companion()->get_total_mass();
	if (get_companion()->remnant()) {
	  magnetic_field  = get_companion()->get_magnetic_field(); 
	  rotation_period = get_companion()->get_rotation_period();
	}
      }
      else if (birth_mass<=0) {
	birth_mass    = cnsts.parameters(kanonical_neutron_star_mass);
      }

      //(SPZ:23March2000)
      // This should have been done in main_sequence::merge_elementes
//      real m_tot    = max(birth_mass, get_total_mass());
//      core_mass     = birth_mass;
//      envelope_mass = m_tot - core_mass;
//      core_radius   = cnsts.parameters(kanonical_neutron_star_radius);

      wind_constant = envelope_mass;
      
// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      // Put star at start of super_giant stage.
      real t_ms = main_sequence_time();
      real t_giant = t_ms + hertzsprung_gap_time()
                          + base_giant_branch_time();
    cerr<<"helium_giant_time not in use anymore, find equivalent"<<endl;  
    //real t_he = helium_giant_time(t_ms, metalicity);
    real t_he = t_giant;
    
    
      relative_age = t_giant + t_he;
      adjust_next_update_age();

      instantaneous_element();
      update();

      post_constructor();
} 

#if 0
void thorne_zytkow::adjust_initial_star() {

  if (relative_age<=0) {
    real t_ms = main_sequence_time();
    real t_giant = t_ms + hertzsprung_gap_time(t_ms)
      + base_giant_branch_time(t_ms);
    real t_he = helium_giant_time(t_ms);
    relative_age = max(t_giant + t_he, relative_age);
  }
}
#endif

void thorne_zytkow::instantaneous_element() {

  real l_g = giant_luminosity();
  real t_ms = main_sequence_time();
  real t_gs = 0.15*t_ms;
  real t_b  = base_giant_time(t_ms);

  luminosity = l_g*pow(t_gs/(next_update_age
			     + t_b - relative_age), 1.17);
  luminosity = min(luminosity, maximum_luminosity());

  radius = (0.25*pow(luminosity, 0.4)
	 + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);
  radius = min(radius, 100.);

  // (SPZ+GN:  1 Aug 2000)
  // coming from previous type the effective readius should 
  // remain the same.
  //    effective_radius = radius;
  
}

// evolve a thorne_zytkow star upto time argument according to
// the super_giant model discribed by Eggleton et al. 1989.
void thorne_zytkow::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {

	// (neutron_star) core accretes from envelope.
         accrete_from_envelope(dt);
         if (core_mass>cnsts.parameters(maximum_neutron_star_mass)) {
	   if (is_binary_component()) 
	     get_binary()->dump("binev.data", false);
	   else
	     dump("binev.data", false);
	   
            star_transformation_story(Black_Hole);
            new black_hole(*this);
            return;
         }

         real l_g = giant_luminosity();
         real t_ms = main_sequence_time();
         real t_gs = 0.15*t_ms;
         real t_b  = base_giant_time(t_ms);

         luminosity = l_g*pow(t_gs/(next_update_age
                    + t_b - relative_age), 1.17);
         luminosity = min(luminosity, maximum_luminosity());
         radius = (0.25*pow(luminosity, 0.4)
                + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);

	 radius = min(radius, 100.);
      }
      else {

	if (is_binary_component()) 
	  get_binary()->dump("binev.data", false);
	else
	  dump("binev.data", false);
	
	 // thorne_zytkow lifetime exceeded.
         if (core_mass>cnsts.parameters(maximum_neutron_star_mass)) {
            star_transformation_story(Black_Hole);
            new black_hole(*this);
            return;
         }
         else {
            star_transformation_story(Neutron_Star);
            new neutron_star(*this);
            return;
         }
      }

      update();
      stellar_wind(dt);
}

void thorne_zytkow::update() {

      detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
      effective_radius = radius;

   }

void thorne_zytkow::accrete_from_envelope(const real dt) {

      real mdot_edd = 1.5e-08*cnsts.parameters(solar_radius)*core_radius*dt;
      if (cnsts.parameters(hyper_critical))
	  mdot_edd *= 1.e+8; 

      real mdot = min(mdot_edd, envelope_mass); 

      core_mass += mdot;
      envelope_mass -= mdot;
}

#if 0
void thorne_zytkow::stellar_wind(const real dt) {

      real wind_mass = 5.5e-7*dt*radius*luminosity/get_total_mass();
      wind_mass = min(wind_mass, envelope_mass);

      if (is_binary_component())
         get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
      else
         reduce_donor_mass(wind_mass);
   }

real thorne_zytkow::helium_core_mass() {

      real m_core = min(cnsts.parameters(kanonical_neutron_star_mass), get_total_mass());

      return m_core;
   }
#endif

//              Mass transfer utilities.
real thorne_zytkow::add_mass_to_accretor(const real mdot) {

      adjust_accretor_age(mdot);

      envelope_mass += mdot;
      if (relative_mass<get_total_mass()) {
	relative_mass = get_total_mass();
	adjust_next_update_age();
      }

      set_spec_type(Accreting);

      return mdot;
}

real thorne_zytkow::add_mass_to_accretor(real mdot, const real dt) {

      mdot = accretion_limit(mdot, dt);

      adjust_accretor_age(mdot);

      envelope_mass += mdot;
      if (relative_mass<get_total_mass()) {
	relative_mass = get_total_mass();
	adjust_next_update_age();
      }

      adjust_accretor_radius(mdot, dt);

      set_spec_type(Accreting);
      
      return mdot;
   }

star* thorne_zytkow::reduce_mass(const real mdot) {

    if (envelope_mass - mdot <=
	cnsts.parameters(neutron_star_accretion_limit)*core_mass) {
 
      envelope_mass = 0;

      if (is_binary_component()) 
	get_binary()->dump("binev.data", false);
      else
	dump("binev.data", false);

      // Envelope mass of neutron star cannot be supported by
      // gravitational pressure. Envelope is blown away by 
      // accretion luminosity of central compact object.
	  
      if (core_mass>cnsts.parameters(maximum_neutron_star_mass)) {
	star_transformation_story(Black_Hole);
	return dynamic_cast(star*, new black_hole(*this));
      }
      else {
	star_transformation_story(Neutron_Star);
	return dynamic_cast(star*, new neutron_star(*this));
      }
    }

    envelope_mass -= mdot;
    return this;
}

star* thorne_zytkow::subtrac_mass_from_donor(const real dt, real& mdot) {

      real mdot_max = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_max);

      if (envelope_mass<mdot) {

	mdot = envelope_mass;
	envelope_mass = 0;

	if (is_binary_component()) 
	  get_binary()->dump("binev.data", false);
	else
	  dump("binev.data", false);

	if (core_mass>cnsts.parameters(maximum_neutron_star_mass)) {
	  star_transformation_story(Black_Hole);
	  return dynamic_cast(star*, new black_hole(*this));
	}
	else {
	  star_transformation_story(Neutron_Star);
	  return dynamic_cast(star*, new neutron_star(*this));
	}
      }

      envelope_mass -= mdot;
      return this;
}

#if 0
void thorne_zytkow::adjust_accretor_age(const real mdot,
					const bool rejuvenate) {

      real m_tot_new = get_total_mass() + mdot;
      real m_rel_new = max(m_tot_new, relative_mass);

      real t_tagb_old = TAGB_time(relative_mass, metalicity);
      real t_du_old = dredge_up_time(relative_mass, metalicity);
      PRC(relative_mass);PRC(t_tagb_old);PRL(t_du_old);

      real z_new = get_metalicity();
      real t_tagb_new = TAGB_time(m_rel_new, z_new);
      real t_du_new = dredge_up_time(m_rel_new, z_new);
      PRC(m_rel_new);PRC(t_tagb_new);PRL(t_du_new);

      real dtime = relative_age - t_tagb_old;

      last_update_age = t_tagb_new;
      relative_age = t_tagb_new
                   + dtime*(t_du_new/t_du_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new);

       relative_age = max(relative_age, 
			  last_update_age + cnsts.safety(minimum_timestep));
      

      // next_update_age should not be reset here
      // next_update_age = t_nuc;
   }
#endif

void thorne_zytkow::adjust_next_update_age() {

      next_update_age = nucleair_evolution_time();
   }

#if 0
real thorne_zytkow::stellar_radius(const real mass, const real age) {

      real t_nuc = nucleair_evolution_time(mass);

      real l_g = giant_luminosity();
      real t_ms = main_sequence_time();
      real t_gs = 0.15*t_ms;
      real t_b  = base_giant_time(t_ms);

      real l_agb = l_g*pow(t_gs/(t_nuc + t_b - age), 1.17);
      l_agb = min(l_agb, maximum_luminosity(mass));

      real r_agb = (0.25*pow(l_agb, 0.4)
             + 0.8*pow(l_agb, 0.67))/pow(mass, 0.27);

      return r_agb;
   }
#endif

real thorne_zytkow::gyration_radius_sq() {

     return cnsts.parameters(convective_star_gyration_radius_sq); 
}
