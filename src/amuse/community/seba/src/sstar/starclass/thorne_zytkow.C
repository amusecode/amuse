//
// thorne_zytkow.C

#include "thorne_zytkow.h"
#include "main_sequence.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.


// (GN Feb 2011)
// Proposal: rewrite TZ as rather simple class in which L and R are determined
// from eddington accretion onto central NS plus wind loss from envelope
// base on TZ papers

thorne_zytkow::thorne_zytkow(main_sequence & m) : single_star(m) {


  PRL(core_mass);
      delete &m;

      // (SPZ: Removed 23March2000)
      //      birth_mass = 0;
      if (is_binary_component()) {

	birth_mass      = get_companion()->get_total_mass();
	core_mass  = birth_mass;
	if (get_companion()->remnant()) {
	  magnetic_field  = get_companion()->get_magnetic_field(); 
	  rotation_period = get_companion()->get_rotation_period();
	}
      }
      else if (birth_mass<=0) {
	birth_mass    = cnsts.parameters(kanonical_neutron_star_mass);
      }

      // (GN+SPZ May  4 1999) last update age is time of previous type change
      //last_update_age = next_update_age;

      // (GN+SilT Feb 2011) TZO: start all times from scratch
      last_update_age = 0.;
      relative_age = 0.;
      // next_update_age initialised to 0 (inst. element has next_update_age += ....x
      next_update_age = 0.;

      PRL(core_mass);
      instantaneous_element();

      post_constructor();
} 


void thorne_zytkow::instantaneous_element() {


  // L is L_Edd onto core
  luminosity = 3.3e4 * core_mass;

  // Radius from L = 4 pi R^2 sig Teff^4
  real Teff = 3500.; //(Wikipedia Teff Betelgeuse)

  radius = sqrt(luminosity * pow(Teff/cnsts.parameters(Tsun), -4));
  effective_radius = radius;

  adjust_next_update_age();

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

	 instantaneous_element();

      }
      else {

	//This should not happen
	cerr << "TZO has reached next_update_age!!! That should not happen" <<endl;
	exit(-1);
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


void thorne_zytkow::stellar_wind(const real dt) {

      // Enhanced Reimers. NOTE: dt in Myr  
      real wind_mass = 5.5e-7*dt*radius*luminosity/get_total_mass();
      wind_mass = min(wind_mass, envelope_mass);

      if (is_binary_component())
         get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
      else
         reduce_mass(wind_mass);
   }

//              Mass transfer utilities.

real thorne_zytkow::add_mass_to_accretor(real mdot, bool, const real dt) {

      mdot = accretion_limit(mdot, dt);

      //adjust_accretor_age(mdot);

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

    mdot = mdot_limit(dt, mdot);
    
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

void thorne_zytkow::adjust_next_update_age() {

      real wind_dMdt = 5.5e-7*radius*luminosity/get_total_mass();
      real dt_wind = 0.1*envelope_mass/wind_dMdt;

      next_update_age += Starlab::min(1.,dt_wind);
   }

real thorne_zytkow::gyration_radius_sq() {

     return cnsts.parameters(convective_star_gyration_radius_sq); 
}
