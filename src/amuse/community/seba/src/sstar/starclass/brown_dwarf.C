//
// brown_dwarf.C
//
//to be based on Dantona, F., Mazzitelli, I., 1985, ApJ 296, 502
// and  2000astro.ph..5557, Chabrier, G.; Baraffe, I.; Allard, F.; 
// Hauschildt, P.

#include "brown_dwarf.h"
#include "main_sequence.h"
#include "proto_star.h"

brown_dwarf::brown_dwarf(proto_star & p) : single_star(p) {
    
      delete &p;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

      last_update_age = 0;
      relative_age = 0;
	  
      instantaneous_element();
      update();

      post_constructor();
}

brown_dwarf::brown_dwarf(main_sequence & m) : single_star(m) {
    
      delete &m;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;
      instantaneous_element();
      update();

      post_constructor();
}

void brown_dwarf::instantaneous_element() {

     luminosity = 1.e-4; 

     core_radius = radius = 0.1;
       
     update();
}

real brown_dwarf::get_evolve_timestep() {
    // (GN+SPZ Apr 28 1999) was a bit too small
    //  return max(next_update_age - relative_age
    //	     -0.5*cnsts.safety(minimum_timestep),
    //	     cnsts.safety(minimum_timestep));
    
    // (GN+SPZ May  5 1999) type change time must be small because of rapid
    // growth of giants at end phase 0.0001 seems to be OK (?)
    // return max(next_update_age - relative_age - (0.5*0.001), 0.001);

    return max(next_update_age - relative_age, 0.0001);
}



void brown_dwarf::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

        next_update_age = relative_age + cnsts.safety(maximum_timestep);

	//Burrows & Libert 1993, J. Rev. Mod. Phys. 65, 301
	luminosity = 938 * pow(relative_mass, 2.64); 
	if(relative_age>1)
	  luminosity = 938 * pow(relative_mass, 2.64) / pow(relative_age, 1.3);

        core_radius = radius = 0.1;
       
        update();
     }

void brown_dwarf::update() {

     detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
     effective_radius = radius;

     }

real brown_dwarf::brown_dwarf_core_mass() {
    
        return 0.01 * get_total_mass();
     }


//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real brown_dwarf::mdot_limit(const real dt){
    real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    return mass_ratio_mdot_limit(mdot);
    
}


star* brown_dwarf::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = mdot_limit(dt);
    
        if (mdot<=envelope_mass)
	  envelope_mass -= mdot;
        else if (mdot>envelope_mass) 
	  envelope_mass = 0;

        return this;
     }


real brown_dwarf::add_mass_to_accretor(real mdot, bool, const real dt) {

        if (mdot<0) {
           cerr << "brown_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   mdot = 0;
        }

        mdot = accretion_limit(mdot, dt);
 
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);
	
        return mdot;
     }

real brown_dwarf::accretion_limit(const real mdot, const real dt) {

  if (dt < 0) return mdot;

        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington)
	  return eddington;

        return mdot;
     }


real brown_dwarf::zeta_thermal() {

     return 0;
}

star* brown_dwarf::merge_elements(star* str) {

     real merger_core = str->get_core_mass();

     add_mass_to_accretor(str->get_envelope_mass(), 
			  cnsts.parameters(spiral_in_time), str->hydrogen_envelope_star());

     if (relative_mass < get_total_mass() + merger_core)
       relative_mass=get_total_mass() + merger_core;
     core_mass += merger_core;

     spec_type[Merger]=Merger;
     instantaneous_element();

     return this;
}

star* brown_dwarf::reduce_mass(const real mdot) {

      if (envelope_mass < mdot)
	envelope_mass = 0;
      else
	envelope_mass -= mdot;

      return this;
}

real brown_dwarf::gyration_radius_sq() {

    return cnsts.parameters(convective_star_gyration_radius_sq); 
}


stellar_type brown_dwarf::get_element_type() {
  
  if (get_total_mass() < 0.1*cnsts.parameters(minimum_main_sequence))
      return Planet;
    else
      return Brown_Dwarf;
}
	
