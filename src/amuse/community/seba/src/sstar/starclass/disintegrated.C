//
// disintegrated.C
//

#include "disintegrated.h"
//#include "hyper_giant.h"
#include "super_giant.h"
//#include "thorne_zytkow.h"
#include "helium_giant.h"


disintegrated::disintegrated(super_giant & g) : single_star(g) {

        delete &g;

        for (int i=Emission; i<no_of_spec_type; i++)
            spec_type[i] = NAC;
	spec_type[Dsntgr] = Dsntgr;

	suddenly_lost_mass = 0;

        super_nova();

        core_mass = envelope_mass = cnsts.safety(minimum_mass_step);
        radius = effective_radius = core_radius = 0;
        luminosity = 1;
        velocity = 0;
        wind_constant = accreted_mass = 0;
	magnetic_field = rotation_period = 0;
	birth_mass=0;

        if (is_binary_component() &&
	    get_binary()->get_bin_type()!=Merged) 
           get_binary()->set_bin_type(Disrupted);

         instantaneous_element();

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

     disintegrated::disintegrated(helium_giant & h) : single_star(h) {

        delete &h;

        for (int i=Emission; i<no_of_spec_type; i++)
            spec_type[i] = NAC;
	spec_type[Dsntgr] = Dsntgr;

	suddenly_lost_mass = 0;

        super_nova();

        core_mass = envelope_mass = cnsts.safety(minimum_mass_step);
        radius = effective_radius = core_radius = 0;
        luminosity = 1;
        velocity = 0;
        wind_constant = accreted_mass = 0;
	magnetic_field = rotation_period = 0;
	birth_mass=0;

        if (is_binary_component() &&
	    get_binary()->get_bin_type()!=Merged) 
           get_binary()->set_bin_type(Disrupted);

         instantaneous_element();

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

     disintegrated::disintegrated(white_dwarf & w) : single_star(w) {


        delete &w;

        for (int i=Emission; i<no_of_spec_type; i++)
            spec_type[i] = NAC;
	spec_type[Dsntgr] = Dsntgr;

	suddenly_lost_mass = 0;

        super_nova();

        core_mass = envelope_mass = cnsts.safety(minimum_mass_step);
        radius = effective_radius = core_radius = 0;
        luminosity = 1;
        velocity = 0;
        wind_constant = accreted_mass = 0;
	magnetic_field = rotation_period = 0;
	birth_mass=0;

        if (is_binary_component() &&
	    get_binary()->get_bin_type()!=Merged)  
           get_binary()->set_bin_type(Disrupted);

         instantaneous_element();

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

void disintegrated::instantaneous_element() {

  evolve_element(relative_age);
}
 
void disintegrated::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

        next_update_age  = relative_age + cnsts.safety(maximum_timestep);

        update();
     }

void disintegrated::update() {

// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

}

bool disintegrated::super_nova() {

	suddenly_lost_mass = get_total_mass();

        bool hit_companion = FALSE;
        spec_type[Dsntgr]=Dsntgr;

	return hit_companion;
}

//              Mass transfer utilities.
real disintegrated::mass_transfer_timescale(mass_transfer_type
                                             &type) {
    
        cerr << "disintegrated::mass_transfer_timescale()" << endl;
        cerr << "Disintegrated cannot be donor star!"<<endl;
        cerr << "ABSOLUTE_DT_MIN (" << cnsts.safety(minimum_timestep)
              << ") returned." << endl;
	
     return cnsts.safety(minimum_timestep);
     }

star* disintegrated::subtrac_mass_from_donor(const real dt, real& mdot) {

        cerr << "disintegrated::subtrac_mass_from_donor(dt="
	     << dt << ", mdot=" << mdot << ")" << endl;
        cerr << "Disintegrated cannot be donor star!"<<endl;

        return this;
     }

real disintegrated::add_mass_to_accretor(const real mdot) {

     cerr << "disintegrated::add_mass_to_accretor(mdot="
          << mdot << ")" <<endl;

     return 0;

}

real disintegrated::add_mass_to_accretor(real mdot, const real dt) {

     cerr << "disintegrated::add_mass_to_accretor(mdot="
          << mdot << ", dt=" << dt << ")" <<endl;

     return 0;

     }

star* disintegrated::merge_elements(star* str) {

        cerr << "void disintegrated::merge_stars()"<<endl;

        spec_type[Merger]=Merger;
	
	return this;
 }


star* disintegrated::reduce_mass(const real mdot) {

      return this;
     }

real disintegrated::temperature() {
      return 1;
   }

real disintegrated::magnitude() {
      return 100;
   }

real disintegrated::bolometric_correction() {
      return 0;
   }

void disintegrated::stellar_wind(const real dt) {
   }

real disintegrated::sudden_mass_loss() {

    real mass_lost = suddenly_lost_mass;
    suddenly_lost_mass = 0;

    return mass_lost;

   }

real disintegrated::gyration_radius_sq() {

  return 0;
}
