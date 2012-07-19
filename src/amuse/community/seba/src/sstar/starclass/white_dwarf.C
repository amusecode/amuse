//
// white_dwarf.C
//

#include "white_dwarf.h"
#include "hertzsprung_gap.h"
#include "sub_giant.h"
#include "super_giant.h"
#include "helium_star.h"
#include "helium_giant.h"

white_dwarf::white_dwarf(super_giant & g) : single_star(g) {
    
      delete &g;

      real m_tot    = get_total_mass();
      core_mass     = min(0.99*cnsts.parameters(kanonical_neutron_star_mass),
			  core_mass); 
      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

      lose_envelope_decent();

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();
}

white_dwarf::white_dwarf(sub_giant & s) : single_star(s) {
      delete &s;

      real m_tot    = get_total_mass();
      core_mass     = min(0.99*cnsts.parameters(Chandrasekar_mass),
			  core_mass); 
      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();

}

white_dwarf::white_dwarf(hertzsprung_gap & s) : single_star(s) {
    
      delete &s;

      real m_tot    = get_total_mass();
      core_mass     = min(0.99*cnsts.parameters(Chandrasekar_mass),
			  core_mass); 
      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();

}


white_dwarf::white_dwarf(helium_star & h) : single_star(h) {
 
        delete &h;

	real m_tot    = get_total_mass();
	core_mass     = min(0.99*cnsts.parameters(Chandrasekar_mass),
			    core_mass); 
	envelope_mass = m_tot - core_mass;
	accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
        last_update_age = next_update_age;

	lose_envelope_decent();
	
    relative_age = 0.;//1 + nucleair_evolution_time();

	instantaneous_element();
	update();

	post_constructor();

}


white_dwarf::white_dwarf(helium_giant & h) :  single_star(h) {

        delete &h;

	real m_tot    = get_total_mass();
	core_mass     = min(0.99*cnsts.parameters(Chandrasekar_mass),
			    core_mass); 
	envelope_mass = m_tot - core_mass;
	accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
        last_update_age = next_update_age;

	lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

	instantaneous_element();
	update();

	post_constructor();

}

#if 0
void white_dwarf::adjust_initial_star() {

  if (core_mass<=0) {
    core_mass = get_total_mass();
    envelope_mass = 0;
    if (core_mass>=cnsts.parameters(Chandrasekar_mass)) {
      envelope_mass = core_mass-0.75;
      core_mass = 0.75;
    }
  }
        
  if(relative_age<=0)
      // nucleair_evolution_time() should not be used anymore
      // should be relative_age = max(current_time, 0.);

    relative_age = max(current_time, nucleair_evolution_time());

      if (relative_mass>cnsts.parameters(super_giant2neutron_star)) {
         cerr<<"White_dwarf with high mass!"
             <<"\n should be turned into a neutron star"
             <<"\nadopt maximum  WD relative mass"<<endl;
         relative_mass = cnsts.parameters(super_giant2neutron_star);
         return;
      }
   }
#endif
      
void white_dwarf::instantaneous_element() {
//cerr << "white_dwarf::instantaneous_element"<<endl;

        next_update_age = relative_age + cnsts.safety(maximum_timestep);

	luminosity = 40;
	// m_rel may not become larger than M_Ch
//	real m_rel = min(0.99, core_mass/cnsts.parameters(Chandrasekar_mass));

	// Nauenberg, M, 1972, Apj 175, 417
	// mu=4.039 is the mean molecular weight.
	real mu = 2.;

	effective_radius =
	core_radius      =
	radius           = white_dwarf_radius(core_mass, relative_age);
	  //               min(0.1, (0.0225/mu)*
	  //               sqrt(1./pow(m_rel, cnsts.mathematics(two_third))
	  //	           - pow(m_rel, cnsts.mathematics(two_third))));
}       


void white_dwarf::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

// done in add_mass_to_accretor
//        accrete_from_envelope(dt);

        if (core_mass > cnsts.parameters(Chandrasekar_mass) ||
	    core_mass <= cnsts.safety(minimum_mass_step)) {
            // || accreted_mass > 0.2) {
	    // (GN Mar 26 1999) test edge lid detanation 
	    // needs more research

	  if (is_binary_component()) 
	    get_binary()->dump("binev.data", false);
	  else
	    dump("binev.data", false);
	  
	  star_transformation_story(Disintegrated);
	  new disintegrated(*this);
	  // new neutron_star(*this);    // AIC
	  return;
        }

//    real t_nuc = nucleair_evolution_time();
        next_update_age = relative_age + cnsts.safety(maximum_timestep);

    if (relative_age <= 0.){
  //      if (t_nuc>=relative_age) {
           luminosity = 40;
        //relative_age = t_nuc;
        relative_age = 0.;
        
        }
        else {
	  // (GN May  4 1999) fit to Driebe et al 1999
	  real fit_mass = min(0.6, max(0.18, get_total_mass()));
	  real l_max = pow(10, (3.83 - 4.77* fit_mass));
	  luminosity = l_max/pow((relative_age), 1.4);
	}

	
//	real m_rel = min(0.99, core_mass/cnsts.parameters(Chandrasekar_mass));

	   // Nauenberg, M, 1972, Apj 175, 417
	   // mu=4.039 is the mean molecular weight.
           real mu = 2.;
	   core_radius =
	     radius           = white_dwarf_radius(core_mass, relative_age);

//             radius = min(0.1, (0.0225/mu)
//	            * sqrt(1./pow(m_rel, cnsts.mathematics(two_third))
//                    - pow(m_rel, cnsts.mathematics(two_third))));

	// (GN+SPZ May  3 1999) critical mass for nova (Livio 1993; Saas-Fee)
	   real m_crit = 1.7e-4*pow(get_total_mass(),-0.7)
	               * pow(69.9*radius,2.8); // radius in 10^9 cm

	if (envelope_mass >= m_crit) 
	   thermo_nucleair_flash(dt);


//	   if (envelope_mass >= core_mass
//	                     * cnsts.parameters(thermo_nuclear_flash))
//	     thermo_nucleair_flash(dt);
       
	   update();
}

void white_dwarf::update() {

        detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
        effective_radius = max(effective_radius, radius);

}



void white_dwarf::accrete_from_envelope(const real dt) {

     if (envelope_mass>0) {
       real mdot = cnsts.parameters(white_dwarf_accretion_limit)
	 * accretion_limit(envelope_mass, dt);

       core_mass += mdot;
       envelope_mass -= mdot;
     }
}


void white_dwarf::thermo_nucleair_flash(const real dt) {

cerr<<"void white_dwarf::thermo_nucleair_flash() mass="
    << envelope_mass<<" "<<get_core_mass()<<endl;
cerr<<" "<<get_total_mass()<<endl;

        real mdot = accretion_limit(envelope_mass, dt);
        if (is_binary_component() && 
	    get_binary()->get_bin_type()!=Merged &&
	    get_binary()->get_bin_type()!=Disrupted) {
//		Spiral in or dwarf nova

	  cerr << "Perform binary action"<<endl;
	  cerr << "       "<<get_total_mass()<<" a = "
                           <<get_binary()->get_semi();

	  if (get_binary()->roche_radius(this)<=1.0) //was 2.5 
	    common_envelope(mdot);
	  else 
	    nova(mdot);

	  cerr<<" -->"<< get_total_mass()<<" a = "
                      << get_binary()->get_semi()<<endl;
        }

      envelope_mass -= mdot;

}

void white_dwarf::nova(const real mdot) {

     real ecc = mdot/(get_binary()->get_total_mass()-mdot);

     get_binary()->set_semi((1+ecc)*get_binary()->get_semi());
}

void white_dwarf::common_envelope(const real mdot) {
   
     if (is_binary_component() &&
	 get_binary()->get_semi()>radius) {

       real semi = get_binary()->get_semi();
       real m_comp = get_companion()->get_total_mass();
       real r_lobe = get_binary()->roche_radius(this)/semi;
       real a_spi = semi*((get_total_mass()-mdot)/get_total_mass())
                  / (1. + (2.*mdot
	          / (cnsts.parameters(common_envelope_efficiency)
	          *  cnsts.parameters(envelope_binding_energy)
		  *  r_lobe*m_comp)));
       
       if (a_spi>=radius) 
	 get_binary()->set_semi(a_spi);
       else
	 get_binary()->set_semi(radius);
     }
}

star* white_dwarf::subtrac_mass_from_donor(const real dt, real& mdot) {

        mdot = relative_mass*dt/get_binary()->get_donor_timescale();
        mdot = mass_ratio_mdot_limit(mdot);

        if (mdot<=envelope_mass)
           envelope_mass -= mdot;
        else {
	  mdot -= envelope_mass;
	  envelope_mass = 0;
	  if (mdot >= core_mass)
	    mdot = core_mass - cnsts.safety(minimum_mass_step);
	  core_mass -= mdot;
	}

        return this;
}

real white_dwarf::add_mass_to_accretor(const real mdot) {

        if (mdot<0) {
           cerr << "white_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   return 0;
        }

        adjust_accretor_age(mdot);
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);
	
        return mdot;

     }

real white_dwarf::add_mass_to_accretor(real mdot, const real dt) {

        if (mdot<0) {
           cerr << "white_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   mdot = 0;
        }

	
	
// (GN+SPZ May  3 1999)
	real mu =1;
	if (is_binary_component() && 
	    !get_companion()->hydrogen_envelope_star()) mu =2;
	

	if (mdot >= eddington_limit(radius, dt, mu)) {
	  effective_radius = min(10., get_binary()->roche_radius(this));
	}


        mdot = accretion_limit(mdot, dt);

// (GN+SPZ May  3 1999) test helium accretion and steady burning accumulate
// helium on core of wd. Steady burning between maximum_steady_burning
// and minimum_steady_burning 
// (Sienkewicz 1980, described in van den Heuvel et al. 1992) 
// If this layer becomes more heavy than 0.2 Msun: boom! see evolve_element

	if (is_binary_component() &&
	    get_companion()->hydrogen_envelope_star() 
	    && mdot < minimum_steady_burning(dt)) {// normal accretion

	  envelope_mass += mdot;
	} 
	else { // helium accretion or steady burning

	  core_mass += mdot;
	  accreted_mass += mdot;
	}
 
        adjust_accretor_age(mdot);
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);

        return mdot;

     }

#if 0
real white_dwarf::accretion_limit(const real mdot, const real dt) {

        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington) return eddington;

        return mdot;
     }
#endif

real  white_dwarf::accretion_limit(const real mdot, const real dt) {

  return min(maximum_steady_burning(dt), mdot);  

}

// (GN Mar 30 1999) steady burning limit 
// (Sienkewicz 1980, described in van den Heuvel et al. 1992)
real  white_dwarf::maximum_steady_burning(const real dt) {

// Fit to Iben & Tutukov 1989
  real steady_burning = 4.e-1*pow(get_total_mass(),1.8)*dt;
  return steady_burning;
}

real  white_dwarf::minimum_steady_burning(const real dt) {

// Fit to Iben & Tutukov 1989
  return 1.8e-7*pow(get_total_mass(),2.7)*dt;

}


void white_dwarf::adjust_accretor_age(const real mdot,
				      const bool rejuvenate) {

        real m_rel_new;
        real m_tot_new = get_total_mass() + mdot;
        if (m_tot_new>relative_mass)
           m_rel_new = m_tot_new;
        else m_rel_new = relative_mass;

    real t_nuc_old = 0.;//nucleair_evolution_time();
	real z_new = get_metalicity();
    real t_nuc_new = 0.;//nucleair_evolution_time(m_rel_new, m_tot_new, z_new);

        real dtime = relative_age - t_nuc_old;

        relative_age = t_nuc_new + dtime
	             * (1-pow(mdot/(get_total_mass()+mdot),
		     	      cnsts.parameters(rejuvenation_exponent)));
	
     }

real white_dwarf::zeta_thermal() {

     // Based on white dwarf mass radius relation.
     // zeta = (m/r)dr/dm;
  
     return -cnsts.mathematics(one_third);
}

real white_dwarf::zeta_adiabatic() {

     // Based on white dwarf mass radius relation.
     // zeta = (m/r)dr/dm;
  
     return -cnsts.mathematics(one_third);
}

star* white_dwarf::merge_elements(star* str) {

     envelope_mass = 0;		//Destroy disc.
	
     real merger_core = str->get_core_mass();

     add_mass_to_accretor(str->get_envelope_mass());

     if (relative_mass<get_total_mass() + merger_core)
       relative_mass=get_total_mass() + merger_core;
     core_mass += merger_core;

     spec_type[Merger]=Merger;
     instantaneous_element();

     return this;
}

star* white_dwarf::reduce_mass(const real mdot) {

      if (envelope_mass<+mdot) {
	real mdot_rest = mdot - envelope_mass;
	envelope_mass = 0;
	if (mdot_rest >= core_mass)
	  mdot_rest = core_mass - cnsts.safety(minimum_mass_step);
	core_mass -= mdot_rest;
      }
      else envelope_mass -= mdot;

      return this;
}


real white_dwarf::gyration_radius_sq() {

  return cnsts.parameters(homogeneous_sphere_gyration_radius_sq); 
}

stellar_type white_dwarf::get_element_type(){

  if (core_mass     < cnsts.parameters(helium_dwarf_mass_limit) &&
      relative_mass < cnsts.parameters(upper_ZAMS_mass_for_degenerate_core))
    return Helium_Dwarf;
  
  else if (relative_mass >= cnsts.parameters(super_giant2neutron_star) &&
           core_mass > cnsts.parameters(carbon_dwarf_mass_limit))
    return Oxygen_Dwarf;
  
  else 
    return Carbon_Dwarf;
}

