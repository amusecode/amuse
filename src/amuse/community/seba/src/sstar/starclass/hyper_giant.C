//
// hyper_giant.C
//

#include "hyper_giant.h"
#include "helium_star.h"
#include "main_sequence.h"

hyper_giant::hyper_giant(main_sequence & m) : single_star(m) {
   
        delete &m;

        real m_tot = get_total_mass();
        core_mass = hyper_giant_core_mass();
        envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

        adjust_next_update_age();

	instantaneous_element();
	update();

	post_constructor();
}

void hyper_giant::instantaneous_element() {

  luminosity = 10e+5;
// (GN+SPZ May 12 1999) WR star grow from ms radius to 1000 Rsun
//  effective_radius = radius = 200;

  // (SPZ+GN:  1 Aug 2000)
  // coming from previous type the effective readius should 
  // remain the same.
  if(radius<=0)
     effective_radius = radius = 5;
  core_radius = 5;


  if (envelope_mass <= 0) 
    radius = core_radius;
  
}

void hyper_giant::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

        if (relative_age>next_update_age) {
           create_remnant();
           return;
        }

	// (SPZ+GN: 8 Oct 1998)
// (GN Apr  6 1999) "Wolf Rayet" stars can become helium stars becuase we use
// wolf rayet for massive giants with Hydrogen envelope
	if (get_total_mass() <= cnsts.parameters(minimum_helium_star) ||
	    envelope_mass < cnsts.safety(minimum_mass_step)) {
	      //star_transformation_story(Helium_Star);
          //    new helium_star(*this);
	      cerr<<"ERROR!!:constructor helium_star(Hyper giant) is commented out "<<endl;
        return;
	}

//(GN Mar 30 1999) quick hack to make WR stars evolve to red supergiants
// de Koter 1999 private communication: stars < 50 Msun become red super-
// giants, only higher mass stars become really WR stars. EFT89 fits not 
// appropriate for stars 25 < M < 50
//	PRL(relative_age);PRL(next_update_age);
	real relative_time = relative_age - last_update_age;
	real t_end = next_update_age - last_update_age;
	radius = 
	  effective_radius = 
	     max(radius, 1000*relative_time/t_end);
//	PRL(radius);

	
        update();
	stellar_wind(dt);
}

real hyper_giant::hyper_giant_core_mass() {

//        real m_core = 0.073*(1 + cnsts.parameters(core_overshoot))
//	            * pow(relative_mass, 1.42);
// (GN+SPZ May  3 1999) Yungelson core
//  real m_core = 0.073*(1 + cnsts.parameters(core_overshoot))
//              * pow(relative_mass, 1.42);

  real m_core_TY = 0.066*pow(relative_mass, 1.54);
  real m_core_max = 20 + 0.27*relative_mass;

  real m_core = min(m_core_TY, m_core_max);

//	// Extra enhanced mass loss for stars with M>80 Msun.
//	// to make low-mass compact objects. (SPZ+GN:24 Sep 1998)
//	if (m_core>=35) {
//	  if (m_core<=55)
//	    m_core = 35 - 1.25 * (m_core -35); 
//	  else
//	    m_core = 10 + 0.75 * (m_core-55);
//	}

       if (m_core>get_total_mass()) m_core = get_total_mass();

        return m_core;
     }

void hyper_giant::create_remnant() {

    if (core_mass>=cnsts.parameters(COcore2black_hole)) {
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

#if 0
real hyper_giant::add_mass_to_accretor(const real mdot) {

        if (mdot<0) {
           cerr << "hyper_giant::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   return 0;
        }

        adjust_accretor_age(mdot);
        if (relative_mass<get_total_mass() + mdot)
           relative_mass = get_total_mass() + mdot;
        envelope_mass += mdot;

	// Wind mass loss for massive stars (SPZ:July 1998)
	if (relative_mass >= cnsts.parameters(massive_star_mass_limit)) {

	  cerr << type_string(get_element_type())
	       << " wind treatment for stars with M >="
	       << cnsts.parameters(massive_star_mass_limit)
	       << " for " << identity << endl
	       << "   M = " << get_total_mass() << " [Msun] "
	       << "   Mdot = " << wind_constant
	       << " [Msun/Myr] "
	       << endl;

	}

	set_spec_type(Accreting);
	
	return mdot;

}

real hyper_giant::add_mass_to_accretor(real mdot, const real dt) {

        real m_tot = get_total_mass();

        if (mdot<0) {
/*
           error << "hyper_giant::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           error << "mdot (" << mdot << ") smaller than zero!" << endl;
           error << "Action: proceed!" << endl;
*/
        }

        mdot = accretion_limit(mdot, dt);
        adjust_accretor_age(mdot);
        if (relative_mass<get_total_mass() + mdot)
           relative_mass = get_total_mass() + mdot;
        envelope_mass += mdot;

	// Wind mass loss for massive stars (SPZ:July 1998)
	if (relative_mass >= cnsts.parameters(massive_star_mass_limit)) {

	  // Mass loss constant in time.
	  wind_constant = envelope_mass/(next_update_age-relative_age);
	  wind_constant = max(wind_constant, 0.);

	  cerr << type_string(get_element_type())
	       << " wind treatment for stars with M >="
	       << cnsts.parameters(massive_star_mass_limit)
	       << " for " << identity << endl
	       << "   M = " << get_total_mass() << " [Msun] "
	       << "   Mdot = " << wind_constant
	       << " [Msun/Myr] "
	       << endl;

	}

	set_spec_type(Accreting);
	
        return mdot;

     }
#endif

real hyper_giant::accretion_limit(const real mdot, const real dt) {

        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington) return eddington;

        return mdot;
     }

#if 0
//Stellar wind routine borrowed from class helium_star 
void hyper_giant::stellar_wind(const real dt) {

  if (relative_mass<cnsts.parameters(massive_star_mass_limit))
    return;

  real wind_mass = wind_constant * dt;

  if (wind_mass>=envelope_mass) {
    wind_mass = envelope_mass;
    radius = core_radius;
  }

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_donor_mass(wind_mass);
  return;
} 

#endif


star* hyper_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot);
	
      if (mdot<=envelope_mass)
	envelope_mass -= mdot;
      else {

	mdot = envelope_mass;
	envelope_mass = 0;
      }

      return this;
}




void hyper_giant::adjust_accretor_age(const real mdot,
				     const bool rejuvenate=true) {

        real m_rel_new;
        real m_tot_new = get_total_mass() + mdot;
        if (m_tot_new>relative_mass)
           m_rel_new = m_tot_new;
        else m_rel_new = relative_mass;

        real t_ms_old  = main_sequence_time();
        real t_nuc_old = nucleair_evolution_time();
        real dt_wr_old = t_nuc_old - t_ms_old;

        real t_ms_new  = main_sequence_time(m_rel_new, metalicity);
	real z_new = get_metalicity();
        real t_nuc_new = nucleair_evolution_time(m_rel_new, m_tot_new, z_new);
        real dt_wr_new = t_nuc_new - t_ms_new;

        real dtime = relative_age - t_ms_old;

        real new_relative_age = t_ms_new 
	                      + dtime*(dt_wr_new/dt_wr_old);
	if (rejuvenate)
	   new_relative_age *= rejuvenation_fraction(mdot/m_tot_new); 

	relative_age = max(relative_age, new_relative_age);
                       

     }


real hyper_giant::zeta_adiabatic() {
  
// See Hjellming and Webbink 1987 ApJ, 318, 804

//  real x = core_mass/get_total_mass();
//  real a = -0.220823;
//  real b = -2.84699;
//  real c = 32.0344;
//  real d = -75.6863;
//  real e = 57.8109;
  
//  real z = a + b*x + c*x*x + d*x*x*x + e*x*x*x*x;

  real z = -cnsts.mathematics(one_third);

  return z;

}



real hyper_giant::zeta_thermal() {
        return 0;
     }

star* hyper_giant::reduce_mass(const real mdot) {

  //  cerr <<"hyper_giant::reduce_mass(mdot="<<mdot<<")"<<endl;
  //  dump(cerr, false);
      if(envelope_mass<=mdot) {
	real mdot_rest = mdot - envelope_mass;
	envelope_mass = 0;
	
	
	if (mdot_rest>=core_mass-cnsts.parameters(minimum_helium_star))
	  core_mass = cnsts.parameters(minimum_helium_star);
	else
	  core_mass -= mdot_rest;
      }
      else envelope_mass -= mdot;
           
      return this;
}

void hyper_giant::adjust_next_update_age() {
        
        real t_ms = relative_age;
	if (relative_age<t_ms)
           relative_age = t_ms;

// (GN+SPZ May 12 1999) heavy stars become helium_stars 
// (actually hyper_giant stars) in wind: timestep must not be too large
	real end_time = nucleair_evolution_time();
	real alpha = cnsts.parameters(massive_star_mass_loss_law);
	real t_wind_stripped = end_time
	                     * pow((envelope_mass/wind_constant + 
				    pow((relative_age/end_time),alpha))
				   , 1/alpha);
	
// (GN+SPZ May 12 1999) + eps because otherwise it becomes a remnant 
	t_wind_stripped += 0.001;

//	PRC(t_wind_stripped);PRL(end_time);
	next_update_age = min(end_time, t_wind_stripped);
     }

real hyper_giant::gyration_radius_sq() {

  return cnsts.parameters(convective_star_gyration_radius_sq); 

}

