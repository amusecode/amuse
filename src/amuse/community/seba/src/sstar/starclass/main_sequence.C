//
// main_sequence.C
//
// derived class of class star.
// class main_sequence describes stellar evolution for a main sequence star.
// Class main_sequence will automatically create the following 
// base classes: star, single star and starbase.
//
// Implementation for metalicity dependence from 
// Hurley, J., Pols, O. & Tout, C., 2000, MNRAS 315, 534--569

#include "main_sequence.h"
#include "brown_dwarf.h"
#include "proto_star.h"

// Default (empty) constructor in main_sequence.h

main_sequence::main_sequence(proto_star & p) : single_star(p) {

       delete &p; 

       last_update_age = 0;
       relative_age = 0;
       relative_mass = envelope_mass + core_mass;
       envelope_mass = relative_mass;// - 0.01;
       core_mass = 0.0; //0.01;
    
       adjust_next_update_age();
       update_wind_constant();

       instantaneous_element();
       update();

       post_constructor();
}

void main_sequence::update() {

  // last update age is set after stellar expansion timescale is set.
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

  detect_spectral_features();

}

void main_sequence::update_wind_constant() {
    cerr<<"enter main_sequence::update_wind_constant"<<endl;
#if 0
    // wind_constant is fraction of envelope lost in nuclear lifetime
    // of stars. Should be updated after mass accretion
    // (SPZ+GN: 1 Oct 1998)
    
  if (relative_mass >= cnsts.parameters(massive_star_mass_limit)) {

    real m_core = 0.073 * (1 + cnsts.parameters(core_overshoot))
                        * pow(relative_mass, 1.42);
    m_core = min(m_core, cnsts.parameters(massive_star_mass_limit));
    
    // extra enhanced mass loss for stars with M>80 Msun.
    // to make low-mass compact objects. (SPZ+GN:24 Sep 1998)
    if (m_core>=35) {
      if (m_core<=55)
	m_core = 35 - 1.25 * (m_core -35); 
      else
	m_core = 10 + 0.75 * (m_core-55);
    }
    
    if (m_core>get_total_mass())
      m_core = get_total_mass();

    wind_constant = (get_relative_mass()-m_core)
                  * cnsts.parameters(massive_star_envelope_fraction_lost);

    cerr << "Main_sequence wind treatment for stars with M >= "
	 << cnsts.parameters(massive_star_mass_limit)
	 << " for " << identity << endl
	 << "   M = " << get_total_mass() << " [Msun] "
	 << "   Mdot = " << wind_constant
	 << " t^" << cnsts.parameters(massive_star_mass_loss_law) 
	 << " [Msun/Myr] "
	 << endl;
  }
  else {

    wind_constant = relative_mass
                  * cnsts.parameters(non_massive_star_envelope_fraction_lost);
  }
#endif

#if 0
    // (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_constant = meader_fit_dm;
    else {// constant
      real final_mass = 43; // final mass after ms
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // no wind for low mass ms stars
    wind_constant = 0;
  }

  wind_constant = max(wind_constant, 0.0);
#endif

    // wind_constant is in solar masses per year
    // Should be updated after mass accretion
    // (ST: 17 Sep 2009)
    
    // Vink 2000, 20001
    // Massive stars, including multi scattering effects
    cerr<< "Vink mass loss: what is the range of validity?"<<endl;
    real dm_v = 0;
//    if (metalicity > cnsts.parameters(solar_metalicity)/30. && metalicity < 3*cnsts.parameters(solar_metalicity)){
//        real sigma;//electron scattering cross section
//        if (temperature() >= 35000){
//            sigma = 0.34;
//        }
//        else if(temperature() < 30000){
//            sigma = 0.31;
//        }
//        else {
//            sigma = 0.32;
//        }
//        real rad_acc = 7.66E-5 * sigma * luminosity / get_total_mass();
//        real log_density = -14.94 + 0.85 * log10(metalicity/cnsts.parameters(solar_metalicity)) +3.2*rad_acc;
//        real Tjump = (61.2 + 2.59*log_density)*1000;
//        real arg_dm_v;
//        if (temperature() >= 12500 && temperature() <= Tjump){
//            arg_dm_v = -6.688 + 2.210*log10(luminosity/1.E5) - 1.339*log10(get_total_mass()/30) 
//            - 1.601*log10(1.3/2.0) + 1.07*log10(temperature()/20000) + 0.85*log10(metalicity/solar_metalicity);
//            dm_v = pow(10, arg_dm_v);
//        }
//        else if(temperature() >= Tjump && temperature() <= 50000){
//            arg_dm_v = -6.697 + 2.194*log10(luminosity/1.e5) - 1.313*log10(get_total_mass()/30) -1.226*log10(2.6/2.0)
//            +0.933*log10(temperature()/40000) -10.92*pow(log10(temperature()/40000),2)
//            + 0.85*log10(metalicity/cnsts.parameters(solar_metalicity));
//            dm_v = pow(10, arg_dm_v);
//        }
//    }    

    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x = min(1.0, (luminosity -4000.0)/500.0);
        cerr<< "x factor de Jager mass loss HG correct?" << endl;
        dm_dj = x * 9.6E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5);
        cerr<< "exponent metallicity dependence de Jager mass loss MS correct?" <<endl<<"ook 9.6->9.6310"<< endl;
}
    
    
    PRC(dm_v);PRL(dm_dj);
    wind_constant = max(max(dm_dj, dm_v), 0.0);
    //wind_constant = 0;//1.0E-5 * get_total_mass();
    
    if (dm_dj > dm_v) cerr<< "MS: de Jager"<<endl;
    else cerr<<"MS: Vink"<<endl;
}


real main_sequence::bolometric_correction()
{
  // temperature() is defined in Kelvin.
  // here we should use old 10^3K implementation 
  // (SPZ+GN: 1 Oct 1998)
  real temp_in_kK = 0.001 * temperature();

  real bc;
  if (temp_in_kK < 4.452)
    bc = 2.5*log10((6.859e-6*pow(temp_in_kK,8) + 9.316e-3)
		   / (1. + 5.975e-10*pow(temp_in_kK,14)));
  else if (temp_in_kK < 10.84)
    bc = 2.5*log10((3.407e-2*pow(temp_in_kK,2.))
		   / (1. + 1.043e-4*pow(temp_in_kK, 4.5)));
  else
    bc = 2.5*log10((2728./pow(temp_in_kK, 3.5) 
		    + 1.878e-2*temp_in_kK)
		   / (1. + 5.362e-5*pow(temp_in_kK,3.5)));

  return bc;
}

// main_sequence stars do not have a compact core.
// for convenience the core is set to a small value.
real main_sequence::main_sequence_core_mass()
{
    real m_core = 0.01;
    m_core = max(core_mass, m_core);
    if (m_core > get_total_mass()) m_core = get_total_mass();
   
    return m_core;
}

real main_sequence::main_sequence_core_radius()
{   cerr<<"Is main_sequence_core_radius still valid with HPT tracks? no ms core!"<<endl;
    return min(0.01, radius);
}

#if 0
// add mass to accretor
// is a separate function (see single_star.C) because of the variable
// stellar wind. The wind_constant must be reset.
real main_sequence::add_mass_to_accretor(const real mdot) {

  if (mdot<=0) {
    cerr << "main_sequence::add_mass_to_accretor(mdot=" << mdot << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
    cerr << "Action: No action!" << endl;

    return 0;
  }
 
    adjust_accretor_age(mdot);
    envelope_mass += mdot;
    accreted_mass += mdot;
    if (relative_mass<get_total_mass()) 
      relative_mass = get_total_mass();

    update_wind_constant();
    
    set_spec_type(Accreting);
  
    return mdot;
}

// Proper accretion routine.
// inpended and consistent.

real main_sequence::add_mass_to_accretor(real mdot, const real dt) {
  //	error checking is just for safety and debugging purposes.
  //cerr<<"real single_star::add_mass_to_accretor(mdot="<<mdot<<", dt="<<dt<<")"<<endl;

  if (mdot<0) {
    cerr << "single_star::add_mass_to_accretor(mdot=" << mdot 
	 << ", dt=" << dt << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
    cerr << "Action: put mdot to zero!" << endl;
    return 0;
  }

  mdot = accretion_limit(mdot, dt);
  adjust_accretor_age(mdot);
  if (relative_mass<get_total_mass() + mdot) 
    relative_mass = get_total_mass() + mdot;

  update_wind_constant();

  envelope_mass += mdot;
  accreted_mass += mdot;

  adjust_accretor_radius(mdot, dt);

  set_spec_type(Accreting);

  //cerr<<relative_age<<" "<<next_update_age<<" "<<main_sequence_time()<<endl;
  //cerr<<"accretor radius adjusted "<<effective_radius<<endl;

  //cerr<<"leave add_mass"<<endl;
  
  return mdot;
}
#endif

// used for RLOF
star* main_sequence::subtrac_mass_from_donor(const real dt, real& mdot)
{
  //  cerr << " main_sequence::subtrac_mass_from_donor"<<endl;

  mdot = relative_mass*dt/get_binary()->get_donor_timescale();

  mdot = mass_ratio_mdot_limit(mdot);

  if (envelope_mass <= mdot) {
    mdot = envelope_mass;
    envelope_mass = 0;
    //star_transformation_story(Helium_Star);
    //return dynamic_cast(star*, new helium_star(*this));
      cerr<<"ERROR!!:constructor helium_star(main_sequence) is commented out"<<endl;
  }

  if (low_mass_star()) { // only when mass transfer timescale = nuc?
    // after mass is subtracted star becomes lower mass star
    // (SPZ+GN:24 Sep 1998)
    adjust_donor_age(mdot);
    update_relative_mass(relative_mass-mdot);
  }

  adjust_donor_radius(mdot);
    
  envelope_mass -= mdot;
  return this;
  
}

star* main_sequence::merge_elements(star* str) {

      star* merged_star = this;

      add_mass_to_core(str->get_core_mass());

      //core_mass += str->get_core_mass();
      //if (relative_mass<get_total_mass())
      //   update_relative_mass(get_total_mass());

      if (str->get_envelope_mass()>0) 
         add_mass_to_accretor(str->get_envelope_mass());

      spec_type[Merger]=Merger;

      switch(str->get_element_type()) {
	 case Hyper_Giant:
         case Hertzsprung_Gap: 	
         case Sub_Giant: 	
         case Horizontal_Branch: 
         case Super_Giant: 
         case Carbon_Star: 
         case Helium_Star: 
         case Helium_Giant: 
         case Carbon_Dwarf: 
  	 case Oxygen_Dwarf:
         case Helium_Dwarf: 
	      if (relative_mass <
		  cnsts.parameters(massive_star_mass_limit)) {
		star_transformation_story(Hertzsprung_Gap);

		// (GN+SPZ May  4 1999) should return now
		//  merged_star = dynamic_cast(star*, 
		//  new hertzsprung_gap(*this));
		//  dump(cerr, false);

		// Chose relative_age to be next update age!
		// otherwise sub_giants become unhappy.
cerr << "Merge MS+wd"<<endl;		
PRC(relative_age);PRC(next_update_age);
//		relative_age = next_update_age;
		
		return dynamic_cast(star*, new hertzsprung_gap(*this));
	      }
	      else {
		star_transformation_story(Hyper_Giant);
		//  merged_star = dynamic_cast(star*, 
		//  new wolf_rayet(*this));
		return dynamic_cast(star*, new hyper_giant(*this));
	      }
         case Thorn_Zytkow :
	 case Xray_Pulsar:
         case Radio_Pulsar:
         case Neutron_Star :
         case Black_Hole   : 
              star_transformation_story(Thorn_Zytkow);
	      // merged_star = dynamic_cast(star*, 
	      // new thorne_zytkow(*this));
	      return dynamic_cast(star*, new thorne_zytkow(*this));
	      default:	   instantaneous_element();
      }
      
      return merged_star;

}
           
// Star is rejuvenated by accretion.
void main_sequence::adjust_accretor_age(const real mdot,
					const bool rejuvenate=true) {
    
  //  cerr << "main_sequence::adjust_accretor_age(mdot = " << mdot<<")"<<endl;
  
      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_ms_old = main_sequence_time();
      real z_new = get_metalicity();
      real t_ms_new = main_sequence_time(m_rel_new, z_new);

      relative_age *= (t_ms_new/t_ms_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
    
      // next_update_age should not be reset here,
      // is done in add_mass_to_accretor, where also relative_mass
      // is updated (SPZ+GN: 1 Oct 1998)
      // next_update_age = t_ms_new; 

   }

void main_sequence::adjust_age_after_wind_mass_loss(const real mdot,
                                        const bool rejuvenate=true) {
    real m_rel_new;
    real m_tot_new = get_total_mass() - mdot;
    //if (m_tot_new<relative_mass)
    //    m_rel_new = m_tot_new;
    //else m_rel_new = relative_mass;
    m_rel_new = m_tot_new;
    
    real t_ms_old = main_sequence_time();
    real z_new = get_metalicity();
    real t_ms_new = main_sequence_time(m_rel_new, z_new);
    
    relative_age *= (t_ms_new/t_ms_old);
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(-1. * mdot/m_tot_new); 
    
    // next_update_age should not be reset here,
    // is done in add_mass_to_accretor, where also relative_mass
    // is updated
    // next_update_age = t_ms_new; 
    
}



// Low-mass main-sequence donor lifetimes are expanded by
// reducing relative_mass
// (SPZ+GN:25 Sep 1998)
void main_sequence::adjust_donor_age(const real mdot) { 

      real m_rel_new = get_relative_mass() - mdot;

      real z_new = get_metalicity();
      relative_age *= main_sequence_time(m_rel_new, z_new)
	            / main_sequence_time();

}


// Adiabatic responce function for main_sequence star.
// Used for determining mass_transfer_timescale.
// Increasing zeta stabilizes binary.
real main_sequence::zeta_adiabatic() {

      real z;

      if (get_relative_mass()<=0.4)         // convective envelope
	z = -cnsts.mathematics(one_third);
      else if(low_mass_star()) {
	z = 2; // was 0.55 but this causes cv's to transfer on a dynamicall
	       // timescale where aml-driven is expected.
      }
      else if(medium_mass_star()) {
	z = 4; // Eggleton's book 
      } 
      else
	z = 4; // somewhare between -100 and 100?

      return z;
   }

// Thermal responce function for main_sequence star.
// Used for determining mass_transfer_timescale.
// (SPZ+GN: 1 Oct 1998)
real main_sequence::zeta_thermal() {

      real z = -1;

      if (get_relative_mass()<=0.4)
         z = 0;                         // Unknown
      else if (low_mass_star())
	z = 0.9;	                // Pols & Marinus 1995
                                        // (GN+SPZ Apr 29 1999) was -0.5
      else 
	z = 0.55; 	                //  (GN+SPZ Apr 29 1999) was -1

      return z;
   }

star* main_sequence::reduce_mass(const real mdot) {

    adjust_age_after_wind_mass_loss(mdot, true);
    envelope_mass -= mdot;
    if (relative_mass > get_total_mass()){
        update_relative_mass(get_total_mass());
    }
    
    if (relative_mass < cnsts.parameters(minimum_main_sequence)) {
        // Main_sequence star will not continue core hydrogen burning.
        star_transformation_story(Brown_Dwarf);
        return dynamic_cast(star*, new brown_dwarf(*this));
    }
    
//   On the MS there is no defined core yet, 
//   no He star can form yet 
//   if (envelope_mass<=mdot) {
//         envelope_mass = 0;
//         //star_transformation_story(Helium_Star);
//         //return dynamic_cast(star*, new helium_star(*this));
//      }

      return this;
   }

void main_sequence::adjust_next_update_age() {

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = 0;
      next_update_age = main_sequence_time();
   }

void main_sequence::detect_spectral_features() {

      single_star::detect_spectral_features();

      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;
      if (get_relative_mass() > turn_off_mass(current_time)
	                   * (1+cnsts.parameters(Blue_straggler_mass_limit)))
	spec_type[Blue_Straggler]=Blue_Straggler;
   }

#if 0
real main_sequence::stellar_radius(const real mass, const real age) {

      real alpha, beta, gamma, delta, kappa, lambda;

      real log_mass = log10(mass);

      if (mass>1.334) {
         alpha = 0.1509 + 0.1709*log_mass;
         beta  = 0.06656 - 0.4805*log_mass;
         gamma = 0.007395 + 0.5083*log_mass;
         delta = (0.7388*pow(mass, 1.679)
               - 1.968*pow(mass, 2.887))
               / (1.0 - 1.821*pow(mass, 2.337));
      }
      else {
         alpha = 0.08353 + 0.0565*log_mass;
         beta  = 0.01291 + 0.2226*log_mass;
         gamma = 0.1151 + 0.06267*log_mass;
         delta = pow(mass, 1.25)
               * (0.1148 + 0.8604*mass*mass)
               / (0.04651 + mass*mass);
      }

      if (mass<1.334) {
         kappa  = 0.2594 + 0.1348*log_mass;
         lambda = 0.144 - 0.8333*log_mass;
      }
      else {
         kappa  = 0.092088 + 0.059338*log_mass;
         lambda = -0.05713 + log_mass*(0.3756 -0.1744);
      }

      real t_ms = main_sequence_time(mass);
      real ff    = age/t_ms;
      real r_ms  = delta*pow(10., (ff*(ff*(ff*gamma + beta) + alpha)));

      return r_ms;
   }
#endif

// Fit to Claret & Gimenez 1990, ApSS 196,215, (SPZ+GN:24 Sep 1998)
real main_sequence::gyration_radius_sq() {

  real m = get_total_mass();

  // gravitational acceleration at surface.
  real g = cnsts.physics(G)
         * m*cnsts.parameters(solar_mass)
         / pow(get_effective_radius() * cnsts.parameters(solar_radius), 2);

  // constant part
  real A = -1.5;
  real B = 0.2;
  
  // linear interpolation
  if (low_mass_star()) {
    A = -3.8 + 1.8*m;
    B = 0.77 - 0.44*m;
  }

  real k = pow(10., (A + B*log10(g)));

  return k*k;

}

real main_sequence::nucleair_evolution_timescale() {
  // t_nuc = 10^10 [years] Msun/Lsun.
  // Assumed that 0.1 Msun is thermalized.

  real fused_mass = 0.1*relative_mass;

  return cnsts.parameters(energy_to_mass_in_internal_units)
       * fused_mass/luminosity;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// New metalicity dependencies from thesis of Hurley, J., 2000
// starts here.
//


// Adjust radius & luminosity at relative_age
void main_sequence::instantaneous_element() {

  luminosity       = main_sequence_luminosity(relative_age, 
					      relative_mass, metalicity);
  radius           = main_sequence_radius(relative_age, 
					  relative_mass, metalicity);
    
  //effective_radius = max(effective_radius, radius);
  effective_radius = radius;
}

// Evolve a main_sequence star upto time argument according to
// the new 2000 models.
void main_sequence::evolve_element(const real end_time) {

     real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;

    if (relative_mass < cnsts.parameters(minimum_main_sequence)) {
      // Main_sequence star will not ignite core hydrogen burning.

	star_transformation_story(Brown_Dwarf);
	new brown_dwarf(*this);
	return;
    }

    if (relative_age <= next_update_age) {

      instantaneous_element();

    } else {

	// Main sequence star's age exceeds hydrogen core burning
	// lifetime.

	//if (relative_mass < cnsts.parameters(massive_star_mass_limit)) {

            star_transformation_story(Hertzsprung_Gap);
            new hertzsprung_gap(*this);
            return;
	//} else {

        //    star_transformation_story(Hyper_Giant);
        //    new hyper_giant(*this);
        //    return;
	//}
    }

    update();
    //stellar_wind(dt);
}



real main_sequence::get_evolve_timestep() {
    
    // (GN+SPZ Apr 28 1999) was a bit too small
    //  return max(next_update_age - relative_age
    //	     -0.5*cnsts.safety(minimum_timestep),
    //	     cnsts.safety(minimum_timestep));
    
    // (GN+SPZ May  5 1999) type change time must be small because of rapid
    // growth of giants at end phase 0.0001 seems to be OK (?)
    //  return max(next_update_age - relative_age - (0.5*0.001), 0.001);
    
    real t_hook = 0.99*next_update_age;
    if (relative_age < t_hook){
        return max(t_hook - relative_age, 0.001);
    }
    else{    
        return max((next_update_age - relative_age), 0.0001);
    }
}




