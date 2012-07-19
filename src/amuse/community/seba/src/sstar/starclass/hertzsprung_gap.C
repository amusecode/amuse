#include "hertzsprung_gap.h"
#include "main_sequence.h"

// This constructor 
// copies the old main_sequence star into the newly created 
// hertzsprung_gap star and destroys the old main_sequence star 
// and finally evolves the hertzsprung_gap in order to determine its
// appearence.
//
// ANSI C++ first creates the base class before the derived classes are
// created. 

     hertzsprung_gap::hertzsprung_gap(main_sequence & m) : single_star(m) {

      delete &m; 	

      //      real m_tot    = get_total_mass();
      //      core_mass     = min(TAMS_helium_core_mass(), m_tot);
      //      envelope_mass = m_tot - core_mass;
      //      core_radius   = helium_core_radius();
                  
         if (relative_mass != get_total_mass()){
             cerr<<"error constructor HG: relative_mass != get_total_mass()"<<endl;   
         }
         
// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      adjust_next_update_age();

      instantaneous_element();
      evolve_core_mass();
      small_envelope_perturbation();   
         
      update();
      post_constructor();
   }


#if 0
void hertzsprung_gap::adjust_initial_star() {

  if(relative_age<=0)
    relative_age = max(main_sequence_time(), relative_age);
}
#endif

star* hertzsprung_gap::reduce_mass(const real mdot) {
    if (envelope_mass<=mdot) {
        envelope_mass = 0;

//        (SPZ+GN: 27 Jul 2000)
//        // non degenerate core < helium_dwarf_mass_limit always(!) become
//        // white dwarfs
//        if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
//            relative_mass < cnsts.parameters(
//                    upper_ZAMS_mass_for_degenerate_core)) {
//            star_transformation_story(Helium_Dwarf);
//            return dynamic_cast(star*, new white_dwarf(*this));
//        } 
//        else {
//            star_transformation_story(Helium_Star);
//            return dynamic_cast(star*, new helium_star(*this));
//        }
        
        real m_HeF = helium_flash_mass(metalicity);
        if (get_total_mass() < m_HeF){
            star_transformation_story(Helium_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this));
        }
        else {
            star_transformation_story(Helium_Star);
            return dynamic_cast(star*, new helium_star(*this));

        }
    }
    real mc_bgb = terminal_hertzsprung_gap_core_mass(get_total_mass()-mdot, metalicity);
    real m_FGB = helium_ignition_mass(metalicity);
    
    if ( core_mass > mc_bgb || relative_mass > m_FGB){
        //the evolution of m_core, luminosity, timescales decouples from M
        // relative mass is no longer kept at same value as total mass
        envelope_mass -= mdot;
    }
    else{
        adjust_age_after_wind_mass_loss(mdot, true);
        envelope_mass -= mdot;
        if (relative_mass > get_total_mass()){
            update_relative_mass(get_total_mass());
        }
    }
    return this;
   }


star* hertzsprung_gap::subtrac_mass_from_donor(const real dt, real& mdot) {
//cerr<<"real hertzsprung_gap::subtrac_mass_from_donor( dt= "<<dt<<")"<<endl;

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();

      mdot = mass_ratio_mdot_limit(mdot);

      if (envelope_mass<=mdot) {
//	  cerr << "Transform hertzsprung_gap star into Heliun stars"<<endl;
         mdot = envelope_mass;
         envelope_mass = 0;

	// (SPZ+GN: 27 Jul 2000)
	// non degenerate core < helium_dwarf_mass_limit always(!) become
	// white dwarfs
	if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
	    relative_mass < cnsts.parameters(
			    upper_ZAMS_mass_for_degenerate_core)) {
	  star_transformation_story(Helium_Dwarf);
	  return dynamic_cast(star*, new white_dwarf(*this));
	} 
	else {
	  star_transformation_story(Helium_Star);
	  return dynamic_cast(star*, new helium_star(*this));
	}
      }

// (GN+SPZ Apr 29 1999) 
      adjust_donor_radius(mdot);

      envelope_mass -= mdot;
      return this;
   }

//  relative age adjusted lineairly with accreted mass.
void hertzsprung_gap::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {

      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_ms_old = main_sequence_time();
      real t_hg_old = hertzsprung_gap_time() - t_ms_old;

      real z_new = get_metalicity();
      real t_ms_new = main_sequence_time(m_rel_new, z_new);

      //For now, we keep metalicity constant (SPZ: 29 May 2001)
      real t_hg_new = hertzsprung_gap_time(m_rel_new, z_new) - t_ms_new;

      real dtime = relative_age - t_ms_old;

// (GN+SPZ May  4 1999) update last_update_age
//      last_update_age = t_ms_new;
    
      relative_age = t_ms_new 
                   + dtime*(t_hg_new/t_hg_old);
    
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new); 

      relative_age = max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep)); 

      // next_update_age should not be reset here
      // next_update_age = t_ms_new + t_hg_new;
}

void hertzsprung_gap::adjust_age_after_wind_mass_loss(const real mdot, const bool rejuvenate=true) {
    
    real m_rel_new;
    real m_tot_new = get_total_mass() - mdot;
//    if (m_tot_new>relative_mass)
//        m_rel_new = m_tot_new;
//    else m_rel_new = relative_mass;
    m_rel_new = m_tot_new; 
    
    real t_ms_old = main_sequence_time();
    real t_hg_old = hertzsprung_gap_time() - t_ms_old;
    
    real z_new = get_metalicity();
    real t_ms_new = main_sequence_time(m_rel_new, z_new);
    //For now, we keep metalicity constant (SPZ: 29 May 2001)
    real t_hg_new = hertzsprung_gap_time(m_rel_new, z_new) - t_ms_new;
    
    real dtime = relative_age - t_ms_old;
    
    // (GN+SPZ May  4 1999) update last_update_age
    //last_update_age = t_ms_new;
    
    relative_age = t_ms_new + dtime*(t_hg_new/t_hg_old);
    
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(-1.*mdot/m_tot_new); 
    
    if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
        cerr<<"relative age updated on HG, but < last_update_age"<<endl;
    }
    relative_age = max(relative_age, 
                       last_update_age + cnsts.safety(minimum_timestep)); 
    
    // next_update_age should not be reset here
    // next_update_age = t_ms_new + t_hg_new;
}


// Adiabatic response function for hertzsprung_gap star.
// Polynomial fit to Hjellming and Webbink 1987 ApJ, 318, 804
real hertzsprung_gap::zeta_adiabatic() {
//cerr<<"real hertzsprung_gap::zeta_adiabatic(): " << endl;

#if 0
      real z;

      real x = core_mass/get_total_mass();
      real A = -0.220823;
      real B = -2.84699;
      real C = 32.0344;
      real D = -75.6863;
      real E = 57.8109;

      if (get_relative_mass()<=0.4)
         z = -cnsts.mathematics(one_third);
      else if (low_mass_star())
         z = A + x*(B + x*(C + x*(D + x*E)));
      else if (medium_mass_star())
         z = 2.25;                 // 15 according to Pols & Marinus 1994
      else                         // We do it differently.
         z = 2.25;                 // lekker puh.
#endif
      real z = 4; // this is neede to prevent Thermal in extreme mass
                     // ratio systems ... 
// (GN+SPZ Apr 29 1999) Pols & Marinus 1994 were maybe right: not!

      return z;
   }

// Thermal response function for hertzsprung_gap star.
real hertzsprung_gap::zeta_thermal() {

      real z;

      if (get_relative_mass()<=0.4)
         z = 0;         // no better estimate present.
      else if (low_mass_star())
         z = -2; 	// -10 according to Pols & Marinus 1994
      else              // Changed to current values
         z = -2;	// by (SPZ+GN: 1 Oct 1998)

      return z;
   }
 
void hertzsprung_gap::adjust_next_update_age() {
   
    real t_ms = main_sequence_time();
    if (relative_age<t_ms) {
        cerr << "WARNING: relative_age < t_ms in Hertzsprung_gap"<<endl;
        relative_age = t_ms;
    }
    real t_bhg = hertzsprung_gap_time();
    if(t_ms<t_bhg) 
        next_update_age = hertzsprung_gap_time();
    else {
        cerr << "WARNING:hertzsprung_gap::adjust_next_update_age()" << endl;
        cerr << "main_sequence time exceeds hertzprung_gap time"<<endl;
        PRC(t_ms);PRL(t_bhg);
        dump(cerr, false);
        cerr << flush;
        exit(-1);
    }


}

void hertzsprung_gap::detect_spectral_features() {

// 		Use standard spectral feature detection.
      single_star::detect_spectral_features();

      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;
   }

#if 0
real hertzsprung_gap::stellar_radius(const real mass, const real age_max) {

      real alpha, beta, gamma, delta;
      real t_ms   = main_sequence_time(mass);
      real age = max(t_ms, age_max);            // Safety
      real t_hg   = hertzsprung_gap_time(t_ms);

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

       real l_g    = giant_luminosity(mass);

       real rt = delta*pow(10., (alpha + beta + gamma));
       real rg = (0.25*pow(l_g, 0.4) + 0.8*pow(l_g, 0.67))
               / pow(mass, 0.27);
       real ff = (age - t_ms)/t_hg;
       real r_hg = rt*pow(rg/rt, ff);

       return r_hg;
   }
#endif


// Stellar Gyration radii squared for detmination of
// angular momentum.
// Implemented by (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

// Helium core at the terminal-age main-sequence
// (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::TAMS_helium_core_mass() {

  // (SPZ+GN: 28 Jul 2000)
  // old Eggleton 2000 (ToBeBook) fit to core mass.
  real mc = (0.11*pow(relative_mass,1.2)
	  + 7.e-5*pow(relative_mass,4))
          / (1+ 2e-4*pow(relative_mass, 3));

  // (SPZ+GN: 28 Jul 2000)
  // produces lower mass core betwee 0 and 20 Msun.
  // is fine for low mass (<5Msun) stars but too small for 
  // high mass (>8Msun) stars.
  //  real mc = (0.05 + 0.08*pow(relative_mass,1.2)
  //          + 7.e-5*pow(relative_mass,4))
  //          / (1+ 2e-4*pow(relative_mass, 3));

  // (SPZ+GN: 31 Jul 2000)
  // TAMS core mass must be smaller than get_total_mass to 
  // assure that it becomes a helium star decently.
  return min(max(mc,core_mass), 
	     get_total_mass()-cnsts.safety(minimum_mass_step));

}

void hertzsprung_gap::update_wind_constant() {
    cerr<<"enter hertzsprung_gap::update_wind_constant"<<endl;

#if 0
// (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense
// wind_constant is fraction of envelope lost in nuclear lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 1 Oct 1998)
    
  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_constant = meader_fit_dm;
    else {// constant
      real final_mass = 30;
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // (GN+SPZ Apr 29 1999) 1% loss on hg

    wind_constant = 0.01*relative_mass;
  }

  wind_constant = max(wind_constant, 0.0);
#endif
    
    // wind_constant is in solar masses per year
    // Should be updated after mass accretion
    // (ST: 17 Sep 2009)
    
    // Vink 2000, 20001
    // Massive stars, including multi scattering effects
    real dm_v = 0;
    if (metalicity > cnsts.parameters(solar_metalicity)/30. && metalicity < 3*cnsts.parameters(solar_metalicity)){
        real sigma;//electron scattering cross section
        if (temperature() >= 35000){
            sigma = 0.34;
        }
        else if(temperature() < 30000){
            sigma = 0.31;
        }
        else {
            sigma = 0.32;
        }
        real rad_acc = 7.66E-5 * sigma * luminosity / get_total_mass();
        real log_density = -14.94 + 0.85 * log10(metalicity/cnsts.parameters(solar_metalicity)) +3.2*rad_acc;
        real Tjump = (61.2 + 2.59*log_density)*1000;
        real arg_dm_v;
        if (temperature() >= 12500 && temperature() <= Tjump){
            arg_dm_v = -6.688 + 2.210*log10(luminosity/1.E5) - 1.339*log10(get_total_mass()/30) 
            - 1.601*log10(1.3/2.0) + 1.07*log10(temperature()/20000) + 0.85*log10(metalicity/solar_metalicity);
            dm_v = pow(10, arg_dm_v);
        }
        else if(temperature() >= Tjump && temperature() <= 50000){
            arg_dm_v = -6.697 + 2.194*log10(luminosity/1.e5) - 1.313*log10(get_total_mass()/30) -1.226*log10(2.6/2.0)
            +0.933*log10(temperature()/40000) -10.92*pow(log10(temperature()/40000),2)
            + 0.85*log10(metalicity/cnsts.parameters(solar_metalicity));
            dm_v = pow(10, arg_dm_v);
        }
    }    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x_dj = min(1.0, (luminosity - 4000.0)/500.0);
        dm_dj = x_dj*9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5);
        cerr<< "exponent metallicity dependence de Jager mass loss HG correct?" << endl;
    }
    
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    cerr <<"Reimers neta correct?"<<endl;
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    
    // Hamann, Koesterke & Wessolowski 1995
    //Reduced WR-like mass loss for small H-envelope mass
    real mu = (get_total_mass()-core_mass)/get_total_mass() * min(5.0,max(1.2, pow(luminosity/7.E4,-0.5)));
    real dm_h = 0;
    if ( mu < 1.){
        dm_h = 1.0E-13 * pow(luminosity, 1.5) * (1.0-mu);
        cerr<<"Hamann mass loss: I'm not convinced about the mu dependence of dm_h"<<endl;
    }
    
    //LBV
    real dm_lbv = 0;
    real x_lbv = 1.0E-5*radius*sqrt(luminosity);
    if(luminosity > 6.0E5 && x_lbv > 1.0) {
        dm_lbv = 0.1 * pow(x_lbv-1.0, 3)*(luminosity/6.0E5-1.0);
    }
    
    PRC(dm_lbv);PRC(dm_v);PRC(dm_dj);PRC(dm_r);PRL(dm_h);
    wind_constant = max(max(max(max(dm_h, dm_dj), dm_v), dm_r), 0.0) + dm_lbv;
    
    if(dm_h > dm_dj && dm_h > dm_v && dm_h > dm_r) cerr<<"HG: WR_like"<<endl;
    else if (dm_dj > dm_v && dm_dj > dm_r) cerr<< "HG: de Jager"<<endl;
    else if (dm_v > dm_r) cerr<<"HG: Vink"<<endl;   
    else cerr<<"HG: Reimers"<<endl;
    
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Metalicity dependency from HPT 200

void hertzsprung_gap::instantaneous_element() {

    luminosity       = hertzsprung_gap_luminosity(relative_age, 
                                                  relative_mass, metalicity);
    radius           = hertzsprung_gap_radius(relative_age, relative_mass, 
                                              get_total_mass(), metalicity);
    
    //      effective_radius = max(effective_radius, radius);
   effective_radius = radius;
}

// Evolve a main_sequence star upto time argument according to
// the new 2000 models.
void hertzsprung_gap::evolve_element(const real end_time) {
      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;
    
      if (relative_age<=next_update_age) {

          instantaneous_element();
          evolve_core_mass();
          small_envelope_perturbation();
      }
      else {
        dump(cerr, false);
        if (relative_mass <= helium_ignition_mass(metalicity) ){
	        star_transformation_story(Sub_Giant);
            new sub_giant(*this);
	    return;
	    }
	    else {
          star_transformation_story(Horizontal_Branch);
          new horizontal_branch(*this);
          return;
	    }
      }
      update();
      //stellar_wind(dt);
}



//Eq.7+
real hertzsprung_gap::terminal_hertzsprung_gap_luminosity(const real mass, 
						      const real z) {
  
  real l_ehg;
  if (mass<helium_ignition_mass(z))
    l_ehg = base_giant_branch_luminosity(mass, z);
  else
    l_ehg = helium_ignition_luminosity(mass, z);

    return  l_ehg;
}

//Eq.7+
real hertzsprung_gap::terminal_hertzsprung_gap_radius(const real mass, 
						      const real mass_tot, const real z) {
  
    real r_ehg;
    if (mass<helium_ignition_mass(z)) {
        real l_bgb = base_giant_branch_luminosity(mass, z);
        r_ehg = giant_branch_radius(l_bgb, mass_tot, z);
    }
    else {
        r_ehg = helium_ignition_radius(mass, mass_tot, z);
        //safety check
        //these lines are not in the HPT2000 article, but they are in the HPT2000 code
        //in case a massive star skips the blue loop phase, 
        // the stellar radius should continue smoothly 
        real l_HeI = helium_ignition_luminosity(mass, z);
        real r_agb =  AGB_radius(l_HeI, mass, mass_tot, z);
        if (r_ehg > r_agb){
            if (mass >= 12.0){
                r_ehg = r_agb;
                cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: R_ehg is set to r_agb"<<endl;
            }
            else {
                cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: R_AGB(L_HeI) < R_mHe, skipping blue loop?"<<endl;
            }
        }
        if (blue_phase_timescale(mass, mass_tot, z) < cnsts.safety(tiny) && abs(r_ehg-r_agb)> cnsts.safety(tiny)) {
            cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: t_bl <0, but r_ehg != r_agb)"<<endl;; 
        }
        
    }  
    return r_ehg;
}

real hertzsprung_gap::terminal_hertzsprung_gap_radius() {
    cerr<<"terminal_hertzsprung_gap_radius() without parameters is used. "<<endl;
  return terminal_hertzsprung_gap_radius(relative_mass, get_total_mass(), metalicity);
}


real hertzsprung_gap::terminal_hertzsprung_gap_luminosity() {

  return terminal_hertzsprung_gap_luminosity(relative_mass, metalicity);
}

//real hertzsprung_gap::base_giant_branch_luminosity() {
//
//  return terminal_hertzsprung_gap_luminosity(relative_mass, metalicity);
//}

//Eq.26
real hertzsprung_gap::hertzsprung_gap_luminosity(const real time,
					     const real mass, 
					     const real z) {

  real t_ms = main_sequence_time(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real tau = (time - t_ms)/(t_bgb - t_ms);
  PRC(time);PRC(t_ms);PRC(t_bgb);PRL(tau);

  real l_ehg = terminal_hertzsprung_gap_luminosity(mass, z);
  real l_tms = terminal_main_sequence_luminosity(mass, z);

  PRC(l_ehg);PRC(l_tms);PRC(tau);

  real l_hg = l_tms * pow(l_ehg/l_tms, tau);
  PRL(l_hg);

  return l_hg;
}

real hertzsprung_gap::hertzsprung_gap_luminosity(const real time) {

  return hertzsprung_gap_luminosity(time, relative_mass, metalicity);
}

real hertzsprung_gap::hertzsprung_gap_luminosity() {

  return hertzsprung_gap_luminosity(relative_age, relative_mass, metalicity);
}


//Eq.27
real hertzsprung_gap::hertzsprung_gap_radius(const real time,
					     const real mass, 
                         const real mass_tot, 
					     const real z) {

  real t_ms = main_sequence_time(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real tau = (time - t_ms)/(t_bgb - t_ms);
    
  real r_tms = terminal_main_sequence_radius(mass, z);
  real r_thg = terminal_hertzsprung_gap_radius(mass, mass_tot, z);
 
  real r_hg = r_tms * pow(r_thg/r_tms, tau);

  return r_hg;
}

real hertzsprung_gap::hertzsprung_gap_radius(const real time) {
 
  return hertzsprung_gap_radius(time, relative_mass, get_total_mass(), metalicity);
}

real hertzsprung_gap::hertzsprung_gap_radius() {

  return hertzsprung_gap_radius(relative_age, relative_mass, get_total_mass(), metalicity);
}

void hertzsprung_gap::evolve_core_mass(const real time,
				       const real mass,
				       const real z) {

  real mc_Hg = hertzsprung_gap_core_mass(time, mass, z);

    if(!update_core_and_envelope_mass(mc_Hg)) {
    cerr << "Update core mass failed in hertzsprung_gap()"<<endl;
  }
}

void hertzsprung_gap::evolve_core_mass() {

  evolve_core_mass(relative_age, relative_mass, metalicity);
}

//Eq.28
real hertzsprung_gap::terminal_hertzsprung_gap_core_mass(const real mass, 
							 const real z) {

  real m_core;
  real m_HeF = helium_flash_mass(z);
  real m_FGB = helium_ignition_mass(z);
  if (mass < m_HeF) {
    real l_bgb = base_giant_branch_luminosity(mass, z);
    m_core = FGB_core_mass_luminosity_relation(l_bgb, mass, z);
    PRL(m_core);
  }
  else if (mass >= m_FGB) {
    real mc_HeI = helium_ignition_core_mass(mass, z);//sect.5.3: Eq.67
    PRL(mc_HeI);
    m_core = mc_HeI;
  }
  else {
    real mc_bgb = base_giant_branch_core_mass(mass, z); //Eq.44    
    PRL(mc_bgb);
    m_core = mc_bgb;

  }
  return m_core;
}

//Eq.30
real hertzsprung_gap::hertzsprung_gap_core_mass(const real time, 
						const real mass,
						const real z) {

  real t_ms = main_sequence_time();
  real t_bgb = base_giant_branch_time(mass, z);
  real tau = (time - t_ms)/(t_bgb - t_ms);
  
  real m5_25 = pow(mass, 5.25);
  real mc_ehg = terminal_hertzsprung_gap_core_mass(mass, z);
  
  real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);
  real m_core = (tau + rho*(1-tau)) *  mc_ehg;
  
  if (m_core < core_mass) {
      cerr<<"On hertzsprung_gap::hertzsprung_gap_core_mass new m_core smaller than core_mass "<<endl;
  }
  m_core = max(m_core, core_mass);
  return m_core;
}

real hertzsprung_gap::helium_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        r_c = white_dwarf_radius(m_core, 0.);
    }
    return r_c;
}

real hertzsprung_gap::helium_core_radius(){
    return helium_core_radius(relative_mass, core_mass, metalicity);
}    

real hertzsprung_gap::small_envelope_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        cerr<<"hg:he_core_rad should return wd radius, new or old one?"<<endl;
        cerr<<"warning core_radius = 5* wd radius, though remnant has radius wd radius"<<endl; 
        r_c = 5.*white_dwarf_radius(m_core, 0.);
    }
    return r_c;
}
  
real hertzsprung_gap::small_envelope_core_radius(){
    return small_envelope_core_radius(relative_mass, core_mass, metalicity);
}    



real hertzsprung_gap::small_envelope_core_luminosity(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real l_c;
    if(mass > m_HeF){
    l_c = helium_star_luminosity_for_solar_metalicity(m_core);
    }
    else{
        l_c = 40.;
    }
    return l_c;
}
  
real hertzsprung_gap::small_envelope_core_luminosity(){
    return small_envelope_core_luminosity(relative_mass, core_mass, metalicity);
}    



