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
 
      // (GN Oct 25 2010) try to fix transition MS -> HG in merger MS + *
//      if (core_mass > 0) add_mass_to_accretor(1e-5, false);


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
        if (relative_mass < m_HeF){
            star_transformation_story(Helium_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
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
        adjust_age_after_mass_loss(mdot, true);
        envelope_mass -= mdot;
        if (relative_mass > get_total_mass()){
            update_relative_mass(get_total_mass());
        }
    }
    return this;
}


//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real hertzsprung_gap::mdot_limit(const real dt){
    real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    return mass_ratio_mdot_limit(mdot);
    
}


star* hertzsprung_gap::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = mdot_limit(dt);
      
      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;

        // (SPZ+GN: 27 Jul 2000)
        // non degenerate core < helium_dwarf_mass_limit always(!) become
        // white dwarfs
        //if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
        //    relative_mass < cnsts.parameters(
        //            upper_ZAMS_mass_for_degenerate_core)) {
        //  star_transformation_story(Helium_Dwarf);
        //  return dynamic_cast(star*, new white_dwarf(*this));
        //} 
        //else {
        //  star_transformation_story(Helium_Star);
        //  return dynamic_cast(star*, new helium_star(*this));
        //}
        
          real m_HeF = helium_flash_mass(metalicity);
          if (relative_mass < m_HeF){
              star_transformation_story(Helium_Dwarf);
              return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
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
        adjust_age_after_mass_loss(mdot, true);
        envelope_mass -= mdot;
        if (relative_mass > get_total_mass()){
            update_relative_mass(get_total_mass());
        }
    }
    
    adjust_donor_radius(mdot);
    return this;    
}



// add mass to accretor
// is a separate function (see single_star.C) because rejuvenation
real hertzsprung_gap::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

    if (mdot<0) {
        cerr << "hertzsprung_gap::add_mass_to_accretor(mdot=" << mdot 
        << ", dt=" << dt << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        cerr << "Action: put mdot to zero!" << endl;
        return 0;
    }
    
    bool update_age = false;

    if(hydrogen){
        //hydrogen accretion
        mdot = accretion_limit(mdot, dt);

	envelope_mass += mdot;
	accreted_mass += mdot;

	if (relative_mass < get_total_mass()){
	  update_age = true;
	  update_relative_mass(get_total_mass());
	}
	
        
        adjust_accretor_radius(mdot, dt);
        
    }
    else{
        //for the moment assume helium accretion
        // for the moment no adjust_accretor_radius
        
        mdot = accretion_limit_eddington(mdot, dt);
        core_mass += mdot;
        update_relative_mass(relative_mass + mdot);
        
	update_age = true;
    }


    if (update_age) {

        //adjust age part
        real mc_ehg = terminal_hertzsprung_gap_core_mass(relative_mass, metalicity);
        real m5_25 = pow(relative_mass, 5.25);            
        real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);
        real tau = (core_mass / mc_ehg - rho ) / (1. - rho);
        real t_ms = main_sequence_time();
        real t_bgb = base_giant_branch_time(relative_mass, metalicity);
        relative_age = t_ms + tau * (t_bgb - t_ms);
        last_update_age = t_ms;

        if (tau < 0.){
            real (single_star::*fptr)(const real, real) = &single_star::initial_hertzsprung_gap_core_mass;                   
            real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
            update_relative_mass(m_rel);
            last_update_age = main_sequence_time(relative_mass, metalicity);
            relative_age = last_update_age;     
            evolve_core_mass();
        }
        if (tau > 1.){
            real (single_star::*fptr)(const real, real) = &single_star::terminal_hertzsprung_gap_core_mass;        
            real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
            update_relative_mass(m_rel);
            last_update_age = main_sequence_time(relative_mass, metalicity);
            relative_age = next_update_age;            
            evolve_core_mass();
        }
    }
    set_spec_type(Accreting);
    return mdot;
}

#if 0
// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
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
      real t_bgb = hertzsprung_gap_time(m_rel_new, z_new); 
      real t_hg_new = t_bgb - t_ms_new;

//      real dtime = relative_age - t_ms_old;

      // (GN+SPZ May  4 1999) update last_update_age
      last_update_age = t_ms_new; 
   
//      relative_age = t_ms_new 
//                   + dtime*(t_hg_new/t_hg_old);
//
//      if (rejuvenate)
//           relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
//


    // rejuvenation based on keeping mc constant and finding relative age for new mass
    real mc_ehg = terminal_hertzsprung_gap_core_mass(relative_mass, metalicity);
    real m5_25 = pow(relative_mass, 5.25);    
    real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);

    real mc_ehg_new = terminal_hertzsprung_gap_core_mass(m_rel_new, z_new);
    real m5_25_new = pow(m_rel_new, 5.25);
    real rho_new = (1.586 + m5_25_new) / (2.434 + 1.02*m5_25_new);

    
    relative_age = t_ms_new + t_hg_new * ((1-rho)/(1-rho_new) * mc_ehg/mc_ehg_new * (relative_age - t_ms_old) / t_hg_old + 
                (rho*mc_ehg - rho_new * mc_ehg_new)/ (1-rho_new) / mc_ehg_new);
                
                
      relative_age = max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep)); 
    
      relative_age = min(relative_age, t_bgb);
      
      
        // next_update_age should not be reset here
        // next_update_age = t_ms_new + t_hg_new;
}
#endif

// Age adjustment especially for (wind) mass loss
// It is part of the single star evolution, 
// so it can include information from tracks
void hertzsprung_gap::adjust_age_after_mass_loss(const real mdot, const bool rejuvenate=true) {

    real m_rel_new;
    real m_tot_new = get_total_mass() - mdot;
    if (m_tot_new<relative_mass)
        m_rel_new = m_tot_new;
    else m_rel_new = relative_mass;
    
    real t_ms_old = main_sequence_time();
    real t_hg_old = hertzsprung_gap_time() - t_ms_old;
    
    real z_new = get_metalicity();
    real t_ms_new = main_sequence_time(m_rel_new, z_new);
    //For now, we keep metalicity constant (SPZ: 29 May 2001)
    real t_bgb = hertzsprung_gap_time(m_rel_new, z_new); 
    real t_hg_new = t_bgb - t_ms_new;
    
    real dtime = relative_age - t_ms_old;
    
    // (GN+SPZ May  4 1999) update last_update_age
    last_update_age = t_ms_new;
    //following HPT tracks only update the relative_age relative to the length of the phase

    relative_age = t_ms_new + dtime*(t_hg_new/t_hg_old);
//    if (rejuvenate){
//        real mdot_fr = -1.*mdot/m_tot_new;
//        real rejuvenation = (1-pow(mdot_fr,
//                                   cnsts.parameters(rejuvenation_exponent)));
//        relative_age *= rejuvenation; 
//    
//    }
//    

    relative_age = max(relative_age, 
                       last_update_age + cnsts.safety(minimum_timestep)); 
    relative_age = min(relative_age, t_bgb);

    // next_update_age should not be reset here
    // next_update_age = t_ms_new + t_hg_new;
}


// Adiabatic response function for hertzsprung_gap star.
// Polynomial fit to Hjellming and Webbink 1987 ApJ, 318, 804
real hertzsprung_gap::zeta_adiabatic() {

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


// Stellar Gyration radii squared for detmination of
// angular momentum.
// Implemented by (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

void hertzsprung_gap::update_wind_constant() {
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
    
    real dm_dj_v = 0;
    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x = min(1.0, (luminosity - 4000.0)/500.0);
        dm_dj = x * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5);
    }
    
    // Vink 2000, 2001
    // Massive stars, including multi scattering effects
    real dm_v = 0;
    if (luminosity > 4000 && metalicity > cnsts.parameters(solar_metalicity)/30. && metalicity < 3*cnsts.parameters(solar_metalicity)){
        real temp = temperature();
        real sigma;//electron scattering cross section
        if (temp >= 35000){
            sigma = 0.34;
        }
        else if(temp < 30000){
            sigma = 0.31;
        }
        else {
            sigma = 0.32;
        }
        real rad_acc = 7.66E-5 * sigma * luminosity / get_total_mass();
        real log_density = -14.94 + 0.85 * log10(metalicity/cnsts.parameters(solar_metalicity)) +3.1857*rad_acc; //Eq.23 Vink 2001
        real Tjump = (61.2 + 2.59*log_density)*1000; //Eq.15 Vink 2001
        real Tjump_low = (100. + 6. * log_density)*1000; //Eq.6 Vink 2000
        real arg_dm_v;
        real T_smooth = 1500.;
        real T_smooth_below = 1000.;
        real cnsts_dm_v[9];
        real cnsts_dm_v_above[] = {2.6, -6.697, 2.194, -1.313, -1.226, 0.933, 0, -10.92, 0.85};
        real cnsts_dm_v_below[] = {1.3, -6.688, 2.210, -1.339, -1.601, 1.07, 1.07, 0, 0.85};
        
        
        if (rad_acc >0.5 || temp > 50000) {
            // vink approaches LBV, stop? transition needed? possible for low metallicities
            dm_v = 0;
            dm_dj_v = dm_dj;
        }
        else {
            if (temp <= Tjump-T_smooth){
                // smooth out second instability jump
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = cnsts_dm_v_below[i_t];
                }
                
                if(temp <= Tjump_low - T_smooth_below){
                    cnsts_dm_v[0] = 0.7;
                    cnsts_dm_v[1] = -5.990;
                }
                else if (temp < Tjump_low + T_smooth_below){
                    real scale_T = (Tjump_low+T_smooth_below - temp)/(2*T_smooth_below);
                    cnsts_dm_v[0] = (1.-scale_T) * cnsts_dm_v_below[0] + scale_T * 0.7;
                    cnsts_dm_v[1] = (1.-scale_T) * cnsts_dm_v_below[1] + scale_T * -5.990;
                }
            }
            else if(temp > Tjump + T_smooth){
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = cnsts_dm_v_above[i_t];
                }
            }
            else {
                //smooth out first instability jump
                real scale_T = (Tjump+T_smooth - temp)/(2*T_smooth);
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = scale_T * cnsts_dm_v_below[i_t] 
                    + (1.-scale_T) * cnsts_dm_v_above[i_t];
                }
            }
            
            arg_dm_v  = cnsts_dm_v[1] 
            +cnsts_dm_v[2]*log10(luminosity/1.E5) 
            +cnsts_dm_v[3]*log10(get_total_mass()/ 30) 
            +cnsts_dm_v[4]*log10(cnsts_dm_v[0]/2.0)
            +cnsts_dm_v[5]*log10(temp/40000)
            +cnsts_dm_v[6]*log10(2.) // (T/20000) below Tjump
            +cnsts_dm_v[7]*pow(log10(temp/40000),2)
            +cnsts_dm_v[8]*log10(metalicity/cnsts.parameters(solar_metalicity));
            
            dm_v = pow(10, arg_dm_v);   
            dm_dj_v = dm_v;
            
            
            if (temp < 8000){
                // line driven winds no longer efficient
                // see Achmad et al 1997
                dm_v = dm_v * 200. / (8200.-temp);
                dm_dj_v = max(max(dm_v, dm_dj), 0.);
            }
        }
    }   
    else
        dm_dj_v = dm_dj;
    
    
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    
//    //Schroder & Cuntz
//    // cool GB like stars
//    real neta_sc = 8.E-14; 
//    real surface_gravity = get_total_mass()/ pow(radius, 2);
//    real dm_sc = neta_sc * 4.E-13 * radius * luminosity / get_total_mass() 
//           * pow(temperature()/4000, 3.5) * (1 + 1./(4300*surface_gravity));
    
    
    //based on Nugis & Lamers
    // eq 8.4 in Gijs' thesis Chapter 8
    //Reduced WR-like mass loss for small H-envelope mass
    real mu = (get_total_mass()-core_mass)/get_total_mass() * min(5.0,max(1.2, pow(luminosity/7.E4,-0.5)));
    real dm_wr = 0;
    if ( mu < 1.){
        //factor (1.-mu) should be checked e.g. with resulting # BH in binaries
        dm_wr = 1.38E-08 * pow(get_total_mass(), 2.87) * (1.-mu);
    }
    
    //LBV
    real dm_lbv = 0;
    real x_lbv = 1.0E-5*radius*sqrt(luminosity);
    if(luminosity > 6.0E5 && x_lbv > 1.0) {
        dm_lbv = 0.1 * pow(x_lbv-1.0, 3)*(luminosity/6.0E5-1.0);
    }
    
    wind_constant = max(max(max(dm_wr, dm_dj_v), dm_r), 0.0) + dm_lbv;
        
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
          
          // if no envelope make transition to remnants
          // just as a procedure: reduce_mass with 1
          if (envelope_mass <= 0){
              reduce_mass(1.);
              return;
          }
      }
      else {
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
      stellar_wind(dt);
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
//        if (blue_phase_timescale(mass, mass_tot, z) < cnsts.safety(tiny) && abs(r_ehg-r_agb)> cnsts.safety(tiny)) {
//            cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: t_bl <0, but r_ehg != r_agb)"<<endl;; 
//        }
        
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

  real l_ehg = terminal_hertzsprung_gap_luminosity(mass, z);
  real l_tms = terminal_main_sequence_luminosity(mass, z);

  real l_hg = l_tms * pow(l_ehg/l_tms, tau);
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
    
  real r_tms = terminal_main_sequence_radius(mass_tot, z);
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
				       const real z, const real m_core_old) {

  real mc_Hg = hertzsprung_gap_core_mass(time, mass, z, m_core_old);

    if(!update_core_and_envelope_mass(mc_Hg)) {
    cerr << "Update core mass failed in hertzsprung_gap()"<<endl;
  }
}


void hertzsprung_gap::evolve_core_mass() {
    
    evolve_core_mass(relative_age, relative_mass, metalicity, core_mass);
}


// Eq.4 (base_giant_branch_time(mass, z);
// Supersedes ::hertzsprung_gap_time(mass, t_ms) in single star
real hertzsprung_gap::hertzsprung_gap_time(const real mass, const real z) {
    
    real t_Hg = base_giant_branch_time(mass, z);
    return t_Hg;
}

real hertzsprung_gap::hertzsprung_gap_time() {
    
    return hertzsprung_gap_time(get_relative_mass(), get_metalicity());
}


//Eq.30
real hertzsprung_gap::hertzsprung_gap_core_mass(const real time, 
						const real mass,
						const real z, const real m_core_old) {

    real t_ms = main_sequence_time(mass, z);
    real t_bgb = base_giant_branch_time(mass, z);
    real tau = (time - t_ms)/(t_bgb - t_ms);

    real mc_ehg = terminal_hertzsprung_gap_core_mass(mass, z);
    real mc_ihg = initial_hertzsprung_gap_core_mass(mass, z);  
    real m_core = mc_ihg + (mc_ehg - mc_ihg) * tau;
        
    // according to HPT this is important in case of mass loss 
    m_core = max(m_core, m_core_old);   
  
  return m_core;
}

real hertzsprung_gap::helium_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        // due to small nucleair burning layer 
        // r_c > white_dwarf_radius
        r_c = 5.* white_dwarf_radius(m_core, 10000.);
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
        r_c = white_dwarf_radius(m_core, 10000.);
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



