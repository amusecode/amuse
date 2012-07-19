//
// sub_giant.C
//

#include "sub_giant.h"
#include "hertzsprung_gap.h"
#include "main_sequence.h"

// ANSI C++ first creates the base class before the dreived classes are
// created. 
 
sub_giant::sub_giant(hertzsprung_gap & h) : single_star(h) {

  delete &h; 

  last_update_age = next_update_age;
  adjust_next_update_age();

  instantaneous_element();
  evolve_core_mass();
  small_envelope_perturbation();   
  update();

  post_constructor();
}

// (GN+ SilT Feb 2011): new constructor for merger products of
// main sequence plus WD/He stars
sub_giant::sub_giant(main_sequence & m) : single_star(m) {

  delete &m; 

  last_update_age = next_update_age;

  // Proper adding of core mass
  if (is_binary_component()) {

    if (get_companion()->get_core_mass() > 0)
      add_mass_to_accretor(get_companion()->get_core_mass(), false);

    // this should not happen....(or hardly) for WD
    // but helium stars have He envelope and CO core....
    if (get_companion()->get_envelope_mass() > 0)
      add_mass_to_accretor(get_companion()->get_envelope_mass(), get_companion()->hydrogen_envelope_star());
  }

  adjust_next_update_age();

  instantaneous_element();
  evolve_core_mass();
  small_envelope_perturbation();   
  update();

  post_constructor();
}



#if 0
void sub_giant::adjust_initial_star() {

  if(relative_age<=0) {
    real t_ms = main_sequence_time();
    relative_age = max(t_ms + hertzsprung_gap_time(t_ms), relative_age);
  }
}
#endif



//		general mass transfer utilities.
// Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real sub_giant::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

    if (mdot<0) {
        cerr << "sub_giant::add_mass_to_accretor(mdot=" << mdot 
        << ", dt=" << dt << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        return 0;
    }
    
    bool update_age = false;
    
    if(hydrogen){
        //hydrogen accretion
        mdot = accretion_limit(mdot, dt);
        
        envelope_mass += mdot;
        accreted_mass += mdot;
        
        // For now, rejuvenation of SG, CHeB, AGB or He giant accretor   
	// only if mtot > relative_mass
	if (relative_mass<get_total_mass())  {

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
        real m_HeF = helium_flash_mass(metalicity);
        real t_bgb = base_giant_branch_time(relative_mass, metalicity);

    	 // (GN Oct 26 2010) for mergers: core_mass can be > max degenerate He core --> jump to next block
        if (relative_mass < m_HeF  && core_mass <  helium_ignition_core_mass(0.5, metalicity)){
            real l_bgb = base_giant_branch_luminosity(relative_mass, metalicity);
            real A_H = sub_giant_Ah_estimator(relative_mass);                    
            relative_age = determine_age(core_mass, relative_mass, metalicity, A_H, t_bgb, l_bgb);
            last_update_age = t_bgb;
            
            if(relative_age < last_update_age){
                real (single_star::*fptr)(const real, real) = &single_star::terminal_hertzsprung_gap_core_mass;        
                
                real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
                update_relative_mass(m_rel);
                last_update_age = base_giant_branch_time(relative_mass, metalicity);
                relative_age = last_update_age;   
                evolve_core_mass();
           }
            if(relative_age > next_update_age){
                real (single_star::*fptr)(const real, real) = &single_star::helium_ignition_core_mass;        
                
                real xmin = cnsts.parameters(minimum_main_sequence);
                real xmax = m_HeF;
                real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity, xmin, xmax);     
                update_relative_mass(m_rel);
                last_update_age = base_giant_branch_time(relative_mass, metalicity);
                relative_age = next_update_age;   
                evolve_core_mass();
            }
        }
        else {

            // (GN Oct 26 2010) see previous block, rel_mass can be < m_HeF which leads to trouble....
            if (relative_mass < m_HeF) relative_mass = m_HeF;

            real mc_bgb = base_giant_branch_core_mass(relative_mass, metalicity);//Eq.44
            real mc_HeI = helium_ignition_core_mass(relative_mass, metalicity);
            real t_HeI = helium_ignition_time(relative_mass, metalicity);
            real tau = (core_mass -mc_bgb) / (mc_HeI - mc_bgb);
            relative_age = t_bgb + tau * (t_HeI - t_bgb);
            last_update_age = t_bgb;

            if (tau < 0.){
                real (single_star::*fptr)(const real, real) = &single_star::terminal_hertzsprung_gap_core_mass;        
                real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
                update_relative_mass(m_rel);
                last_update_age = base_giant_branch_time(relative_mass, metalicity);		
                relative_age = last_update_age;
                evolve_core_mass();

            }
            if (tau > 1.){
                    real m_FGB = helium_ignition_mass(metalicity);
                    real mcore_max =  helium_ignition_core_mass(m_FGB, metalicity); 
                    real m_rel = m_FGB;
                    // (GN Oct 27 2010) for mergers core_mass can be larger than max core mass of class
                    // dirty fix: reduce core mass and conserve total mass
		    //PRC(relative_mass);PRC(core_mass);PRL(mcore_max);
                    if (core_mass > mcore_max) {                       
                        real m_tot = get_total_mass();
                        core_mass = mcore_max;
                        envelope_mass = m_tot - core_mass;
                    } else {
		    
		      real (single_star::*fptr)(const real, real) = &single_star::helium_ignition_core_mass;        
		      real xmin = m_HeF;
		      real xmax = helium_ignition_mass(metalicity);// m_FGB < m_rel not possible for gb star
		      m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity, xmin, xmax);     

		    }

                update_relative_mass(m_rel);
                last_update_age = base_giant_branch_time(relative_mass, metalicity);
                relative_age = next_update_age;
                evolve_core_mass();
            }
        } 
    }
    set_spec_type(Accreting);
    return mdot;
}



star* sub_giant::reduce_mass(const real mdot) {

    if (envelope_mass<=mdot) {
        envelope_mass = 0;
	
	 
//        // (SPZ+GN: 27 Jul 2000)
//        // non degenerate core < helium_dwarf_mass_limit always(!) become
//        // white dwarfs
//        if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
//             relative_mass < cnsts.parameters(
//		             upper_ZAMS_mass_for_degenerate_core)) {
//            star_transformation_story(Helium_Dwarf);
//            return dynamic_cast(star*, new white_dwarf(*this));
//        } 
//        else {
//           star_transformation_story(Helium_Star);
//           return dynamic_cast(star*, new helium_star(*this));
//        }
//        
        real m_HeF = helium_flash_mass(metalicity);
        if (relative_mass< m_HeF || core_mass < cnsts.parameters(minimum_helium_star)){
            star_transformation_story(Helium_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
        }
        else {
            star_transformation_story(Helium_Star);
            return dynamic_cast(star*, new helium_star(*this));
            
        }
    }
    else
        envelope_mass -= mdot;
    return this;
}


//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real sub_giant::mdot_limit(const real dt){
    real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    return mass_ratio_mdot_limit(mdot);
    
}

star* sub_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = mdot_limit(dt);
      
      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;

        // (SPZ+GN: 27 Jul 2000)
        // non degenerate core < helium_dwarf_mass_limit always(!) become
        // white dwarfs
        // if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
        //     relative_mass < cnsts.parameters(
        //                 upper_ZAMS_mass_for_degenerate_core)) {
        //   star_transformation_story(Helium_Dwarf);
        //   return dynamic_cast(star*, new white_dwarf(*this));
        // } 
        // else {
        //   star_transformation_story(Helium_Star);
        //   return dynamic_cast(star*, new helium_star(*this));
        // }

          real m_HeF = helium_flash_mass(metalicity);
          if (relative_mass< m_HeF || core_mass < cnsts.parameters(minimum_helium_star)){
              star_transformation_story(Helium_Dwarf);
              return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
          }
          else {
              star_transformation_story(Helium_Star);
              return dynamic_cast(star*, new helium_star(*this));
          }
      }

      else{// (GN+SPZ Apr 29 1999)
          adjust_donor_radius(mdot);

          envelope_mass -= mdot;
      }
      return this;
   }


// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void sub_giant::adjust_accretor_age(const real mdot,
				    const bool rejuvenate=true) {
     cerr<<"sub_giant::adjust_accretor_age is currently not used"<<endl;
     real tend_hg_old = hertzsprung_gap_time();
     real dt_bgb_old = helium_ignition_time() - tend_hg_old ;

     real m_tot_new = get_total_mass() + mdot;
     real m_rel_new = max(m_tot_new, relative_mass);

     //For now, we keep metalicity constant (SPZ: 29 May 2001)
     real z_new = metalicity;
     real tend_hg_new = hertzsprung_gap_time(m_rel_new, z_new);
     real t_HeI = helium_ignition_time(m_rel_new, z_new);
     real dt_bgb_new = t_HeI - tend_hg_new;

     real dtime = relative_age - tend_hg_old;

     // For relative_mass > helium_ignition_mass(z) ~ 13Msolar
     // sub_giants can not exist. (SPZ+GN:10 Oct 1998, SilT:7 Oct 2009)
     
     if (dt_bgb_new>0) {

       // (GN+SPZ May  4 1999) update last_update_age
       last_update_age = tend_hg_new;
       relative_age = tend_hg_new 
                    + dtime*(dt_bgb_new/dt_bgb_old);

       if (rejuvenate) {
           relative_age *= rejuvenation_fraction(mdot/m_tot_new);
       }

       if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
             cerr<<"In sub_giant::adjust_accretor_age relative age updated on SG, but < last_update_age"<<endl;
       }
       relative_age = max(relative_age, 
			  last_update_age  + cnsts.safety(minimum_timestep));
       relative_age = min(relative_age, t_HeI);
     }
     else {
       // Relative_age should be set to the predicted next update age.
       // Instead use tend_hg_new, which is end point of HG.
       //       relative_age = next_update_age;
       relative_age = tend_hg_new;
     }

     // next_update_age should not be reset here
     // next_update_age = tend_hg_new+t_bgb_new;

}



void sub_giant::adjust_next_update_age() {

  real t_bgb = base_giant_branch_time(relative_mass, metalicity);

//  if( relative_age < t_bgb - cnsts.safety(tiny)) {
//    cerr << "WARNING: relative_age != t_Hg in sub_giant"<<endl;
//    relative_age = t_bgb;
//  }

  real t_HeI = helium_ignition_time();
  next_update_age = t_HeI;
}

void sub_giant::detect_spectral_features() {

      single_star::detect_spectral_features();


      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;

}

real sub_giant::gyration_radius_sq() {

  return cnsts.parameters(convective_star_gyration_radius_sq); 
}


real sub_giant::zeta_adiabatic() {

// (GN+SPZ Apr 28 1999) fit from Lev Yungelson private communication
// for giants with not too convective envelope = radiative envelope

  real r_dconv = 2.4*pow(relative_mass,1.56);
  if (relative_mass > 10 )
    r_dconv = 5.24*pow(relative_mass,1.32);
  else if (relative_mass > 5)
    r_dconv = 1.33*pow(relative_mass,1.93);
    
    //(SilT Sep 1 2010) Need factor 1.5 with new HPT tracks in order to get
   // stable mass transfer on early giant branch  
    r_dconv = 1.5* r_dconv;
  if (radius < r_dconv) {

    return 4;
  }
  else {
//   		Hjellming and Webbink 1987 ApJ, 318, 804
    real x = core_mass/get_total_mass();
    real A = -0.220823;
    real B = -2.84699;
    real C = 32.0344;
    real D = -75.6863;
    real E = 57.8109;

    return A + x*(B + x*(C + x*(D + x*E)));

  }
}

// Values of zeta are changed (SPZ+GN:28 Sep 1998)
real sub_giant::zeta_thermal() {

  real z;
  if (low_mass_star())
    z = 0;
  else 
    z = 0; // (GN+SPZ Apr 28 1999) radius determined by core only (was -1) 

  return z;
}

void sub_giant::update_wind_constant() {
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
  else { // no wind for low mass ms stars
// (GN Apr 16 1999) 0.2 Msun loss for degenerate core stars
//    real t_ms = main_sequence_time();

    // (SPZ+GN: 26 Jul 2000) see Nelemans, YPZV 2000 (A&A Submitted)
    // wind_constant = 0.2; is a slight improvement on
    // Though a minimum of 0.0 is somewhat on the low side.
    wind_constant = max(0., (2.5 - relative_mass)/7.5);

// (GN+SPZ May  4 1999) not needed: single_star::stellar_wind changed
//                  /(1 - 
//		   pow((t_ms + hertzsprung_gap_time(t_ms))
//                      /next_update_age,
//                     cnsts.parameters(massive_star_mass_loss_law)));
  }
  
  wind_constant = max(wind_constant, 0.0);
#endif
    
    // wind_constant is in solar masses per year
    // Should be updated after mass accretion
    // (ST: 17 Sep 2009)
    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x_dj = min(1.0, (luminosity - 4000.0) / 500.0);
        dm_dj = x_dj * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5);
    }
    
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();

//    //Schroder & Cuntz
//    // cool GB like star
//    real neta_sc = 8.E-14; 
//    real surface_gravity = pow(radius, 2) / get_total_mass();
//    real dm_sc = neta_sc * 4.E-13 * radius * luminosity / get_total_mass() 
//    * pow(temperature()/4000, 3.5) * (1 + 1./(4300*surface_gravity));
    
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
    
    wind_constant = max(max(max(dm_wr, dm_dj), dm_r), 0.0)+dm_lbv;
    
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Metalicity dependencies: Hurley, Pols & tout, 2000, MNRAS 315543 +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void sub_giant::instantaneous_element() {
  luminosity =  FGB_luminosity_core_mass_relation(relative_age,
						  relative_mass, 
						  metalicity);
  radius = giant_branch_radius(luminosity, get_total_mass(), metalicity);

}

// Evolve a main_sequence star upto time argument according to
// the new 2000 models.
void sub_giant::evolve_element(const real end_time) {

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
          //		sub_giant lifetime exceeded. Transform star into
          //		horizontal branch star.
          star_transformation_story(Horizontal_Branch);
          new horizontal_branch(*this);
          return;
      }

      update();
      stellar_wind(dt);
}


real sub_giant::get_evolve_timestep() {
    
    real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                        next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   
        
    //extra safety measure
    // when L and R increase rapidly, so will mdot
    real A_H = sub_giant_Ah_estimator(relative_mass);
    real l_x = FGB_x_luminosity(relative_mass, metalicity);
    real l_bgb = base_giant_branch_luminosity(relative_mass, metalicity);
    real r_bgb = giant_branch_radius(l_bgb, relative_mass, metalicity);

    real dt_mdot = timestep;
    if (relative_mass < helium_flash_mass(metalicity)){
        if (luminosity < l_x){
            real p = sub_giant_p_parameter(relative_mass, metalicity);
            dt_mdot = core_mass / ( p * luminosity * A_H) * r_bgb/ radius * 1.0;
        }
        else{
            real q = sub_giant_q_parameter(relative_mass, metalicity);
            dt_mdot = core_mass / ( q * luminosity * A_H)  * r_bgb /radius * 1.0;
        }
    }
        
    return max(min(timestep, dt_mdot), cnsts.safety(minimum_timestep));
    
}



void sub_giant::evolve_core_mass(const real time,
				 const real mass,
				 const real z) {

  
  real mc_sg = sub_giant_core_mass(time, mass, z);
  if(!update_core_and_envelope_mass(mc_sg)) {
      cerr << "Update core mass failed in sub_giant()"<<endl;
  }
}



void sub_giant::evolve_core_mass() {

  evolve_core_mass(relative_age, relative_mass, metalicity);
}

real sub_giant::sub_giant_core_mass(const real time,
                                    const real mass,
                                    const real z) {
    
    real m_core;
    real t_bgb = base_giant_branch_time(mass, z);

    if (mass <= helium_flash_mass(z)){
        real l_bgb = base_giant_branch_luminosity(mass, z);
        real A_H = sub_giant_Ah_estimator(mass);
  
        m_core = determine_core_mass(time, mass, z, 
 				     A_H, t_bgb, l_bgb);
    }
    else{
        real mc_bgb = base_giant_branch_core_mass(mass, z);//Eq.44
        real mc_HeI = helium_ignition_core_mass(mass, z);
        real t_HeI = helium_ignition_time(mass,z);
        real tau = (time- t_bgb)/(t_HeI-t_bgb);

	if (tau > 1) tau = 1.; // Safety

        m_core = mc_bgb + (mc_HeI - mc_bgb)* tau;
    }
    return m_core;
}

real sub_giant::helium_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        // due to small nucleair burning layer 
        // r_c > white_dwarf_radius
        r_c = 5.*white_dwarf_radius(m_core, 10000.);
    }
    return r_c;
}
real sub_giant::helium_core_radius(){
    return helium_core_radius(relative_mass, core_mass, metalicity);
}

real sub_giant::small_envelope_core_radius(const real mass, const real m_core, const real z){
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
real sub_giant::small_envelope_core_radius(){
    return small_envelope_core_radius(relative_mass, core_mass, metalicity);
}


real sub_giant::small_envelope_core_luminosity(const real mass, const real m_core, const real z){
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
  
real sub_giant::small_envelope_core_luminosity(){
    return small_envelope_core_luminosity(relative_mass, core_mass, metalicity);
}    





