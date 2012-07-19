//
// horizontal_branch.C
//

#include "horizontal_branch.h"
#include "sub_giant.h"
#include "hertzsprung_gap.h"
#include "main_sequence.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.
//

horizontal_branch::horizontal_branch(sub_giant & g) : single_star(g) {

    delete &g;
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    last_update_age = next_update_age;

    adjust_next_update_age();
    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();   
    
    update();
    post_constructor();
    
}

horizontal_branch::horizontal_branch(hertzsprung_gap & h) : single_star(h) {

    delete &h;

// (GN+SPZ May  4 1999) last update age is time of previous type change
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
horizontal_branch::horizontal_branch(main_sequence & m) : single_star(m) {

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

// possible track hydrogen accreting helium star can turn into horizontal branch star
//not implemented currently
//horizontal_branch::horizontal_branch(helium_star & h) : single_star(h) {
//
//    delete &h;
//    
//    real t_hems = helium_main_sequence_time_for_solar_metalicity(core_mass);
//    relative_age = relative_age / t_hems;  
//    PRL(relative_age);
//    
//   
//    real m_rel_min = max(helium_flash_mass(metalicity), base_AGB_relative_mass(core_mass, metalicity));
//    real m_rel_max = cnsts.parameters(maximum_main_sequence);//different than HPT, that calculate
//    //for which m_rel holds core_mass = helium_ignition_core_mass(relative_mass, metalicity)
//    
//    PRL(helium_ignition_core_mass(helium_flash_mass(metalicity), metalicity));
//    
//    
//    real (single_star::*fptr)(const real, real) = &single_star::core_helium_burning_core_mass_from_helium_star;        
//    real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity, m_rel_min, m_rel_max);             
//    update_relative_mass(m_rel);
//    
//  
//    real t_HeI = helium_ignition_time(relative_mass, metalicity);
//    real t_eagb = t_HeI + core_helium_burning_timescale(relative_mass, metalicity); 
//
//    relative_age *=  t_eagb;  
//    adjust_next_update_age();
//    last_update_age = t_HeI;
//    
//    instantaneous_element();
//    evolve_core_mass();
//    small_envelope_perturbation();   
//
//    update();
//    post_constructor();
//      
//}


//		general mass transfer utilities.
// Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real horizontal_branch::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

    if (mdot<0) {
        cerr << "horizontal_branch::add_mass_to_accretor(mdot=" << mdot 
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
        real mc_bagb = base_AGB_core_mass(relative_mass, metalicity);
        real mc_HeI = helium_ignition_core_mass(relative_mass, metalicity);
        
        real t_HeI = helium_ignition_time(relative_mass, metalicity); 
        real t_He  = core_helium_burning_timescale(relative_mass, metalicity);
        real tau = (core_mass - mc_HeI) / (mc_bagb - mc_HeI);
        relative_age = t_HeI + tau * t_He;
        last_update_age = t_HeI;

        if (tau < 0.){            
            real xmin, xmax;
            real m_HeF = helium_flash_mass(metalicity);
            if (relative_mass < m_HeF){
                xmin = cnsts.parameters(minimum_main_sequence);
                xmax = m_HeF;
            }
            else{ 
                xmin = m_HeF;
                xmax = cnsts.parameters(maximum_main_sequence);
            }
            
            real (single_star::*fptr)(const real, real) = &single_star::helium_ignition_core_mass; 
            real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity, xmin, xmax);     
            update_relative_mass(m_rel);
            last_update_age = helium_ignition_time(relative_mass, metalicity); 
            relative_age = last_update_age;     
            evolve_core_mass();
        }
        if (tau > 1.){
            // in principle this is not correct because the helium core has reached it's max in this phase
            // because of accretion, but what happens with the co core?
            update_relative_mass(base_AGB_relative_mass(core_mass, metalicity));
            relative_age = next_update_age;
            last_update_age = helium_ignition_time(relative_mass, metalicity); 
            evolve_core_mass();
        }
    }
    set_spec_type(Accreting);
    return mdot;
}



star* horizontal_branch::reduce_mass(const real mdot) {

      if (envelope_mass<=mdot) {
         envelope_mass = 0;
          
//         if (get_total_mass() < cnsts.parameters(minimum_helium_star)) {
//              // horizontal branch star will not continue core helium burning.
//              cerr<<"Warning: not homogeneous WD"<<endl;
//              star_transformation_story(Helium_Dwarf);
//              return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
//          }
//         else{
             star_transformation_story(Helium_Star);
             return dynamic_cast(star*, new helium_star(*this));
//         }
      }

      else 
          envelope_mass -= mdot;
      
    return this;
   }
   
   
//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real horizontal_branch::mdot_limit(const real dt){
    real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    return mass_ratio_mdot_limit(mdot);
    
}
   

star* horizontal_branch::subtrac_mass_from_donor(const real dt,
						 real& mdot) {

      mdot = mdot_limit(dt);
      
      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;
//         if (relative_mass < cnsts.parameters(minimum_helium_star)) {
//              // horizontal branch star will not continue core helium burning.
//              cerr<<"Warning: not homogeneous WD"<<endl;
//              star_transformation_story(Helium_Dwarf);
//              return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
//         }
//         else{   
             star_transformation_story(Helium_Star);
             return dynamic_cast(star*, new helium_star(*this));
//         }
      }
      else
          envelope_mass -= mdot;
      return this;
   }

// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void horizontal_branch::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {
     cerr<<"sub_giant::adjust_accretor_age is currently not used"<<endl;

      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_HeI_old = helium_ignition_time();
      real dt_cHe_old = core_helium_burning_timescale();

      real z_new = metalicity;
      real t_HeI_new = helium_ignition_time(m_rel_new, z_new);
      real dt_cHe_new = core_helium_burning_timescale(m_rel_new, z_new);

      real dtime = relative_age - t_HeI_old;

      last_update_age = t_HeI_new;
      relative_age = t_HeI_new
                   + dtime*(dt_cHe_new/dt_cHe_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new);
     
      if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
         cerr<<"In horizontal_branch::adjust_accretor_age relative age updated on HB, but < last_update_age"<<endl;
      }
    
      relative_age = max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep));
      relative_age = min(relative_age, t_HeI_new + dt_cHe_new);

}

real horizontal_branch::zeta_adiabatic() {	
//    return 15;	

    // (SilT 25 October 2010) new tracks require new zeta
    // definition of horizontal branch changed
    real m_FGB = helium_ignition_mass(metalicity);
    real m_HeF = helium_flash_mass(metalicity);
    if (relative_mass < m_HeF){ //low mass stars   
        return 4;    
    }
    else if (relative_mass < m_FGB){ //intermediate mass stars
        real t_HeI = helium_ignition_time(relative_mass, metalicity); //#eq.43    
        real t_He = core_helium_burning_timescale(relative_mass, metalicity);
        real tau = (relative_age - t_HeI)/t_He;
        real tau_x = relative_age_at_start_of_blue_phase(relative_mass, metalicity);
        
        if (tau < tau_x){ //decent along GB
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
            if (radius < r_dconv)
                return 4;
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
        else //blue loop phase 
          return 4;            
    }
    else{// high mass stars
        real t_HeI = helium_ignition_time(relative_mass, metalicity); //#eq.43    
        real t_He = core_helium_burning_timescale(relative_mass, metalicity);
        real tau = (relative_age - t_HeI)/t_He;
        real tau_y = relative_age_at_end_of_blue_phase(relative_mass, metalicity);

        if (tau < tau_y){ //blue phase before reaching the giant branch
            return 4;
        }
        else {// red (super)giant phase
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
}

real horizontal_branch::zeta_thermal() {

//      return 15;

    // (SilT 25 October 2010) new tracks require new zeta
    // definition of horizontal branch changed
    real m_FGB = helium_ignition_mass(metalicity);
    real m_HeF = helium_flash_mass(metalicity);
    if (relative_mass < m_HeF){ // low mass stars   
        return 4;
    }
    else if (relative_mass < m_FGB){//intermediate mass stars
        real t_HeI = helium_ignition_time(relative_mass, metalicity); //#eq.43    
        real t_He = core_helium_burning_timescale(relative_mass, metalicity);
        real tau = (relative_age - t_HeI)/t_He;
        real tau_x = relative_age_at_start_of_blue_phase(relative_mass, metalicity);
        
        if (tau < tau_x) //decent along GB
          return 0;            
        else //blue loop phase 
            return 4;
    }
    else{//high mass stars
        real t_HeI = helium_ignition_time(relative_mass, metalicity); //#eq.43    
        real t_He = core_helium_burning_timescale(relative_mass, metalicity);
        real tau = (relative_age - t_HeI)/t_He;
        real tau_y = relative_age_at_end_of_blue_phase(relative_mass, metalicity);

        if (tau < tau_y){ //blue phase before reaching the giant branch
            return -2;
        }
        else {// red (super)giant phase
            return 0;
        }        
    }
}

void horizontal_branch::adjust_next_update_age() {
  real t_HeI = helium_ignition_time();
  real dt_cHe = core_helium_burning_timescale();
    
//  if(relative_age < t_HeI - cnsts.safety(tiny)) {
//    cerr << "WARNING: relative_age != t_HeI in horizontal_branch"<<endl;
//    relative_age = t_HeI;
//  }
  next_update_age = t_HeI + dt_cHe;
}

real horizontal_branch::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

void horizontal_branch::update_wind_constant() {
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
  else { // 5% of initial envelope

    wind_constant = 0.05*(relative_mass - final_core_mass());
  }
  
  wind_constant = max(wind_constant, 0.0);
#endif
    // wind_constant is in solar masses per year
    // Should be updated after mass accretion
    // (ST: 17 Sep 2009)
    
    
    
    real dm_dj_v = 0;
    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    // devided by two to get reasonable single wolf rayet stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x = min(1.0, (luminosity - 4000.0)/500.0);
        dm_dj = x * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5)/2.;
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
            //vink approaches LBV, stop? transition needed? possible for low metallicities
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
    
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    
//    //Schroder & Cuntz
//    // cool GB like stars
//    real neta_sc = 8.E-14; 
//    real surface_gravity = pow(radius, 2) / get_total_mass();
//    real dm_sc = neta_sc * 4.E-13 * radius * luminosity / get_total_mass() 
//    * pow(temperature()/4000, 3.5) * (1 + 1./(4300*surface_gravity));
    
    
    //based on Nugis & Lamers
    // eq 8.4 in Gijs' thesis Chapter 8
    //Reduced WR-like mass loss for small H-envelope mass
    //real mu = (get_total_mass()-core_mass)/get_total_mass() * min(5.0,max(1.2, pow(luminosity/7.E4,-0.5)));
    //PRL(mu);
    real dm_wr = 0;
    //if ( mu < 1.){
    //    //factor (1.-mu) should be checked e.g. with resulting # BH in binaries
    //    dm_wr = 1.38E-08 * pow(get_total_mass(), 2.87) * (1.-mu);
    //}
    
    
    
    //LBV
    real dm_lbv = 0;
    real x_lbv = 1.0E-5*radius*sqrt(luminosity);
    if(luminosity > 6.0E5 && x_lbv > 1.0) {
        dm_lbv = 0.1 * pow(x_lbv-1.0, 3)*(luminosity/6.0E5-1.0);
    }

    wind_constant = max(max(max(dm_wr, dm_dj_v), dm_r), 0.0) + dm_lbv;

    //PRC(luminosity);PRC(dm_wr);PRC(dm_dj_v);PRC(dm_r);PRC(dm_lbv);PRL(wind_constant);
    
}


void horizontal_branch::instantaneous_element() {

  luminosity = core_helium_burning_luminosity(relative_age,
					      relative_mass,
                          get_total_mass(),
					      metalicity);
  radius = core_helium_burning_radius(relative_age,
				      relative_mass,
                      get_total_mass(),
				      metalicity,
                      luminosity);
    // don't do:
    //effective_radius = radius
    //because of small_envelope_perturbation

}


// Post-evolution stellar update.
// Makes sure age and radius are updated.
void horizontal_branch::update() {

  // New core mass determination occurs in ::evolve_element.
  // (SPZ+GN:09/1998)
  // real m_tot = get_total_mass();
  // core_mass = helium_core_mass();
  // envelope_mass = m_tot - core_mass;

  core_radius = helium_core_radius();

  // (GN+SPZ Apr 28 1999)
  // effective_radius can be larger than  radius
  // (SPZ:  8 Jul 2001) 
  // except for horzontal branch stars
  //effective_radius = max(radius,effective_radius);
  effective_radius = radius;

  // last update age is set after stellar expansion timescale is set.
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

  detect_spectral_features();
}

// Evolve a horizontal branch star upto time argument according to
// the new 2000 models.
void horizontal_branch::evolve_element(const real end_time) {
    
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
        // Core mass remains unchanged: second dredge-up    
        // accounted for in Super_Giant (SPZ+GN: 27 Jul 2000)
        star_transformation_story(Super_Giant);
        new super_giant(*this);
         return;
      }
      update();
      stellar_wind(dt);
    
}

void horizontal_branch::evolve_core_mass(const real time,
					 const real mass,
					 const real z) {


  real mc_hb = core_helium_burning_core_mass(time, mass, z); //Eq.67
  if(!update_core_and_envelope_mass(mc_hb)) {
    cerr << "Update core mass failed in horizontal_branch()"<<endl;
  }
}

void horizontal_branch::evolve_core_mass() {

  evolve_core_mass(relative_age, relative_mass, metalicity);
}


real horizontal_branch::helium_core_radius(const real time, 
                                                   const real mass, const real m_core, const real z){
    real t_HeI = helium_ignition_time(mass, z); 
    real t_He = core_helium_burning_timescale(mass, z);
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(m_core);
    real time_scaled = (time-t_HeI) / t_He * t_Hems;
    
    real r_c = helium_main_sequence_radius(time_scaled, m_core, m_core);
    return r_c;
}
real horizontal_branch::helium_core_radius(){
    return helium_core_radius(relative_age, relative_mass, core_mass, metalicity);
}



real horizontal_branch::small_envelope_core_radius(const real time, 
                                           const real mass, const real m_core, const real z){
    real t_HeI = helium_ignition_time(mass, z); 
    real t_He = core_helium_burning_timescale(mass, z);
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(m_core);
    real time_scaled = (time-t_HeI) / t_He * t_Hems;
    
    real r_c = helium_main_sequence_radius(time_scaled, m_core, m_core);
    return r_c;
}

real horizontal_branch::small_envelope_core_radius(){
    //small_envelope_core_radius equals the helium_core_radius
    //for horizontal branch stars
    return small_envelope_core_radius(relative_age, relative_mass, core_mass, metalicity);
}
    


real horizontal_branch::small_envelope_core_luminosity(const real time, 
                                const real mass, const real m_core, const real z){
    real t_HeI = helium_ignition_time(mass, z); 
    real t_He = core_helium_burning_timescale(mass, z);
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(m_core);
    real time_scaled = (time-t_HeI) / t_He * t_Hems;
    real l_c = helium_main_sequence_luminosity(time_scaled, m_core);
    return l_c;
}
real horizontal_branch::small_envelope_core_luminosity(){
    return small_envelope_core_luminosity(relative_age, relative_mass, core_mass, metalicity);
}


// Eq.61
real horizontal_branch::core_helium_burning_luminosity(const real time, 
                                                 const real mass, 
                                                 const real mass_tot,
                                                 const real z) {
    
    real t_HeI = helium_ignition_time(mass, z); //#eq.43
    real t_He = core_helium_burning_timescale(mass, z);
    real tau = (time - t_HeI)/t_He;
    real tau_x = relative_age_at_start_of_blue_phase(mass, z);
    real l_x = helper_x_luminosity(mass, z);
    real l_cHe;
    
    //tau can be slightly larger than one 
    //as t_He is precise to 1e-17 and t_HeI to 1e-16 
    if (tau > 1.0 && tau < 1.0 + cnsts.safety(tiny)) {
        tau=1.0;
    }
    
    if(tau>=tau_x && tau<=1.0) {
        real r_mHe = minimum_blue_loop_radius(mass, mass_tot, z);
        real r_x = helper_x_radius(mass, mass_tot,z);
        real m_FGB = helium_ignition_mass(z);
        
        
        if (mass >= max(12.0, m_FGB) && abs(r_x-r_mHe) >cnsts.safety(tiny)){
            cerr<<"warning in single_star::l_cheb"<<endl;
            cerr<<"erase this function"<<endl;
        }
        
        real xi = min(2.5, max(0.4, r_mHe/r_x));
        real lambda = pow((tau - tau_x)/(1 - tau_x), xi);
        real l_bagb = base_AGB_luminosity(mass, z);
        l_cHe = l_x*pow(l_bagb/l_x, lambda);     
    }
    else if (tau>=0.0 && tau<tau_x) {
        real lambda = pow((tau_x - tau)/tau_x, 3);
        real l_HeI = helium_ignition_luminosity(mass, z);
        l_cHe = l_x*pow(l_HeI/l_x, lambda); 
    }
    else {
        cerr << "WARNING: tau not within allowed interval [0, 1]"<<endl;
        cerr << "in single_star::core_helium_burning_luminosity(...)"<<endl;
        cerr.precision(HIGH_PRECISION);
        PRI(4);PRC(time);PRC(tau_x);PRL(tau);
        PRC(t_HeI);PRL(t_He);
        cerr.precision(STD_PRECISION);
        cerr << flush;
        exit(-1);
    }
    
    return l_cHe;
}

// Eq.64
real horizontal_branch::core_helium_burning_radius(const real time, 
                                             const real mass, 
                                             const real mass_tot,                    
                                             const real z,
                                             const real lum) {
    real t_HeI = helium_ignition_time(mass, z); //#eq.43    
    real t_He = core_helium_burning_timescale(mass, z);
    real tau = (time - t_HeI)/t_He;
    real tau_x = relative_age_at_start_of_blue_phase(mass, z);
    real tau_y = relative_age_at_end_of_blue_phase(mass, z);
    
    //tau can be slightly larger than one 
    //as t_He is precise to 1e-17 and t_HeI to 1e-16 
    if (tau >1.0 && tau < 1.0+cnsts.safety(tiny)){
        tau = 1.0;
    }
    
    real r_cHe;
    if (tau>=0 && tau<tau_x) {
        r_cHe = giant_branch_radius(lum, mass_tot, z);
    } 
    else if(tau>=tau_y && tau<=1.0) {
        r_cHe = AGB_radius(lum, mass, mass_tot, z);
    }
    else if (tau>=tau_x && tau<tau_y) {
        real r_mHe = minimum_blue_loop_radius(mass, mass_tot, z);
        real r_x = helper_x_radius(mass, mass_tot, z);
        real r_y = helper_y_radius(mass, mass_tot, z);
        real dt = tau_y-tau_x;
        real r_min, rho;
                
        if (r_mHe < r_x){
            r_min = r_mHe;
            rho = pow(log(r_y/r_min), ONE_THIRD) * (tau - tau_x)/dt
            - pow(log(r_x/r_min), ONE_THIRD) * (tau_y - tau)/dt;
        }
        else {
            r_min = r_x;
            rho = pow(log(r_y/r_min), ONE_THIRD) * (tau - tau_x)/dt;   
        }
        r_cHe = r_min*exp(pow(abs(rho), 3));
    }
    else {
        cerr << "WARNING: tau not within allowed interval [0, 1]"<<endl;
        cerr << "in single_star::core_helium_burning_radius(...)"<<endl;
        PRI(4);PRC(time);PRC(tau_x);PRC(tau_y);PRL(tau);
        cerr << flush;
        exit(-1);
    }
    
    return r_cHe;
}


real horizontal_branch::helper_x_luminosity(const real mass, 
                                      const real z) {
    
    real l_x;
    if (mass < helium_flash_mass(z)) {
        l_x = base_horizontal_branch_luminosity(mass, z);
    }
    else if (mass < helium_ignition_mass(z)) {
        l_x = minimum_horizontal_branch_luminosity(mass, z);
    }
    else {
        l_x = helium_ignition_luminosity(mass, z);
    }
    return l_x;
}

// Eq.60
real horizontal_branch::helper_x_radius(const real mass, 
                                  const real mass_tot, const real z) {
    
    real r_x;
    if (mass < helium_flash_mass(z)) {
        r_x = base_horizontal_branch_radius(mass, mass_tot, z);
    }
    else if (mass < helium_ignition_mass(z)) {
        real l_mHe = minimum_horizontal_branch_luminosity(mass, z);
        r_x = giant_branch_radius(l_mHe, mass_tot, z);
    }
    else {
        r_x = helium_ignition_radius(mass, mass_tot, z);
    }
    return r_x;
}

real horizontal_branch::helper_y_luminosity(const real mass, 
                                      const real mass_tot, const real z) {
    
    real l_y;
    real m_FGB = helium_ignition_mass(z);
    real l_bagb = base_AGB_luminosity(mass, z);
    
    if(mass<m_FGB){
        l_y = l_bagb;
    }
    else {       
        real tau_y = relative_age_at_end_of_blue_phase(mass, z);
        real tau_x = relative_age_at_start_of_blue_phase(mass, z);
        
        real r_mHe = minimum_blue_loop_radius(mass, mass_tot, z);
        real r_x = helper_x_radius(mass, mass_tot, z);
        real xi = min(2.5, max(0.4, r_mHe/r_x));
        real lambda = pow((tau_y - tau_x)/(1 - tau_x), xi);
        real l_x = helper_x_luminosity(mass, z);
        l_y = l_x*pow(l_bagb/l_x, lambda); 
    }
    return l_y;
}

real horizontal_branch::helper_y_radius(const real mass, 
                                  const real mass_tot, const real z) {
    
    real l_y = helper_y_luminosity(mass, mass_tot, z);
    real r_y = AGB_radius(l_y, mass, mass_tot, z);
    
    return r_y;
}


// Eq.67
real horizontal_branch::core_helium_burning_core_mass(const real time,
                                                const real mass, 
                                                const real z) {
    
    real mc_bagb = base_AGB_core_mass(mass, z);
    real mc_HeI = helium_ignition_core_mass(mass, z);
    
    real t_HeI = helium_ignition_time(mass, z); 
    real t_He  = core_helium_burning_timescale(mass, z);
    real tau = (time - t_HeI)/t_He;    
    real m_core = (1-tau)*mc_HeI + tau*mc_bagb;

    return m_core;
    
}


//Eq. Tau_x
real horizontal_branch::relative_age_at_start_of_blue_phase(const real mass, 
                                                      const real z) {
    
    real tau_x = 0;
    if( mass >= helium_flash_mass(z) && mass <= helium_ignition_mass(z))
        tau_x = 1 - blue_phase_timescale(mass, z);
    
    return tau_x;
}    

//Eq. Tau_y
real horizontal_branch::relative_age_at_end_of_blue_phase(const real mass, 
                                                    const real z) {
    
    real tau_y = 1;
    if (mass > helium_ignition_mass(z))
        tau_y = blue_phase_timescale(mass, z);
    
    return tau_y;
}


// Eq.58
real horizontal_branch::blue_phase_timescale(const real mass, const real z) {
    
    real b45 = smc.b(45,z);
    
    real t_bl = 1;
    real m_HeF = helium_flash_mass(z);
    real m_FGB = helium_ignition_mass(z);
    real m_HeF_FGB = m_HeF/m_FGB;
    real b46 = -1*smc.b(46,z) * log10(m_HeF_FGB);
    
    if(mass>=m_HeF && mass<=m_FGB) {
        real c_bl = (1 - b45*pow(m_HeF_FGB, 0.414));
        t_bl = b45*pow(mass/m_FGB, 0.414) + c_bl*pow(log10(mass/m_FGB)/log10(m_HeF_FGB), b46);
    }
    else if(mass>m_FGB) {
        t_bl = (1-smc.b(47,z))*f_bl(mass, z)/f_bl(m_FGB, z);
    }
    
//    if (t_bl > 1.0){
//        cerr<<"WARNING in single_star::blue_phase_timescale: t_bl > 1"<<endl;
//    }
//    else if (t_bl < 0.0){
//        cerr<<"WARNING in single_star::blue_phase_timescale: t_bl < 0"<<endl;
//    }
    
    t_bl = min(1.0, max(0.0, t_bl));  
    return t_bl;
}


real horizontal_branch::f_bl(const real mass, const real z) {
    real b48 = smc.b(48,z);
    real b49 = smc.b(49,z);
    
    // r_mHe & r_agb are not functions of mass_tot (get_total_mass()) 
    // in accordance with HPT, 
    // otherwise things go wrong for stars that loose much mass
    real r_mHe = minimum_blue_loop_radius(mass, mass, z);
    real l_HeI = helium_ignition_luminosity(mass, z);
    real r_agb = AGB_radius(l_HeI, mass, mass, z);
    real ratio = 1- r_mHe / r_agb;
    if (ratio < 0) {
        ratio = 0;
        cerr<<"WARNING in single_star::f_bl: ratio 1-r_mHe / r_agb< 0"<<endl;
    }
    real f = pow(mass, b48) * pow(ratio, b49);    
    return f;
}
