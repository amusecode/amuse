//
// horizontal_branch.C
//

#include "horizontal_branch.h"
#include "sub_giant.h"
#include "hertzsprung_gap.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.
//

horizontal_branch::horizontal_branch(sub_giant & g) : single_star(g) {

    delete &g;
    real m_HeF = helium_flash_mass(metalicity);
    if ( relative_mass< m_HeF){
        relative_mass = get_total_mass();
        relative_age = helium_ignition_time(relative_mass, metalicity);
    }    
    real t_HeI = helium_ignition_time(relative_mass, metalicity); 
    PRC(relative_mass);PRC(get_total_mass());PRC(relative_age);PRL(t_HeI);

    
// (GN+SPZ May  4 1999) last update age is time of previous type change
    last_update_age = next_update_age;

    adjust_next_update_age();
    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();   
    
    update();
    post_constructor();
    t_HeI = helium_ignition_time(relative_mass, metalicity); 
    PRC(relative_mass);PRC(get_total_mass());PRC(relative_age);PRL(t_HeI);

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

star* horizontal_branch::reduce_mass(const real mdot) {

      if (envelope_mass<=mdot) {
         envelope_mass = 0;
         star_transformation_story(Helium_Star);
         return dynamic_cast(star*, new helium_star(*this));
      }

      else 
          envelope_mass -= mdot;
      
    return this;
   }

star* horizontal_branch::subtrac_mass_from_donor(const real dt,
						 real& mdot) {

      real mdot_temp = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_temp);

      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;
         star_transformation_story(Helium_Star);
         return dynamic_cast(star*, new helium_star(*this));
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
    cerr<<"hb::zeta_adiabatic is used?"<<endl;

    return 15;	
   }

real horizontal_branch::zeta_thermal() {
    cerr<<"hb::zeta_thermal is used?"<<endl;
      return 15;
   }

void horizontal_branch::adjust_next_update_age() {
  real t_HeI = helium_ignition_time();
  real dt_cHe = core_helium_burning_timescale();

  if(relative_age!=t_HeI) {
    cerr << "WARNING: relative_age != t_HeI in horizontal_branch"<<endl;
    cerr.precision(HIGH_PRECISION);
    PRC(relative_age);PRL(t_HeI);
    cerr.precision(STD_PRECISION);
    relative_age = t_HeI;
  }
  next_update_age = t_HeI + dt_cHe;
}

real horizontal_branch::gyration_radius_sq() {
    cerr<<"hb::gyration_radius_sq is used?"<<endl;

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
        real x_dj = min(1.0, (luminosity -4000.0)/500.0);
        dm_dj = x_dj * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
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
  
    if(dm_h > dm_dj && dm_h > dm_v && dm_h > dm_r) cerr<<"CHeB: WR_like"<<endl;
    else if (dm_dj > dm_v && dm_dj > dm_r) cerr<< "CHeB: de Jager"<<endl;
    else if (dm_v > dm_r) cerr<<"CHeB: Vink"<<endl;   
    else cerr<<"CHeB: Reimers"<<endl;
    
    
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
          
          if (envelope_mass == 0){
              star_transformation_story(Helium_Star);
              new helium_star(*this);
              //return dynamic_cast(star*, new helium_star(*this));
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
      //stellar_wind(dt);
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
    if (intermediate_mass_star(mass, z)){
        l_x = minimum_horizontal_branch_luminosity(mass, z);
    }
    else if(high_mass_star(mass, z)){
        l_x = helium_ignition_luminosity(mass, z);
    }
    else{
        l_x = base_horizontal_branch_luminosity(mass, z);
    }
    
    return l_x;
}

// Eq.60
real horizontal_branch::helper_x_radius(const real mass, 
                                  const real mass_tot, const real z) {
    
    real r_x;
    if (low_mass_star(mass, z)) {
        r_x = base_horizontal_branch_radius(mass, mass_tot, z);
    }
    else if (intermediate_mass_star(mass, z)) {
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
    if(intermediate_mass_star(mass, z))
        tau_x = 1 - blue_phase_timescale(mass, z);
    
    return tau_x;
}    

//Eq. Tau_y
real horizontal_branch::relative_age_at_end_of_blue_phase(const real mass, 
                                                    const real z) {
    
    real tau_y = 1;
    if(high_mass_star(mass, z))
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
    
    if (t_bl > 1.0){
        cerr<<"WARNING in single_star::blue_phase_timescale: t_bl > 1"<<endl;
        PRC(mass);PRC(m_FGB);PRC(m_HeF);PRL(t_bl);
        exit(-1);
    }
    else if (t_bl < 0.0){
        cerr<<"WARNING in single_star::blue_phase_timescale: t_bl < 0"<<endl;
        PRC(mass);PRC(m_FGB);PRL(t_bl);
        exit(-1);
    }
    
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
