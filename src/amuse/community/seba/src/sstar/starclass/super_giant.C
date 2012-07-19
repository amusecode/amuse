//
// super_giant.C
//

#include "super_giant.h"
#include "horizontal_branch.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.

super_giant::super_giant(horizontal_branch & h) : single_star(h) {

  delete &h;

  // (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;

  adjust_next_update_age();

  //Initialize EAGB helium and co core mass 
  McL_core_mass = 0;
  core_mass = base_AGB_core_mass(relative_mass, metalicity);  
  evolve_core_mass(relative_age, relative_mass, metalicity);

  effective_radius = 0;
  instantaneous_element();
  small_envelope_perturbation();   
 
  update();
  post_constructor();
}

star* super_giant::reduce_mass(const real mdot) {
    //if (envelope_mass<=mdot) {
//        envelope_mass = 0;
//        
//        // (GN+SPZ Apr 29 1999) stripped super_giants become 
//        // white_dwarf or helium stars
//        if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
//           core_mass     >= cnsts.parameters(helium2neutron_star)) {
//            
//            // (SPZ+GN: 27 Jul 2000)
//            // Initialize core_mass as CO core mass and envelope mass
//            // as helium core mass to make super giant ready to become
//            // a helium giant.
//            real m_tot = core_mass;
//            core_mass = COcore_mass;
//            envelope_mass = m_tot - core_mass;
//            
//            star_transformation_story(Helium_Giant);
//            return dynamic_cast(star*, new helium_giant(*this));
//        }
//        else {
//            
//            if(relative_age <= dredge_up_time(relative_mass, metalicity)) {
//                
//                real m_tot = core_mass;
//                core_mass = COcore_mass;
//                envelope_mass = m_tot - core_mass;
//                
//                star_transformation_story(Helium_Giant);
//                return dynamic_cast(star*, new helium_giant(*this));
//            }
//            else {
//                star_transformation_story(Carbon_Dwarf);	   
//                return dynamic_cast(star*, new white_dwarf(*this));
//            }
//        }
//    }
//    
  if (envelope_mass<=mdot) {
      envelope_mass = 0;

      real t_du = dredge_up_time(relative_mass, metalicity);    
      if (relative_age < t_du){    
           real m_tot = core_mass;
           core_mass = COcore_mass;
           envelope_mass = m_tot - core_mass;

           star_transformation_story(Helium_Giant);
           return dynamic_cast(star*, new helium_giant(*this));
      }
	  else {
          real mc_bagb = base_AGB_core_mass(relative_mass, metalicity);
          if(mc_bagb >= 1.6) {
             
               star_transformation_story(Oxygen_Dwarf);
               return dynamic_cast(star*, new white_dwarf(*this));
           }
          else {
               star_transformation_story(Carbon_Dwarf);	   
               return dynamic_cast(star*, new white_dwarf(*this));
          }     
      }
  }

  envelope_mass -= mdot;
  return this;
}

star* super_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      real mdot_temp = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_temp);

//      if (envelope_mass<=mdot) {
//         mdot = envelope_mass;
//         envelope_mass = 0;
//	 
//          // (GN+SPZ Apr 29 1999) stripped super_giants become 
//          // white_dwarf or helium stars
//         if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
//            core_mass     >= cnsts.parameters(helium2neutron_star)) {
//
//           real m_tot = core_mass;
//           core_mass = COcore_mass;
//           envelope_mass = m_tot - core_mass;
//
//           star_transformation_story(Helium_Giant);
//           return dynamic_cast(star*, new helium_giant(*this));
//         }
//         else {
//	   
//            // (SPZ+GN: 27 Jul 2000)
//            if(relative_age <= dredge_up_time(relative_mass, metalicity)) {
//
//            real m_tot = core_mass;
//            core_mass = COcore_mass;
//            envelope_mass = m_tot - core_mass;
//
//             star_transformation_story(Helium_Giant);
//             return dynamic_cast(star*, new helium_giant(*this));
//            }
//            else {
//             star_transformation_story(Carbon_Dwarf);	   
//             return dynamic_cast(star*, new white_dwarf(*this));
//            }
//         }
//      }

    if (envelope_mass<=mdot) {
        mdot = envelope_mass;
        envelope_mass = 0;
        
        real t_du = dredge_up_time(relative_mass, metalicity);    
        if (relative_age < t_du){    
            real m_tot = core_mass;
            core_mass = COcore_mass;
            envelope_mass = m_tot - core_mass;
            
            star_transformation_story(Helium_Giant);
            return dynamic_cast(star*, new helium_giant(*this));
        }
        else {
            real mc_bagb = base_AGB_core_mass(relative_mass, metalicity);
            if(mc_bagb >= 1.6) {
                
                star_transformation_story(Oxygen_Dwarf);
                return dynamic_cast(star*, new white_dwarf(*this));
            }
            else {
                star_transformation_story(Carbon_Dwarf);	   
                return dynamic_cast(star*, new white_dwarf(*this));
            }     
        }
    }
    
      // (GN+SPZ Apr 29 1999)
      adjust_donor_radius(mdot);

      envelope_mass -= mdot;
      return this;
}

// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void super_giant::adjust_accretor_age(const real mdot, 
				      const bool rejuvenate) {
    cerr<<"sub_giant::adjust_accretor_age is currently not used"<<endl;

    
    cerr<<"Maybe super_giant::adjust_accretor_age should be changed into two parts."<<endl;
    cerr<<"The relative_age should be changed with limits:"<<endl;
    cerr<<" time < t_du: t_bagb & t_du "<<endl;
    cerr<<"time > t_du: t_du & t_tagb"<<endl;
      real m_tot_new = get_total_mass() + mdot;
      real m_rel_new = max(m_tot_new, relative_mass);

      real t_bagb_old = base_AGB_time(relative_mass, metalicity);  
      real dt_tagb_old = TAGB_time(relative_mass, metalicity)-t_bagb_old;

      real z_new = get_metalicity();
      real t_bagb_new = base_AGB_time(m_rel_new, z_new); 
      real t_tagb_new = TAGB_time(m_rel_new, z_new);
      real dt_tagb_new = t_tagb_new - t_bagb_new; 
  
      real dtime = relative_age - t_bagb_old;

      last_update_age = t_bagb_new;
      relative_age = t_bagb_new
                   + dtime*(dt_tagb_new/dt_tagb_old);
    if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new);

      if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
         cerr<<"In super_giant::adjust_accretor_age relative age updated on AGB, but < last_update_age"<<endl;
      }
    
       relative_age = max(relative_age, 
			  last_update_age + cnsts.safety(minimum_timestep));
       relative_age = min(relative_age, t_tagb_new);


      // next_update_age should not be reset here
      // next_update_age = t_nuc;
   }

void super_giant::adjust_next_update_age() {
    
  real t_tagb = TAGB_time(relative_mass, metalicity);

  if(relative_age>t_tagb) {
    cerr << "WARNING: relative_age > t_tagb in super_giant"<<endl;
    relative_age = t_tagb;
  }
  next_update_age = t_tagb;
}


real super_giant::zeta_thermal() {
    cerr<<"supg::zeta_thermal is used?"<<endl;

  real z = 0.; // (GN+SPZ Apr 29 1999) was -10 in 1998; 
               // was -0.64 somewhere in the past (~1992 or so).

      return z;
   }

real super_giant::gyration_radius_sq() {
    cerr<<"supg::gyration_radius_sq is used?"<<endl;

  return cnsts.parameters(convective_star_gyration_radius_sq); 
}

void super_giant::update_wind_constant() {
#if 0  
// (GN Apr  1 1999) fit for massive stars to Maeder (but not so much wind loss)
// (GN Apr 16 1999) low mass stars need separate treatment

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    wind_constant = meader_fit_dm;

  }
  else {
// (GN+SPZ May  4 1999) nor neede: see single_star::stelar_wind
//    real factor = 1- pow(relative_age/next_update_age,cnsts.parameters(
//                                      massive_star_mass_loss_law));
//    wind_constant = 0.2*(get_total_mass() - final_core_mass())/factor;

    // (SPZ+GN: 27 Jul 2000) 0.8 of inital envelope is lost on AGB
    // see Nelemans YPZV 2000
    wind_constant = 0.8*(relative_mass - final_core_mass());
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
    
    // Vassiliadis & Wood 1993
    // AGB (including and leading up to superwind phase)
    cerr<<"Vassiliadis & Wood only valid in small range, HPT uses it everywhere on the AGB"<<endl;
    real P = pow(10, -2.07 + 1.94 * log10(radius) - 0.9 * log10(get_total_mass()));//Pulsation period
    real dm_vw;
    if (get_total_mass() > 2.5) {
        dm_vw = -11.4 + 0.0125*(P-100*(get_total_mass()-2.5));
    }
    else{
        dm_vw = -11.4 + 0.0123* P;
    }
    dm_vw = pow(10, dm_vw);
    
    real v_exp = -13.5 + 0.056*P;
    v_exp = min(15.0, max(3.0, v_exp));
    dm_vw = min(dm_vw, luminosity / v_exp * 2.0589E-8);
    
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
    wind_constant = max(max(max(max(max(dm_h, dm_dj), dm_v), dm_r), dm_vw), 0.0) +dm_lbv;

    if(dm_h > dm_dj && dm_h > dm_v && dm_h > dm_r && dm_h > dm_vw) cerr<<"AGB: WR_like"<<endl;
    else if (dm_dj > dm_v && dm_dj > dm_r && dm_dj > dm_vw) cerr<< "AGB: de Jager"<<endl;
    else if (dm_v > dm_r && dm_v > dm_vw) cerr<<"AGB: Vink"<<endl;   
    else if (dm_r >dm_vw) cerr<<"AGB: Reimers"<<endl;
    else cerr<< "AGB:Vassiliadis & Wood" <<endl;
    
    
}

void super_giant::instantaneous_element() {
  luminosity = AGB_luminosity(McL_core_mass,
			      relative_mass,
			      metalicity);
  radius = AGB_radius(luminosity, relative_mass, get_total_mass(), metalicity);
}

// Evolve a super_giant upto time argument according to
// the new 2000 models.
void super_giant::evolve_element(const real end_time) {
      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;
    
        if (relative_age<=next_update_age) {
            evolve_core_mass(relative_age, relative_mass, metalicity);
            instantaneous_element();
            small_envelope_perturbation();   
        }
        else {
            create_remnant(relative_mass, metalicity);
            return;
        }

      update();
      //stellar_wind(dt);
}


void super_giant::create_remnant(const real mass, const real z) {

    if (is_binary_component()) 
        get_binary()->dump("binev.data", false);
    else
        dump("binev.data", false);

    stellar_type type;
    real mc_bagb = base_AGB_core_mass(mass, z);
    real mc_SN = maximum_AGB_core_mass(mass, z);
    real mc_max = min(get_total_mass(), mc_SN);

    // if mc_max equals get_total_mass()
    // core mass reaches outside of star, no envelope anymore
    if (mc_max < cnsts.parameters(Chandrasekar_mass)){
        if (mc_bagb < 1.6)
            type = Carbon_Dwarf;
        else if (mc_bagb <= 2.25)
            type = Oxygen_Dwarf;
        else {
            cerr<<"Warning: in super_giant::create_remnant: "
                <<"mc_bagb > 2.25 has lost envelope"<<endl;
        }
    }
    else {    
        if (mc_bagb < 1.6) 
            type = Disintegrated;
        else {
            if (mc_SN <= 7.)
                type = Neutron_Star;
            else
                type = Black_Hole;
        }
    }
    
    switch (type) {
        case Black_Hole : star_transformation_story(Black_Hole);
        new black_hole(*this); 
        return;
        case Neutron_Star : star_transformation_story(Neutron_Star);
        new neutron_star(*this);
        return;
        case Disintegrated : star_transformation_story(Disintegrated);
        new disintegrated(*this);
        return;
        case Carbon_Dwarf : star_transformation_story(Carbon_Dwarf);
        new white_dwarf(*this);
        return;
        case Oxygen_Dwarf : star_transformation_story(Oxygen_Dwarf);
        new white_dwarf(*this);
        return;
        default :   cerr << "super_giant::create_remnant(mass, z)" <<endl;
        cerr << "star_type not recognized." << endl;
        exit(-1);
    }

}

void super_giant::evolve_core_mass(const real time,
				   const real mass,
				   const real z) {
    
    real t_du = dredge_up_time(mass, z);    
    real mco;
    if (time <= t_du) {
        
        if(abs(core_mass-base_AGB_core_mass(mass,z)) >cnsts.safety(tiny)){
            cerr<<"Core_mass needs update on EAGB, this should not be the case"<<endl;
            if(!update_core_and_envelope_mass(COcore_mass)) {
                cerr << "Update core mass failed in super_giant()"<<endl;
            }
        }
        
        real A_He = AGB_A_He_estimator();
        real t_bagb = base_AGB_time(mass, z);
        real l_bagb = base_AGB_luminosity(mass, z);
        mco = determine_core_mass(time, mass, z, 
				      A_He, t_bagb, l_bagb);
  
        if(mco >= COcore_mass && mco <= core_mass) {
            McL_core_mass = COcore_mass = mco;
        }
        else {
            cerr << "WARNING: in void super_giant::evolve_core_mass(...)"<<endl;
            cerr << "New COcore_mass < current CO core mass"<<endl;
            PRC(mco);PRC(COcore_mass);PRL(core_mass);
            dump(cerr, false);
        }
    }
    else {
        // TPAGB
        // Take the growth of the core mass into account for the 
        // third dredge-up
        cerr<<"Discontinuity in Mc,He"<<endl; 
        real L_du = dredge_up_luminosity(mass, z);
        real L_x = FGB_x_luminosity(mass, z);
        real AH_He = TPAGB_AH_He_estimator();
    
        real m_core;
        if (L_du <= L_x){
            m_core = determine_core_mass(time, mass, z, 
                           AH_He, t_du, L_du);
        }
        else{
            real B = sub_giant_B_factor(mass);
            real q = sub_giant_q_parameter(mass, z);

            real t_inf2 = specific_time_limit(AH_He, t_du,
                                          B, L_du, q);
            m_core = pow((q-1)*AH_He*B*(t_inf2-time), 1./(1-q));
        
            //safety
            real D = sub_giant_D_factor(mass, z); 
            real p = sub_giant_p_parameter(mass, z);
            real t_x = specific_time_boundary(mass, AH_He, t_du, L_du, D, p, L_x);
            if(time<=t_x) {
                cerr<<"ERROR in super_giant::evolve_core_mass"<<endl;
                cerr<<"time <= t_x need tinf1"<<endl;
            }
        }
    
        real mc_du = dredge_up_core_mass(mass, z);
        real lambda =  min(0.9, 0.3+0.001*pow(mass, 5)); // Eq.73
        mco = mc_du + (1-lambda)*(m_core - mc_du);
        
        
        if(mco >= COcore_mass) {
            COcore_mass = mco;
            McL_core_mass = m_core;
            real mc_bagb = base_AGB_core_mass(mass, z);
            if (mc_bagb <= 0.8){
                if(!update_core_and_envelope_mass(COcore_mass)) {
                    cerr << "Update core mass failed in super_giant()"<<endl;
                }
            }
            else if (mc_bagb < 2.25){
                if(!update_core_and_envelope_mass_TPAGB(COcore_mass)) {
                    cerr << "Update core mass failed in super_giant()"<<endl;
                }
            }
            else {
                cerr<<"WARNING: in void super_giant::evolve_core_mass"<<endl;
                cerr<<"mc_bagb >2.25 no TPAGB phase"<<endl;
            }
        }
        else {
            cerr << "WARNING: void super_giant::evolve_core_mass(...)"<<endl;
            cerr << "New COcore_mass < current CO core mass"<<endl;
            PRC(mco);PRC(COcore_mass);PRL(core_mass);
            dump(cerr, false);
        }
    }
}


real super_giant::helium_core_radius(const real time, const real mass, const real mass_tot, const real m_core, const real z){
    real r_c, l_c;
    real t_du = dredge_up_time(mass, z);    
    if(time < t_du){
        //EAGB
        real t_tagb = TAGB_time(mass, z);
        real t_bagb = base_AGB_time(mass, z);
        real tau = 3.*(time-t_bagb) / (t_tagb-t_bagb);    
        l_c = helium_giant_luminosity_from_core_mass(COcore_mass, core_mass, z);if (tau < 1.){
            real l_x = terminal_helium_main_sequence_luminosity(m_core);
            l_c = l_x * pow(l_c/l_x, tau);
        }
        r_c = helium_giant_radius(l_c, m_core, m_core, z);
    }
    else{    
        //TPAGB
        // due to small nucleair burning layer 
        // r_c > white_dwarf_radius
        r_c = 5.*white_dwarf_radius(m_core, 0.); 
    }
    return r_c;
}
real super_giant::helium_core_radius(){
    return helium_core_radius(relative_age, relative_mass, get_total_mass(), core_mass, metalicity);
}

real super_giant::small_envelope_core_radius(const real time, const real mass, const real mass_tot, const real m_core, const real z){
    real r_c, l_c;
    real t_du = dredge_up_time(mass, z);    
    if(time < t_du){
        //EAGB
        real t_tagb = TAGB_time(mass, z);
        real t_bagb = base_AGB_time(mass, z);
        real tau = 3.*(time-t_bagb) / (t_tagb-t_bagb);    
        l_c = helium_giant_luminosity_from_core_mass(COcore_mass, core_mass, z);
        if (tau < 1.){
            real l_x = terminal_helium_main_sequence_luminosity(m_core);
            l_c = l_x * pow(l_c/l_x, tau);
        }
        r_c = helium_giant_radius(l_c, m_core, m_core, z);
    }
    else{    
        //TPAGB
        r_c = white_dwarf_radius(m_core, 0.); 
    }
    return r_c;
}
real super_giant::small_envelope_core_radius(){
    return small_envelope_core_radius(relative_age, relative_mass, get_total_mass(), core_mass, metalicity);
}

    
real super_giant::small_envelope_core_luminosity(const real time, const real mass, const real mass_tot, const real m_core, const real z){
    real l_c;
    real t_du = dredge_up_time(mass, z);    
    if(time < t_du){
        //EAGB
        real t_tagb = TAGB_time(mass, z);
        real t_bagb = base_AGB_time(mass, z);
        real tau = 3.*(time-t_bagb) / (t_tagb-t_bagb);    
        l_c = helium_giant_luminosity_from_core_mass(COcore_mass, core_mass, z);
        if (tau < 1.){
            real l_x = terminal_helium_main_sequence_luminosity(m_core);
            l_c = l_x * pow(l_c/l_x, tau);
        }   
    }
    else{    
        //TPAGB
        l_c = 40.;
    }
    return l_c;
}
real super_giant::small_envelope_core_luminosity(){
    return small_envelope_core_luminosity(relative_age, relative_mass, get_total_mass(), core_mass, metalicity);
}
    

