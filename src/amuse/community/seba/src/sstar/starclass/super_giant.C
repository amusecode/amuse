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

    //First initialize EAGB helium and co core mass 
    McL_core_mass = 0;
    evolve_core_mass(relative_age, relative_mass, metalicity);

    // then calculate luminosity and radius
    instantaneous_element();
    small_envelope_perturbation();   
    update();
    post_constructor();
}



//		general mass transfer utilities.
// Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real super_giant::add_mass_to_accretor(const real mdot, bool hydrogen) {
    
    if (mdot<=0) {
        cerr << "super_giant::add_mass_to_accretor(mdot=" << mdot << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        
        set_spec_type(Accreting, false);
        
        return 0;
    }
    else {
        
        if(hydrogen){
            // hydrogen accretion
            // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
            //adjust_accretor_age(mdot);
            envelope_mass += mdot;
            accreted_mass += mdot;
            
            // only neccessary for AGB & He giant accretor as  
            // next_update_age is a function of total mass
            // as the maximal core mass can be the total mass
            // when total mass < chandrasekhar mass      
            adjust_next_update_age();  
            
            //if (relative_mass<get_total_mass()) 
            //  update_relative_mass(get_total_mass());
        }
        else{
            //for the moment assume helium accretion will not change a star from being a EAGB or TPAGB
            real t_du = dredge_up_time(relative_mass, metalicity);    
            real t_bagb = base_AGB_time(relative_mass, metalicity);

            if (relative_age < t_du){
                //EAGB
                core_mass += mdot;
                accreted_mass += mdot;
                
                //part to adjust relative_mass
                real b36 = smc.b(36, metalicity);
                real b37 = smc.b(37, metalicity);
                real b38 = smc.b(38, metalicity);     
                real new_relative_mass =  pow((pow(core_mass, 4.) - b38) / b36, 1./b37);
                update_relative_mass(new_relative_mass);
                
                // only neccessary for AGB & He giant accretor as  
                // next_update_age is a function of total mass
                // as the maximal core mass can be the total mass
                // when total mass < chandrasekhar mass      
                adjust_next_update_age();  
                
                //part to adjust age using the co_core_mass
                real A_He = AGB_A_He_estimator();
                real l_bagb = base_AGB_luminosity(relative_mass, metalicity);    
                relative_age = determine_age(COcore_mass, relative_mass, metalicity, A_He, t_bagb, l_bagb);
                last_update_age = t_bagb;
                
                if(relative_age < last_update_age){
                    relative_age = last_update_age;
                    real mco = determine_core_mass(relative_age, relative_mass, metalicity, 
                                              A_He, t_bagb, l_bagb);
                    
                    if(mco >= COcore_mass && mco <= core_mass) {
                        McL_core_mass = mco;
                        if(!update_COcore_mass(mco)) {
                            cerr << "Update COcore mass failed in super_giant()"<<endl;
                        }
                        
                    }
                    exit(-1);
                }
                if(relative_age > t_du){
                    //this should not be possible
                    cerr<<"EAGB helium accretion add_mass_to_accretor mc_co > mc_du ?"<<endl;
                    exit(-1);
                }
            }
            else{
                //TPAGB
                core_mass += mdot;
                COcore_mass += mdot;
                accreted_mass += mdot;
                update_relative_mass(relative_mass + mdot);
                if(core_mass != COcore_mass){
                    cerr<<"on TPAGB add_mass_to_accretor core_mass not equal to co_core_mass"<<endl;
                    exit(-1);
                }
                
                // only neccessary for AGB & He giant accretor as  
                // next_update_age is a function of total mass
                // as the maximal core mass can be the total mass
                // when total mass < chandrasekhar mass      
                adjust_next_update_age();  
                
                //part to adjust age using the co_core_mass
                real mc_du = dredge_up_core_mass(relative_mass, metalicity);
                real lambda =  min(0.9, 0.3+0.001*pow(relative_mass, 5)); // Eq.73
                McL_core_mass = (core_mass - mc_du) / (1 - lambda) + mc_du;
                
                real L_du = dredge_up_luminosity(relative_mass, metalicity);
                real L_x = FGB_x_luminosity(relative_mass, metalicity);
                real AH_He = TPAGB_AH_He_estimator();
 
                if (L_du <= L_x){
                    relative_age = determine_age(McL_core_mass, relative_mass, metalicity, AH_He, t_du, L_du);
                    last_update_age = t_bagb; 
                }
                else{
                    
                    real q = sub_giant_q_parameter(relative_mass, metalicity);
                    real B = sub_giant_B_factor(relative_mass);
                    real t_inf2 = specific_time_limit(AH_He, t_du,
                                                      B, L_du, q);
                    relative_age = t_inf2 - pow(McL_core_mass, 1.-q)/AH_He/B/(q-1.);
                    last_update_age = t_bagb;
                    
                }                    
                if(relative_age < t_du){
                    //this should not be possible
                    cerr<<"TPAGB helium accretion add_mass_to_accretor mc_co < mc_du ?"<<endl;
                    exit(-1);
                }
                if(relative_age > next_update_age){
                    //original mrel < 2.25
                    // updating mrel only usefull if new mrel > 2.25
                    
//                    real mc_max = maximum_AGB_core_mass(relative_mass, metalicity);
//                    //core should minimally grow 5% on the AGB
//                    real A_He = AGB_A_He_estimator();
//                    real t_bagb = base_AGB_time(relative_mass, metalicity);
//                    real l_bagb = base_AGB_luminosity(relative_mass, metalicity);
//                    real mc_bagb = determine_core_mass(t_bagb, relative_mass, metalicity, 
//                                                       A_He, t_bagb, l_bagb); //co core mass
//                    mc_max = max(mc_max, 1.05*mc_bagb);
//                    mc_max = min(mc_max, get_total_mass());
//
//                    
//                    real dm = core_mass - mc_max;
//                    core_mass -= dm;
//                    COcore_mass -= dm;
//                    accreted_mass -= dm;
//                    update_relative_mass(relative_mass - dm);
//                    if(core_mass != COcore_mass){
//                        cerr<<"on TPAGB add_mass_to_accretor core_mass not equal to co_core_mass"<<endl;
//                        exit(-1);
//                    }
//                    
//                    // only neccessary for AGB & He giant accretor as  
//                    // next_update_age is a function of total mass
//                    // as the maximal core mass can be the total mass
//                    // when total mass < chandrasekhar mass      
//                    adjust_next_update_age();  

                    create_remnant(relative_mass, get_total_mass(), core_mass, metalicity);
                    exit(-1);
                }
            }                
                
            cerr<<"agb add-mass_to-accretor helium"<<endl;
            exit(-1);
            
        }
    }
    
    set_spec_type(Accreting);    
    return mdot;
}

real super_giant::add_mass_to_accretor(real mdot, const real dt, bool hydrogen) {
    if (mdot<0) {
        cerr << "super_giant::add_mass_to_accretor(mdot=" << mdot 
        << ", dt=" << dt << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        return 0;
    }
    
    
    if(hydrogen){
        //hydrogen accretion
        mdot = accretion_limit(mdot, dt);
        
        // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
        // adjust_accretor_age(mdot);
        envelope_mass += mdot;
        accreted_mass += mdot;
        
        
        // only neccessary for AGB & He giant accretor as  
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass      
        adjust_next_update_age();  
        
        //  if (relative_mass<get_total_mass()) 
        //    update_relative_mass(get_total_mass());
        
        adjust_accretor_radius(mdot, dt);
        
    }
    else{
        //for the moment assume helium accretion, will not change a star from being a EAGB or TPAGB
        // for the moment no helium_accretion_limit and adjust_accretor_radius

        real t_du = dredge_up_time(relative_mass, metalicity);    
        real t_bagb = base_AGB_time(relative_mass, metalicity);
        
        if (relative_age < t_du){
            //EAGB
            core_mass += mdot;
            accreted_mass += mdot;
            
            //part to adjust relative_mass
            real b36 = smc.b(36, metalicity);
            real b37 = smc.b(37, metalicity);
            real b38 = smc.b(38, metalicity);     
            real new_relative_mass =  pow((pow(core_mass, 4.) - b38) / b36, 1./b37);
            update_relative_mass(new_relative_mass);
            
            // only neccessary for AGB & He giant accretor as  
            // next_update_age is a function of total mass
            // as the maximal core mass can be the total mass
            // when total mass < chandrasekhar mass      
            adjust_next_update_age();  
            
            //part to adjust age using the co_core_mass
            real A_He = AGB_A_He_estimator();
            real l_bagb = base_AGB_luminosity(relative_mass, metalicity);    
            relative_age = determine_age(COcore_mass, relative_mass, metalicity, A_He, t_bagb, l_bagb);
            last_update_age = t_bagb;
            
            if(relative_age < last_update_age){
                relative_age = last_update_age;
                real mco = determine_core_mass(relative_age, relative_mass, metalicity, 
                                               A_He, t_bagb, l_bagb);
                
                if(mco >= COcore_mass && mco <= core_mass) {
                    McL_core_mass = mco;
                    if(!update_COcore_mass(mco)) {
                        cerr << "Update COcore mass failed in super_giant()"<<endl;
                    }
                    
                }
                exit(-1);
            }
            if(relative_age > t_du){
                //this should not be possible
                cerr<<"EAGB helium accretion add_mass_to_accretor mc_co > mc_du ?"<<endl;
                exit(-1);
            }
        }
        else{
            //TPAGB
            core_mass += mdot;
            COcore_mass += mdot;
            accreted_mass += mdot;
            update_relative_mass(relative_mass + mdot);
            if(core_mass != COcore_mass){
                cerr<<"on TPAGB add_mass_to_accretor core_mass not equal to co_core_mass"<<endl;
                exit(-1);
            }
            
            // only neccessary for AGB & He giant accretor as  
            // next_update_age is a function of total mass
            // as the maximal core mass can be the total mass
            // when total mass < chandrasekhar mass      
            adjust_next_update_age();  
            
            //part to adjust age using the co_core_mass
            real mc_du = dredge_up_core_mass(relative_mass, metalicity);
            real lambda =  min(0.9, 0.3+0.001*pow(relative_mass, 5)); // Eq.73
            McL_core_mass = (core_mass - mc_du) / (1 - lambda) + mc_du;
            
            real L_du = dredge_up_luminosity(relative_mass, metalicity);
            real L_x = FGB_x_luminosity(relative_mass, metalicity);
            real AH_He = TPAGB_AH_He_estimator();
            
            if (L_du <= L_x){
                relative_age = determine_age(McL_core_mass, relative_mass, metalicity, AH_He, t_du, L_du);
                last_update_age = t_bagb; 
            }
            else{
                
                real q = sub_giant_q_parameter(relative_mass, metalicity);
                real B = sub_giant_B_factor(relative_mass);
                real t_inf2 = specific_time_limit(AH_He, t_du,
                                                  B, L_du, q);
                relative_age = t_inf2 - pow(McL_core_mass, 1.-q)/AH_He/B/(q-1.);
                last_update_age = t_bagb;
                
            }                    
            if(relative_age < t_du){
                //this should not be possible
                cerr<<"TPAGB helium accretion add_mass_to_accretor mc_co < mc_du ?"<<endl;
                exit(-1);
            }
            if(relative_age > next_update_age){
                //original mrel < 2.25
                // updating mrel only usefull if new mrel > 2.25
                
//                real mc_max = maximum_AGB_core_mass(relative_mass, metalicity);
//                //core should minimally grow 5% on the AGB
//                real A_He = AGB_A_He_estimator();
//                real t_bagb = base_AGB_time(relative_mass, metalicity);
//                real l_bagb = base_AGB_luminosity(relative_mass, metalicity);
//                real mc_bagb = determine_core_mass(t_bagb, relative_mass, metalicity, 
//                                                   A_He, t_bagb, l_bagb); //co core mass
//                mc_max = max(mc_max, 1.05*mc_bagb);
//                mc_max = min(mc_max, get_total_mass());
//                
//                real dm = core_mass - mc_max;
//                core_mass -= dm;
//                COcore_mass -= dm;
//                accreted_mass -= dm;
//                update_relative_mass(relative_mass - dm);
//                if(core_mass != COcore_mass){
//                    cerr<<"on TPAGB add_mass_to_accretor core_mass not equal to co_core_mass"<<endl;
//                    exit(-1);
//                }
//                
//                // only neccessary for AGB & He giant accretor as  
//                // next_update_age is a function of total mass
//                // as the maximal core mass can be the total mass
//                // when total mass < chandrasekhar mass      
//                adjust_next_update_age();  

                create_remnant(relative_mass, get_total_mass(), core_mass, metalicity);
                exit(-1);
            }
        }                
        
        cerr<<"agb add-mass_to-accretor helium"<<endl;
        exit(-1);
        
    }
    set_spec_type(Accreting);
    return mdot;
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
               return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
           }
          else {
               star_transformation_story(Carbon_Dwarf);	   
               return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
          }     
      }
  }
  else {
      envelope_mass -= mdot;
      // next_update_age is a function of total mass
      // as the maximal core mass can be the total mass
      // when total mass < chandrasekhar mass
      adjust_next_update_age();        
  }

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
                return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
            }
            else {
                star_transformation_story(Carbon_Dwarf);	   
                return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
            }     
        }
    }
    else{ 
        envelope_mass -= mdot;
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass
        adjust_next_update_age();        
    
        // (GN+SPZ Apr 29 1999)
        adjust_donor_radius(mdot);
    }
    
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
      real dt_tagb_old = TAGB_time(relative_mass, get_total_mass(), metalicity)-t_bagb_old;

      real z_new = get_metalicity();
      real t_bagb_new = base_AGB_time(m_rel_new, z_new); 
      real t_tagb_new = TAGB_time(m_rel_new, m_tot_new, z_new);
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
    
  real t_tagb = TAGB_time(relative_mass, get_total_mass(), metalicity);

  if(relative_age>t_tagb) {
    cerr << "WARNING: relative_age > t_tagb in super_giant"<<endl;
    exit(-1);
    //relative_age = t_tagb;
  }
    
  next_update_age = t_tagb;
}


real super_giant::zeta_thermal() {

  real z = 0.; // (GN+SPZ Apr 29 1999) was -10 in 1998; 
               // was -0.64 somewhere in the past (~1992 or so).

      return z;
   }

real super_giant::gyration_radius_sq() {

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
    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x_dj = min(1.0, (luminosity -4000.0)/500.0);
        dm_dj = x_dj * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.5);
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
    
    
    // Vassiliadis & Wood 1993
    // AGB (including and leading up to superwind phase)
    //Vassiliadis & Wood only valid in small range, HPT uses it everywhere on the AGB
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
        
    wind_constant = max(max(max(max(dm_wr, dm_dj), dm_r), dm_vw), 0.0) +dm_lbv;
    
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
            //first calculate EAGB helium and co core mass             
            evolve_core_mass(relative_age, relative_mass, metalicity);
            // then calculate luminosity and radius
            instantaneous_element();
            small_envelope_perturbation();   
            // if no envelope make transition to remnants
            // just as a procedure: reduce_mass with 1
            if (envelope_mass <= 0){
                reduce_mass(1.);
                return;    
            }
        }
        else {
            create_remnant(relative_mass, get_total_mass(), core_mass, metalicity);
            return;
        }

      update();
      stellar_wind(dt);
}



real super_giant::get_evolve_timestep() {
    
    real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                        next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   
    
    //extra safety measure
    // when L and R increase rapidly, so will mdot
    real l_du = dredge_up_luminosity(relative_mass, metalicity);
    real A;
    if (luminosity < l_du){
        A = AGB_A_He_estimator();
    }
    else{
        A = TPAGB_AH_He_estimator();
    }
    
    real l_x = FGB_x_luminosity(relative_mass, metalicity);
    // radius should be a function of get_total_mass, but M_rel is the best approximation for
    // M_tot_bgb
    real l_bagb = base_AGB_luminosity(relative_mass, metalicity);
    real r_bagb = AGB_radius(l_bagb, relative_mass, get_total_mass(), metalicity);

    real dt_mdot = timestep;
    if (luminosity < l_x){
        real p = sub_giant_p_parameter(relative_mass, metalicity);
        dt_mdot = McL_core_mass / ( p * luminosity * A) * 0.1 * r_bagb/ radius;
    }
    else{
        real q = sub_giant_q_parameter(relative_mass, metalicity);
        dt_mdot = McL_core_mass / ( q * luminosity * A) * 0.1 * r_bagb /radius;
    }
    
    return max(min(timestep, dt_mdot), cnsts.safety(minimum_timestep));
    
}


void super_giant::create_remnant(const real mass, const real mass_tot, const real mc_core, const real z) {

    if (is_binary_component()) 
        get_binary()->dump("binev.data", false);
    else
        dump("binev.data", false);

    stellar_type type;
    real mc_bagb = base_AGB_core_mass(mass, z);
    real mc_SN = maximum_AGB_core_mass(mass, z);

    // if mc_core equals get_total_mass()
    // core mass reaches outside of star, no envelope anymore
    if (mc_core < cnsts.parameters(Chandrasekar_mass)){
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
        new white_dwarf(*this, Carbon_Dwarf);
        return;
        case Oxygen_Dwarf : star_transformation_story(Oxygen_Dwarf);
        new white_dwarf(*this, Oxygen_Dwarf);
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
    real mc_bagb = base_AGB_core_mass(mass, z);
    real mco;
    if (time <= t_du) {
        if(!update_core_and_envelope_mass(mc_bagb)) {
            cerr << "Update core mass failed in super_giant()"<<endl;
        }
        
//        if(abs(core_mass-base_AGB_core_mass(mass,z)) >cnsts.safety(tiny)){
//            cerr<<"Core_mass needs update on EAGB, this should not be the case"<<endl;
//        }
        
        real A_He = AGB_A_He_estimator();
        real t_bagb = base_AGB_time(mass, z);
        real l_bagb = base_AGB_luminosity(mass, z);
        mco = determine_core_mass(time, mass, z, 
				      A_He, t_bagb, l_bagb);
  
        if(mco >= COcore_mass && mco <= core_mass) {
            McL_core_mass = mco;
            if(!update_COcore_mass(mco)) {
                cerr << "Update COcore mass failed in super_giant()"<<endl;
            }
            
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
        
        
        if(!update_COcore_mass(mco)) {
            cerr << "Update COcore mass failed in helium_giant()"<<endl;
        }
        
        McL_core_mass = m_core;
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
}


real super_giant::helium_core_radius(const real time, const real mass, const real mass_tot, const real m_core, const real z){
    real r_c, l_c;
    real t_du = dredge_up_time(mass, z);    
    if(time < t_du){
        //EAGB
        real t_tagb = TAGB_time(mass, mass_tot, z);
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
        real t_tagb = TAGB_time(mass, mass_tot, z);
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
        real t_tagb = TAGB_time(mass, mass_tot, z);
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
    

