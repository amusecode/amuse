//
// helium_giant.C
//

#include "super_giant.h"
#include "helium_star.h"
#include "helium_giant.h"


// (GN+SPZ May  3 1999) only for stars more massive than 
// super_giant2neutron_star, below become white_dwarf
// to be done (SPZ+GN: 27 Jul 2000)
// (ST: 10 Sep 2000)  it's different now with the HPT tracks 
helium_giant::helium_giant(super_giant & g) : single_star(g) {

    delete &g;
//    real second_dredge_up_time = next_update_age 
//                          * min(1., relative_mass
//			  / cnsts.parameters(super_giant2neutron_star));
//
//    real remaining_time = second_dredge_up_time - relative_age;
//    next_update_age = helium_time();
//    relative_age = next_update_age - remaining_time;
    
    update_relative_mass(get_total_mass());
    relative_age = helium_giant_age_core_mass_relation(core_mass, relative_mass);
    last_update_age = helium_main_sequence_time_for_solar_metalicity(relative_mass); 
    
    if (relative_age < last_update_age){
        cerr<<"in helium giant constructor from super giant age_rel < last_update_age"<<endl;
        real m_rel= get_relative_mass_from_core_mass("Mc_HeMS.txt", core_mass, relative_mass, metalicity);
        update_relative_mass(m_rel);
        last_update_age = helium_main_sequence_time_for_solar_metalicity(relative_mass);
        relative_age = last_update_age;
        //for orderliness core_mass should be updated,
        //though it also happens at end of constructor
        evolve_core_mass();
        //helium_giant_core_mass(relative_age, relative_mass); // not to have the safetycheck on mc_max
    }
    
    adjust_next_update_age();  
    // (SPZ+GN: 27 Jul 2000)
    // core mass has been set in super_giant before constructor is called.
    
    if (relative_age >= next_update_age) {
        // this can happen when mc_max = 1.45*get_total_mass()-0.31        
        // or mc_max = 0.773*m_rel - 0.35 --> in this case should m_rel be updated?
        
        // happens in evolve_element
        // create_remnant(relative_mass, get_total_mass(), core_mass);  
        return;
    }

    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();       
    update();
    post_constructor();
}


helium_giant::helium_giant(helium_star & h) : single_star(h) {
    delete &h;
    
// (GN+SPZ May  4 1999) last update age is time of previous type change
    last_update_age = next_update_age;    
    adjust_next_update_age();

    update_wind_constant();    
    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();   
    
    update();
    post_constructor();

}

// Adjust radius & luminosity at relative_age
void helium_giant::instantaneous_element() {

    luminosity       = helium_giant_luminosity_core_mass_relation(relative_age, relative_mass, metalicity);
    radius           = helium_giant_radius(luminosity, relative_mass, get_total_mass(), metalicity);
    //effective_radius = max(effective_radius, radius);
    effective_radius = radius;
}

// Evolve a helium_giant upto time argument according to
// the new 2000 models.
void helium_giant::evolve_element(const real end_time) {
    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;
    
    if (relative_age<=next_update_age) {
        instantaneous_element();
        evolve_core_mass();
        small_envelope_perturbation();
        
        // if no envelope make transition to remnants
        // just as a procedure: reduce_mass with 1
        if(envelope_mass <= 0){
            reduce_mass(1.);
            return;    
        }       
    }
    else {
        create_remnant(relative_mass, get_total_mass(), core_mass);
        return;
    }
    
    update();
    stellar_wind(dt);
}

real helium_giant::nucleair_evolution_timescale() {        
    return helium_giant_end_time(relative_mass, get_total_mass());
}


void helium_giant::update() {

//  real m_tot = get_total_mass();
//  core_mass = COcore_mass = CO_core_mass();
//  envelope_mass = m_tot - core_mass;

  // removed (SPZ+GN:10 Nov 1998)
  //core_radius = helium_core_radius();
  //added (ST: 15 Sept 2009)
    core_radius = co_core_radius();
    
  detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
  effective_radius = radius;
}

// should only be used in constructors.
// (SPZ+GN:30 Sep 1998)
// ST:8 Sep 2009 can also be used in other functions.
void helium_giant::adjust_next_update_age() {

    //next_update_age /= cnsts.parameters(helium_star_lifetime_fraction);
    
    real t_Heg = helium_giant_end_time(relative_mass, get_total_mass());
    
//    if(relative_age>t_Heg) {
//        cerr << "WARNING: relative_age > t_Heg in helium_giant"<<endl;
//        relative_age = t_Heg;
//    }

    next_update_age = t_Heg;
    
}

real helium_giant::helium_giant_end_time(const real mass, const real mass_tot) {
    
    //stop evolving star if the core mass reaches
    // max(m_Ch = cnsts.parameters(Chandrasekar_mass),
    // mc_ignite_CO = 0.773*relative_mass - 0.35) or
    // min(get_total_mass(), 1.45*get_total_mass()-0.31)
    real mc_max = min(mass_tot, 1.45*mass_tot-0.31);
    mc_max = min(maximum_helium_giant_core_mass(mass), mc_max);
    
    
    //core should minimally grow 5% as a helium giant
    real t_hems = helium_main_sequence_time_for_solar_metalicity(mass);   
    real mco_hems = helium_giant_core_mass(t_hems, mass);  //co core mass
    mc_max = max(mc_max, 1.05*mco_hems);
    
    real t_Heg = helium_giant_age_core_mass_relation(mc_max, mass);
    real t_Hems  = helium_main_sequence_time_for_solar_metalicity(mass);
    if (t_Hems > t_Heg){
        cerr<<"In helium_giant_end_time t_Hems > t_Heg"<<endl;
        exit(-1);
    }
    
    return t_Heg;
}


void helium_giant::create_remnant(const real mass, const real mass_tot, const real mc_core) {
#if 0
     if (is_binary_component()) 
       get_binary()->dump("binev.data", false);
     else
       dump("binev.data", false);

        stellar_type type = NAS;
     //if (get_total_mass() >= cnsts.parameters(helium2black_hole)) 
//           type = Black_Hole;

     if (core_mass >= cnsts.parameters(COcore2black_hole)) 
       type = Black_Hole;
     else if(core_mass >= cnsts.parameters(Chandrasekar_mass))
       type = Neutron_Star;
     else
       type = Carbon_Dwarf;
#endif
#if 0 //(SPZ+GN: 27 Jul 2000)
     // (GN+SPZ May  3 1999) core mass > 1.4 or total mass > 2.2 
        else if(core_mass >= cnsts. //(GN+SPZ May  3 1999) was 1.3 * cnsts
		             parameters(Chandrasekar_mass) ||
		get_total_mass() >= cnsts.parameters(helium2neutron_star))
	  type = Neutron_Star;
	
	// else if(get_total_mass()>=
	//         1.3*cnsts.parameters(kanonical_neutron_star_mass) &&
	//         relative_mass>=8) 
	// type = Disintegrated;
        else
           type = Carbon_Dwarf;
#endif // (SPZ+GN: 27 Jul 2000)

    stellar_type type;
    real mc_SN = maximum_helium_giant_core_mass(mass);
    
    if (mc_core < cnsts.parameters(Chandrasekar_mass)){
        // if mc_core equals 1.45*get_total_mass()-0.31
        // shell burning stops before whole envelope is converted into C and O
        if (mc_core < mass_tot){
            cerr<<"Warning: not homogeneous WD"<<endl;
            type = Carbon_Dwarf;
        }
        // mc_core equals total mass
        // core mass reaches outside of star, no envelope anymore
        else if (relative_mass < 1.6)
            type = Carbon_Dwarf;
        else if (relative_mass <= 2.25)
            type = Oxygen_Dwarf;
        else {
             type = Oxygen_Dwarf;
        }
    }
    else {    
        if (relative_mass < 1.6) 
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
                break;
        case Neutron_Star : star_transformation_story(Neutron_Star);
                new neutron_star(*this);
                break;
        case Disintegrated : star_transformation_story(Disintegrated);
			    new disintegrated(*this);
			    break;
        case Carbon_Dwarf : star_transformation_story(Carbon_Dwarf);
	            new white_dwarf(*this, Carbon_Dwarf);
                break;
        case Oxygen_Dwarf : star_transformation_story(Oxygen_Dwarf);
                new white_dwarf(*this, Oxygen_Dwarf);
                break;
            
        default :   cerr << "helium_giant::create_remnant()" <<endl;
                    cerr << "star_type not recognized." << endl;
                    exit(-1);
       }

}


star* helium_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot);
      
    if (mdot<envelope_mass){
        envelope_mass -= mdot;
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass
        adjust_next_update_age();
    }
      else {
        mdot = envelope_mass;
        envelope_mass = 0.;	
        // if (core_mass <= cnsts.parameters(minimum_helium_star)) {
        //    star_transformation_story(Carbon_Dwarf);
        //    return dynamic_cast(star*, new white_dwarf(*this));
        // }
          
        if (core_mass >= cnsts.parameters(Chandrasekar_mass)){
            real mc_SN = maximum_helium_giant_core_mass(relative_mass);
            if (mc_SN <= 7.){
                star_transformation_story(Neutron_Star);
                return dynamic_cast(star*, new neutron_star(*this));
            }
            else{
                star_transformation_story(Black_Hole);
                return dynamic_cast(star*, new black_hole(*this));
            }
        }         
        else if(relative_mass >= 1.6) {      
              star_transformation_story(Oxygen_Dwarf);
              return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
        }
        else {
              star_transformation_story(Carbon_Dwarf);	   
              return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
        }          
      }
      return this;
}

star* helium_giant::reduce_mass(const real mdot) {   
    if (mdot < envelope_mass){
        envelope_mass -= mdot;
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass
        adjust_next_update_age();        
    }
    else {
 //       real rest_mdot = mdot;
//        rest_mdot -= envelope_mass;
//        envelope_mass = 0.;
//	   
//        if (core_mass>rest_mdot) {
//            if (core_mass-rest_mdot<=
//                cnsts.parameters(minimum_helium_star)) {
//                core_mass -= rest_mdot;
//                COcore_mass = core_mass;
//            }
//            else {
//                core_mass -= rest_mdot;
//                COcore_mass = core_mass;
//            }
//        }
//        else {
//            cerr << "ERROR!:"<<endl;
//            cerr << "void helium_giant::reduce_mass(mdot="
//                << rest_mdot<<")"<<endl;
//            cerr << "mdot exceeds helium core mass ("<<core_mass
//                << ")"<<endl;
//            cerr << "Decision: Disintegrate helium star!"<<endl;
//	  
//            star_transformation_story(Disintegrated);
//            return dynamic_cast(star*, new disintegrated(*this));
//        }
        envelope_mass = 0.;
        
        if (core_mass >= cnsts.parameters(Chandrasekar_mass)){
            real mc_SN = maximum_helium_giant_core_mass(relative_mass);
            if (mc_SN <= 7.){
                star_transformation_story(Neutron_Star);
                return dynamic_cast(star*, new neutron_star(*this));
            }
            else{
                star_transformation_story(Black_Hole);
                return dynamic_cast(star*, new black_hole(*this));
            }
        }
        else if(relative_mass >= 1.6) {   
            star_transformation_story(Oxygen_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
        }
        else {
            star_transformation_story(Carbon_Dwarf);	   
            return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
            
        }     
    }
    return this;
}

#if 0
// bool hydrogen not used currently
real helium_giant::add_mass_to_accretor(const real mdot, bool hydrogen) {
    cerr<<"He_g::add_mass_to_accretor is used?"<<endl;

        if (mdot<0) {

           cerr << "helium_star::add_mass_to_accretor(mdot=" 
		<< mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   return 0;

        }

        adjust_accretor_age(mdot);
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	// next_update_age should nog be altered here (SPZ+GN: 3 Oct 1998)

// (GN+SPZ May  5 1999) wind_constant must be zero-age envelope mass, 
// not current
//	update_wind_constant();
	wind_constant += mdot;

	set_spec_type(Accreting);
	
        return mdot;

     }

// bool hydrogen not used currently
real helium_giant::add_mass_to_accretor(real mdot, const real dt, bool hydrogen) {
    cerr<<"He_g::add_mass_to_accretor is used?"<<endl;

        if (mdot<0) {

	  cerr << "helium_giant::add_mass_to_accretor(mdot=" 
	       << mdot << ")" <<endl;
	  cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        }

        mdot = accretion_limit(mdot, dt);
        adjust_accretor_age(mdot);
        envelope_mass += mdot;
    
// (GN+SPZ May  5 1999) wind_constant must be zero-age envelope mass, 
// not current
//	update_wind_constant();
	wind_constant += mdot;
	relative_mass = max(relative_mass, get_total_mass());
    
	set_spec_type(Accreting);
	
        return mdot;

     }

#endif 


//		general mass transfer utilities.
// Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real helium_giant::add_mass_to_accretor(const real mdot, bool hydrogen) {
    
    if (mdot<=0) {
        cerr << "helium_giant::add_mass_to_accretor(mdot=" << mdot << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        
        set_spec_type(Accreting, false);
        
        return 0;
    }
        
    if(hydrogen){
        // hydrogen accretion
        cerr<<"hydrogen accretion on helium star! adjust routine"<<endl;

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
        //for the moment assume helium accretion

        // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
        //adjust_accretor_age(mdot);
        envelope_mass += mdot;
        accreted_mass += mdot;
        
        // only neccessary for AGB & He giant accretor as  
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass      
        adjust_next_update_age();  
            
    }
    
    set_spec_type(Accreting);    
    return mdot;
}

real helium_giant::add_mass_to_accretor(real mdot, const real dt, bool hydrogen) {
    if (mdot<0) {
        cerr << "helium_giant::add_mass_to_accretor(mdot=" << mdot 
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
        //for the moment assume helium accretion
        
        // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
        //adjust_accretor_age(mdot);
        envelope_mass += mdot;
        accreted_mass += mdot;
        
        // only neccessary for AGB & He giant accretor as  
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass      
        adjust_next_update_age();  
        
    }
    set_spec_type(Accreting);
    return mdot;
}



real helium_giant::accretion_limit(const real mdot, const real dt) {

     real mdot_limit = mdot;

     real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;
     if (mdot>=eddington)
       mdot_limit =  eddington;

     if (envelope_mass<=0)
       mdot_limit = 0.;

     return mdot_limit;

}

// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void helium_giant::adjust_accretor_age(const real mdot,
				       const bool rejuvenate) {
    
    real m_rel_new;
    real m_tot_new = get_total_mass() + mdot;
    if (m_tot_new>relative_mass)
        m_rel_new = m_tot_new;
    else m_rel_new = relative_mass;
    
    real t_hems_old = helium_main_sequence_time_for_solar_metalicity(relative_mass); 
    real dt_Heg_old = helium_giant_end_time(relative_mass, get_total_mass()) - t_hems_old;

    real t_hems_new = helium_main_sequence_time_for_solar_metalicity(m_rel_new);    
    //He star+giant not a function of metalicity for now
    real t_Heg_new = helium_giant_end_time(m_rel_new, m_tot_new);
    real dt_Heg_new = t_Heg_new - t_hems_new; 
    
    real dtime = relative_age - t_hems_old;
    
    // (GN+SPZ May  4 1999) update last_update_age
    last_update_age = t_hems_new;
    relative_age = t_hems_new 
    + dtime*(dt_Heg_new/dt_Heg_old);
    
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
    
    if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
        cerr<<"In helium_giant::adjust_accretor_age relative age updated on HeGiant, but < last_update_age"<<endl;
    }
    relative_age = max(relative_age, 
                       last_update_age + cnsts.safety(minimum_timestep)); 
    relative_age = min(relative_age, t_Heg_new);

    // next_update_age should not be reset here
    // next_update_age = t_ms_new + t_hg_new;
    
//    real frac = (1-pow(mdot/(get_total_mass()+mdot),
//			     cnsts.parameters(rejuvenation_exponent)));
//      last_update_age *= frac;
//      relative_age *= frac;
}

// see helium_star.C
real helium_giant::zeta_adiabatic() {
    cerr<<"helium_giant::zeta_adiabatic: include difference for He HG and He SG star"<<endl;

     real z = 0.;
//      Hjellming and Webbink 1987 ApJ, 318, 804
     real x = core_mass/get_total_mass();
     real A = -0.220823;
     real B = -2.84699;
     real C = 32.0344;
     real D = -75.6863;
     real E = 57.8109;

     if (get_total_mass()<=0.4)
        z = -cnsts.mathematics(one_third);
     else
        z = A + x*(B + x*(C + x*(D + x*E)));

     return z;

   }

real helium_giant::zeta_thermal() {
    cerr<<"helium_giant::zeta_thermal: include difference for He HG and He SG star"<<endl;

    real z = -2;

    if (get_core_mass()<=0.4)  // like a white dwarf
	z = -cnsts.mathematics(one_third);

    return z;

}

#if 0
real helium_giant::CO_core_mass() {
      // C/O core of helium star grows linearly with time
      // (SPZ+GN:26 Sep 1998)

// (GN+SPZ May  4 1999) core groth totally in helium_star phase: core constant
//  real m_core = get_total_mass()*relative_age/next_update_age;
//  m_core = max(core_mass,m_core);
//    return min(m_core, get_total_mass());

  return min(core_mass, get_total_mass());
}
#endif
#if 0
void helium_giant::stellar_wind(const real dt) {

//  PRL(last_update_age);
//  PRL(next_update_age);
//  PRL(relative_age);
//  PRL(previous.relative_age);
//  PRL(dt);
// (GN+SPZ Apr 28 1999) wind for low mass stars per phase
    real end_time = next_update_age - last_update_age;
    real prev_rel_time = max(0.,previous.relative_age - last_update_age);
    real relative_time = min(relative_age - last_update_age, end_time);

//    PRL(end_time);
//    PRL(relative_time);
//    PRL(prev_rel_time);
//    PRL(wind_constant);
    
    real wind_mass = wind_constant 
                   * (pow(relative_time/end_time,
			cnsts.parameters(massive_star_mass_loss_law))
	           -  pow((relative_time-dt)/end_time,
			cnsts.parameters(massive_star_mass_loss_law)));

// (GN+SPZ May  6 1999) try low wind: 
//    wind_mass = 0.;

//    PRL(wind_mass);
//    PRL(envelope_mass);

  if (wind_mass>=envelope_mass) {
    wind_mass = envelope_mass;
    effective_radius = radius = core_radius;
  }

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_mass(wind_mass);
  return;
}
#endif

real helium_giant::gyration_radius_sq() {
    cerr<<"He_s::gyration_radius_sq is used?"<<endl;

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

// (GN+SPZ May  3 1999) helium_giants loose complete envelope in wind
void helium_giant::update_wind_constant() {

//  wind_constant = (1 - cnsts.parameters(helium_star_final_core_fraction))
//                * get_total_mass(); 

// (GN+SPZ May  7 1999) envelope is about 30% of total mass,
// we loose 10% of total mass ....
//  wind_constant = 0.3*envelope_mass;
    
// (SilT Nov 24 2009) follow eq 8.4 in Gijs' thesis Chapter 8
    //based on Nugis & Lamers
 wind_constant = 1.38E-08 * pow(get_total_mass(), 2.87); 


}

stellar_type helium_giant::get_element_type() {

     stellar_type type = Helium_Giant;
//     if (envelope_mass <= 0)
//	 type = Carbon_Star;

     return type;
     }

real helium_giant::temperature() {

  real T_eff = cnsts.parameters(Tsun)
             * sqrt(sqrt(luminosity)/effective_radius);
  
  return T_eff;

  //return sqrt(33.45*sqrt(luminosity)/effective_radius);
}


bool helium_giant::giant_star() {
    
    return TRUE;
}

bool helium_giant::remnant() {

    return FALSE;
}


//Eq. 84
real helium_giant::helium_giant_luminosity_core_mass_relation(const real time, const real mass, const real z){
    real A_He = AGB_A_He_estimator();
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real l_x = helium_giant_x_luminosity(mass);
    real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
    real l_He;
    if (time  <= t_x){
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        real arg = (p-1)*A_He*D*(t_inf1-time);   
        l_He = D * pow(arg, p/(1-p));      
    }
    else {        
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real t_inf2 = specific_time_limit(A_He, t_x,
                                          B, l_x, q);
        real arg = (q-1)*A_He*B*(t_inf2-time);   
        l_He = B * pow(arg, q/(1-q));
    }  
    return l_He; 
}


real helium_giant::helium_giant_age_core_mass_relation(const real m_core, const real mass){
    real age;  
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real A_He    = AGB_A_He_estimator();
    real t_Hems  = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    
    real mx = helium_giant_x_mass(mass);
    if (m_core <= mx){        
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        age = t_inf1 - pow(m_core, 1.-p)/A_He/D/(p-1.);
    }
    else{
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real l_x = helium_giant_x_luminosity(mass);
        real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
        real t_inf2 = specific_time_limit(A_He, t_x,
                                           B, l_x, q);
        age = t_inf2 - pow(m_core, 1.-q)/A_He/B/(q-1.);
    }
    return age;
}


void helium_giant::evolve_core_mass(const real time,
                    const real mass, const real mass_tot) {
    real mco_core = helium_giant_core_mass(time, mass);    
    
//    real mc_SN = maximum_helium_giant_core_mass(mass);
//    real mc_max = min(mass_tot, 1.45*mass_tot-0.31);
//    mc_max = min(mc_SN, mc_max);
//    
//    if (mco_core > mc_max + cnsts.safety(tiny)){
//        cerr<<"helium_giant::helium_giant_core_mass m_core > mc_max"<<endl;
//        cerr.precision(HIGH_PRECISION);
//        PRC(mco_core);PRC(mc_max);PRC(mass_tot);PRL(relative_mass);
//        PRL(time);
//        cerr.precision(STD_PRECISION);
//        exit(-1);
//    }
//    
//    mco_core = min(mco_core, mc_max);
    
    
    if(!update_COcore_mass(mco_core)) {
        cerr << "Update COcore mass failed in helium_giant()"<<endl;
    }
    if(!update_core_and_envelope_mass(mco_core)) {
        cerr << "Update core mass failed in helium_giant()"<<endl;
    }
}

void helium_giant::evolve_core_mass() {
    evolve_core_mass(relative_age, relative_mass, get_total_mass());
}

//related to Eq. 84
real helium_giant::helium_giant_core_mass(const real time, const real mass){
    
    
    real A_He = AGB_A_He_estimator();
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real l_x = helium_giant_x_luminosity(mass);
    real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
    real m_core;
    
    if (time  <= t_x){
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        real arg = (p-1)*A_He*D*(t_inf1-time);   
        m_core = pow(arg, 1.0/(1-p));     
    }
    else {
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real t_inf2 = specific_time_limit(A_He, t_x,
                                          B, l_x, q);
        real arg = (q-1)*A_He*B*(t_inf2-time);   
        m_core = pow(arg, 1.0/(1-q));      
    }  
    return m_core;     
}


// Eq.75
real helium_giant::maximum_helium_giant_core_mass(const real mass) {
    real m_Ch = cnsts.parameters(Chandrasekar_mass);
    real mc_max = max(m_Ch, 0.773*mass - 0.35);
    return mc_max;
}


void helium_giant::small_envelope_perturbation(){
    real mu = small_envelope_mu(luminosity, get_total_mass(), core_mass);
    if(mu < 1.){
        real lum_c = small_envelope_core_luminosity();
        luminosity = perturb_luminosity(luminosity, lum_c, get_total_mass(), core_mass, mu);
        real rad_c = small_envelope_core_radius();
        if(rad_c < radius){
            radius = perturb_radius(radius, rad_c, get_total_mass(), core_mass, mu);
        }
    }
}


//Eq.98
real helium_giant::small_envelope_mu(const real lum, const real mass_tot, const real m_core){
    // the function parameter lum is not needed here, 
    // but needed in single_star::small_envelope_mu
    // in this way the correct function is called automatically
    
    real mc_max = min(mass_tot, 1.45*mass_tot-0.31);    
    real mu = 5.*(mc_max-m_core) / mc_max;
    return mu;
}

real helium_giant::co_core_radius(const real m_core){
    // due to small nucleair burning layer 
    // r_c > white_dwarf_radius
    real r_c = 5.*white_dwarf_radius(m_core, 0.);
    return r_c;
}

real helium_giant::co_core_radius(){
    return co_core_radius(core_mass);
}


real helium_giant::small_envelope_core_radius(const real m_core){
    real r_c = white_dwarf_radius(m_core, 0. );
    return r_c;
}

real helium_giant::small_envelope_core_radius(){
    return small_envelope_core_radius(core_mass);
}

real helium_giant::small_envelope_core_luminosity(){
    real l_c = 40.;
    return l_c;
}    





