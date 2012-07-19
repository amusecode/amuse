//
// helium_star.C
//

#include "helium_star.h"
#include "main_sequence.h"
#include "hyper_giant.h"
#include "sub_giant.h"
#include "horizontal_branch.h"
#include "super_giant.h"
#include "hertzsprung_gap.h"

#if 0
helium_star::helium_star(main_sequence & m) : single_star(m) {

    delete &m;

    lose_envelope_decent();

    adjust_next_update_age();
// (GN+SPZ May  4 1999) last update age is time of previous type change
// in NEW relative_age units!!!
    last_update_age = 0;
    relative_age = 0;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = 0;
    final_core_mass = final_CO_core_mass(m_tot);
    core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();

    post_constructor();
    
}
#endif
#if 0
helium_star::helium_star(hyper_giant & w) : single_star(w) {

    delete &w;

    lose_envelope_decent();


//  PRC(relative_age);PRC(main_sequence_time());PRL(next_update_age);
// fraction of Hyper_Giant phase spend before stripped
    real t_frac = (relative_age - main_sequence_time())
               / (nucleair_evolution_time()-main_sequence_time());

//    PRC(nucleair_evolution_time());PRL(t_frac);
    adjust_next_update_age();
    relative_age = t_frac*next_update_age;
//    PRC(relative_age);PRL(next_update_age);
    
// (GN+SPZ May  4 1999) last update age is time of previous type change
// in NEW relative_age units!!
    last_update_age = relative_age;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = 0;
    final_core_mass = final_CO_core_mass(m_tot);
    core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();
    
    post_constructor();

}
#endif

helium_star::helium_star(hertzsprung_gap & h) : single_star(h) {

    delete &h;
    lose_envelope_decent();
    relative_mass = get_total_mass();
    
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age units
    last_update_age = 0.;
    adjust_next_update_age();
    relative_age = 0;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();
    
    instantaneous_element();
    update();
    
    post_constructor();

}

helium_star::helium_star(sub_giant & g) : single_star(g) {

    delete &g;
    lose_envelope_decent();
    relative_mass = get_total_mass();
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age units
    last_update_age = 0.;
    adjust_next_update_age();
    relative_age = 0;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();
    
    post_constructor();

}

helium_star::helium_star(horizontal_branch & h) : single_star(h) {

    delete &h;
    lose_envelope_decent();
    relative_mass = get_total_mass();
    

//    real t_ms = main_sequence_time();
//    real t_giant = t_ms + hertzsprung_gap_time()
//                 + base_giant_branch_time();
//    real t_he = helium_giant_time(t_ms, metalicity);
//    real t_frac = min(0.9,(relative_age - t_giant)/t_he);

    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age uits
    adjust_next_update_age();
    //relative_age = t_frac* next_update_age;
    relative_age = (relative_age - helium_ignition_time(relative_mass, metalicity)) / core_helium_burning_timescale(relative_mass, metalicity)
    * helium_main_sequence_time_for_solar_metalicity(relative_mass);
    last_update_age = relative_age;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();
    
    post_constructor();

 }

#if 0
void helium_star::adjust_initial_star() {

  update_wind_constant();

  if(relative_age<=0)
    relative_age = max(current_time, 0.0);
}
#endif

// Use only in constructors.
void helium_star::adjust_next_update_age() {

  //next_update_age = cnsts.parameters(helium_star_lifetime_fraction)
  //                * helium_time();
  
    next_update_age = helium_main_sequence_time_for_solar_metalicity(relative_mass);
}

#if 0
void helium_star::instantaneous_element() {

  
    real temp;
    real m_tot = get_total_mass();

    temp = log10(m_tot)
      * (0.4509 - 0.1085*log10(m_tot)) + 4.7143;
    luminosity  = 8.33*temp - 36.8;
    temp = 1.e-3*pow(10., temp);
    luminosity  = pow(10., luminosity);
    effective_radius = radius = 33.45*sqrt(luminosity)/(temp*temp);
    core_radius  = 0.2*radius;

    if (envelope_mass <= 0) {
      effective_radius = radius  = core_radius;
    }
}
#endif 
#if 0
void helium_star::evolve_element(const real end_time)
{
    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;
    
    real m_tot = get_total_mass();
    if (relative_age>next_update_age) {
      
      //core_mass = COcore_mass = CO_core_mass();
      envelope_mass = m_tot - core_mass;
	
      stellar_wind(dt);

      star_transformation_story(Helium_Giant);
      new helium_giant(*this);
      return;
    }
    else if (m_tot < cnsts.parameters(helium_dwarf_mass_limit) &&
	     relative_mass < cnsts.parameters(
			     upper_ZAMS_mass_for_degenerate_core)) {

      // low-mass degenerate core has no envelope (SPZ+GN: 4 Oct 1998)
      core_mass = get_total_mass();
      envelope_mass = 0;

      if (is_binary_component()) 
	get_binary()->dump("binev.data", false);
      else
	dump("binev.data", false);
      
      star_transformation_story(Helium_Dwarf);
      new white_dwarf(*this);
      return;
    }
    else { // normal helium star evolution
	real tmp = log10(m_tot)
	* (0.4509 - 0.1085*log10(m_tot)) + 4.7143;
	luminosity  = pow(10., 8.33*tmp - 36.8);
	tmp = pow(1.e-3*pow(10., tmp), 2);
	radius = 33.45*sqrt(luminosity)/tmp;
	core_radius  = 0.2*radius;

// (GN+SPZ May  3 1999) according to Savonije et al 1986
// semi-degenerate helium_stars; have r = 0.029 m^-0.19 and live forever
	if (m_tot < 0.2) { 
	  radius = 0.029 * pow(m_tot, -0.19);
	  next_update_age = relative_age + cnsts.safety(maximum_timestep);
	}

	if (envelope_mass <= 0) {
	    radius  = core_radius;
	}
    }

    //core_mass = COcore_mass = CO_core_mass();
    
    envelope_mass = m_tot - core_mass;
    //core_radius = helium_core_radius();
    
    update();
    stellar_wind(dt);

}
#endif

// Adjust radius & luminosity at relative_age
void helium_star::instantaneous_element() {
    luminosity  = helium_main_sequence_luminosity(relative_age, relative_mass);
    
    radius      = helium_main_sequence_radius(relative_age, relative_mass, get_total_mass());

    //effective_radius = max(effective_radius, radius);
    effective_radius = radius;
}

// Evolve a helium_star upto time argument according to
// the new 2000 models.
void helium_star::evolve_element(const real end_time) {
    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;
        
    if (relative_mass < cnsts.parameters(minimum_helium_star)) {
        // Helium Main_sequence star will not ignite core helium burning.
        star_transformation_story(Helium_Dwarf);
        new white_dwarf(*this);
        return;
    }
    
    if (relative_age <= next_update_age) {
        instantaneous_element();
    } 
    else {
        star_transformation_story(Helium_Giant);
        new helium_giant(*this);
        return;
    }
    
    update();
    //stellar_wind(dt);
}


void helium_star::update() {
    
    detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
    effective_radius = radius;

}

#if 0
void helium_star::create_remnant() {

/*
//cout<<"void helium_star::create_remnant() called"<<endl;
       // core_mass = get_total_mass();
        stellar_type type = NAS;
        if (get_total_mass()>=cnsts.parameters(helium2black_hole))
           type = Black_Hole;
//			was 1.3*M_NS
        else if((core_mass>=1.2*cnsts.mass.neutron_star     &&
                 relative_mass>=10) ||
                get_total_mass()>=cnsts.mass.helium2neutron_star)
           type = Neutron_Star;
//        else if(get_total_mass()>=1.3*M_NS && relative_mass>=8) 
//           type = Disintegrated;
        else
           type = White_Dwarf;

//        base_element * s = make(type, this);
       switch (type) {
          case Black_Hole : if (is_binary_component())
	                       get_binary()->dump("sn.dat");
                            else dump("sn.dat");
			    
			    star_transformation_story(this,Black_Hole);
                            new black_hole(*this); 
                            break;
          case Neutron_Star : 
                           if (is_binary_component()) get_binary()->dump("sn.dat");
                           else dump("sn.dat");
			   star_transformation_story(this,Neutron_Star);
			   new neutron_star(*this);
                           break;
          case Disintegrated : if (is_binary_component()) get_binary()->dump("sn.dat");
                               else dump("sn.dat");
			       star_transformation_story(this,Disintegrated);
                               new disintegrated(*this);
                               break;
         case White_Dwarf : star_transformation_story(this,White_Dwarf);
	                    new white_dwarf(*this);
			    break;
          default : ;//error << "helium_star::create_remnant()" <<endl;
                    //error << "star_type not recognized." << endl;
                    //error << "Action: Proceed normally." << endl;
       }

*/
     }
#endif

/*
real helium_star::helium_time() {

//              Helium Star lifetime from
     real t_ms = main_sequence_time();
     real t_he;
     real m_tot = get_total_mass();

     if (exp_HE_lifetime) {
//              Iben I., 1967 ??; en Iben and Tutukov, 1985, ApJSS, 58, 661.
        t_he =  0.56*t_ms/pow(relative_mass, 0.52);
     }
     else {
//              Pols O., 1993, proefschrift, P13, Eq. 2.15
     if (m_tot<=0.7) t_he = 10.76*pow(m_tot, -3.75);
     else if(m_tot<=1.6) t_he = 17.1*pow(m_tot, -2.45);
     else if(m_tot<=4.8) t_he = 11.48*pow(m_tot, -1.6);
     else                    t_he = 2.37*pow(m_tot, -0.6);
     }

     if(t_he>t_ms) t_he=t_ms;
     return t_he;
}
*/

star* helium_star::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = relative_mass*dt/get_binary()->get_donor_timescale();

    mdot = mass_ratio_mdot_limit(mdot);

    if (envelope_mass >= mdot)
      envelope_mass -= mdot;
    else {
      mdot = envelope_mass;
      envelope_mass = 0;
      
      //	   star_transformation_story(this,White_Dwarf);
      //return dynamic_cast(star*, new white_dwarf(*this));
    }

    return this;
}

star* helium_star::reduce_mass(const real mdot) {

    adjust_age_after_wind_mass_loss(mdot, true);
    envelope_mass -= mdot;
    if (relative_mass > get_total_mass()){
        update_relative_mass(get_total_mass());
    }
    
    if (relative_mass < cnsts.parameters(minimum_helium_star)) {
        // Helium Main_sequence star will not continue core helium burning.
        star_transformation_story(Helium_Dwarf);
        return dynamic_cast(star*, new white_dwarf(*this));
    }
    
    //   On the He-MS there is no defined core yet, 
    //   no CO star/WD can form yet 
    //   if (envelope_mass<=mdot) {
    //         envelope_mass = 0;
    //         //star_transformation_story(Helium_Star);
    //         //return dynamic_cast(star*, new helium_star(*this));
    //      }
    
//    if (envelope_mass >= mdot)
//        adjust_age_after_wind_mass_loss(mdot, true);
//        envelope_mass -= mdot;
//        if (relative_mass > get_total_mass()){
//            update_relative_mass(get_total_mass());
//        }
//    else {
//      cerr<<"In helium_star::reduce_mass: what does this part do? only for constructing he star?"<<endl;
//      real mass_reduced = mdot;
//      mass_reduced -= envelope_mass;
//      envelope_mass = 0;
//      if (core_mass>mass_reduced) {
//        if (core_mass-mass_reduced<=cnsts.parameters(minimum_helium_star)) {
//            core_mass -= mass_reduced;
//    //		 star_transformation_story(this,White_Dwarf);
//    //                 new white_dwarf(*this);
//    //                 return;
//	    }
//        else {
//          core_mass -= mass_reduced;
//          COcore_mass = core_mass;
//        }
//      }
//      else {
//        cerr<<"ERROR!:"<<endl;
//        cerr<<"void helium_star::reduce_mass(mdot="
//            <<mass_reduced<<")"<<endl;
//        cerr<<"mdot exceeds helium core mass ("<<core_mass
//            <<")"<<endl;
//        cerr<<"Decision: Disintegrate helium star!"<<endl;
//        //new disintegrated(*this);
//        //return;
//       }
//    }
    
    return this;
}


void helium_star::adjust_age_after_wind_mass_loss(const real mdot,
                                                    const bool rejuvenate=true) {
    real m_tot_new = get_total_mass() - mdot;
    real m_rel_new = m_tot_new;
    real t_hems_old = helium_main_sequence_time_for_solar_metalicity(relative_mass); 
    real t_hems_new = helium_main_sequence_time_for_solar_metalicity(m_rel_new);
    
    relative_age *= (t_hems_new/t_hems_old);
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(-1.0 * mdot/m_tot_new); 
    
    // next_update_age should not be reset here,
    // is done in add_mass_to_accretor, where also relative_mass
    // is updated
    // next_update_age = t_ms_new; 
}



real helium_star::add_mass_to_accretor(const real mdot) {

        if (mdot<0) {
           cerr << "helium_star::add_mass_to_accretor(mdot=" 
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
	   return 0;
        }

        adjust_accretor_age(mdot);
        if (relative_mass<get_total_mass() + mdot)
	  relative_mass = get_total_mass() + mdot;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//	update_wind_constant();
	
        envelope_mass += mdot;
	
	set_spec_type(Accreting);
	
        return mdot;

}

real helium_star::add_mass_to_accretor(real mdot, const real dt) {

        if (mdot<0) {
           cerr << "helium_star::add_mass_to_accretor(mdot=" 
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
	   return 0;
        }

        mdot = accretion_limit(mdot, dt);
        adjust_accretor_age(mdot);
        if (relative_mass<get_total_mass() + mdot)
           relative_mass = get_total_mass() + mdot;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//	update_wind_constant();

        envelope_mass += mdot;

	set_spec_type(Accreting);

        return mdot;

     }

real helium_star::accretion_limit(const real mdot, const real dt) {
       
        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington) return eddington;

        return mdot;
}



void helium_star::adjust_accretor_age(const real mdot,
				      const bool rejuvenate) {

    relative_age *= (1-pow(mdot/(get_total_mass()+mdot),
		    cnsts.parameters(rejuvenation_exponent)));

}

// Helium stars have radiative envelopes.
// this requires more research (SPZ+GN: 3 Oct 1998)
real helium_star::zeta_adiabatic() {

  real z = 15;
  
  return z;

}

real helium_star::zeta_thermal() {

// (GN+SPZ May  3 1999) see evolve_element; mass - radius relation
  if (get_total_mass() < 0.2 ) {
    return -0.19;
  } else {
    return 1;
  }
}
#if 0
real helium_star::final_CO_core_mass(const real initial_mass) {

  // Implementation of Nelemans YPZV 2000 (A&A submitted)
  // Itroduces discontinuity at relative_mass = 2.2
  // bases on Habets 1986 & IT85
  // (SPZ+GN: 27 Jul 2000)
  real final_coremass_fraction;
  if(initial_mass <= 0.8) 
    final_coremass_fraction = 1;
  else if(initial_mass >= cnsts.parameters(helium2neutron_star)) 
    final_coremass_fraction = 0.65;
  else 
    final_coremass_fraction = 1 - 0.32 * (initial_mass - 0.8);

  return final_coremass_fraction*initial_mass;
}
#endif
#if 0
real helium_star::CO_core_mass() {

  real m_core = final_core_mass * relative_age/next_update_age;
  m_core = max(core_mass, m_core);

  return min(m_core, get_total_mass());
}
#endif

void helium_star::stellar_wind(const real dt) {

//cerr<<"void helium_star::stellar_wind(dt="<<dt<<") = ";

	    real kappa = pow(get_total_mass(),2.5);

//	real kappa = wind_constant;
//        real wmlc = wind_constant;
//        wind_constant = kappa;

//		WR Wind
//		Langer, N., 1989, AA 220, 135.
        if(get_total_mass()>core_mass) {
//           kappa = 0.5*(kappa + wmlc);
//           real wind_mass = 0.005*dt*kappa;
// (GN+SPZ May  3 1999) Lager wind (again)
// Langer: M_dot = 5-10 10^-8 M^2.5

	  real wind_mass = 0.05*dt*kappa;
// (GN+SPZ May 12 1999) proper integration
	  real m = get_total_mass();
	  real constant = 0.05;

	  real m_next = pow((pow(m,-1.5) + 1.5*constant*dt),-1/1.5);
	  wind_mass = m - m_next;
//	  PRC(m);PRC(m_next);PRL(wind_mass);
//           real wind_mass = 0.0025*dt*kappa;
// (GN+SPZ May  3 1999) low mass helium stars have no wind
// helium_giant looses envlope
	  if (get_total_mass() < 2.5 ) wind_mass = 0.;

           if (wind_mass>=envelope_mass) {
              wind_mass = envelope_mass;
              radius = core_radius;
           }

           if (is_binary_component())
              get_binary()->adjust_binary_after_wind_loss(
                          this, wind_mass, dt);
           else
              reduce_mass(wind_mass);
              return;
        }
   } 


real helium_star::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}


// (GN+SPZ May  3 1999) NOT USED!, but called by post_constructor
void helium_star::update_wind_constant() {

  wind_constant = (1 - cnsts.parameters(helium_star_final_core_fraction))
                      * get_total_mass();
}

stellar_type helium_star::get_element_type() {
  if (envelope_mass <= 0) 
    return Carbon_Star;
  else
    return Helium_Star;
}

real helium_star::temperature() {
  real T_eff = cnsts.parameters(Tsun)
             * sqrt(sqrt(luminosity)/effective_radius);
  return T_eff;
  //return sqrt(33.45*sqrt(luminosity)/effective_radius);
}
