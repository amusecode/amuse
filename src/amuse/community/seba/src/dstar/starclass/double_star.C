//
//
// double_star.C
//
// test if I can submit to Amuse svn

#include "double_star.h"

#define REPORT_BINARY_EVOLUTION    false
#define REPORT_RECURSIVE_EVOLUTION false
#define REPORT_FUNCTION_NAMES      false
#define REPORT_TRANFER_STABILITY   false

// (SPZ+GN: 28 Jul 2000) Obsolete
//#define MAXIMUM_BINARY_UPDATE_TIMESTEP cnsts.star_to_dyn(binary_update_time_fraction) 
// GIJS: If you want to increase the timestep for population synthesis,
//       please use the the parameter: MAXIMUM_BINARY_UPDATE_TIMESTEP
//       For example by defining:
//       #define MAXIMUM_BINARY_UPDATE_TIMESTEP 0.9
//       Which should be identical to your previous results.

double_star * new_double_star(node* n, real sma, real ecc, 
			      real start_time, int id,	// Defaults specified
                              binary_type type)		// in double_star.h
    {
     double_star* element = new double_star(n);
     element->initialize(type, sma, ecc, start_time, id);

     element->set_donor_timescale(element->get_primary(), true);
     
     element->determine_minimal_timestep();

     element->evolve_element(start_time);

     element->post_constructor();

     return element;
}

double_star::double_star(node* n) : star(n) {

           semi=eccentricity=binary_age=minimal_timestep=velocity=0;
           donor_timescale=0;
           identity = donor_identity=0;
           donor_type=NAS;
           first_contact=FALSE;
}

void double_star::initialize(binary_type type, real sma, 
  			real ecc, real age, int id) {
   
      semi         = sma;
      eccentricity = ecc;
      bin_type     = type;
      binary_age   = age;
      identity     = id;

      initial.semi         = semi;
      initial.q            = mass_ratio();
      initial.mass_prim    = get_primary()->get_total_mass();
      initial.eccentricity = eccentricity;

      if(get_primary()->get_identity()<0)
	  get_primary()->set_identity(0);
      if(get_secondary()->get_identity()<0)
	  get_secondary()->set_identity(1);
      
      refresh_memory();       //implemented 25-1-95
}

real double_star::mass_ratio()           {
     real mp = get_primary()->get_total_mass();
     real ms = get_secondary()->get_total_mass();

     real q = ms/mp;
     if (q>1) q = mp/ms;

     return q;
  }

real double_star::get_evolve_timestep() {

  real delta_t = cnsts.safety(maximum_timestep); //donor_timescale;
  // changed (again) spz:1 Jan 2006 
  //  real delta_t = get_donor_timescale();

  if (bin_type != Merged) {
    if (get_primary()->get_evolve_timestep() <
	get_secondary()->get_evolve_timestep())
      return Starlab::min(delta_t, get_primary()->get_evolve_timestep());
    else
      return Starlab::min(delta_t, get_secondary()->get_evolve_timestep());
  }
  else {
    return Starlab::min(delta_t, get_primary()->get_evolve_timestep());
  }

}

// Note Common-envelope and spiral-in are not detected.
binary_type double_star::obtain_binary_type() {

  binary_type type_of_binary = Unknown_Binary_Type;
  
  real rl_p = roche_radius(get_primary());
  real rl_s = roche_radius(get_secondary());

  real rp = get_primary()->get_effective_radius();
  real rs = get_secondary()->get_effective_radius();

  //  cerr << "obtain_binary_type: " << identity
  //       << " rp:" <<log10(rl_p/rp)
  //       << " rs:" << log10(rl_s/rs) 
  //       << " a/r* = " << semi/max(rp, rs);
  // Just for memory and output routines.
  // Primary fills Roche-lobe.
  if (semi <= 0 && eccentricity <= 0)
    type_of_binary = Merged;
  else if (eccentricity >= 1)
    type_of_binary = Disrupted;
  else if (rp >= rl_p || rs > rl_s)
    type_of_binary = Semi_Detached;
  else if (rp >= rl_p && rs >= rl_s)
    type_of_binary = Contact;
  else if (semi <= rp * cnsts.parameters(tidal_circularization_radius) ||
           semi <= rs * cnsts.parameters(tidal_circularization_radius)) 
    type_of_binary = Synchronized;
  else
    type_of_binary = Detached;

  //  cerr << " type= " << type_string(type_of_binary) << endl;
  
  return type_of_binary;
}

real double_star::roche_radius(star * str) {
  // Eggleton PP., ApJ, 1983, 268, 368.

  // (SPZ+GN: 28 Jul 2000)
  // semi-major axis van merged and disrupted is now zero (for convenience)
  // So, we require this safety...
  if (bin_type == Merged || bin_type == Disrupted)
    return VERY_LARGE_NUMBER;

  real mr = str->get_total_mass()
    / str->get_companion()->get_total_mass();
  real q1_3 = pow(mr, cnsts.mathematics(one_third));
  real q2_3 = pow(q1_3, 2);                  //pow(mr, TWO_THIRD);

  real Rl = semi*0.49*q2_3/(0.6*q2_3 + log(1 + q1_3));
  return Rl;
}

real double_star::roche_radius(const real a, const real m1, const real m2) {

  real q = m1/m2;
  real q1_3 = pow(q, cnsts.mathematics(one_third));
  real q2_3 = pow(q1_3, 2);   //pow(mr, TWO_THIRD);
  
  return a*0.49*q2_3/(0.6*q2_3 + log(1 + q1_3));
       }

void double_star::put_element() {
//              should be made starlab compatible


     ofstream outfile("binev.data", ios::app|ios::out);
     if(!outfile) cerr << "error: couldn't create file binev.data"<<endl;

          outfile << identity << "\t" << binary_age << endl;
          outfile << identity << endl;
          get_primary()->put_element();
	  get_secondary()->put_element();

     initial.put_element();
     outfile << bin_type << " " << eccentricity << " "
             << get_period() << " " << semi << " " << velocity << "\n" << endl;

     outfile.close();
        }

void double_star::put_hrd(ostream & s) {

        s << "2\n  ";
        get_primary()->put_hrd(s);
        s << "  "; 
        get_secondary()->put_hrd(s);
   }

//	Pretty print.
void double_star::print_status() {
        cout << "status at time = " << binary_age << endl;

        if (get_primary()->no_of_elements()>1) 
           get_primary()->print_status();
        if (get_secondary()->no_of_elements()>1)   
           get_secondary()->print_status();

        cout << type_string(bin_type);
        switch (bin_type) {
           case Merged: 
              cout << " {" << type_string(get_primary()->get_element_type())
                   << "}" << endl;
              cout << "M = "<<get_primary()->get_total_mass()
                   << "\t R = "<<get_primary()->get_effective_radius()
                   << "\t V = "<<get_primary()->get_velocity() << endl;
              break;
           case Disrupted:
              cout << " )" << type_string(get_primary()->get_element_type())
                   << ", " << type_string(get_secondary()->get_element_type())
                   << "(" << endl;
              cout << "M = "<<get_primary()->get_total_mass()
                   << "\t R = "<<get_primary()->get_effective_radius()
                   << "\t V = "<<get_primary()->get_velocity() << endl;
              cout << "m = "<<get_secondary()->get_total_mass()
                   << "\t r = "<<get_secondary()->get_effective_radius()
                   << "\t v = "<<get_secondary()->get_velocity() << endl;
              break;
           default:
              if (get_primary()->get_spec_type(Rl_filling))
                 cout << " [";
              else
                 cout << " (";
              cout << type_string(get_primary()->get_element_type())
                   << ", " << type_string(get_secondary()->get_element_type());
              if (get_secondary()->get_spec_type(Rl_filling))
                 cout << "]" << endl;
              else
                 cout << ")" << endl;
              cout << "M = "<<get_primary()->get_total_mass()
                   << "\t R = "<<get_primary()->get_effective_radius() << endl;
              cout << "m = "<<get_secondary()->get_total_mass()
                   << "\t r = "<<get_secondary()->get_effective_radius() << endl;
              cout << "a = "<<semi<<"\t e = "<<eccentricity
                   << "\t v = "<<get_velocity() << endl;
           }
     }

void double_star::put_state() {
  
  star *a, *b;
  if (get_use_hdyn()) {
   a = get_primary();
   b = get_secondary();
  }
  else {
    a = get_initial_primary();
    b = get_initial_secondary();
  }

        switch (bin_type) {
           case Merged:
              cerr << " {"; 
              get_primary()->put_state();
              cerr << "}" << endl;
              break;
           case Disrupted:
              cerr << " )";
              a->put_state();
              cerr << ", "; 
              b->put_state();
              cerr << "(" << endl;
              break;
           default:
              if (a->get_spec_type(Rl_filling))
                 cerr << " [";
              else
                 cerr << " (";
              a->put_state();
              cerr << ", ";
              b->put_state();
	      if (b->get_spec_type(Rl_filling))
		cerr << "]";
	      else
		cerr << ")";

	     real pday = get_period();
	      if (pday<0.04) 
		cerr << "   Porb = " << pday*1440 << " minutes"<<endl;
	      else if (pday<1) 
		cerr << "   Porb = " << pday*24 << " hours"<<endl;
	      else if (pday<365.25) 
		cerr << "   Porb = " << pday << " days"<<endl;
	      else
		cerr << "   Porb = " << pday/365.25 << " years"<<endl;	      
           }
     }

void double_star::print_roche() {
//	Gives compatible output for EPJ Roche program.
  
     if (bin_type!=Merged) {
       get_initial_primary()->print_roche();
       get_initial_secondary()->print_roche();
     }
     else get_primary()->print_roche();

     cerr<<get_period();
   }

void double_star::dump(ostream & s, bool brief) {
//	Zero IQ dump program.

  if (brief) {
    star* stp = get_initial_primary();
    star* sts = get_initial_secondary();

    real m1 = stp->get_total_mass();
    real m2 = sts->get_total_mass();

    real primary_roche_lobe, secondary_roche_lobe;
    if (bin_type==Merged || bin_type==Disrupted) {
	primary_roche_lobe   = 2 * stp->get_effective_radius();
	secondary_roche_lobe = 2 * sts->get_effective_radius();
    }
    else {
	primary_roche_lobe   = roche_radius(semi, m1, m2);
	secondary_roche_lobe = roche_radius(semi, m2, m1);

    }
    
    //    if(identity>0)
    //	s << identity << " ";
    //    else
    //	s << get_node()->format_label() << " ";

    s << identity << " "
      << bin_type << " "
      << current_mass_transfer_type << " "
      << binary_age << " "
      << semi << " "
      << eccentricity << "     ";



	s << stp->get_identity() << " "	  
	  << stp->get_element_type() << " "
	  << m1 << " "
	  << stp->get_effective_radius() << " "
	  << log10(stp->temperature()) << " "
	  << stp->get_core_mass() << " "
//      << stp->get_effective_radius()/primary_roche_lobe  << "    "
//      << primary_roche_lobe << "  "
//      << stp->get_radius()<< "      "

	  << sts->get_identity() << " "	  
	  << sts->get_element_type() << " "
	  << m2 << " "
	  << sts->get_effective_radius() << " "
	  << log10(sts->temperature()) << " "
	  << sts->get_core_mass()
//      << " " << sts->get_effective_radius()/secondary_roche_lobe
//      << secondary_roche_lobe << "  "
//      << sts->get_radius()
	  << endl;

   
      //get_initial_primary()->dump(s, brief);
      //      s << "    ";
      //    get_initial_secondary()->dump(s, brief);
      //    s << endl;
    
  }
  else {

      if(bin_type!=Merged) {

      int n_elements = no_of_elements();

      s << n_elements
        << "\n " << identity
	<< " " << bin_type 
        << " " << binary_age
        << " " << bin_type
        << " " << eccentricity
        << " " << get_period()
        << " " << semi
        << " " << velocity
        << endl;
      initial.dump(s);

      s << "  ";
      get_primary()->dump(s, brief);
      s << "  ";
      get_secondary()->dump(s, brief);
      }
      else get_primary()->dump(s, brief);
  }
}

void double_star::dump(char * filename, bool brief) {


  ofstream s(filename, ios::app|ios::out);
  if (!s) cerr << "error: couldn't create file " << filename <<endl;

  dump(s, brief);
  s.close();
}

void double_star::adjust_binary_after_wind_loss(star * donor,
                                                const real md_dot,
                                                const real dt) {
// Binary separation is adjusterd upon mass loss by
// stellar wind of one of the companions. 
// Fraction of the mass lost by the donor is dumped on its
// companion.
// Actually this routine is only applicable to
// non-dynamically evolving systems.
// In a N-body system the dynamics should take care of the
// stellar wind mass loss.


//		check if object is still bound.
     if (bin_type!=Merged && bin_type!=Disrupted) {
        star* accretor = donor->get_companion();

        real M_old = get_total_mass();
        real old_donor_mass = donor->get_total_mass();
        real old_accretor_mass = accretor->get_total_mass();

        real ma_dot = accretor->accrete_from_stellar_wind(md_dot, dt);
        real M_new = M_old - md_dot + ma_dot;
        real new_donor_mass = old_donor_mass - md_dot;
        real new_accretor_mass = old_accretor_mass + ma_dot; 
        
        if (md_dot>0 && ma_dot>=0) {
	  // real alpha = 1 - ma_dot/md_dot;
          // semi *= pow(pow(new_donor_mass/old_donor_mass, 1-alpha)
          //      *  new_accretor_mass/old_accretor_mass, -2)
          //      *  M_old/M_new;
	  // Changed alpha to eta (SPZ+GN:09/1998)
	   real eta =  ma_dot/md_dot;
           semi *= pow(pow(new_donor_mass/old_donor_mass, eta)
                *  new_accretor_mass/old_accretor_mass, -2)
                *  M_old/M_new; 
        }
        donor = donor->reduce_mass(md_dot);
     }
     else {
       donor = donor->reduce_mass(md_dot);
       donor->get_companion()->set_spec_type(Accreting, false);
     }

}

real double_star::angular_momentum() {
//	Determine the angular momentum of the binary.

// 	Shore, S.N., Livio, M., Heuvel, EPJ, 1992, in Interacting Binaries
// 	(edt. N Nussbaumer and A. Orr) Spirnger verlag p.\ 389.
        real omega = TWO_PI/(get_period()*cnsts.physics(days));
        real mp = get_primary()->get_total_mass()*cnsts.parameters(solar_mass);
        real ms = get_secondary()->get_total_mass()*cnsts.parameters(solar_mass);
        real a =  semi*cnsts.parameters(solar_radius);
        // (SilT: 6 Jan 2010) square root over (1-eccentricy*eccentricy) was missing
        return omega*mp*ms*a*a*sqrt(1-eccentricity*eccentricity)/(mp+ms);         
     }

void double_star::circularize() {
//	Tidal circularization with 
//	conservation of angular momentum.
//	Binary is not nesecerely circularized but 
//	periostron adapted to current TIDAL_RANGE.
//	This method is more realistic than instantanious 
//	curcularization.

     real pericenter = semi*(1-eccentricity);

     if((bin_type != Merged && bin_type != Disrupted)          &&
        (pericenter<= cnsts.parameters(tidal_circularization_radius)
	            * get_primary()->get_effective_radius()         ||
         pericenter<= cnsts.parameters(tidal_circularization_radius)
	            * get_secondary()->get_effective_radius())      &&
        (eccentricity>=0 && eccentricity<1.)) /*safety*/         {
            real peri_new = cnsts.parameters(tidal_circularization_radius)
	                  * Starlab::max(get_primary()->get_effective_radius(),
				    get_secondary()->get_effective_radius());
            real circ_semi = semi*(1-eccentricity*eccentricity);
            real new_ecc = circ_semi/peri_new - 1;
//cerr<<"old: "<<" "<<semi<<" "<<eccentricity<<" "<<pericenter<<endl;
//cerr<<"new: "<<circ_semi<<" "<<new_ecc<<" "<<peri_new<<endl;
            if(new_ecc<=0) { 
              semi=circ_semi;
              eccentricity = 0;
            }
            else { 
              semi = circ_semi/(1 - new_ecc*new_ecc);
              eccentricity = new_ecc;
            }
//cerr<<"final: "<<semi<<" "<<eccentricity<<endl;
            //semi *= (1-eccentricity*eccentricity);
            //eccentricity = 0;
	    
// (SilT 26 October 2010) Rotation_period should not be set every timestep to synchronisation                    	    
//	    real rotation_period = cnsts.physics(days)*get_period();
//            get_primary()->set_rotation_period(rotation_period);
//            get_secondary()->set_rotation_period(rotation_period);
        }
//force_circularization();	//	For model E.
}

//	Makes binaries immediately circular.
void double_star::force_circularization() {


     if((bin_type != Merged && bin_type != Disrupted) &&
        (eccentricity>0 && eccentricity<1.)) /*safety*/          {
            semi *= (1-eccentricity*eccentricity);
            eccentricity = 0;
	    real rotation_period = cnsts.physics(days)*get_period();
            get_primary()->set_rotation_period(rotation_period);
            get_secondary()->set_rotation_period(rotation_period);
        }
}

//	Minimum timestep determination is needed for
//	preventing infinit integration upon dynamical timescale
//	mass-transfer.
void double_star::determine_minimal_timestep() {

     minimal_timestep = cnsts.safety(timestep_factor)*donor_timescale;

     if(minimal_timestep<cnsts.safety(minimum_timestep)) {
//        cerr << " minimal_timestep < safety minimum ";
//	cerr << " in double_star::determine_minimal_timestep: ";
//        cerr << minimal_timestep << endl;
        minimal_timestep = cnsts.safety(minimum_timestep);
     }
}

// (SPZ+GN: 28 Jul 2000)
// Limits timestep with safety which may be different for 
// dynamics or pop. synth.
real double_star::internal_time_step(const real evolve_timestep) {

  real dt;
  if (get_use_hdyn()) {

    dt = cnsts.star_to_dyn(binary_update_time_fraction) * evolve_timestep;
  } 
  else {

    dt = cnsts.safety(maximum_binary_update_time_fraction) 
       * evolve_timestep;
  }

  return dt;
}

real double_star::determine_dt(const real ageint, const real time_done) {
//	Determined the evolution timestep.
//	This is the hard of the program.
//	DO NOT CHANGE THIS ROUTINE!

     real dtime = ageint - time_done;

     // (SPZ+GN: 28 Jul 2000)
     // Introduced local internal_time_step(...)
     real dtime_ev = internal_time_step(get_primary()->get_evolve_timestep());

     if (bin_type != Merged) {
       dtime_ev = Starlab::min(dtime_ev, internal_time_step(
		                get_secondary()->get_evolve_timestep()));

       dtime = Starlab::min(dtime, dtime_ev);

       if(bin_type != Disrupted) {

	 real dt_orbit = orbital_timescale();

	 dtime = Starlab::min(dtime, dt_orbit);
	 
//	 if(dtime>minimal_timestep && bin_type==Semi_Detached) 
//	   dtime = minimal_timestep;
       }
     }
     else 
       dtime = Starlab::min(dtime, dtime_ev);

     if(dtime>cnsts.safety(maximum_timestep))
       dtime = cnsts.safety(maximum_timestep);
     if(dtime<cnsts.safety(minimum_timestep)) 
       dtime = cnsts.safety(minimum_timestep);

     return dtime;
}

// First time a dstar is initialized nucleair_evolution_timescale
// is set. (SPZ+GN:24 Sep 1998)
void double_star::set_donor_timescale(star* donor,
				      bool first) { // default = false
  

  if (first) {
//    donor_timescale = donor->nucleair_evolution_timescale();
//    current_mass_transfer_type = Nuclear;
    // (SPZ+GN: 28 Jul 2000)
    // First step may be shaky and it is not clear what is best
    // Therefore, to be conservative, we use
    // the thermal timescale of the ACCRETOR!
    // So, Mdot is computed for conservative mass trasnfer, but
    // on the shortest possible timestep.
    //donor_timescale = donor->get_companion()->kelvin_helmholds_timescale();
    current_mass_transfer_type = Unknown;
  }
  else{
    donor_timescale =
      donor->mass_transfer_timescale(current_mass_transfer_type); 
  }
  donor_identity = donor->get_identity();
  donor_type = donor->get_element_type();
}

#if 0
// Not used any more due to replacement of safeties.
void double_star::binary_in_contact(const real dt) {
//	Search which star is actually filling its Roche lobe.

     real rl_p = roche_radius(get_primary());
     real rl_s = roche_radius(get_secondary());

     if (get_primary()->get_effective_radius() >= rl_p) {
        if (get_secondary()->get_effective_radius() < rl_s) {
           if (!ready_for_mass_transfer(get_primary())) return;
           semi_detached(get_primary(), get_secondary(), dt);
        }
        else
	  contact_binary(dt);   // used to be called: common_envelope(dt);
     }
     else if (get_secondary()->get_effective_radius() >= rl_s) {
        if (!ready_for_mass_transfer(get_secondary())) return;
        semi_detached(get_secondary(), get_primary(), dt);
     }

}

bool double_star::ready_for_mass_transfer(star* donor) {
//	Make sure that the star that fills its Roche lobe is the same as
//	the one that provided its mass-transfer timescale.

     real old_dt_min = minimal_timestep;
     bool mass_transfer_allowed = TRUE;

     mass_transfer_type type = Unknown;
     
     if (donor->get_identity()!=donor_identity) {
        donor_identity    = donor->get_identity();
        donor_type        = donor->get_element_type();
        donor_timescale   = donor->mass_transfer_timescale(type);

        determine_minimal_timestep();

        if(old_dt_min<=minimal_timestep) 
           mass_transfer_allowed = TRUE;
        else
           mass_transfer_allowed = FALSE;
     }

     return mass_transfer_allowed;
}

#endif
       

void double_star::semi_detached(star* donor,
                                star* accretor,
                                real dt) {
//	Binary is semi-detached.
//	Donor should lose mass while the accretor should accrete.

  if (!stable(donor)) {
    cerr << "semi_deatched not stable => ::common_envelope" << endl; 
    
    //    dynamic_mass_transfer();
    // (GN+SilT Mar  2 2011)
    current_mass_transfer_type = Darwin;
    tidal_instability();
  }
  else if (get_current_mass_transfer_type()==Dynamic) {
    cerr << "dynamic mass transfer => ::dynamic_mass_transfer" << endl; 
    
    dynamic_mass_transfer();
// in common_envelope dynamic_mass_trasfer is called
//    dynamic_mass_transfer(donor,accretor);
  }
  else {

    if (REPORT_TRANFER_STABILITY) {
      cerr << "stable mass transfer => ::perform_mass_transfer" << endl; 
    }

   // dump(cerr, false);
    perform_mass_transfer(dt, donor, accretor);
  }

  //get_seba_counters()->semi_detached++;

}

// Is not used but provides interesting options.
void double_star::contact(star* donor,
                          star* accretor,
                          real dt) {
//	Both stars in contact.
//	W Uma (stable) system or spiral-in may occur.

//      Mass loss according to AM accretor.
     real M_old = get_total_mass();
     real old_donor_mass = donor->get_total_mass();
     real old_accretor_mass = accretor->get_total_mass();
     //real q_old = old_accretor_mass/old_donor_mass;

//	Old mechanism
//     real md_dot = 0;
//     md_dot = donor->subtrac_mass_from_donor(dt, md_dot);

//	New mechanism
     real md_dot=0;
     donor = donor->subtrac_mass_from_donor(dt, md_dot);
/*
 *    try { 
 *       donor->subtrac_mass_from_donor(dt, md_dot);
 *      }
 *    catch (star * element) {
 *       donor = element;
 *    }
 */
     if (md_dot>0) {
        real ma_dot = accretor->add_mass_to_accretor(md_dot, donor->hydrogen_envelope_star(), dt);

        real M_new = get_total_mass();
        real new_donor_mass = donor->get_total_mass();
        real new_accretor_mass = accretor->get_total_mass();
        //real q_new = new_accretor_mass/new_donor_mass;

        real a_fr;
        if (!accretor->remnant()) {
	  //real beta = 3;
           a_fr  = pow(old_donor_mass*old_accretor_mass
                 / (new_donor_mass*new_accretor_mass), 2);
           semi *= pow(M_new/M_old,
		       2*cnsts.parameters(specific_angular_momentum_loss) + 1)*a_fr;
        }
        else {
           real alpha = 1 - ma_dot/md_dot;
           if (alpha<1) {
              a_fr  = (new_donor_mass/old_donor_mass)
                    * pow(new_accretor_mass/old_accretor_mass, 1/(1-alpha));
              semi *= (M_old/M_new)/pow(a_fr, 2);
           }
           else {
              a_fr  = exp(2*(new_donor_mass-old_donor_mass)/new_accretor_mass);
              semi *= (M_old/M_new)*a_fr/pow(new_donor_mass/old_donor_mass, 2);
           }
        }

      // redundant due to implementation of new force donor timescale
      //                                         (SPZ+GN:23 Sep 1998)
#if 0	
       real t_donor = donor->mass_transfer_timescale();
       if (((q_old<=1 && q_new>=1) ||
          (t_donor<0.1*donor_timescale || t_donor>10*donor_timescale)) &&
          (bin_type!=Disrupted && bin_type!=Merged)) {
          set_donor_timescale(t_donor,
                              donor->get_element_type(),
                              donor->get_identity());
          determine_minimal_timestep();
       }
#endif
       
     }

     //get_seba_counters()->contact++;

     //     if (bin_type!=Merged && bin_type!=Disrupted)
     //        update_binary_parameters();
}

bool  double_star::stable(star* st) {	// default = NULL
//	Stability test for mass transfer.
//	Single stellar angular momentum should be smaller than
//	1/3-th the binary angular momentum.

    // do not allow mass transfer between planets
//    if(st->get_total_mass() < 2*cnsts.safety(minimum_mass_step)) {
//	return false;
//    }

  real J_star;
  if(st) {
    J_star = st->angular_momentum();
  }
  else {
    real J_prim = get_primary()->angular_momentum();
    real J_sec  = get_secondary()->angular_momentum();

    J_star = Starlab::max(J_prim, J_sec);

  }
  
  real J_bin  = angular_momentum();
  
  if(REPORT_TRANFER_STABILITY) {
    PRC(J_star);PRL(J_bin);
  }

  if (J_star >= J_bin *
      cnsts.parameters(Darwin_Riemann_instability_factor)) {

    cerr << "Spiral-in: Darwin Riemann instability"<<endl;
    
    return FALSE;
  }
 
  return TRUE;
}

#if 0
// is not used any more (SPZ+GN:28 Sep 1998)
void double_star::update_binary_parameters() {
//	Update the evolved binary.

     real primary_mass = get_primary()->get_total_mass();
     real secondary_mass = get_secondary()->get_total_mass();

     get_primary()->set_previous_radius(get_primary()->get_radius());
     get_secondary()->set_previous_radius(get_secondary()->get_radius());
}
#endif

void double_star::perform_mass_transfer(const real dt,
					star * donor,
					star * accretor) {
//	Adjust binary parameter according to the model
//	Take account of mass lost and transfered from the donor to the
//	accretor.

// (GN+SilT Mar  2 2011) new dumping regime
    if (!first_contact) {
    
        dump("SeBa.data", true);
        first_contact = true;
    }
    
    
    if (REPORT_TRANFER_STABILITY) {
        cerr<<"enter perform_mass_transfer("<<dt<<", "
            <<donor->get_element_type()<<", "<<accretor->get_element_type()
        <<") semi= "<< semi<<endl;
    }
    
//      Mass loss according to AM accretor.
    real M_old = get_total_mass();
    real old_donor_mass = donor->get_total_mass();
    real old_accretor_mass = accretor->get_total_mass();
    
    // (SilT 13 Apr 2012) Only calculate md_dot, subtract the mass later 
    // so that in the SeBa output the correct separation is given
    real md_dot = donor->mdot_limit(dt);
    md_dot = min(md_dot, donor->get_envelope_mass());

    if (md_dot>0) {

        real ma_dot = accretor->add_mass_to_accretor(md_dot, donor->hydrogen_envelope_star(), dt);
    
        real M_new = M_old - md_dot + ma_dot;
        real new_donor_mass = old_donor_mass - md_dot;
        real new_accretor_mass = old_accretor_mass + ma_dot;
        //PRC(M_new);PRC(new_donor_mass);PRC(new_accretor_mass);
    
        real a_fr;
        if (!accretor->remnant()) {
            // General case: semi-conservative mass transfer.
            a_fr  = pow(old_donor_mass*old_accretor_mass
                    / (new_donor_mass*new_accretor_mass), 2);
            semi *= pow(M_new/M_old,
                2*cnsts.parameters(specific_angular_momentum_loss)
                + 1)*a_fr;	
        }
        else {
            // Special case: mass transfer to compact object as accretor.
            //               Two possibilities:
            //               1) eta>0: mass lost as wind from accretor.
            //               2) eta==0: exponential spiral-in.
            real eta = ma_dot/md_dot; 
            if (eta>0) {
                a_fr  = (new_donor_mass/old_donor_mass)
                    * pow(new_accretor_mass/old_accretor_mass, 1/eta);
                semi *= (M_old/M_new)/pow(a_fr, 2); 
            }
            else {
                a_fr  = exp(2*(new_donor_mass-old_donor_mass)
            	 /new_accretor_mass); 
                semi *= (M_old/M_new)*a_fr
                    / pow(new_donor_mass/old_donor_mass, 2);
            }
        
            //   Let outer component accrete from mass lost of inner binary.
            //	Routine block for triple evolution.
            //	Tricky but fun!
            //   if (is_binary_component()) 
            if (is_star_in_binary()) {
                real mdot_triple = M_old - M_new;
                if (mdot_triple>0) {
                    get_binary()->adjust_triple_after_wind_loss(this,
            					   mdot_triple, dt);
                }
                else if (is_binary_component()) {
                    if(mdot_triple<0) {
                        cerr << "enter perform_mass_transfer(" << dt << ", "
                        << donor->get_element_type()<<", "
                        <<accretor->get_element_type()
                        <<")"<<endl;
                        cerr<<"Mass lost during non-conservative mass transfer is negative."
                        <<"mlost ="<<mdot_triple<<endl;
                    }
                }
            else {
                cerr<<"Presumed binary component is not a binary root, "
                << "nor a hierarchical root"<<endl;
            }
        }
    }
    // md_dot2 is equal to md_dot
    real md_dot2 = 0 ;
    donor = donor->subtrac_mass_from_donor(dt, md_dot2);
//    PRC(md_dot);PRC(md_dot2);
}

#if 0
  switch(current_mass_transfer_type) {
      case Nuclear:           get_seba_counters()->nuclear++;
           break;     
      case AML_driven:        get_seba_counters()->aml_driven++;
           break;     
      case Thermal:           get_seba_counters()->thermal++;
           break;     
      case Dynamic:           get_seba_counters()->dynamic++;
           break;     
  }
#endif
}

// final part of ::perform_mass_transfer()... 
#if 0
	// Since donor timescale is set within recursive with
	// the appropriate donor this part is redundant.
	// (SPZ+GN:23 Sep 1998)
//	Adjust mass-transfer timescale if needed.
        real t_donor = donor->mass_transfer_timescale();
        if (((q_old<=1 && q_new>=1) ||
	    (t_donor<0.1*donor_timescale || t_donor>10*donor_timescale)) && 
           (bin_type!=Disrupted && bin_type!=Merged)) {
              set_donor_timescale(t_donor, 
                           donor->get_element_type(), 
                           donor->get_identity());
              determine_minimal_timestep();
        }
#endif
	

void double_star::post_sn_companion_velocity(star* companion,
					     const real mdot) {
//		Result of impact on companion velocity!
//		Temorary removed!
        return ;

//		Statement never reached.
//		Removed 16_01_94!
//cerr<<"void double_star::post_sn_companion_velocity()"<<endl;
//cerr<<"pre sn velocity: "<< companion->get_velocity();
     real new_total_binary_mass = get_total_mass()-mdot;
     real vel = abs((1 - 2*new_total_binary_mass/get_total_mass())
                           *pow(get_total_mass()
                           /companion->get_total_mass(), 2));
     //if(vel<=-1) vel = -1;
     vel = companion->get_velocity()*sqrt(1 + vel);
     companion->set_velocity(vel);
//cerr<<"post sn velocity: "<< companion->get_velocity()<<endl;
}

void double_star::instantaneous_element() {

  try_zero_timestep();
}


// evolve_element 
// evolves a binary to the new time which is close to end_time 
// NOTE: the new time of the binary is generally != end_time
//       but slightly smaller.
//       The BINARY_UPDATE_TIMESTEP must therefore be small.
void double_star::evolve_element(const real end_time) 
{
  //  cerr<<"Evolve the binary"<<flush <<endl;
  //    real current_time = binary_age
  //			  + BINARY_UPDATE_TIMESTEP*get_evolve_timestep();
  //  real current_time = binary_age
  //                    + get_evolve_timestep();
    
    // Add try_zero_timestep for the case that current_time==end_time.
    // Might be completely wrong, anyway, numerical precision prevents
    // such a constructiona anyway.
    // Implemented by (SPZ:2/1998)
    //if (current_time<=end_time)
 
    real current_time = binary_age;

    //PRC(binary_age);PRL(end_time);
    if (binary_age<end_time)
	do {
	  if (REPORT_BINARY_EVOLUTION) {
	    cerr << "converge binary to system age at t = "
		 << current_time << "  binary age = " << binary_age
		 << "  time step is " << get_evolve_timestep() << endl;
	  }

//	  current_time = min(end_time, 
//			     current_time+BINARY_UPDATE_TIMESTEP *
//			     get_evolve_timestep();
//	  current_time += BINARY_UPDATE_TIMESTEP * get_evolve_timestep();

	  // Attempt a small timestep for consistency reasons.
	  // but this timestep may be a bit too small.
          // current_time += cnsts.safety(binary_update_time_fraction) 
          //               * get_evolve_timestep();

	  // A bigger timestep may be allowed, but the choise of the 
	  // timestep must be consistemt, 
	  // in steps of binary_update_time_fraction
	  current_time += internal_time_step(get_evolve_timestep());
	    
	  // safety: Never evolve passed end_time.
	  // removed SPZ, 7 Febr 2003
	  //if (current_time>=end_time)
	  //  break;
	  // The previous if-statment has been removed as it prevents 
	  // a binary to be regularly updated if small time steps are 
	  // required. This may not have been a problem in the passed.
	  // Low mass binaries ware likely not updated, even though they
	  // may have undergone severe dynamical evolution, which again 
	  // could have affected the binary evolution.
	  // for high mass binaries (R136) this is no problem, as the 
	  // evolution timescale for these massive binaries are short
	  // to begin with.
	  // However, we may have ignored some aspects of binary evolution 
	  // in the passed.

	  if (current_time>=end_time){
	    evolve_the_binary(end_time);
	  }
	  else {
	    evolve_the_binary(current_time);
	  }

	} while (binary_age<end_time);
	
    else if (end_time == 0){
      try_zero_timestep();
    }
    else {// (SPZ: 21 Mar 2001): added the circularizaton in case....
     circularize();
    }
}


void double_star::try_zero_timestep() {

  if (bin_type!=Merged && bin_type!=Disrupted) {

     star *p = get_primary();
     star *s = get_secondary();
     if (p) 
	 p->evolve_element(binary_age);

     if(s) {
	 s->evolve_element(binary_age);
     }

     circularize();

    real rl_p = roche_radius(get_primary());
    real rl_s = roche_radius(get_secondary());

    real rp = get_primary()->get_effective_radius();
    real rs = get_secondary()->get_effective_radius();

    // Just for memory and output routines.
    // Primary fills Roche-lobe.
       if (rp >= rl_p) {
	 //get_primary()->set_effective_radius(rl_p);
          get_primary()->set_spec_type(Rl_filling);
	  get_primary()->set_spec_type(Accreting, false);
          get_secondary()->set_spec_type(Accreting);
       }
       else
	 get_primary()->set_spec_type(Rl_filling, false);

       if (rs >= rl_s) {

	 //get_secondary()->set_effective_radius(rl_s);
	 get_secondary()->set_spec_type(Rl_filling);
	 get_secondary()->set_spec_type(Accreting, false);
	 get_primary()->set_spec_type(Accreting);

	  // One could check here for contact binary to cause
	  // the secondary to be the donor.
       }
       else
	 get_secondary()->set_spec_type(Rl_filling, false);
  }

}

void double_star::evolve_the_binary(const real end_time) {
//	Evolve both stars synchonously.

     real dt, old_binary_age;
     real time_done = 0;
     real ageint = end_time - binary_age;

     do {
       
        determine_minimal_timestep();
        dt = determine_dt(ageint, time_done);

        recursive_counter = 0;             
        old_binary_age = binary_age;

        //redundant line 	
        get_primary()->set_previous_radius(get_primary()
		     ->get_effective_radius());
        get_secondary()->set_previous_radius(get_secondary()
		       ->get_effective_radius());

    	// Stellar wind mass loss is taken care of by the stars themselves.
    	// line removed at 20/08/97 by SPZ
    	// perform_wind_loss(dt);

        refresh_memory();

       	recursive_binary_evolution(dt, binary_age+dt);
        calculate_velocities();
        time_done += binary_age-old_binary_age;
     }
     while(time_done<ageint);
}

void double_star::recursive_binary_evolution(real dt,
                                const real end_time) {

//	Evolve two stars synchronously.
//	If one star start filling its Roche-lobe
//	return one time-step and search for the moment
//	of Roche lobe contact.
//	Then start to transfer mass.

    recursive_counter++;

    if (REPORT_RECURSIVE_EVOLUTION) {
      cerr << "recursive_binary_evolution(dt " << dt 
	   << ", end_time " << end_time << "): " <<bin_type<< endl;
      cerr<<dt<<" "<<binary_age<<" "<<end_time<< " "<<recursive_counter<<endl;
      cerr << "recursive # " << recursive_counter << endl;
      dump(cerr, false);
    }
      
    if(recursive_counter >= cnsts.safety(maximum_recursive_calls)) {
	cerr << "WARNING: ";
	cerr << "recursive_binary_evolution(dt " << dt 
	    << ", end_time " << end_time << "): " <<bin_type<< endl;
	cerr<<dt<<" "<<binary_age<<" "<<end_time<<endl;
	cerr << "recursive_counter == " << recursive_counter << endl;
	dump(cerr);

	return;
    }

    if (binary_age>=end_time) return;
    
    dt = Starlab::min(dt, determine_dt(end_time, binary_age));
    binary_age += dt;

    // Angular momentum loss must be computed here: binary can merge!
    // However must be before evolve_element of stars, since the first step
    // of the recursive binev reaches to far, causing the a giant to become
    // huge: -> merger due to magnetic braking!
    angular_momentum_loss(dt);
       
    star *p = get_primary();
    star *s = get_secondary();
    if (p) p->evolve_element(binary_age);
    if (bin_type!=Merged && s)  s->evolve_element(binary_age);
    p = s = NULL;

    if (bin_type!=Merged && bin_type!=Disrupted) {
      circularize();

      real rl_p = roche_radius(get_primary());
      real rl_s = roche_radius(get_secondary());
      
      real rp = get_primary()->get_radius();
      real rs = get_secondary()->get_radius();
      
      star* donor    = get_primary();
      star* accretor = get_secondary();
      
      // Primary fills Roche-lobe?
      if (rp >= rl_p) {
	get_primary()->set_spec_type(Rl_filling);
	get_primary()->set_spec_type(Accreting, false);
	get_secondary()->set_spec_type(Accreting);
      }
      else
	get_primary()->set_spec_type(Rl_filling, false);

      if (rs >= rl_s) {

	// secondary is donor
	donor  = get_secondary();
	accretor = get_primary();

	get_secondary()->set_spec_type(Rl_filling);
	get_secondary()->set_spec_type(Accreting, false);
	get_primary()->set_spec_type(Accreting);
	
	// One could check here for contact binary to cause
	// the secondary to be the donor.
      }
      else
	get_secondary()->set_spec_type(Rl_filling, false);

	 
//  Determines if we really have to do with a mass
//  transferring binary.
      if (rp >= rl_p    ||
	  rs >= rl_s) {
	if (REPORT_RECURSIVE_EVOLUTION) {
	  cerr << " (donor, accretor) = (" << donor->get_identity()
	       << ", " << accretor->get_identity()
	       << ")";

	}
	
	// set the mass transfer timescale to the
	// timescale of the newly determined donor
	// implemented (SPZ+GN:23 Sep 1998)
	mass_transfer_type prev_mt = current_mass_transfer_type;

	set_donor_timescale(donor);
	if (prev_mt > 0 && prev_mt < 4 && current_mass_transfer_type != prev_mt) first_contact = false;

	determine_minimal_timestep();

	// Check for Contact
	if (rp >= rl_p &&  rs >= rl_s){       

	  if (!first_contact) {
	  // (SilT Nov 25 2012)
   	  // Goes wrong when first stable mass transfer. Need to set in ::contact_binary...
      // bin_type = Contact;
	  // first_contact=true;

	    if (REPORT_RECURSIVE_EVOLUTION) 
	      cerr << "\tFirst contact" << endl;

	    cerr << "First Roche-lobe contact for: ";
	    put_state();
	    cerr << endl;
	  // (SilT Nov 25 2012)
	  // Goes wrong when first stable mass transfer. Need to set in ::contact_binary...	    
	  //dump("SeBa.data", true);
	  }

	  if (REPORT_RECURSIVE_EVOLUTION) {
	    cerr << "\tContact" << endl;
	    put_state();
	    cerr << endl;
	    dump(cerr, false);
	  }
	  
	  contact_binary(dt);   

	  star *p = get_primary();
	  star *s = get_secondary();
	  if (p) p->evolve_element(binary_age);
	  if (bin_type!=Merged && s)  s->evolve_element(binary_age);
	  p = s = NULL;
	  refresh_memory();

	  if (binary_age>=end_time)
	    return;
	  else {
	    real max_dt = end_time - binary_age;

	    recursive_binary_evolution(max_dt, end_time);
	    return;
	  }
	} // end Contact
	
	else if(dt <= minimal_timestep) { // Ready to go!

	  if (REPORT_RECURSIVE_EVOLUTION) {
	    cerr << "\tMass transfer" << endl;
	    put_state();
	    cerr << endl;
	    dump(cerr, false);
	  }
	  
	  bin_type = Semi_Detached;
	  
	  if (!first_contact) {
	    if (REPORT_RECURSIVE_EVOLUTION)
	      cerr << "\tFirst contact" << endl;
	    
	    donor->first_roche_lobe_contact_story(accretor->get_element_type());

	    cerr << "First Roche-lobe contact for: ";
        put_state();
	    cerr << endl;
	    
	    // (GN+SilT Mar  2 2011)
	    // do dump later, so that we know exactly what type of MT
	    // (i.e. Spiral_In, Common_Envelope etc...
	    //dump("SeBa.data", true);
	  }
	    
	  semi_detached(donor, accretor, dt);
	  
	  // (GN+SilT Mar  2 2011)
	  //first_contact=true; 
	  // (GN+SilT Apr 19 2011)
	  // Goes wrong with Common-envelope etc. Need to set in ::perform_mass_transfer...

	  star *p = get_primary();
	  star *s = get_secondary();
	  if (p) p->evolve_element(binary_age);
	  if (bin_type!=Merged && s)  s->evolve_element(binary_age);
	  p = s = NULL;	      
	  refresh_memory();

	  if (binary_age>=end_time) return;

	  dt = minimal_timestep;

	  if (binary_age + dt > end_time ||
	      (bin_type==Merged || bin_type == Disrupted)) 
	    dt = end_time - binary_age;
	  
	  recursive_binary_evolution(dt, end_time);
	  return;
	} // end semi_detached and ready to go
	else {// semi detached and dt >= minimal_timestep --> Total recall
	
	  if (REPORT_RECURSIVE_EVOLUTION) {
	    cerr << "\tTotal recall" << endl;
	    put_state();
	    cerr << endl;
	    dump(cerr, false);
	  }
	  
	  bin_type = Semi_Detached;
	  recall_memory();
	  dt *= 0.5;
	  recursive_binary_evolution(dt, end_time);
	  return;
	}
      } // end Total recall
      else if(bin_type == Semi_Detached || bin_type == Contact) { // restep

	bin_type = Detached;
	get_primary()->set_spec_type(Rl_filling, false);
	get_secondary()->set_spec_type(Rl_filling, false);
	
	if (REPORT_RECURSIVE_EVOLUTION) {
	  cerr << "\tRestep" << endl;
	  put_state();
	  cerr << endl;
	  dump(cerr, false);
	}

	// (GN+SilT Mar  2 2011)
	// NEEDS TO BE CHECKED WITH SPZ
	//if (dt>minimal_timestep) {

	current_mass_transfer_type = Unknown;
	refresh_memory();
	  //}

	recursive_binary_evolution(dt, end_time);

	return;
      }
      else { // Detached evolution

	if (REPORT_RECURSIVE_EVOLUTION) {
	  cerr << "\tDetached" << endl;
	  put_state();
	  cerr << endl;
	  dump(cerr, false);
	}

	if(bin_type != Detached) {
	  cerr<< "ERROR in ::recursive bin_type != Detached "<<endl;
	  cerr<< "where should be Detached"<<endl;
	  dump(cerr, false);
	  bin_type = Detached;
	}
	
	// (GN+SPZ May  3 1999) effective_radius 
	get_primary()->
	  set_effective_radius(get_primary()->get_radius());
	get_secondary()->
	  set_effective_radius(get_secondary()->get_radius());
	
	current_mass_transfer_type = Unknown;

	refresh_memory();

	real max_dt = end_time - binary_age;
	if (max_dt < end_time - binary_age) {
	  cerr << "Time step reduction in double_star::recursive..."
	       << endl;
	  cerr << "     max_dt (" << max_dt
	       << ") < end_time - binary_age ("
	       << end_time - binary_age << ")." << endl;
	}
	
	//if (max_dt == 0 && dt > 2 * minimum_timestep) first_contact = false;

	recursive_binary_evolution(max_dt, end_time);

	return;
      } // end Detached
    } // end evolution surviving binaries 
    else { // Mergers and disrupted systems

      if (REPORT_RECURSIVE_EVOLUTION) {
	cerr << "\tTerminated" << endl;
	put_state();
	cerr << endl;
	dump(cerr, false);
      }
      get_primary()->set_spec_type(Rl_filling, false);
      get_secondary()->set_spec_type(Rl_filling, false);
      // (SPZ+GN: 28 Jul 2000) was a 'Gedoogde bug' but now ok.
      // merged object should not accrete.
      get_primary()->set_spec_type(Accreting, false);
      get_secondary()->set_spec_type(Accreting, false);
  
      current_mass_transfer_type = Unknown;

      refresh_memory();
      real max_dt = end_time - binary_age;
      recursive_binary_evolution(max_dt, end_time);
    } 
}


// redundant due to implementation of new force donor timescale
//                                         (SPZ+GN:23 Sep 1998)
#if 0
void double_star::set_donor_timescale(star * donor) {

//cerr <<"void double_star::set_donor_timescale(donor="
//     << donor->get_element_type()<<")"<<endl;

  
     if (bin_type!=Merged && bin_type!=Disrupted) {
       int id          = donor->get_identity();
       stellar_type st = donor->get_element_type();
       mass_transfer_type type = Unknown;
       real t          = donor->mass_transfer_timescale(type);
       current_mass_transfer_type = type;
	
	 if (st!=NAS) {
           if (donor_identity!=id) {
              donor_identity = id;
              donor_type = st;
              donor_timescale = t;
           }
           else if (donor_timescale<=0 || donor_type!=st) {
	     donor_identity = id;
	     donor_type = st;
	     donor_timescale = t;
           }
	 }
     }
}


void double_star::set_donor_timescale(const real t, const stellar_type st,
                                      const int id) {

     if (bin_type!=Merged && bin_type!=Disrupted) {
        if (st!=NAS) {
           if (donor_identity==id) {
              donor_type = st;
              donor_timescale = t;
           }
           else if(donor_timescale<=0) {
              donor_identity = id;
              donor_type = st;
              donor_timescale = t;
           }
        }
     }

}
#endif

// Determine in which way common_envelope systems spiral-in
void double_star::tidal_instability() {

  //cerr << "tidal_instability"<<endl;

  if (bin_type!=Merged && bin_type!=Disrupted) {
       
    if (get_primary()->giant_star()) {
      if (get_secondary()->giant_star()) 
	double_spiral_in();
      else if(cnsts.parameters(use_angular_momentum_tidal))
	angular_momentum_envelope_ejection(get_primary(), get_secondary());
      else
	spiral_in(get_primary(), get_secondary());
    }
    else if (get_secondary()->giant_star()) {
      // (SPZ+GN:2002Dec4)
      // based on: Nelemans 2003 (to be written)
      if(cnsts.parameters(use_angular_momentum_tidal)) 
	angular_momentum_envelope_ejection(get_secondary(), get_primary());
      else
	spiral_in(get_secondary(), get_primary());
    }
    else {
      cerr << "Merger double_star::tidal_instability" << endl;
      dump(cerr, false);
	   
      // (GN+PimvO Jan 18 2013)
      // if primary is remnant, secondary should be consumer
      if (get_primary()->remnant()) 
	  merge_elements(get_secondary(), get_primary());
      else 
	  merge_elements(get_primary(), get_secondary());

    }
    //get_seba_counters()->common_envelope++;
  }
}

// Determine in which way common_envelope systems spiral-in
void double_star::dynamic_mass_transfer() {

  //cerr << "dynamic_mass_transfer"<<endl;

  if (bin_type!=Merged && bin_type!=Disrupted) {
       
    if (get_primary()->giant_star()) {
      if (get_secondary()->giant_star()) 
	double_spiral_in();
      // (SPZ+GN:2002Dec4)
      // based on: Nelemans 2003 (to be written)
      // (SPZ+GN:2007June29)
      else if (get_secondary()->remnant()) 
	if(cnsts.parameters(use_common_envelope_gamma_gamma))
	  angular_momentum_envelope_ejection(get_primary(), get_secondary());
	else
	  spiral_in(get_primary(), get_secondary());
      else {
	if(cnsts.parameters(use_common_envelope_alpha_alpha))
	  spiral_in(get_primary(), get_secondary());
	else
	  angular_momentum_envelope_ejection(get_primary(), get_secondary());
      }
    }
    else if (get_secondary()->giant_star()) {
      // (SPZ+GN:2002Dec4)
      // based on: Nelemans 2003 (to be written)
      // (SPZ+GN:2007June29)
      // based on: Nelemans 2003 (to be written)
      if (get_primary()->remnant()) {
	if(cnsts.parameters(use_common_envelope_gamma_gamma))
	  angular_momentum_envelope_ejection(get_secondary(), get_primary());
	else
	  spiral_in(get_secondary(), get_primary());
      }
      else {
	if(cnsts.parameters(use_common_envelope_alpha_alpha))
	  spiral_in(get_secondary(), get_primary());
	else
	  angular_momentum_envelope_ejection(get_secondary(), get_primary());
      }
    }
    else {
      cerr << "Merger double_star::dynamic_mass_transfer" << endl;
      dump(cerr, false);
	   
      // (GN+PimvO Jan 18 2013)
      // if primary is remnant, secondary should be consumer
      if (get_primary()->remnant()) 
	  merge_elements(get_secondary(), get_primary());
      else 
	  merge_elements(get_primary(), get_secondary());

    }
    //get_seba_counters()->common_envelope++;
  }
}

// common-envelope for two giant stars in contact.
// results in spiral in of both cores in both envelopes.
// final separation is computed from ejection of both giant envelopes.
// if one of the cores fills its Roche-lobe a merger occurs.
// The merger product loses some mass.
void double_star::double_spiral_in() {

       // (GN+SilT Mar  2 2011) new dumping regime
       bin_type = Double_Spiral_In;
       dump("SeBa.data", true);
       current_mass_transfer_type = Unknown;


       star *p = get_primary();
       star *s = get_secondary();
       
       real mcore_p = p->get_core_mass();
       real menv_p = p->get_envelope_mass();
       real mtot_p = mcore_p + menv_p;
       real mcore_s = s->get_core_mass();
       real menv_s = s->get_envelope_mass();
       real mtot_s = mcore_s + menv_s;

       real r_p = p->get_radius();
       real r_s = s->get_radius();

        //real r_l_p = roche_radius(p);
        //real r_l_s = roche_radius(s);
	//PRC(r_l_p);PRC(r_p);PRC(r_l_s);PRL(r_s); 	
		
       real alpha_lambda = cnsts.parameters(common_envelope_efficiency)
	                 * cnsts.parameters(envelope_binding_energy);

       real post_ce_orbital_energy = mtot_p*menv_p/(alpha_lambda*r_p)
	                           + mtot_s*menv_s/(alpha_lambda*r_s)
	                           + mtot_p*mtot_s/(2 * semi);

       real post_ce_semi = mcore_p*mcore_s/(2*post_ce_orbital_energy);
       
       real rl_p = roche_radius(post_ce_semi, mcore_p, mcore_s);
       real rl_s = roche_radius(post_ce_semi, mcore_s, mcore_p);

       real prc = p->get_core_radius();
       real src = s->get_core_radius();

       if (prc >= rl_p    ||
	   src >= rl_s) {

#if 0       
         // see double_star::spiral_in()
	  
         real semi_p = p->get_core_radius()
	             / roche_radius(1, p->get_core_mass(),
			               s->get_core_mass());
	 real semi_s = s->get_core_radius()
	             / roche_radius(1, s->get_core_mass(),
			               p->get_core_mass());

	  real post_ce_semi = max(semi_p, semi_s);

	  // the mass lost before the two stars merge is computed from
	  // change in orbital energy assuming that
	  // the final separation << initial separation.
	  
	  real mass_loss_fr = post_ce_semi
	                    / (alpha_lambda * mtot_p * mtot_s)
	                    * (pow(mtot_p, 2)/r_p + pow(mtot_s, 2)/r_s);

	  real mlf_limit = (menv_p + menv_s)/(mtot_p + mtot_s);
#endif
	 
	 
         // see double_star::spiral_in()

	 real semi_p = prc                         // was p->get_core_radius()
		     / roche_radius(1, mcore_p, mcore_s);
	 real semi_s = src                         // was s->get_core_radius()
		     / roche_radius(1, mcore_s, mcore_p);
	 
	 real post_ce_semi = Starlab::max(semi_p, semi_s);

	 // the mass lost before the two stars merge is computed from
	 // change in orbital energy assuming that
	 // the final separation << initial separation.
	 // (GN Oct 27 1998) was wrong, we assumed f = (1-f) !
	 //real mass_loss_fr = post_ce_semi
	 //                  / (alpha_lambda * mtot_p * mtot_s)
	 //                  * (pow(mtot_p, 2)/r_p + pow(mtot_s, 2)/r_s);


	 // assume mass left after merger = k * m_tot
	 // ( (1-k) [ M^2/R + m^2/r] ) / alpha_lambda =
	 //                              k^2 M m /(2 af) -  M m /(2 ai)
	 // is quadratic equation in k with
	 // a = M m / (2 af)
	 // b = [ M^2/R + m^2/r] / alpha_lambda
	 // c = -b - M m / (2 ai)
	 //
         // (GN Oct 27 1998) solve k with abc-equation

         real a = (mtot_p * mtot_s)/(2 * post_ce_semi);
         real b = (pow(mtot_p, 2)/r_p + pow(mtot_s, 2)/r_s)
                /  alpha_lambda;
         real c = -b -(mtot_p * mtot_s)/(2 * semi);

         real mass_not_lost_fraction = (-b + sqrt(pow(b, 2) - 4 * a * c))
                                     / (2 * a);

         real mass_loss_fr = 1 - mass_not_lost_fraction;
         real mlf_limit = (menv_p + menv_s)/(mtot_p + mtot_s);

	 if (mass_loss_fr >= 1) {
	      cerr << "ERROR: in double_star::double_spiral_in()" << endl;
	      cerr << "         Fraction of mass lost (" << mass_loss_fr 
	           << ") greater then 1" << endl;
	    
	    dump(cerr, false);
	    mass_loss_fr = 0.99;
	  }
	  else if (mass_loss_fr > mlf_limit) {
	      cerr << "WARNING: in double_star::double_spiral_in()" << endl;
	      cerr << "         Fraction of mass lost (" << mass_loss_fr 
	           << ") greater then limit (" << mlf_limit << ")" << endl;
	    dump(cerr, false);
	  }

	  // Only envelope mass may be removed to prevent
	  // helium star formation.
	  real mlost_p = Starlab::min(mass_loss_fr*mtot_p, 0.99*menv_p);
	  p->reduce_mass(mlost_p);
	  real mlost_s = Starlab::min(mass_loss_fr*mtot_s, 0.99*menv_s);
	  s->reduce_mass(mlost_s);

	  cerr << "Merger double_star::double_spiral_in()"<<endl;
	  dump(cerr, false);

	  merge_elements(p, s);

        } 
        else {

	  cerr << "Survival double_star::double_spiral_in()"<<endl;
	  dump(cerr, false);
	  
	  semi = post_ce_semi;
	  p = p->reduce_mass(menv_p);
	  //I don't know why but sometimes the envelope_mass of the secondary changes, so:
      menv_s = s->get_envelope_mass();
	  s = s->reduce_mass(menv_s);
        }

       //get_seba_counters()->double_spiral_in++;
}


void double_star::spiral_in(star* larger,
                            star* smaller) {
//	Compare total binary energy with the binding energy of the
//	donor's envelope.
//	If two stars are in contact after the common envelope,
//	stars are merged.

      // (GN+SilT Mar  2 2011) new dumping regime
      bin_type = Spiral_In;
      dump("SeBa.data", true);
      current_mass_transfer_type = Unknown;



     if (bin_type!=Merged && bin_type!=Disrupted) {
        real mcore_l = larger->get_core_mass();
        real menv_l = larger->get_envelope_mass();
        real mtot_l = mcore_l + menv_l;

	// What is accreted by the accretor should be reduce in the donor.
	// (SPZ:  2 Jun 1999)
	//real m_acc = smaller->accretion_limit(m_env_l, 
        //             cnsts.parameters(spiral_in_time));
        real m_acc = 0;
        real mtot_s = smaller->get_total_mass() + m_acc;
	menv_l -= m_acc;

        real r_lobe = roche_radius(larger)/semi;

	real alpha_lambda = cnsts.parameters(common_envelope_efficiency)
	                  * cnsts.parameters(envelope_binding_energy);
	real a_spi = semi*(mcore_l/mtot_l)/(1. + (2.*menv_l
                   / (alpha_lambda *r_lobe*mtot_s)));
       

        real rl_l = roche_radius(a_spi, mcore_l, mtot_s);
        real rl_s = roche_radius(a_spi, mtot_s, mcore_l);

       real lrc = larger->get_core_radius();
       real sr = smaller->get_radius();
       real ser = smaller->get_effective_radius();

       // SPZ July 2002: why is there a difference bewteen sr and ser?
	// Merger
        if (lrc>=rl_l || sr>=rl_s) 	{

	  // Determine the minimum separation before either the
	  // core of the larger star or
	  // its companion (the smaller star) fills its Roche-lobe.
	  // substituted for: 
	  // real da = semi*(1 - r_lobe) - smaller->get_effective_radius();
	  // (SPZ+GN: 1 Oct 1998)
	  
	  real sma_larger = lrc / roche_radius(1, larger->get_core_mass(),
					       smaller->get_total_mass());
	  real sma_smaller = ser / roche_radius(1, smaller->get_total_mass(),
						larger->get_core_mass());
	  real sma_post_ce = Starlab::max(sma_larger, sma_smaller);

      //SilT + GN October 2009 in new tracks sma_post_ce can be larger
      // than semi!!
      sma_post_ce = Starlab::min(sma_post_ce, semi);

	  real mass_lost = (mtot_l*mtot_s/(2*sma_post_ce) 
                             - mtot_l*mtot_s/(2*semi))
			 / (mtot_l/(alpha_lambda*r_lobe*semi) 
                             + mtot_s/(2*sma_post_ce));

	  //real da = semi - sma_post_ce;
	  //real semi_fr = 1-da/semi;
	  //	  real mass_lost = (1 - semi_fr)
	  //           / (1./mtot_l 
	  //           + 2*semi_fr/(mtot_s
	  //	 *   cnsts.parameters(common_envelope_efficiency)
	  //         *   cnsts.parameters(envelope_binding_energy)
	  //	 *   r_lobe));
           if (mass_lost>=menv_l) {
              mass_lost = 0.99*menv_l;
           }
	   
	   cerr << "Merger double_star::spiral_in()"<<endl;
	   dump(cerr, false);

	   //           larger->set_core_mass(mcl_old);
	   //           larger->set_envelope_mass(mel_old);
	   larger = larger->reduce_mass(mass_lost);
	   
	   // (GN Feb 2011) Always have larger be consumer
	   // else problems with MS + giant mergers
	   //if (mtot_l>mtot_s) 
	   merge_elements(larger, smaller);
           //else 
	   //  merge_elements(smaller, larger);
        } 
        else {
	  semi = a_spi;

	  // note that the accreted mass to the smaller 
	  // has not affected the spiral in. (SPZ:  2 Jun 1999)
	  menv_l -= smaller->add_mass_to_accretor(
                             larger->get_envelope_mass(),
                             larger->hydrogen_envelope_star(),
                             cnsts.parameters(spiral_in_time));
	  smaller->set_effective_radius(smaller->get_radius());
	  larger = larger->reduce_mass(larger->get_envelope_mass());

	  // redundant due to implementation of new force donor timescale
	  //                                         (SPZ+GN:23 Sep 1998)
	  //set_donor_timescale(smaller);
	  //determine_minimal_timestep();  
        }
     }

     //get_seba_counters()->spiral_in++;

//cerr<<"leave spiral_in: semi = "<<semi<<endl;
}

void double_star::merge_elements(star* consumer,
                                 star* dinner) {

//cerr<<"void double_star::merge_elements()"<<endl;
//cerr<<"bin: "<<identity<<endl;
  bin_type = Merged;
//  dump(cerr, false);
  dump("SeBa.data", true);

  if (!get_use_hdyn()) {
    consumer->set_velocity(velocity); 
    dinner->set_velocity(velocity);
    
// (GN+SPZ Apr 28 1999) star type can change (see single_star::merge_elements)
    consumer = consumer->merge_elements(dinner);
    semi = 0;
    dinner->set_envelope_mass(0);
    dinner->set_core_mass(0);

  }
  else {
    // solved in dstar_to_dyn merge_elements
  }

//  dump("binev.data", false);
cerr << "Merger is: "<<endl;
  dump(cerr, false);
  dump("SeBa.data", true);

  //get_seba_counters()->mergers++;
}

// This function is not used at the moment.
// But should be used.
// Mass loss must be performed for both

// stars simultaniously.
void double_star::perform_wind_loss(const real dt) {

  star* p = get_primary();
  star* s = get_secondary();
  p->stellar_wind(dt);
  s->stellar_wind(dt);
  p = s = NULL;

  //cerr<<"aMm'="<<semi<<" "<<get_primary()->get_total_mass()<<" "
  //                        <<get_secondary()->get_total_mass()<<endl;

}

void double_star::angular_momentum_loss(const real dt) {
//	Lose angular momentum.
//	First due to magnetic breaking then by gravitational 
//	wave radiation.

  if (bin_type != Merged && bin_type != Disrupted) {
      magnetic_stellar_wind(dt);
      // gravrad is called from within magnetic_stellar_wind
      // because binary can merge or spiral-in
  }
}

//	Lose angular momentum due to magnetic breaking.
//	Rappaport, Verbunt and Joss 1983.
//
//	Magnetic breaking takes place if the binary eccentricity
//	becomes smaller than COROTATION_E_LIMIT due to tidal circularization.
void double_star::magnetic_stellar_wind(const real dt) {
//cerr<<"void double_star::magnetic_stellar_wind(dt="<<dt<<)"<<endl;

    real magnetic_braking_aml = mb_angular_momentum_loss();

    real a_dot = 2*dt*magnetic_braking_aml;
    real semi_new = semi*(1 + a_dot);

    if (semi_new<=0) {

      cerr << "Merger::magnetic_stellar_wind"<<endl;
      dump(cerr, false);

      if (get_primary()->get_effective_radius() >
	  get_secondary()->get_effective_radius())
	merge_elements(get_primary(), get_secondary());
      else
	merge_elements(get_secondary(), get_primary());

      //get_seba_counters()->aml_mergers++;

    }
    else if(semi_new <= get_primary()->get_core_radius() +
			     get_secondary()->get_core_radius()) {
      semi = semi_new;
cerr << "magnetic stellar wind => ::dynamics_mas_transfer" << endl; 
      dynamic_mass_transfer();
      //      if (get_primary()->get_effective_radius() >
      //  get_secondary()->get_effective_radius())
      //	spiral_in(get_secondary(), get_primary());
      //      else
      //spiral_in(get_primary(), get_secondary());
    }
    else { 
      semi = semi_new;
      gravrad(dt);
    }
}

//	eccentricity variation due to the radiation
//	of gravitational waves.
//	Peters. P.C., 1964, {\em Phys. Rev.}, {\bf 136}, 1224.
real double_star::de_dt_gwr(const real e) {

    real de_dt = (-305/15)*e*(1 + (121/304)*e*e)/pow(1-e*e, 2.5);
    de_dt *= get_primary()->get_total_mass()
           * get_secondary()->get_total_mass()
           *  get_total_mass()/pow(semi, 4.);

        return de_dt;
}

//	Equations from:
//	Peters P.C., 1964, Phys. Ref., 136, 4B, B1224.
void double_star::gravrad(const real dt) {

     real G3M3_C5R4 = pow(cnsts.physics(G)
		    * cnsts.parameters(solar_mass), 3.)
                    / (pow(cnsts.physics(C), 5.)
		    * pow(cnsts.parameters(solar_radius), 4.));

//	Only appicable to small separation.
     if(semi/sqrt(get_total_mass())<=6.) {
        real e_new = eccentricity 
                   + G3M3_C5R4*de_dt_gwr(eccentricity)
	           * dt*cnsts.physics(Myear);

        real a_new;
        if (e_new > 0 && eccentricity > 0) {
           real c0 = (semi*(1-pow(eccentricity,2))/pow(eccentricity, 12./19.))
                   / pow(1 + 121*pow(eccentricity, 2.)/304, 870./2299.);
       
           a_new = (c0*pow(e_new, 12./19.)/(1-pow(e_new, 2.)))
                 * pow(1 + 121*pow(e_new, 2.)/304, 870./2299.);
        }
        else {
           e_new = 0; 		// Circularize
           semi *= (1-eccentricity*eccentricity);
           real G3M3_C5 = pow(cnsts.physics(G)
			      *cnsts.parameters(solar_mass), 3.)
	                / pow(cnsts.physics(C), 5.);
           real c0 = G3M3_C5* 256
                   * get_primary()->get_total_mass()
                   * get_secondary()->get_total_mass()
                   * get_total_mass()/5.;
           c0 = pow(semi*cnsts.parameters(solar_radius), 4.)
	      - c0*dt*cnsts.physics(Myear);
           c0 = Starlab::max(c0, 0.);
           a_new = pow(c0, 0.25)/cnsts.parameters(solar_radius);
        }
    
        if (a_new<=0) {

	  cerr << "Merger::gravrad"<<endl;
	  dump(cerr, false);

	  if (get_primary()->get_total_mass() >
	      get_secondary()->get_total_mass())
	    merge_elements(get_primary(), get_secondary());
	  else
	    merge_elements(get_secondary(), get_primary());

	  //get_seba_counters()->gwr_mergers++;
        }
        else if(a_new <= get_primary()->get_core_radius() +
			      get_secondary()->get_core_radius()) {
	  semi=a_new;
	  eccentricity = e_new;
	  cerr << "gravrad => ::dynamics_mass_transfer" << endl; 
	  dynamic_mass_transfer();
	  //	  if (get_primary()->get_effective_radius() >
	  //  get_secondary()->get_effective_radius())
	  //spiral_in(get_secondary(), get_primary());
	  //	  else
	  //spiral_in(get_primary(), get_secondary());
        }
        else {
	  semi=a_new;
	  eccentricity = e_new;
        }
     }
}


//      Mass lost from donor due to loss of angular momentum by
//      magnetic braking and gravitational wave-radiation.
// New implementations (SilT+GN: April 2012)
// New implementations (SilT+GN: Dec 2012)
real double_star::mdot_according_to_roche_radius_change(star* donor,
							star* accretor) {

  real mtr = 0;
  if (bin_type != Merged && bin_type != Disrupted) {

    real m_don_rel = donor->get_relative_mass();
    real m_acc = accretor->get_total_mass();   
    real m_don = donor->get_total_mass();

    real J_mb = mb_angular_momentum_loss();
    real J_gwr = gwr_angular_momentum_loss(m_don, m_acc, semi);
    real mtt = 1./abs(J_mb + J_gwr);
    real mtt_nuc = donor->nucleair_evolution_timescale();
    
    if (J_mb+J_gwr == 0 | mtt_nuc < mtt){
        return m_don / mtt_nuc / 2.; // so that donor_timescale = nuclear_timescale        
    
    }
    
    real md_dot = Starlab::min(donor->get_total_mass(), 
    			 cnsts.safety(minimum_mass_step));            
    real dt = md_dot * mtt/donor->get_total_mass();
    real ma_dot = accretor->accretion_limit(md_dot, dt);
    real eta = ma_dot/md_dot; 
    real z_star = donor->zeta_adiabatic();
    real zeta_l = zeta_scaled(donor, accretor, mtt);

    // Note: dm must be negative, because we talk about the mass losing star
    if (!accretor->remnant()) {

        // General case: semi-conservative mass transfer.
       real c_Jloss = cnsts.parameters(specific_angular_momentum_loss);
       mtr = m_don*abs(2*(J_mb+J_gwr))
    	 / (2 - 2*eta*m_don/m_acc - (1-eta)*(1+2*c_Jloss)*m_don/(m_don+m_acc) + z_star - zeta_l);
    }
    else {
        
        // Special case: mass transfer to compact object as accretor.
        // Only taken into account the posibility
        // 1) eta>0: mass lost as wind from accretor.
       mtr = m_don*abs(2*(J_mb+J_gwr))
    	 / (2 - 2*m_don/m_acc - (1-eta)*m_don/(m_don+m_acc) + z_star - zeta_l);
    }
  }
  return mtr;
}


//	Compute radial orbital velocities a binary components.
//	neglecting the eccenticity.
void double_star::calculate_velocities() {

     if(bin_type != Merged && bin_type != Disrupted) {
        real mu = get_primary()->get_total_mass()/get_total_mass();
        real mean_r = semi*(1+0.5*eccentricity*eccentricity);
        real v_rel = sqrt(cnsts.physics(G)*get_total_mass()
			  *cnsts.parameters(solar_mass)
			  *((2/(mean_r*cnsts.parameters(solar_radius)))
			 - (1/(semi*cnsts.parameters(solar_radius)))));
        v_rel /= cnsts.physics(km_per_s);
        get_primary()->set_velocity((1-mu)*v_rel);
        get_secondary()->set_velocity(mu*v_rel);
     }
}

//	Alternative of the double_star::contact();
//	No mass transfer occurs during stable case.
//void double_star::common_envelope(real dt) {
void double_star::contact_binary(real dt) {

    if (REPORT_FUNCTION_NAMES)
      cerr << "::contact_binary(dt = "<<dt << ")" << endl;

// (SilT Nov 25 2012) new dumping regime
    if (!first_contact) {
    
        bin_type = Contact;
        dump("SeBa.data", true);
        first_contact = true;
    }



// for the moment: only [ms,ms] stable, but what about [bd,ms] and [he,ms] etc.
    if (get_primary()->get_element_type() == Main_Sequence && 
	get_secondary()->get_element_type() == Main_Sequence) {

      // 'stable' contact binaries:
      // double_star::contact transferres matter from secondary to
      // primary: results in flip-flop, which rejuvenates both stars
      // (This might be not so bad, PPE private communication)
      // For now: do nothing
	
//        if (primary->get_effective_radius()>
//	      secondary->get_effective_radius()) {
//           if (ready_for_mass_transfer(secondary))
//              contact(secondary, primary, dt);
//        }
//        else {
//           if (ready_for_mass_transfer(primary)) 
//              contact(primary, secondary, dt);
//        }

	//get_seba_counters()->contact++;

    } else {
      current_mass_transfer_type = Dynamic;
      dynamic_mass_transfer();
    }

}

// (SPZ+GN: 28 Jul 2000)
// Dynamic mass transfer as in Eq. 5 of Nelemans VYPZ 2000 A&A in press

void double_star::angular_momentum_envelope_ejection(star* larger, 
						     star* smaller) {
  
//  cerr << "double_star::angular_momentum_envelope_ejection()"<<endl;
  // (GN+SilT Mar  2 2011) new dumping regime
  bin_type = Common_Envelope;
  dump("SeBa.data", true);
  current_mass_transfer_type = Unknown;


  if (bin_type!=Merged && bin_type!=Disrupted) {
 
	real M_core = larger->get_core_mass();
	real M_env  = larger->get_envelope_mass();
	real M_tot  = M_env + M_core;

	// At the moment no mass is accreted dureing dynamic mass tansfer
	// May be changed easily.
	real m_acc = 0;
	real m_new  = smaller->get_total_mass() - m_acc;
	M_env -= m_acc;
	
// (GN Sep 22 1999) ang.mom balance: \Delta J = \gamma * J * \Delta M / M
// See Eq. 5 of Nelemans VYPZ 2000 A&A
        real gamma = cnsts.parameters(dynamic_mass_transfer_gamma);

	real J_i = sqrt(semi) * (M_tot  * m_new)/sqrt(M_tot  + m_new);
	real J_f_over_sqrt_af = (M_core * m_new)/sqrt(M_core + m_new);
	real J_lost           = gamma * M_env * J_i/ (M_tot  + m_new);
	real sqrt_a_f         = Starlab::max(0.,(J_i - J_lost)/J_f_over_sqrt_af);
	real a_f              = pow(sqrt_a_f, 2);

	if(REPORT_BINARY_EVOLUTION) {

	  PRL(pow(M_tot/M_core, 2));
	  PRL((M_tot + m_new)/(M_core + m_new));
	  PRL(exp(-2. * M_env/m_new));
	  PRL(a_f);

	}

	  real rl_l = roche_radius(a_f, M_core, m_new);
	  real rl_s = roche_radius(a_f, m_new, M_core);


	// Merger
	// Note that we hare want the equilibrium radius of the accretor
	// The effective radius of the accretor is already used to
	// determine how conservative mass transfer is. (SPZ: 2 Jun 1999)
	  if (larger->get_core_radius()>=rl_l 	||
	      smaller->get_radius()>=rl_s) 	{

	    //	    PRC(smaller->get_radius());PRC(rl_s);PRL(a_f);

	    cerr << "Merger double_star::angular_momentum_envelope_ejection(star* larger, star* smaller)" << endl;
	    dump(cerr, false);
	    merge_elements(larger,smaller);
	  }
	  else {


	    if(REPORT_BINARY_EVOLUTION) {

	      cerr << "semi before = " << semi << endl;
	      cerr << "semi after = " << a_f << endl;
	    
	      real Eorb1 = (M_tot * m_new)/(2 * semi);
	      real Ebind = (M_tot * (M_tot - M_core))/larger->get_radius();
	      real Eorb2 = (M_core * m_new) / (2*a_f);

	      PRC(Eorb1);PRC(Eorb2);PRL(Ebind);
	      cerr << "(Eorb2 - Eorb1)/Ebind = "
		   << (Eorb2 - Eorb1)/Ebind << endl;
	    }

	    semi = a_f; 
	    smaller->add_mass_to_accretor(larger->get_envelope_mass(),
					  larger->hydrogen_envelope_star(), 
					  cnsts.parameters(spiral_in_time));
	    smaller->set_effective_radius(smaller->get_radius());

	    larger = larger->reduce_mass(larger->get_envelope_mass());

	  }

	//get_seba_counters()->dynamic++;

  }

//  cerr <<"Leave ::dynamics mass transfer"<<endl;
}

// (SPZ+GN: 28 Jul 2000) 
// Old 'beta' mechanism, not used any more.
#if 0
// (GN Apr 22 1999) artifical experiment to make giant - main_sequence
// contact system much less dramatic spiral-in (no separation change!)
void double_star::dynamic_mass_transfer(star* larger, star* smaller) {
  
//  cerr << "double_star::dynamic_mass_transfer()" << endl;

  if (bin_type!=Merged && bin_type!=Disrupted) {
       
// beta is cnst or specific ang. mom of donor/accretor
	real beta = cnsts.parameters(specific_angular_momentum_loss);
// donor
//	real beta = mtot_s/mtot_l;
// accretor
//	real beta = mtot_l/mtot_s;

//	  real J_i = sqrt(semi)*(mtot_l*mtot_s)/sqrt(mtot_l + mtot_s);
//	  real J_f_over_sqrt_af = (mcore_l*mtot_s)/sqrt(mcore_l + mtot_s);
//	  real J_lost = beta * menv_l * J_i / (mtot_l + mtot_s);
//	  real sqrt_a_f = max(0.,(J_i - J_lost)/J_f_over_sqrt_af);

//	cerr <<"larger " << endl;
//	larger->dump(cerr, false);
//	cerr<< endl;
//	cerr <<"smaller " << endl;
//	smaller->dump(cerr, false);
//	  cerr << endl;

 
	real M_core = larger->get_core_mass();
	real M_env  = larger->get_envelope_mass();
	real M_tot  = M_env + M_core;

	// What is accreted by the accretor should be reduce in the donor.
	//m_acc -= smaller->accretion_limit(M_env, 
        //                  cnsts.parameters(spiral_in_time));
	real m_acc = 0;
	real m_new  = smaller->get_total_mass() - m_acc;
	M_env -= m_acc;
	
	  real a_f = semi*pow(M_tot/M_core,2)
                         *pow((M_core + m_new)/(M_tot + m_new), 
			      2 * beta + 1);
  
	  real Eorb1 = (M_tot * m_new)/(2 * semi);
	  real Ebind = (M_tot * (M_tot - M_core))/larger->get_radius();

	  real rl_l = roche_radius(a_f, M_core, m_new);
	  real rl_s = roche_radius(a_f, m_new, M_core);


	// Merger
	// Note that we hare want the equilibrium radius of the accretor
	// The effective radius of the accretor is already used to
	// determine how conservative mass transfer is. (SPZ: 2 Jun 1999)
	  if (larger->get_core_radius()>=rl_l 	||
	      smaller->get_radius()>=rl_s) 	{

	    PRC(smaller->get_radius());PRC(rl_s);PRL(a_f);

	    cerr << "Merger double_star::dynamic_mass_transfer" << endl;
	    dump(cerr, false);
	    merge_elements(larger,smaller);
	  }
	  else {

	    cerr << "semi before = " << semi << endl;
	    semi = a_f; 
	    cerr << "semi after = " << semi << endl;
	    
	    real Eorb2 = (M_core * m_new) / (2*semi);

	    PRC(Eorb1);PRC(Eorb2);PRL(Ebind);
	    cerr << "(Eorb2 - Eorb1)/Ebind = " << (Eorb2 - Eorb1)/Ebind << endl;
	    smaller->add_mass_to_accretor(larger->get_envelope_mass(),
					  larger->hydrogen_envelope_star(), 
					  cnsts.parameters(spiral_in_time));
	    smaller->set_effective_radius(smaller->get_radius());

	    larger = larger->reduce_mass(larger->get_envelope_mass());
            // larger->set_effective_radius(larger->get_radius());
	  }

	//get_seba_counters()->dynamic++;

  }

//  cerr <<"Leave ::dynamics mass transfer"<<endl;
}
#endif

bool double_star::low_mass_star() {
  
    return (get_total_mass()<=cnsts.parameters(low_mass_star_mass_limit))
           ?true: false;
}

bool double_star::medium_mass_star() {
  
  return (!low_mass_star() && !high_mass_star())?true:false;
}

bool double_star::high_mass_star() {
  
  return (get_total_mass()>cnsts.parameters(medium_mass_star_mass_limit))
         ?true: false;
}


void double_star::refresh_memory() {

     star *p = get_primary();
     star *s = get_secondary();
     get_primary()->refresh_memory();
     get_secondary()->refresh_memory();
     p=s=NULL;

     previous.binary_age      = binary_age;
     previous.semi            = semi;
     previous.eccentricity    = eccentricity;

     previous.donor_timescale = donor_timescale;

}

void double_star::recall_memory() {

     star *p = get_primary();
     star *s = get_secondary();
     p->recall_memory();
     s->recall_memory();
     p=s=NULL;

     binary_age      = previous.binary_age;
     semi            = previous.semi;
     eccentricity    = previous.eccentricity;

     donor_timescale = previous.donor_timescale;
}


// (SilT+GN Feb 2011) new version
// Call of zeta explicitly gives the donor timescale as parameter
// does not include aml losses -> otherwise donor_time_scale becomes thermal instead of aml
real double_star::zeta(star * donor, 
                       star * accretor,
		              const real md_timescale) {

//	Find zeta Roche-lobe.
     if (bin_type!=Merged && bin_type!=Disrupted) {


       // Use total mass here as md_dot is only used to determine zeta
       real md_dot = Starlab::min(donor->get_total_mass(), 
       			 cnsts.safety(minimum_mass_step));
       //       real md_dot = min(donor->get_envelope_mass(), 
       //			 cnsts.safety(minimum_mass_step));
              
       real dt = md_dot * md_timescale/donor->get_relative_mass();
       real ma_dot = accretor->accretion_limit(md_dot, dt);

       real M_old = get_total_mass();
       real old_donor_mass = donor->get_total_mass();
       real old_accretor_mass = accretor->get_total_mass();
       real q_old = old_accretor_mass/old_donor_mass;

       real M_new = get_total_mass() - md_dot +  ma_dot;
       real new_donor_mass = donor->get_total_mass() - md_dot;
       real new_accretor_mass = accretor->get_total_mass() + ma_dot;
       real q_new = new_accretor_mass/new_donor_mass;

       real a_fr, new_semi;
       real beta = cnsts.parameters(specific_angular_momentum_loss);

       if (!accretor->remnant()) {

	 // General case: semi-conservative mass transfer.
	 a_fr  = pow(old_donor_mass*old_accretor_mass
		  / (new_donor_mass*new_accretor_mass), 2);
	 new_semi = semi * pow(M_new/M_old, 2*beta + 1) * a_fr;	
       }
       else {

	 // Special case: mass transfer to compact object as accretor.
	 //               Two possibilities:
	 //               1) eta>0: mass lost as wind from accretor.
	 //               2) eta==0: exponential spiral-in.
	 
	 real eta = ma_dot/md_dot; 

	 if (eta>0) {

	   a_fr  = (new_donor_mass/old_donor_mass)
	         * pow(new_accretor_mass/old_accretor_mass, 1/eta);
	   new_semi = semi * (M_old/M_new)/pow(a_fr, 2); 
	 }
	 else {

	   a_fr  = exp(2*(new_donor_mass-old_donor_mass)
		       /new_accretor_mass); 
	   new_semi = semi * (M_old/M_new)*a_fr
	            / pow(new_donor_mass/old_donor_mass, 2);
	 }
       }
       
//        real a_dot=0;
//
//        real magnetic_braking_aml = mb_angular_momentum_loss();
//        real grav_rad_aml =  gwr_angular_momentum_loss(get_primary()
//						      ->get_total_mass(),
//						      get_secondary()
//						      ->get_total_mass(),
//						      semi);
//
//        //       PRC(magnetic_braking_aml);PRL(grav_rad_aml);
//        //       PRC(dt);PRC(semi);
//        a_dot = 2*dt*semi*(magnetic_braking_aml+grav_rad_aml);
//        //       PRC(a_dot);
//        new_semi += a_dot;
//        //PRL(a_dot);
       
       real rl = roche_radius(donor);

       real rl_d = roche_radius(new_semi,
				new_donor_mass,
				new_accretor_mass);

//       PRC(new_semi);PRC(rl_d);PRL(rl);
       real d_lnr = (rl_d - rl)/rl;
       real d_lnm = (new_donor_mass - donor->get_total_mass()) 
	          /  donor->get_total_mass();
       
       real zeta;
       if(d_lnm==0) {
	 cerr << "WARNING: d_lnm (= " << d_lnm << ") has an illegal value"
	      << endl;
	 zeta = 0;
	 dump(cerr, true);
       }
       else {
	 zeta = d_lnr/d_lnm;
       }


       return zeta;
     }

     // Merged or disrupted.
     return 1;
}

// (SilT 22 Mar 2012) 
// return dln(Rl/a)/dlnM in stead of dlnRl/dlnM
// used in mdot_according_to_roche_radius_change
// includes aml losses
real double_star::zeta_scaled(star * donor, 
                       star * accretor,
		              const real md_timescale) {


//	Find zeta Roche-lobe.
     if (bin_type!=Merged && bin_type!=Disrupted) {


       // Use total mass here as md_dot is only used to determine zeta
       real md_dot = Starlab::min(donor->get_total_mass(), 
       			 cnsts.safety(minimum_mass_step));
       //       real md_dot = min(donor->get_envelope_mass(), 
       //			 cnsts.safety(minimum_mass_step));
              
       real dt = md_dot * md_timescale/donor->get_relative_mass();
       real ma_dot = accretor->accretion_limit(md_dot, dt);

       real M_old = get_total_mass();
       real old_donor_mass = donor->get_total_mass();
       real old_accretor_mass = accretor->get_total_mass();
       real q_old = old_accretor_mass/old_donor_mass;

       real M_new = get_total_mass() - md_dot +  ma_dot;
       real new_donor_mass = donor->get_total_mass() - md_dot;
       real new_accretor_mass = accretor->get_total_mass() + ma_dot;
       real q_new = new_accretor_mass/new_donor_mass;

       real a_fr, new_semi;
       real beta = cnsts.parameters(specific_angular_momentum_loss);

       if (!accretor->remnant()) {

	 // General case: semi-conservative mass transfer.
	 a_fr  = pow(old_donor_mass*old_accretor_mass
		  / (new_donor_mass*new_accretor_mass), 2);
	 new_semi = semi * pow(M_new/M_old, 2*beta + 1) * a_fr;	
       }
       else {

	 // Special case: mass transfer to compact object as accretor.
	 //               Two possibilities:
	 //               1) eta>0: mass lost as wind from accretor.
	 //               2) eta==0: exponential spiral-in.
	 
	 real eta = ma_dot/md_dot; 

	 if (eta>0) {

	   a_fr  = (new_donor_mass/old_donor_mass)
	         * pow(new_accretor_mass/old_accretor_mass, 1/eta);
	   new_semi = semi * (M_old/M_new)/pow(a_fr, 2); 
	 }
	 else {

	   a_fr  = exp(2*(new_donor_mass-old_donor_mass)
		       /new_accretor_mass); 
	   new_semi = semi * (M_old/M_new)*a_fr
	            / pow(new_donor_mass/old_donor_mass, 2);
	 }
       }
       
        real a_dot=0;

        real magnetic_braking_aml = mb_angular_momentum_loss();
        real grav_rad_aml =  gwr_angular_momentum_loss(get_primary()
						      ->get_total_mass(),
						      get_secondary()
						      ->get_total_mass(),
						      semi);

        //       PRC(magnetic_braking_aml);PRL(grav_rad_aml);
        //       PRC(dt);PRC(semi);
        a_dot = 2*dt*semi*(magnetic_braking_aml+grav_rad_aml);
        //       PRC(a_dot);
        new_semi += a_dot;
        //PRL(a_dot);
       
       real rl = roche_radius(donor);

       real rl_d = roche_radius(new_semi,
				new_donor_mass,
				new_accretor_mass);

//       PRC(new_semi);PRC(rl_d);PRL(rl);
//       real d_lnr = (rl_d - rl)/rl;
       real d_lnr_scaled = (rl_d/new_semi - rl/semi)/(rl/semi);
       real d_lnm = (new_donor_mass - donor->get_total_mass()) 
	          /  donor->get_total_mass();
       
       real zeta;
       if(d_lnm==0) {
	 cerr << "WARNING: d_lnm (= " << d_lnm << ") has an illegal value"
	      << endl;
	 zeta = 0;
	 dump(cerr, true);
       }
       else {
//	 zeta = d_lnr/d_lnm;
    zeta = d_lnr_scaled/d_lnm;
       }


       return zeta;
     }

     // Merged or disrupted.
     return 1;

}

     
void double_star::enhance_cluster_profile(cluster_profile& c_prof) {

      double_profile binary;
      //make_profile(identity, initial.start_time, binary, initial);
      binary.enhance_double_profile(this);
      c_prof.enhance_cluster_profile(binary);      
   }

// Used by determine_dt
real double_star::orbital_timescale() {

  real magnetic_braking_aml = mb_angular_momentum_loss();
  real grav_rad_aml =  gwr_angular_momentum_loss(get_primary()
						 ->get_total_mass(),
						 get_secondary()
						 ->get_total_mass(),
						 semi);
  real dt = cnsts.safety(maximum_timestep);

  real aml_loss = abs(magnetic_braking_aml)
                + abs(grav_rad_aml);
  // about 10 steps to lose ang. momentum.

  if (aml_loss> 0) 
    dt=.1/aml_loss;

  return dt;
}

// Peters & Mathews, 1963/4
real double_star::gwr_angular_momentum_loss(const real m_prim,
					    const real m_sec,
					    const real sma) {

  real m_tot = m_prim + m_sec;
  real J_gwr = 0;
  
  if(sma/sqrt(m_tot)<=6.) {
  
    real c_gwr = -32*pow(cnsts.parameters(solar_mass)*
			cnsts.physics(G), 3)
               / (5*pow(cnsts.physics(C), 5)*
		 pow(cnsts.parameters(solar_radius), 4));

    J_gwr = c_gwr*m_prim*m_sec*m_tot/pow(sma, 4);  

  }
  
  return J_gwr*cnsts.physics(Myear);
}

// Rappaport, Verbunt & Joss, 1983
real double_star::mb_angular_momentum_loss() {
    
  real c_mb = -3.8e-30*cnsts.parameters(solar_mass)*cnsts.physics(G)
            /          cnsts.parameters(solar_radius);

  real J_mb_prim = 0;
  if (get_primary()->magnetic() && eccentricity<=
      cnsts.parameters(corotation_eccentricity))
    J_mb_prim = c_mb * pow(get_total_mass(), 2)
              * pow(get_primary()->get_effective_radius(),
	            cnsts.parameters(magnetic_braking_exponent))
              / (get_secondary()->get_total_mass()*pow(semi, 5));
  
  real J_mb_sec = 0;
  if (get_secondary()->magnetic() && eccentricity<=
      cnsts.parameters(corotation_eccentricity))
        J_mb_sec = c_mb * pow(get_total_mass(), 2)
                 * pow(get_secondary()->get_effective_radius(),
		       cnsts.parameters(magnetic_braking_exponent))
                 / (get_primary()->get_total_mass()*pow(semi, 5));

  real J_mb=(J_mb_prim+J_mb_sec)*cnsts.physics(Myear);

  return J_mb;

}

//		The circularization radius for the accretion stream 
//		for star with mass m1 is fiven by 
//		(Shore, S., Livio, M., vanden Heuvel EPJ., saas-Fee
//		Advanced Course 22 for Astronomy and Astrophysics, 
//		Interacting Binaries p145):
real double_star::circularization_radius(const real m1, const real m2) {
  
     real q = m2/m1;
     real r_c = semi*(1+q)*pow(0.5 - 0.227 * log10(q), 4);
     
     return r_c;
}


real double_star::get_period() {

     real p = 2*PI*sqrt(pow(semi*cnsts.parameters(solar_radius), 3.)
			/ (cnsts.physics(G)*get_total_mass()
			   *cnsts.parameters(solar_mass)))
            /  cnsts.physics(days);
	  
     if(p<=0) p=1;
	  
     return p; 
}


real double_star::potential_energy() {
  
     real GM2_R = cnsts.physics(G)*pow(cnsts.parameters(solar_mass), 2)
                / cnsts.parameters(solar_radius);
     real u = GM2_R*get_primary()->get_total_mass()
            * get_secondary()->get_total_mass()/semi;
     
     return -u;
}
real double_star::kinetic_energy() {
  
     real M_km_s = cnsts.parameters(solar_mass)
                 * pow(cnsts.physics(km_per_s), 2);
     real k = 0.5*M_km_s*get_total_mass()*velocity*velocity;
     
     return k;
}

real double_star::total_energy() {
  return kinetic_energy() + potential_energy();
}



